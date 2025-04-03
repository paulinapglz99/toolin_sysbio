# Cargar las librerías necesarias
library(vroom)
library(ggplot2)
library(dplyr)
library(minpack.lm) # Para la función nlsLM de ajuste no lineal

# Cargar el archivo CSV
csv1 <- vroom("~/ToolingSystemsBiology/RawData/Demengeot/2023-01-11 plate_1.csv")

# Definir la función de Michaelis-Menten, con el parámetro c (ruido) sin el +1
michaelis_menten <- function(dilution, Vmax, Km, c) {
  c + (Vmax * dilution) / (Km + dilution)  # Sumar solo c, sin el +1
}

# Definir la función de ajuste con arcTangente, con el parámetro c (ruido) +1
arctangent <- function(dilution, a, b, c, d) {
  (c-1) + (a/2) * atan((dilution - b)/d)  # Agregar c - 1 al principio
}

# Ajuste no lineal de Michaelis-Menten, agregando el parámetro c
fit_mm <- nlsLM(csv1$Reference1 ~ michaelis_menten(dilution, Vmax, Km, c),
                data = csv1,
                start = list(Vmax = max(csv1$Reference1), Km = median(csv1$dilution), c = 2.77))

# Ajuste no lineal con arcTangente, agregando el parámetro c
fit_atan <- nlsLM(csv1$Reference1 ~ arctangent(dilution, a, b, c),
                  data = csv1,
                  start = list(a = max(csv1$Reference1), b = 1 / median(csv1$dilution), c = 2.77))

# Calcular valores ajustados (predichos) para cada modelo
csv1 <- csv1 %>%
  mutate(pred_mm = predict(fit_mm),     # Valores predichos por Michaelis-Menten
         pred_atan = predict(fit_atan), # Valores predichos por arcTangente
         resid_mm = Reference1 - pred_mm,    # Residuos para Michaelis-Menten
         resid_atan = Reference1 - pred_atan, # Residuos para arcTangente
         log_dilution = log10(dilution)) # Escala logarítmica de dilución

# Calcular R^2 para ambos modelos
r_squared_mm <- 1 - sum((csv1$Reference1 - csv1$pred_mm)^2) / sum((csv1$Reference1 - mean(csv1$Reference1))^2)
r_squared_atan <- 1 - sum((csv1$Reference1 - csv1$pred_atan)^2) / sum((csv1$Reference1 - mean(csv1$Reference1))^2)

# Crear la gráfica con los dos ajustes (Michaelis-Menten y arcTangente), incluyendo el R^2
ggplot(csv1, aes(x = dilution, y = Reference1)) +
  geom_point(color = "blue") + 
  scale_x_log10() +
  labs(title = "Michaelis-Menten y ArcTangente Ajuste con Ruido",
       x = "Log10(Dilution)",
       y = "Reference01 (Shifted by Noise)") +
  expand_limits(x = 0, y = 0) +
  theme_minimal() +
  # Curva de Michaelis-Menten con ruido y R^2
  stat_function(fun = function(x) coef(fit_mm)["c"] + (coef(fit_mm)["Vmax"] * x) / (coef(fit_mm)["Km"] + x), 
                color = "red", linetype = "solid") +
  annotate("text", x = max(csv1$dilution), y = max(csv1$Reference1), 
           label = paste("R² Michaelis-Menten =", round(r_squared_mm, 3)), color = "red", hjust = 1) +
  # Curva de arcTangente con ruido y R^2
  stat_function(fun = function(x) (coef(fit_atan)["c"] - 1) + coef(fit_atan)["a"] * atan(coef(fit_atan)["b"] * x), 
                color = "green", linetype = "dashed") +
  annotate("text", x = max(csv1$dilution), y = max(csv1$Reference1) * 0.9, 
           label = paste("R² ArcTangente =", round(r_squared_atan, 3)), color = "green", hjust = 1)

# Crear el gráfico de los residuos en función del log de la dilución
ggplot(csv1, aes(x = log_dilution)) +
  geom_point(aes(y = resid_mm, color = "Residuos Michaelis-Menten"), alpha = 0.6) +
  geom_point(aes(y = resid_atan, color = "Residuos ArcTangente"), alpha = 0.6) +
  labs(title = "Residuos de los Ajustes de Michaelis-Menten y ArcTangente",
       x = "Log10(Dilution)",
       y = "Residuos") +
  theme_minimal() +
  scale_color_manual(values = c("Residuos Michaelis-Menten" = "red", "Residuos ArcTangente" = "green")) +
  # Línea de regresión para los residuos de Michaelis-Menten con R^2
  geom_smooth(aes(y = resid_mm), method = "lm", color = "red", se = FALSE, linetype = "solid") +
  annotate("text", x = max(csv1$log_dilution), y = max(csv1$resid_mm), 
           label = paste("R² Residuos MM =", round(summary(lm(resid_mm ~ log_dilution, data = csv1))$r.squared, 3)), color = "red", hjust = 1) +
  # Línea de regresión para los residuos de arcTangente con R^2
  geom_smooth(aes(y = resid_atan), method = "lm", color = "green", se = FALSE, linetype = "dashed") +
  annotate("text", x = max(csv1$log_dilution), y = max(csv1$resid_atan) * 0.9, 
           label = paste("R² Residuos ArcTangente =", round(summary(lm(resid_atan ~ log_dilution, data = csv1))$r.squared, 3)), color = "green", hjust = 1)
