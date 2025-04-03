# Cargar las librerías necesarias
library(vroom)
library(ggplot2)
library(dplyr)
library(minpack.lm) # Para la función nlsLM de ajuste no lineal

# Cargar el archivo CSV
csv1 <- vroom("~/tooling_up_systems_bio/ToolingSystemsBiology/RawData/Demengeot/2023-01-11 plate_1.csv")

# Definir la función de Michaelis-Menten
michaelis_menten <- function(dilution, Vmax, Km, c) {
  ((Vmax * dilution) / (Km + dilution)) + c
}

# Definir la función de ajuste con arcTangente
arctangent <- function(dilution, a, b, c) {
 (c-1) + (a/2) * atan((x - b /a))
}

# Ajuste no lineal de Michaelis-Menten
fit_mm <- nlsLM(csv1$Reference1 ~ michaelis_menten(dilution, Vmax, Km),
                data = csv1,
                start = list(Vmax = max(csv1$Reference1), Km = median(csv1$dilution)))

# Ajuste no lineal con arcTangente
fit_atan <- nlsLM(csv1$Reference1 ~ arctangent(dilution, a, b),
                  data = csv1,
                  start = list(a = max(csv1$Reference1), b = 1 / median(csv1$dilution)))

# Extraer los parámetros ajustados
Vmax <- coef(fit_mm)["Vmax"]
Km <- coef(fit_mm)["Km"]
a <- coef(fit_atan)["a"]
b <- coef(fit_atan)["b"]

# Crear la gráfica con los dos ajustes
ggplot(csv1, aes(x = dilution, y = Reference1)) +
  geom_point(color = "blue") + 
  scale_x_log10() +
  labs(title = "Michaelis-Menten and ArcTangente Fit for Plate1$Reference1",
       x = "Log10(Dilution)",
       y = "Reference01") +
  expand_limits(x = 0, y = 0) +
  theme_minimal() +
  stat_function(fun = function(x) (Vmax * x) / (Km + x), color = "red", linetype = "solid") +  # Curva de Michaelis-Menten
  stat_function(fun = function(x) a * atan(b * x), color = "green", linetype = "dashed")       # Curva de ajuste arcTangente

