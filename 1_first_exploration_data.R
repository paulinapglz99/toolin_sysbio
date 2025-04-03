#Exploration of data
setwd("~/tooling_up_systems_bio/RawData/")
dirnames <- dir()

dirnames[1]

#
x<- vroom::vroom("/home/paulinapg/tooling_up_systems_bio/RawData/HowardNew/Fri Apr 19 2024 (P1).csv")

# Pivot data
x_long <- x %>%
  pivot_longer(cols = -dilution, names_to = "sample", values_to = "value")

#signal vs dilution
# Plot
ggplot(x_long, aes(x = dilution, y = value, color = sample)) +
  geom_line() +
  geom_point() +
  labs(x = "Dilution", y = "Signal", title = "Dilution vs Signal") +
  theme_minimal()

# Plot
ggplot(x_long, aes(x = dilution, y = value, color = sample)) +
  geom_line() +
  geom_point() +
  scale_x_log10() + #logaritmic
  labs(x = "Dilution", y = "Signal", title = "Dilution vs Signal") +
  theme_minimal()


# Definir colores para cada grupo
group_colors <- c(
  Control = "blue",
  Infected = "red",
  Reference = "green"
)

# Asignar colores según el prefijo de cada muestra
x_long <- x_long %>%
  mutate(group = case_when(
    grepl("^Control", sample) ~ "Control",
    grepl("^Infected", sample) ~ "Infected",
    grepl("^Reference", sample) ~ "Reference"
  ))

# Graficar
howard_new.p <- ggplot(x_long, aes(x = dilution, y = value, color = group, group = sample)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_color_manual(values = group_colors) +
  labs(x = "log10Dilution", y = "Signal", title = "Dilution vs Signal by Sample Group - Howard new", subtitle = "Fri Apr 19 2024 (P1)") +
  theme_minimal()
howard_new.p

#Howard old lol

p <-  vroom::vroom("/home/paulinapg/tooling_up_systems_bio/RawData/Howard/Wed Sep 13 2023 (P1).csv")

# Pivot data
p_long <- p %>%
  pivot_longer(cols = -dilution, names_to = "sample", values_to = "value")

p_long <- p_long %>%
  mutate(group = case_when(
    grepl("^Control", sample) ~ "Control",
    grepl("^Infected", sample) ~ "Infected",
    grepl("^Reference", sample) ~ "Reference"
  ))


# Graficar
howard.p <- ggplot(p_long, aes(x = dilution, y = value, color = group, group = sample)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_color_manual(values = group_colors) +
  labs(x = "log10Dilution", y = "Signal", title = "Dilution vs Signal by Sample Group- Howard", subtitle = "Wed Sep 13 2023 (P1)") +
  theme_minimal()
howard.p

png("howard_and_howardnew_dilutionsgrid.png", width = 20, height = 20, units = "cm", res = 300)
gridExtra::grid.arrange(howard.p, howard_new.p, top = "")
dev.off()

# Establecer el dispositivo gráfico PNG con la resolución deseada

# Generar y guardar el arreglo de gráficos
gridExtra::grid.arrange(
  howard.p, howard_new.p,
  top = textGrob("sg", gp = gpar(fontsize = 15, fontface = "bold"))
)

# Cerrar el dispositivo gráfico
dev.off()


#Scatterplot

x %>% ggplot(aes(x = dilution, y = Reference1)) +
  geom_point() +
  scale_x_log10() +
  ggtitle("Reference1") +
  xlab("Dilution") + ylab("Signal")

#Calculate standard deviation, mean, and plot 

x.x <- x %>% summarise(
  across(everything(), list(
    mean = ~ mean(.),
    sd = ~ sd(.),
    mean_plus_1sd = ~ mean(.) + sd(.),
    mean_plus_2sd = ~ mean(.) + 2 * sd(.),
    mean_plus_3sd = ~ mean(.) + 3 * sd(.)
  ))
) 

# Cargar librerías
library(tidyverse)

# Leer los datos (asegúrate de cambiar la ruta del archivo)
datos <- read.csv("ruta/del/archivo/datos.csv")

# Seleccionar solo las columnas deseadas para calcular estadísticas
datos_seleccion <- datos %>% select(Reference1, Reference2, Control001)

# Calcular las estadísticas
estadisticas <- x %>%
  summarise(
    across(everything(), list(
      mean = ~ mean(.),
      sd = ~ sd(.),
      mean_plus_1sd = ~ mean(.) + sd(.),
      mean_plus_2sd = ~ mean(.) + 2 * sd(.),
      mean_plus_3sd = ~ mean(.) + 3 * sd(.)
    ))
  ) %>%
  pivot_longer(cols = everything(), names_to = c("Variable", "Statistic"), names_sep = "_") %>%
  pivot_wider(names_from = "Statistic", values_from = "value")

# Transformar los datos al formato largo para el gráfico
x_long <- x %>%
  pivot_longer(cols = -dilution, names_to = "sample", values_to = "value")

# Graficar
ggplot(x_long, aes(x = dilution, y = value, color = sample)) +
  geom_line() +
  geom_point() +
  scale_x_log10() + # Opcional, si quieres que el eje x sea logarítmico
  labs(x = "Dilution", y = "Value", title = "Plot de todas las muestras vs Dilution") +
  theme_minimal() +
  # Agregar líneas horizontales para cada estadística
  geom_hline(data = estadisticas %>% filter(Variable == "Reference1"),
             aes(yintercept = mean, linetype = "Reference1 Mean"), color = "blue") +
  geom_hline(data = estadisticas %>% filter(Variable == "Reference1"),
             aes(yintercept = mean_plus_1sd, linetype = "Reference1 Mean + 1 SD"), color = "blue", linetype = "dashed") +
  geom_hline(data = estadisticas %>% filter(Variable == "Reference1"),
             aes(yintercept = mean_plus_2sd, linetype = "Reference1 Mean + 2 SD"), color = "blue", linetype = "dotted")
             