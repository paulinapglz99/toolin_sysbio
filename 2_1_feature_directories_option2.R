setwd("~/Tooling_up/RawData/")

normalized_dictionary <- vroom::vroom(file = "~/ToolingSystemsBiology/normalized_dictionary.csv")

###
library(dplyr)

# Función para leer el archivo, extraer la última fila y calcular el promedio y la desviación estándar
get_last_row_stats <- function(file_path) {
  # Leer el archivo CSV
  data <- read.csv(file_path)
  
  # Extraer la última fila
  last_row <- tail(data, 1)
  
  # Ignorar la primera columna si es de identificación (como "dilution" o similar)
  values <- unlist(last_row[,-9])  # Suponiendo que la primera columna es un identificador
  
  # Calcular el promedio y la desviación estándar
  avg <- mean(values, na.rm = TRUE)
  sd <- sd(values, na.rm = TRUE)
  
  return(list(average = avg, std_dev = sd))
}

# Aplicar la función a cada archivo en la tabla normalized_dictionary
normalized_dictionary <- normalized_dictionary %>%
  rowwise() %>%  # Procesar fila por fila
  mutate(
    file_path = paste0("~/Tooling_up/RawData/", Lab, "/", Title, ".csv"),  # Construir la ruta completa del archivo
    stats = list(get_last_row_stats(file_path))  # Obtener las estadísticas
  ) %>%
  mutate(
    average = stats$average,
    std_dev = stats$std_dev
  ) %>%
  select(-stats, -file_path)  # Eliminar las columnas temporales


# Agregar una nueva columna con la suma de (average + 2 * std_dev)
normalized_dictionary <- normalized_dictionary %>%
  mutate(threshold = average + 2 * std_dev)

# Ver los primeros resultados

write.csv(normalized_dictionary, "~/ToolingSystemsBiology/data_option2.csv")


