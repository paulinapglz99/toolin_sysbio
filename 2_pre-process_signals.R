#Script for extract features of the csv

setwd("~/tooling_up_systems_bio/ToolingSystemsBiology/")

#Get directory in a string --- --- 
dir <- getwd()

#Libraries --- ---
pacman::p_load("dplyr", "purrr", "stringr")

#Functions --- --- 

# Function to extract information from each file

extract_info <- function(filepath, lab_name) {
  file_name <- basename(filepath) #Get file_name
  file_name <- tools::file_path_sans_ext(file_name) # Extract the file identifier (file name without the extension)
  date <- str_extract(file_name, "\\d{4}-\\d{2}-\\d{2}")  # Get dates in YYYY-MM-DD
  if (is.na(date)) {
    date <- str_extract(file_name, "[A-Za-z]{3}\\s[A-Za-z]{3}\\s\\d{2}\\s\\d{4}")  # formato como 'Fri Aug 18 2023'
  }
  # Extract plate names (P1, P2, plate_1, plate_2, etc.)
  plate <- str_extract(file_name, "P\\d+|plate_\\d+")
  # Create a dataframe
  return(data.frame(
    file_name = file_name,
    Date = date,
    Lab = lab_name,
    Plate = plate,
    #FilePath = filepath,
    stringsAsFactors = FALSE
  ))
}

# Function to calculate mean and standard deviation of the last row (excluding “dilution”)
calculate_last_row_stats <- function(filepath) {
  data <- read.csv(filepath) #Get data
  #Get last row without the dilution column
  last_row <- data[nrow(data), -ncol(data)]
  #Calculate mean and sd
  mean_value <- mean(as.numeric(last_row), na.rm = TRUE)
  sd_value <- sd(as.numeric(last_row), na.rm = TRUE)
  
  return(list(mean = mean_value, sd = sd_value))
}


#Function to obtain the value and dilution associated with each column exceeding the threshold.
extract_titer <- function(filepath, threshold) {
  # Leer el archivo CSV
  data <- read.csv(filepath)
  
  # Crear una lista para almacenar los resultados por columna
  result <- list()
  
  # Recorremos cada columna
  for (col_name in names(data)[1:8]) {
    # Extraer la columna de valores y la de diluciones
    values <- data[[col_name]]
    dilutions <- data$dilution
    
    # Identificar el primer valor que esté por encima del threshold
    above_threshold_idx <- max(which(values > threshold))
    
    if (!is.na(above_threshold_idx)) {
      # Guardamos el valor y su dilución asociada
      result[[paste0(col_name, "_value")]] <- values[above_threshold_idx]
      result[[paste0(col_name, "_dilution")]] <- dilutions[above_threshold_idx]
    } else {
      # Si no hay ningún valor superior al threshold, almacenamos NA
      result[[paste0(col_name, "_value")]] <- NA
      result[[paste0(col_name, "_dilution")]] <- NA
    }
    
  }
  result <- as.data.frame(result)
  new_cols <- c(
    "Reference1_val", "Reference1_titre",
    "Reference2_val", "Reference2_titre",
    "Control1_val", "Control1_titre",
    "Control2_val", "Control2_titre",
    "Infected1_val", "Infected1_titre",
    "Infected2_val", "Infected2_titre",
    "Infected3_val", "Infected3_titre",
    "Infected4_val", "Infected4_titre"
  )
  
  colnames(result) <- new_cols
  return(result)
}

# Process directories function
process_directory <- function(dir_path, lab_name) {
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  data <- do.call(rbind, lapply(files, extract_info, lab_name = lab_name))
  return(data)
}

#Pre-process data --- ---

# Define directories
lab_directories <- list(
  Demengeot = paste0(dir, "/RawData/Demengeot"),
  Howard = paste0(dir,"/RawData/Howard"),
  HowardNew = paste0(dir,"/RawData/HowardNew"),
  Vilanova = paste0(dir,"/RawData/Vilanova")
)

#Process directorise and combine results
combined_data <- do.call(rbind, lapply(names(lab_directories), function(lab) {
  process_directory(lab_directories[[lab]], lab_name = lab)
}))

# Check plates and dates 
unique(combined_data$Plate)
unique(combined_data$Date)

# Date and plate normalisation 
combined_data$Plate <- ifelse(grepl("^P[0-9]+$", combined_data$Plate),
                              paste0("plate_", sub("^P", "", combined_data$Plate)),
                              combined_data$Plate)

unique(combined_data$Plate)

# Convert dates to “YYYYY-MM-DD” format using different formats as appropriate.
combined_data$Date <- ifelse(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", combined_data$Date),
                             as.Date(combined_data$Date, format = "%Y-%m-%d"),
                             as.Date(combined_data$Date, format = "%a %b %d %Y"))

# Convert to date
combined_data$Date <- as.Date(combined_data$Date, origin = "1970-01-01")

#Check result
unique(combined_data$Date)

#Unique identifier 

combined_data$date_lab_plate <- paste(combined_data$Date, combined_data$Lab, combined_data$Plate, sep = "_")

combined_data <- combined_data %>% dplyr::select("file_name", "Date", "Lab", "Plate")

#Calculate mean and sd 
combined_data <- combined_data %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_last_row_stats(paste0(dir, "/RawData/", Lab, "/", file_name, ".csv"))),
    background_mean = stats$mean,
    background_SD = stats$sd
  ) %>%
  select(-stats) %>%
  ungroup()

# Add titre `threshold` column by calculating the value of 2 standard deviations over the background_mean

combined_data <- combined_data %>%
  mutate(threshold = background_mean + (2 * background_SD))

#Create a list of complete filepaths using the `Lab` column in `combined_data`

filepaths <- map2(
  combined_data$Lab, 
  combined_data$file_name, 
  ~ file.path(lab_directories[[.x]], paste0(.y, ".csv"))
)

#Apply the `extract_titer` function to each file with a corresponding threshold
results <- map2(filepaths, combined_data$threshold, ~ extract_titer(.x, .y))

#Combine the results in a single dataframe, adding the columns `file_name` and `threshold` of combined_data
final_data <- bind_rows(results, .id = "file_name") %>%
  bind_cols(combined_data)

#Use map2_dfr to apply the function and merge the results into a single dataframe
final_data <- map2_dfr(filepaths, combined_data$threshold, ~ {
  result <- extract_titer(.x, .y)
  result$file_name <- basename(.x) # Agregar file_name como columna
  result$threshold <- .y # Agregar threshold como columna
  return(result)
})

#Add additional combined_data columns
final_data <- final_data %>%
  left_join(combined_data, by = "threshold")

colnames(final_data)

final_data <- final_data %>%
  select(
    file_name.x, Date, Lab, Plate, background_mean, background_SD,
    threshold,
    everything()
  )

# Invert

#final_data <- final_data %>%
  # mutate(
  #   across(
  #     ends_with("_titre"),
  #     ~ 1 / .,  # calculates the inverse of each value
  #     .names = "inv_{.col}"  # Crea una nueva columna con el prefijo "inv_"
  #   )
  # )

#and make columns of average of references

# final_data <- final_data %>%
#   mutate(
#     Reference_average = rowMeans(select(., Reference1_titre, Reference2_titre), na.rm = TRUE)
#   )

#final_data$Reference_average <- (final_data$Reference1_titre + final_data$Reference2_titre) / 2

#Normalize data by dividing all inverted values by the average of references

# final_data <- final_data %>%
#   mutate(
#     across(
#       starts_with("inv_"),
#       ~ . / Reference_average,  # Divide each column inv_ by Reference_average
#       .names = "{.col}_norm"  # Create new columns with the suffix “_norm”.
#     )
#   )

# Identificar las columnas que empiezan con "inv_"
inv_columns <- grep("^inv_", names(final_data), value = TRUE)

# Crear las nuevas columnas dividiendo por Reference_average
for (col in inv_columns) {
  final_data[[paste0(col, "_norm")]] <- final_data[[col]] / final_data$Reference_average
}

######CUACK
# Paso 1: Pivotamos columnas terminadas en `_val` y `_titre` para obtener `val` y `titre`
final_data_long <- final_data %>%
  pivot_longer(
    cols = ends_with("_val") | ends_with("_titre"),  # Selecciona columnas que terminan en _val o _titre
    names_to = c("sample", ".value"),                  # "type" indica el tipo (Reference1, Control1, etc.)
    names_sep = "_"
  )

# Paso 3: Seleccionamos columnas en el orden deseado
final_data_long <- final_data_long %>%
  select(file_name.x, Date, Lab, Plate, background_mean, background_SD, threshold, 
         sample, val, titre)

final_data_long <- final_data_long %>% 
  filter(!is.na(val))

final_data_long$inv <- 1/final_data_long$titre


final_data_long <- final_data_long %>%
  mutate(group = rep(1:(n()/8), each = 8)) # Asignar un número de grupo cada 8 filas

# Calcular el promedio de `inv` para `Reference1` y `Reference2` dentro de cada grupo
reference_avg <- final_data_long %>%
  filter(sample %in% c("Reference1", "Reference2")) %>% # Filtrar solo las referencias
  group_by(group) %>%                                   # Agrupar por cada grupo de 8 filas
  summarise(reference_average = mean(inv, na.rm = TRUE)) # Calcular el promedio de `inv` en cada grupo

# Unir el promedio de referencia con el dataframe original
final_data_long <- final_data_long %>%
  left_join(reference_avg, by = "group") %>% # Unir usando la columna "group"
  select(-group)                             # Opcional: eliminar la columna "group" si ya no es necesaria

final_data_long <- final_data_long %>%
  mutate(norm = inv / reference_average)

#save table

vroom::vroom_write(final_data, "final_data.csv")

#END