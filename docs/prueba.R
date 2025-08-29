# ==========================
# ANÁLISIS DETALLADO DE "OTRAS COMORBILIDADES"
# Script para identificar las enfermedades más prevalentes
# que NO son relacionadas con estilo de vida ni cánceres
# ==========================

# Asumiendo que ya tienes cargado data_2023_cancer del script anterior
# Si no, ejecuta primero el script principal

# ==========================
# 1) Extraer todos los códigos de "otras comorbilidades"
# ==========================

# Verificar que existen las variables necesarias
if(!exists("data_2023_cancer") | !exists("c3_cols") | !exists("all_lifestyle_codes")){
  stop("Primero debe ejecutar el script principal para cargar data_2023_cancer")
}

cat("========================================\n")
cat("ANÁLISIS DETALLADO DE OTRAS COMORBILIDADES\n")
cat("(No relacionadas con estilo de vida, excluyendo cánceres)\n")
cat("========================================\n\n")

# Recopilar TODOS los códigos de otras comorbilidades
other_codes_detailed <- list()
patient_codes <- list()  # Para guardar códigos por paciente

for(i in 1:nrow(data_2023_cancer)){
  row_codes <- unlist(data_2023_cancer[i, c3_cols])
  valid_codes <- row_codes[!is.na(row_codes) & row_codes != ""]
  
  # Filtrar códigos que NO son de estilo de vida Y NO son cánceres
  other_codes <- valid_codes[!valid_codes %in% all_lifestyle_codes & !str_starts(valid_codes, "C")]
  
  if(length(other_codes) > 0){
    other_codes_detailed[[length(other_codes_detailed) + 1]] <- other_codes
    patient_codes[[i]] <- other_codes
  } else {
    patient_codes[[i]] <- character(0)
  }
}

# ==========================
# 2) Análisis de frecuencias
# ==========================

all_other_codes <- unlist(other_codes_detailed)

if(length(all_other_codes) > 0){
  
  # Crear tabla de frecuencias
  other_codes_table <- table(all_other_codes) |>
    as.data.frame() |>
    rename(code = all_other_codes, frequency = Freq) |>
    mutate(
      description = icd3_label(as.character(code)),
      n_patients = frequency,
      prevalence_pct = round(frequency / nrow(data_2023_cancer) * 100, 2),
      pct_of_those_with_other = round(frequency / sum(data_2023_cancer$other_comorbidities) * 100, 2)
    ) |>
    arrange(desc(frequency))
  
  # ==========================
  # 3) TOP 50 códigos más frecuentes
  # ==========================
  cat("TOP 50 'OTRAS COMORBILIDADES' MÁS FRECUENTES\n")
  cat("================================================\n\n")
  
  top50_other <- other_codes_table |> head(50)
  
  # Imprimir con formato mejorado
  for(i in 1:nrow(top50_other)){
    cat(sprintf("%2d. %s - %s\n", 
                i, 
                top50_other$code[i], 
                top50_other$description[i]))
    cat(sprintf("    Pacientes: %d (%.1f%% del total, %.1f%% de los que tienen otras comorb.)\n\n",
                top50_other$n_patients[i],
                top50_other$prevalence_pct[i],
                top50_other$pct_of_those_with_other[i]))
  }
  
  # ==========================
  # 4) Agrupar por categorías CIE-10
  # ==========================
  cat("\n================================================\n")
  cat("AGRUPACIÓN POR CATEGORÍAS CIE-10 PRINCIPALES\n")
  cat("================================================\n\n")
  
  # Función para obtener la categoría principal del código CIE-10
  get_icd_category <- function(code){
    first_letter <- substr(code, 1, 1)
    code_num <- as.numeric(substr(code, 2, 3))
    
    category <- case_when(
      first_letter == "A" | first_letter == "B" ~ "A00-B99: Infectious diseases",
      first_letter == "D" & code_num >= 50 & code_num <= 89 ~ "D50-D89: Blood and immune disorders",
      first_letter == "D" & code_num < 50 ~ "D00-D49: Neoplasms (benign/uncertain)",
      first_letter == "E" & code_num >= 0 & code_num <= 35 ~ "E00-E35: Endocrine disorders",
      first_letter == "E" & code_num >= 70 & code_num <= 90 ~ "E70-E90: Metabolic disorders",
      first_letter == "F" ~ "F00-F99: Mental disorders",
      first_letter == "G" ~ "G00-G99: Nervous system diseases",
      first_letter == "H" & code_num <= 59 ~ "H00-H59: Eye diseases",
      first_letter == "H" & code_num >= 60 ~ "H60-H95: Ear diseases",
      first_letter == "I" ~ "I00-I99: Circulatory diseases",
      first_letter == "J" ~ "J00-J99: Respiratory diseases",
      first_letter == "K" ~ "K00-K93: Digestive diseases",
      first_letter == "L" ~ "L00-L99: Skin diseases",
      first_letter == "M" ~ "M00-M99: Musculoskeletal diseases",
      first_letter == "N" ~ "N00-N99: Genitourinary diseases",
      first_letter == "O" ~ "O00-O99: Pregnancy/childbirth",
      first_letter == "P" ~ "P00-P96: Perinatal conditions",
      first_letter == "Q" ~ "Q00-Q99: Congenital malformations",
      first_letter == "R" ~ "R00-R99: Symptoms and signs",
      first_letter == "S" | first_letter == "T" ~ "S00-T98: Injury and poisoning",
      first_letter == "V" | first_letter == "W" | first_letter == "X" | first_letter == "Y" ~ "V01-Y98: External causes",
      first_letter == "Z" ~ "Z00-Z99: Factors influencing health status",
      TRUE ~ "Other"
    )
    return(category)
  }
  
  # Agregar categoría a la tabla
  other_codes_table$icd_category <- sapply(other_codes_table$code, get_icd_category)
  
  # Resumen por categoría
  category_summary <- other_codes_table |>
    group_by(icd_category) |>
    summarise(
      n_unique_codes = n(),
      total_occurrences = sum(frequency),
      pct_of_all_other = round(sum(frequency) / sum(other_codes_table$frequency) * 100, 2)
    ) |>
    arrange(desc(total_occurrences))
  
  print(category_summary)
  
  # ==========================
  # 5) Top 10 códigos por categoría principal
  # ==========================
  cat("\n================================================\n")
  cat("TOP 5 CÓDIGOS MÁS FRECUENTES POR CATEGORÍA\n")
  cat("================================================\n\n")
  
  # Mostrar solo las categorías más relevantes
  main_categories <- category_summary$icd_category[1:min(8, nrow(category_summary))]
  
  for(cat_name in main_categories){
    cat(paste0("\n", cat_name, ":\n"))
    cat(paste0(rep("-", nchar(cat_name) + 1), collapse = ""), "\n")
    
    top_in_category <- other_codes_table |>
      filter(icd_category == cat_name) |>
      head(5)
    
    for(j in 1:nrow(top_in_category)){
      cat(sprintf("  %s - %s (%.1f%%)\n", 
                  top_in_category$code[j],
                  top_in_category$description[j],
                  top_in_category$prevalence_pct[j]))
    }
  }
  
  # ==========================
  # 6) Análisis por tipo de cáncer
  # ==========================
  cat("\n\n================================================\n")
  cat("PREVALENCIA POR TIPO DE CÁNCER (TOP 10 COMORBILIDADES)\n")
  cat("================================================\n\n")
  
  # Crear matriz de prevalencia por tipo de cáncer
  cancer_types <- unique(data_2023_cancer$cancer_site)
  top10_codes <- head(other_codes_table$code, 10)
  
  prevalence_by_cancer <- matrix(0, 
                                 nrow = length(cancer_types), 
                                 ncol = length(top10_codes))
  rownames(prevalence_by_cancer) <- cancer_types
  colnames(prevalence_by_cancer) <- top10_codes
  
  for(ct in cancer_types){
    cancer_subset <- data_2023_cancer[data_2023_cancer$cancer_site == ct, ]
    n_cancer <- nrow(cancer_subset)
    
    for(code in top10_codes){
      # Contar cuántos pacientes con este tipo de cáncer tienen este código
      count_with_code <- 0
      for(i in 1:nrow(cancer_subset)){
        row_codes <- unlist(cancer_subset[i, c3_cols])
        if(code %in% row_codes) count_with_code <- count_with_code + 1
      }
      prevalence_by_cancer[ct, code] <- round(count_with_code / n_cancer * 100, 1)
    }
  }
  
  # Imprimir la matriz
  prevalence_df <- as.data.frame(prevalence_by_cancer)
  
  # Agregar descripciones cortas
  cat("Códigos:\n")
  for(code in top10_codes){
    desc <- other_codes_table$description[other_codes_table$code == code][1]
    # Truncar descripción si es muy larga
    desc_short <- substr(desc, 1, 40)
    if(nchar(desc) > 40) desc_short <- paste0(desc_short, "...")
    cat(sprintf("%s: %s\n", code, desc_short))
  }
  
  cat("\nPrevalencia (%) por tipo de cáncer:\n")
  print(prevalence_df)
  
  # ==========================
  # 7) Análisis de co-ocurrencia
  # ==========================
  cat("\n\n================================================\n")
  cat("CO-OCURRENCIA DE LAS 10 COMORBILIDADES MÁS FRECUENTES\n")
  cat("================================================\n\n")
  
  # Crear matriz de co-ocurrencia
  cooc_matrix <- matrix(0, 
                        nrow = length(top10_codes), 
                        ncol = length(top10_codes))
  rownames(cooc_matrix) <- top10_codes
  colnames(cooc_matrix) <- top10_codes
  
  # Contar co-ocurrencias
  for(i in 1:length(patient_codes)){
    codes_patient <- patient_codes[[i]]
    codes_in_top10 <- intersect(codes_patient, top10_codes)
    
    if(length(codes_in_top10) > 1){
      for(j in 1:(length(codes_in_top10)-1)){
        for(k in (j+1):length(codes_in_top10)){
          code1 <- codes_in_top10[j]
          code2 <- codes_in_top10[k]
          cooc_matrix[code1, code2] <- cooc_matrix[code1, code2] + 1
          cooc_matrix[code2, code1] <- cooc_matrix[code2, code1] + 1
        }
      }
    }
  }
  
  # Convertir a porcentajes
  diag(cooc_matrix) <- other_codes_table$frequency[match(top10_codes, other_codes_table$code)]
  
  cat("Número de pacientes con ambas condiciones:\n")
  print(cooc_matrix)
  
  # ==========================
  # 8) Guardar resultados
  # ==========================
  
  # Crear directorio si no existe
  if(!dir.exists(file.path(base_path, "output"))) {
    dir.create(file.path(base_path, "output"))
  }
  
  # Guardar tabla completa
  write_csv(other_codes_table, 
            file.path(base_path, "output/other_comorbidities_detailed_analysis.csv"))
  
  # Guardar resumen por categorías
  write_csv(category_summary,
            file.path(base_path, "output/other_comorbidities_by_category.csv"))
  
  # Guardar prevalencia por tipo de cáncer
  prevalence_df$cancer_type <- rownames(prevalence_df)
  write_csv(prevalence_df,
            file.path(base_path, "output/other_comorbidities_by_cancer_type.csv"))
  
  cat("\n\n================================================\n")
  cat("Análisis completado. Archivos guardados en output/\n")
  cat("================================================\n")
  
} else {
  cat("No se encontraron otras comorbilidades en los datos.\n")
}