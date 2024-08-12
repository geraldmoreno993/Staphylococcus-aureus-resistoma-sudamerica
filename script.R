#merge de 3 bases de datos: amrfinderplus + staphopia + metadata de Patric
setwd("/home/gerald/Documentos/analisis/SA_sudamerica/merge")

library("readr")
library("dplyr")
################################################################################3
#data1
df_comma <- as.data.frame(read_delim("BVBRC_genome_2.csv", delim = ","))
class(df_comma)
barplot(table(df_comma$`Isolation Country`))
df_comma$`Assembly Accession`
# Cambiar nombres de columnas específicas
colnames(df_comma)[colnames(df_comma) == "Assembly Accession"] <- "accession"
df_comma$accession

#############################################################################
#data2
df_tab <- as.data.frame(read_delim("presence_absence.tsv", delim = "\t"))

# Transponer el data.frame
df_transposed <- as.data.frame(t(df_tab))

# La primera fila del transpuesto se convierte en nombres de columna
colnames(df_transposed) <- df_transposed[1, ]

# Eliminar la primera fila (ahora redundante)
df_transposed <- df_transposed[-1, ]

# Crear la nueva columna con la parte hasta el segundo guion bajo
df_transposed$new_column <- sub("^([^_]*_[^_]*)_.*", "\\1", rownames(df_transposed))
str(df_transposed)


# Cambiar nombres de columnas específicas
colnames(df_transposed)[colnames(df_transposed) == "new_column"] <- "accession"
df_transposed$accession
str(df_transposed)
# Convertir las columnas de "0" y "1" a numéricos 0 y 1
df_transposed <- df_transposed %>%
  mutate_if(~ all(. %in% c("0", "1")), ~ as.numeric(.))

# Verificar el resultado
print(df_transposed)
str(df_transposed)

#############################################################################
#data3

df_tab_2 <- as.data.frame(read_delim("output_final.tsv", delim = "\t"))
# Crear la nueva columna con la parte hasta el segundo guion bajo
df_tab_2$new_column <- sub("^([^_]*_[^_]*)_.*", "\\1", df_tab_2$sample)
# Cambiar nombres de columnas específicas
colnames(df_tab_2)[colnames(df_tab_2) == "new_column"] <- "accession"
df_tab_2$accession

str(df_tab_2)

# Convertir columnas lógicas a numéricas y reemplazar TRUE por 1 y FALSE por 0
df_tab_2 <- df_tab_2 %>%
  mutate_if(is.logical, ~ as.numeric(.))

# Verificar el resultado
print(df_tab_2)
str(df_tab_2)

#Merge de 3 objetos(databases) 
# Unir df1 y df2 por 'accession'
class(df_transposed$accession)
class(df_comma$accession)
df_merged <- full_join(df_comma, df_transposed, by = "accession")

# Unir el resultado con df3 por 'accession'
df_merged <- full_join(df_merged, df_tab_2, by = "accession")

# Verificar el resultado
print(df_merged)
df_merged$accession

#ver duplicados
barplot(table(df_merged$accession)) #No hay duplicados

#Guardando dataframe
write.csv(df_merged, file = "merge_csv.csv", row.names = FALSE)

#Eliminando observaciones con demasiados missing data (probablemente mala calidad u otras especies)
# Eliminar las filas 3 a 5 (por ejemplo)
df_merged_2 <- df_merged[-c(522:558), ]
df_merged_2 <- df_merged[-c(403:407), ]

#Reiniciando numeración de filas
rownames(df_merged_2) <- NULL

#Guardar tabla final y limpia
write.csv(df_merged_2, file = "merge_csv_2.csv", row.names = FALSE)


#Frecuencias de genes de resistencia
# Asegúrate de que el paquete dplyr esté instalado y cargado

library(dplyr)
library(ggplot2)

# Supongamos que tienes un dataframe llamado df

# Especifica el rango de columnas por su nombre
df_merged_2 <- df_merged_2 %>%
  mutate(suma_genes = rowSums(across(c(bleO:`lsa(A)`))))


plot(df_merged_2$suma_genes, df_merged_2$`GC Content`)
table(df_merged_2$'Collection Year')
plot(df_merged_2$'Collection Year', df_merged_2$suma_genes)

# Filtrar por el rango de años
df_filtered_year<- df_merged_2 %>%
  filter(`Collection Year` >= 2000 & `Collection Year` <= 2021)

# Crear el plot con los datos filtrados
plot(df_filtered_year$`Collection Year`, df_filtered_year$suma_genes,
     xlab = "Año de Colección",
     ylab = "Número Total de Genes de Resistencia",
     main = "Relación entre Año de Colección y Número de Genes de Resistencia",
     pch = 16, col = "blue")


# Filtrar por mecA = 1
df_filtered_meca <- df_merged_2 %>%
  filter(meca == 1)

# Crear el plot con los datos filtrados por mecA = 1
plot(df_filtered_meca$`Collection Year`, df_filtered_meca$suma_genes,
     xlab = "Año de Colección",
     ylab = "Número Total de Genes de Resistencia",
     main = "Relación entre Año de Colección y Número de Genes de Resistencia (mecA = 1)",
     pch = 16, col = "blue")



library(ggplot2)

df_filtered <- df_merged_2 %>%
  filter(!is.na(mecA) & !is.na(suma_genes))


ggplot(df_filtered, aes(x = mecA, y = suma_genes)) +
  geom_boxplot() +
  facet_wrap(~ mecA) +
  labs(title = "Boxplot de suma_genes por cada categoría de mecA (sin NA)",
       x = "mecA",
       y = "suma_genes") +
  theme_minimal()

# Verificar los niveles de la variable mecA
levels(df_filtered$mecA)

# Si tiene solo dos niveles, realizamos la prueba t
sd.test(suma_genes ~ mecA, data = df_filtered)
t_test_result <- t.test(suma_genes ~ mecA, data = df_filtered)
print(t_test_result)


#Hacemos comparaciones entre distintos tipos de mecA

library(tidyr)
library(dplyr)

# Seleccionar el rango de columnas
categorical_columns <- c("Ia","IIa","IIb","IIIa","IVa","IVb","IVc","IVd","IVg","IVh")

for (i in categorical_columns) {
  gg <- ggplot(df_filtered, aes(x = get(i), y = suma_genes)) +
    geom_boxplot() +
    labs(title = paste("Boxplot de suma_genes por", i, "(sin NA)"),
         x = i,
         y = "suma_genes") +
    theme_minimal()
  print(gg)
}



# Seleccionar el rango de columnas
categorical_columns <- c("I","II","III","IV","V","IV","VII","VIII","IX")

for (i in categorical_columns) {
  print(table(df_filtered[[i]]))
}


for (i in categorical_columns) {
  # Convertir la variable categórica a factor
  df_filtered[[i]] <- as.factor(df_filtered[[i]])
  
  gg <- ggplot(df_filtered, aes(x = get(i), y = suma_genes)) +
    geom_boxplot() +
    labs(title = paste("Boxplot de suma_genes por", i, "(sin NA)"),
         x = i,
         y = "suma_genes") +
    theme_minimal()
  
  print(gg)
}



# Convertir la variable mecA en factor (si es necesario)
df_merged_2$mecA <- as.factor(df_merged_2$mecA)

# Seleccionar un rango de variables para incluir en el modelo
# Supongamos que las columnas que te interesan están entre las columnas 10 y 30
# Cambia los números de columna según tu dataset
rango_columnas <- df_merged_2[, 88:210]

# Combinar la variable dependiente mecA con el rango de columnas seleccionadas
data_seleccionada <- cbind(df_merged_2$mecA, rango_columnas)

# Limpiar los nombres de las columnas
colnames(data_seleccionada) <- make.names(colnames(data_seleccionada))

# Crear la fórmula inicial para el modelo de regresión
formula_inicial <- as.formula(paste("mecA ~", paste(colnames(data_seleccionada)[-1], collapse = " + ")))

# Ejecutar el modelo inicial
modelo_inicial <- glm(formula_inicial, data = data_seleccionada, family = binomial)

# Realizar la selección stepwise
modelo_stepwise <- step(modelo_inicial, direction = "both")

# Resumen del modelo final
summary(modelo_stepwise)

# Predicción y evaluación del modelo final
predicciones <- predict(modelo_stepwise, type = "response")
data_seleccionada$predicted_mecA <- ifelse(predicciones > 0.5, 1, 0)

# Matriz de confusión
conf_matrix <- table(Predicted = data_seleccionada$predicted_mecA, Actual = data_seleccionada$mecA)
print(conf_matrix)

# Calcular la precisión del modelo
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Precisión del modelo:", round(accuracy, 4)))



##Frecuencia de todos los genes de resistencia:
colnames(data_seleccionada)

# Definir la variable categórica fija
var_fija <- "mecA"
df <- data_seleccionada
# Crear tablas de contingencia con las demás variables
tablas_contingencia <- lapply(names(df)[names(df) != var_fija], function(var) {
  tabla <- table(df[[var_fija]], df[[var]])
  names(dimnames(tabla)) <- c(var_fija, var) # Nombrar las dimensiones
  return(tabla)
})

print(tablas_contingencia)



# Crear tablas de proporciones generales para cada columna
tablas_proporciones <- lapply(df, function(x) {
  prop.table(table(x))
})

# Imprimir cada tabla de proporciones generales
for (i in seq_along(tablas_proporciones)) {
  cat("Tabla de proporciones de la variable", names(df)[i], ":\n")
  print(tablas_proporciones[[i]])
  cat("\n")
}


##GPAS

# Asegúrate de que las librerías necesarias están instaladas y cargadas
library(dplyr)

# El dataframe `df` tiene la primera columna como el outcome y el resto como genes
outcome <- data_seleccionada[[1]]  # Primera columna como outcome
gene_matrix <- data_seleccionada[, -1]  # Todas las demás columnas como genes

# Función para realizar GPAS
gpas_analysis <- function(outcome, gene_matrix) {
  p_values <- apply(gene_matrix, 2, function(gene) {
    model <- glm(outcome ~ gene, family = binomial)
    summary(model)$coefficients[2, 4]  # p-valor del coeficiente del gen
  })
  return(p_values)
}

# Realizar el análisis de asociación
p_values <- gpas_analysis(outcome, gene_matrix)

# Crear un dataframe con los resultados
results <- data.frame(Gene = colnames(gene_matrix), P_value = p_values)

# Ordenar los resultados por el p-valor
results <- results %>% arrange(P_value)

# Mostrar los resultados
print(results)

# Aplicar la corrección de p-valores utilizando el método de Benjamini-Hochberg
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Crear un dataframe con los resultados
results <- data.frame(Gene = colnames(gene_matrix), P_value = p_values, Adjusted_P_value = adjusted_p_values)

# Mostrar los resultados
print(results)

# Ordenar los resultados por el p-valor ajustado
results <- results %>% arrange(Adjusted_P_value)


# Filtrar los genes con un p-valor ajustado menor a 0.001
significant_results <- results %>% filter(Adjusted_P_value < 0.0001)

# Mostrar los resultados significativos
print(significant_results)

##Aplicando Bon ferroni
# Aplicar la corrección de p-valores utilizando el método de Bonferroni
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")

# Crear un dataframe con los resultados
results <- data.frame(Gene = colnames(gene_matrix), P_value = p_values, Adjusted_P_value = adjusted_p_values)

# Ordenar los resultados por el p-valor ajustado
results <- results %>% arrange(Adjusted_P_value)

# Filtrar los genes con un p-valor ajustado menor a 0.001
significant_results <- results %>% filter(Adjusted_P_value < 0.0001)

# Mostrar los resultados significativos
print(significant_results)



##ORs
# Asegúrate de que las librerías necesarias están instaladas y cargadas
library(dplyr)

# El dataframe `df` tiene la primera columna como el outcome y el resto como genes
outcome <- data_seleccionada[[1]]  # Primera columna como outcome
gene_matrix <- data_seleccionada[, -1]  # Todas las demás columnas como genes

# Función para realizar GPAS y obtener OR
gpas_analysis_or <- function(outcome, gene_matrix) {
  results <- apply(gene_matrix, 2, function(gene) {
    model <- glm(outcome ~ gene, family = binomial)
    coef <- summary(model)$coefficients[2, 1]  # Coeficiente del gen
    p_value <- summary(model)$coefficients[2, 4]  # p-valor del coeficiente del gen
    or <- exp(coef)  # Odds ratio
    c(OR = or, P_value = p_value)
  })
  return(as.data.frame(t(results)))
}

# Realizar el análisis de asociación y obtener OR
results <- gpas_analysis_or(outcome, gene_matrix)

# Aplicar la corrección de p-valores utilizando el método de Benjamini-Hochberg
results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")

# Ordenar los resultados por el p-valor ajustado
results <- results %>% arrange(Adjusted_P_value)

# Filtrar los genes con un p-valor ajustado menor a 0.001
significant_results <- results %>% filter(Adjusted_P_value < 0.0001)

# Mostrar los resultados significativos con OR
print(significant_results)

# Aplicar la corrección de p-valores utilizando el método de Bonferroni
results$Bonferroni_P_value <- p.adjust(results$P_value, method = "bonferroni")

# Ordenar los resultados por el p-valor de Bonferroni
results <- results %>% arrange(Bonferroni_P_value)

# Filtrar los genes con un p-valor ajustado menor a 0.001 según Bonferroni
significant_results_bonferroni <- results %>% filter(Bonferroni_P_value < 0.0001)

# Mostrar los resultados significativos con OR después de Bonferroni
print(significant_results_bonferroni)
