load("C:/Users/carme/Downloads/summarized_experiment.Rda")

# Primero extraemos los datos metabolómicos y el tipo de tratamiento de cada muestra
assay_data <- assay(se)
assay_data <- assay_data[complete.cases(assay_data), ]
treatment <- as.character(colData(se)$Treatment)

# Partimos nuestros datos según el tratamiento
baseline_data <- assay_data[, treatment == "Baseline"]
apple_data <- assay_data[, treatment == "Apple"]
cranberry_data <- assay_data[, treatment == "Cranberry"]

# Creamos un dataframe con la media, mediana y varianza de cada metabolito, y ordenamos de mayor a menor según la media.

# Comenzamos con los datos basales
baseline_stats<-data.frame(
  Metabolite = rownames(baseline_data),
  Mean = rowMeans(baseline_data),
  Median = apply(baseline_data, 1, median),
  Variance = apply(baseline_data, 1, var)
)
baseline_stats$Metabolite_Name <- rowData(se)$names[match(baseline_stats$Metabolite, rownames(rowData(se)))]
baseline_stats <- baseline_stats[order(-baseline_stats$Mean), ]
head(baseline_stats)

# Realizamos el mismo procedimiento con el grupo de zumo de manzana.
apple_stats <- data.frame(
  Metabolite = rownames(apple_data),
  Mean = rowMeans(apple_data),
  Median = apply(apple_data, 1, median),
  Variance = apply(apple_data, 1, var)
)
apple_stats$Metabolite_Name <- rowData(se)$names[match(apple_stats$Metabolite, rownames(rowData(se)))]
apple_stats <- apple_stats[order(-apple_stats$Mean), ]
head(apple_stats)

# Realizamos el mismo procedimiento con el grupo de zumo de arándano rojo.
cranberry_stats <- data.frame(
  Metabolite = rownames(cranberry_data),
  Mean = rowMeans(cranberry_data),
  Median = apply(cranberry_data, 1, median),
  Variance = apply(cranberry_data, 1, var)
)
cranberry_stats$Metabolite_Name <- rowData(se)$names[match(cranberry_stats$Metabolite, rownames(rowData(se)))]
cranberry_stats <- cranberry_stats[order(-cranberry_stats$Mean), ]
head(cranberry_stats)
# Transformación logarítmica de los datos
log_assay_data <- log2(assay_data + 1)

baseline_values <- as.vector(log_assay_data[, treatment == "Baseline"])
apple_values <- as.vector(log_assay_data[, treatment == "Apple"])
cranberry_values <- as.vector(log_assay_data[, treatment == "Cranberry"])


baseline_density <- density(baseline_values, na.rm = TRUE)
apple_density <- density(apple_values, na.rm = TRUE)
cranberry_density <- density(cranberry_values, na.rm = TRUE)

# Realizamos el gráfico de densidad
plot(baseline_density, main = "Grafico de densidad de la expresion de metabolitos por tratamiento",
     xlab = "Log2(Nivel de expresion)", ylab = "Densidad", col = "blue", lwd = 2, 
     ylim = range(0, max(baseline_density$y, apple_density$y, cranberry_density$y)))
lines(apple_density, col = "green", lwd = 2)
lines(cranberry_density, col = "red", lwd = 2)
legend("topright", legend = c("Baseline", "Apple", "Cranberry"),
       col = c("blue", "green", "red"), lwd = 2)

# Escalamos los datos 
# Transponemos los datos ya que scale() trabaja por defecto por columnas, luego volvemos a transponer
assay_data_scaled <- t(scale(t(log_assay_data)))

#Realizamos el PCA 
pca_results <- prcomp(t(assay_data_scaled), scale. = TRUE)
groups <- colData(se)$Treatment
pca_data <- as.data.frame(pca_results$x)
pca_data$sample <- colnames(assay_data)
pca_data$group  <- groups  

# Seleccionamos los 2 componentes principales para hacer el gráfico
library(ggplot2)

ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "PCA por grupo de tratamiento", x = "PC1", y = "PC2") +
  theme_minimal() +
  scale_color_manual(values = c("green", "blue", "red"))

#https://aspteaching.github.io/Analisis_de_datos_omicos-Ejemplo_0-Microarrays/ExploreArrays.html
group_baseline <- which(treatment == "Baseline")
group_apple <- which(treatment == "Apple")
group_cranberry <- which(treatment == "Cranberry")

# Creamos nuestra función para realizar el t-test entre dos grupos
ttest <- function(x, group1, group2) {
  tt <- t.test(x[group1], x[group2])
  return(c(tt$statistic, tt$p.value, tt$estimate[1] - tt$estimate[2]))
}
# Probamos nuestra función ttest con el grupo baseline y apple.
ans_apple_baseline <- apply(assay_data, 1, function(x) ttest(x, group_apple, group_baseline))
ts_apple_baseline <- ans_apple_baseline[1, ]
pvals_apple_baseline <- ans_apple_baseline[2, ]
fc_apple_baseline <- ans_apple_baseline[3, ]

# Creamos un histograma de los estadísticos t
hist(ts_apple_baseline, breaks=100)
#https://www.htgmolecular.com/blog/2022-08-25/understanding-volcano-plots
# Realizamos el volcano plot
plot(fc_apple_baseline, log(pvals_apple_baseline),
     xlab = "Fold Change (Apple vs Baseline)",
     ylab = "-log10(p-value)")

# Establecemos el umbrak como 0.01 para centrarnos en menos metabolitos, adoptamos un enfoque más conservador
pval_threshold <- 0.01
# Filtramos los índices de los metabolitos con los pvalores que buscamos
significant_indices <- which(pvals_apple_baseline < pval_threshold)

# Obtenemos los nombres de los metabolitos y sus p-valores
significant_metabolite_names <- rowData(se)$names[significant_indices]
significant_pvalues <- pvals_apple_baseline[significant_indices]


significant_metabolites_df <- data.frame(
  Metabolite = significant_metabolite_names,
  P_Value = significant_pvalues
)
significant_metabolites_df

#Probamos nuestra función ttest con el grupo baseline y cranberry

ans_cranberry_baseline <- apply(assay_data, 1, function(x) ttest(x, group_cranberry, group_baseline))
ts_cranberry_baseline <- ans_cranberry_baseline[1, ]
pvals_cranberry_baseline <- ans_cranberry_baseline[2, ]
fc_cranberry_baseline <- ans_cranberry_baseline[3, ]
hist(ts_cranberry_baseline, breaks=100)

plot(fc_apple_baseline, log(pvals_cranberry_baseline),
     xlab = "Fold Change (Cranberry vs Baseline)",
     ylab = "-log10(p-value)")

#Creamos el dataframe con los metabolitos con p-valores significativos en cranberry vs baseline
significant_indices <- which(pvals_cranberry_baseline < pval_threshold)


significant_metabolite_names <- rowData(se)$names[significant_indices]
significant_pvalues <- pvals_cranberry_baseline[significant_indices]


significant_metabolites_df <- data.frame(
  Metabolite = significant_metabolite_names,
  P_Value = significant_pvalues
)
significant_metabolites_df

#Probamos nuestra función ttest con el grupo apple y cranberry.
ans_apple_cranberry <- apply(assay_data, 1, function(x) ttest(x, group_apple, group_cranberry))
ts_apple_cranberry <- ans_apple_cranberry[1, ]
pvals_apple_cranberry <- ans_apple_cranberry[2, ]
fc_apple_cranberry <- ans_apple_cranberry[3, ]
hist(ts_apple_cranberry, breaks=100)


plot(fc_apple_baseline, log(pvals_apple_cranberry),
     xlab = "Fold Change (Apple vs Cranberry)",
     ylab = "-log10(p-value)")

#Creamos el dataframe con los metabolitos con p-valores significativos en cranberry vs apple
significant_indices <- which(pvals_apple_cranberry < pval_threshold)

significant_metabolite_names <- rowData(se)$names[significant_indices]
significant_pvalues <- pvals_apple_cranberry[significant_indices]

significant_metabolites_df <- data.frame(
  Metabolite = significant_metabolite_names,
  P_Value = significant_pvalues
)
significant_metabolites_df

# Seleccionamos los datos con expresión diferencial significativa del paso anterior (apples vs cranberries)
assay_data_significant <- assay_data[significant_indices, ]

# Seleccionamos las columnas que comienzan con "a" (apple) o "c" (cranberry)
app_cran_samples <- grep("^(a|c)", colnames(assay_data_significant), value = TRUE)
assay_data_cranberry_apple <- assay_data_significant[, app_cran_samples ]

library(pheatmap)
pheatmap(
  assay_data_cranberry_apple,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",  
  clustering_method = "complete",          
  show_rownames = TRUE,                    
  show_colnames = TRUE
)

