# ==============================================================================
# Pipeline de Análisis de Expresión Diferencial (Limma Trend)
# ==============================================================================

# 1. Cargar librerías
library(limma)
library(edgeR)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)
library(pheatmap)

# 2. Cargar funciones
source("dea_functions.R") 

# 3. Configuración
sample_col <- "sample"
atributo <- "TCGA_GTEX_main_category"
metadata_path <- file.path("../filtered_datasets/filtered_metadata.txt")
rnaseq_path <- file.path("../filtered_datasets/filtered_counts.txt")

if (!dir.exists("DEA_output")) { dir.create("DEA_output") } 

# 4. Carga y verificacion de Datos
cat("Cargando datos...\n")
# check.names = FALSE es vital para que no rompa tus IDs si tienen guiones
metadata <- read.delim(metadata_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
expression_data <- process_rnaseq_data(rnaseq_path)

cat("Alineando orden de muestras (Metadata vs Counts)...\n")

# Aseguramos que solo usamos las muestras que están en ambos archivos
common_samples <- intersect(metadata[[sample_col]], colnames(expression_data))

# Reordenamos las columnas de la matriz
expression_data <- expression_data[, common_samples]

# Reordenamos las filas de la metadata para que coincidan con la matriz
metadata <- metadata[match(common_samples, metadata[[sample_col]]), ]

# Sanity Check Final (Si esto falla, el script se detiene)
if(!all(colnames(expression_data) == metadata[[sample_col]])) {
  stop("ERROR CRÍTICO: Los IDs de muestras estan desordenados. Revisa los nombres de columnas.")
}

cat("Orden corregido. Muestras listas:", length(common_samples), "\n")

# 5. Normalización (Between Arrays)
# Asumimos que el orden Metadata vs Matriz es idéntico.
cat("Normalizando...\n")
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

# Ejecutamos el análisis MDS/PCA
create_pca_mds(
    expression_data = expression_data,
    metadata = metadata,
    group_col = atributo,
    comparison_title = "TCGA_vs_GTEx",
    output_dir = "DEA_output"
)

# NOTA: En este punto, el usuario debe abrir el plot generado 
# en DEA_output/QC_MDS_plot_TCGA_vs_GTEx.png y verificar que:
# 1. Los puntos Tumor y Normal formen dos nubes distintas.
# 2. Ninguna muestra se haya salido volando (Outlier).
# Forzamos limpieza de memoria por si quedaron temporales
# gc(verbose = FALSE)

# 6. Filtro (Sobrescribiendo)
cat("Filtrando genes...\n")
# Seguimos sobrescribiendo la misma variable
expression_data <- filter_lowly_expressed_genes(expression_data)

# 7. Ejecución
cat("Ejecutando Análisis Diferencial...\n")
# Pasamos la variable única
all_results <- perform_differential_expression(expression_data, metadata, atributo)
nombres_comparaciones <- names(all_results)

# 8. Guardado
cat("Guardando resultados...\n")
for (nombre in nombres_comparaciones) {
  df_res <- all_results[[nombre]]
  
  ordered_results <- order_results_by_pvalue(df_res)
  top_50 <- get_top_50_by_pvalue(ordered_results)
  
  write.csv(ordered_results, file = paste0("DEA_output/", nombre, "_results.csv"), row.names = TRUE)
  write.csv(top_50, file = paste0("DEA_output/", nombre, "_top50_results.csv"), row.names = TRUE)
  
  results_plot <- create_volcano_plot(
    ordered_results,
    title = paste0("Volcano: ", nombre)
  )
  
  ggsave(filename = paste0("DEA_output/", nombre, "_volcano_plot.png"),
         plot = results_plot, width = 10, height = 7)

  create_heatmap(
           full_matrix = expression_data,
           results_df = df_res,
           metadata = metadata,
           group_col = atributo,
           comparison_name = nombre,
           top_n = 50,
           output_dir = "DEA_output" 
       )
}
# gc(verbose = FALSE)
cat("Análisis completado. Resultados en carpeta DEA_output/\n")