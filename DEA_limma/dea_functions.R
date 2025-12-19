# ==============================================================================
# Archivo: dea_functions.R
# Descripción: Funciones auxiliares para el análisis de expresión diferencial
# ==============================================================================

# --- Procesamiento y Carga ---
# En dea_functions.R

process_rnaseq_data <- function(rnaseq_path) {
  # 1. Cargar datos
  cat("Leyendo archivo...\n")
  rnaseq_data <- read.table(rnaseq_path, header = TRUE, sep = "\t", check.names = FALSE)

  # Asegurar nombre de columna gen
  if (!"gen" %in% colnames(rnaseq_data)) {
    colnames(rnaseq_data)[1] <- "gen"
  }

  # 2. Separar genes y datos numéricos
  genes <- rnaseq_data$gen
  # Convertimos a matriz numérica directamente, excluyendo la columna 'gen'
  # Esto es más ligero que pasar por data.frames de dplyr
  matriz_num <- as.matrix(rnaseq_data[, -1])
  
  # Liberamos memoria del objeto original inmediatamente
  rm(rnaseq_data) 
  gc() # Garbage Collector: fuerza a R a limpiar la RAM

  # 3. Promediar duplicados usando LIMMA (Mucho más eficiente)
  if(anyDuplicated(genes)) {
    cat("Promediando genes duplicados con limma::avereps...\n")
    # avereps hace exactamente lo que hacíamos con dplyr pero optimizado para matrices
    matriz_final <- avereps(matriz_num, ID = genes)
  } else {
    matriz_final <- matriz_num
    rownames(matriz_final) <- genes
  }

  return(matriz_final)
}

# --- Filtros ---
filter_lowly_expressed_genes <- function(log_data, threshold_percentile = 0.15) {
  avg_expression <- rowMeans(log_data, na.rm = TRUE)
  threshold <- quantile(avg_expression, threshold_percentile)
  
  genes_to_keep <- names(avg_expression[avg_expression >= threshold])
  log_data_filtered <- log_data[genes_to_keep, ]
  
  cat("Genes eliminados por baja expresión:", nrow(log_data) - length(genes_to_keep), "\n")
  return(log_data_filtered)
}

# --- Análisis Diferencial (Limma) ---
perform_differential_expression <- function(log_data_filtered, metadata, clinical_attribute) {
  
  variables_count <- length(unique(metadata[[clinical_attribute]]))
  if (variables_count < 2) stop("Se necesitan al menos dos categorías para el análisis.")
  
  group <- factor(metadata[[clinical_attribute]])
  levels(group) <- make.names(levels(group)) 

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(log_data_filtered, design)

  grupos <- levels(group)
  pares <- combn(grupos, 2, simplify = FALSE)
  
  contrastes <- sapply(pares, function(par) paste0(par[2], " - ", par[1]))
  cat("Contrastes generados:\n")
  print(contrastes)
  
  contrast_matrix <- makeContrasts(contrasts = contrastes, levels = design)

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  results <- lapply(seq_len(ncol(contrast_matrix)), function(i) {
    topTable(fit2, coef = i, number = Inf, adjust.method = "BH")
  })
  names(results) <- colnames(contrast_matrix)

  return(results)
}

# --- Helpers de Ordenamiento ---
order_results_by_pvalue <- function(results) { 
  results[order(results$adj.P.Val), ] 
}

get_top_50_by_pvalue <- function(ordered_results) { 
  head(ordered_results, 50) 
}

# --- Visualización ---
# En dea_functions.R

create_pca_mds <- function(expression_data, metadata, group_col, comparison_title, output_dir, top_n_genes = 1000) {
  cat(paste("Calculando MDS/PCA usando los", top_n_genes, "genes más variables...\n"))

  # 1. Identificar los genes más variables
  # Calculamos la varianza de cada gen (fila)
  variances <- matrixStats::rowVars(expression_data)
  
  # Ordenamos y seleccionamos los N primeros
  top_genes_idx <- order(variances, decreasing = TRUE)[1:top_n_genes]
  
  # 2. Creamos el subset de la matriz
  subset_data <- expression_data[top_genes_idx, ]

  # 3. Aplicamos plotMDS solo a la matriz reducida
  mds <- plotMDS(subset_data, plot = FALSE)
  
  # Creamos el dataframe para ggplot
  df_mds <- data.frame(
    Dim1 = mds$x, 
    Dim2 = mds$y, 
    Group = metadata[[group_col]]
  )
  
  plot_mds <- ggplot(df_mds, aes(x = Dim1, y = Dim2, color = Group)) +
    geom_point(size = 3) +
    labs(title = paste("MDS Plot - Control de Calidad:", comparison_title),
         subtitle = paste0("Muestras visualizadas usando el Top ", top_n_genes, " genes más variables."),
         x = paste("Leading LogFC Dimension 1 (Var:", round(mds$var.explained[1]*100, 1), "%)"),
         y = paste("Leading LogFC Dimension 2 (Var:", round(mds$var.explained[2]*100, 1), "%)")) +
    theme_bw()
  
  # Guardar el plot
  ggsave(filename = paste0(output_dir, "/QC_MDS_plot_", comparison_title, ".png"),
         plot = plot_mds, width = 8, height = 7)
  
  cat("Gráfico QC (MDS/PCA) guardado. Por favor, revísalo.\n")

  
  return(plot_mds)
}

create_volcano_plot <- function(results, p_value_threshold = 0.05, logFC_threshold = 1, title = "Volcano Plot", top_genes_label = 20) {
  
  # 1. Clasificación para Colores (Igual que antes)
  results$diff_expressed <- "NO"
  results$diff_expressed[results$adj.P.Val < p_value_threshold & abs(results$logFC) > logFC_threshold] <- "YES"
  
  # 2. Selección INTELIGENTE de etiquetas
  # En lugar de pasar TODOS los significativos a ggrepel, seleccionamos solo el Top N
  # Ordenamos por p-valor (los más significativos primero)
  top_genes <- results %>%
    filter(diff_expressed == "YES") %>%
    arrange(adj.P.Val) %>%
    head(top_genes_label) # Nos quedamos solo con los N primeros (ej. 20)
  
  # 3. Generar el Plot
  volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = diff_expressed)) +
    
    # Puntos de fondo (Todos los genes)
    geom_point(size = 1.5, alpha = 0.6) + 
    
    scale_color_manual(values = c("NO" = "grey60", "YES" = "firebrick"), 
                       labels = c("No Significativo", "Significativo")) +
    
    # Líneas de corte
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "black", alpha=0.5) + 
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black", alpha=0.5) + 
    
    # ETIQUETAS: Usamos el subset 'top_genes' en lugar de 'significant_genes'
    geom_text_repel(data = top_genes, 
                    aes(label = rownames(top_genes)), 
                    size = 3.5, 
                    box.padding = 0.5,
                    max.overlaps = Inf, # Permitimos overlaps infinitos porque son pocos genes # nolint: line_length_linter.
                    show.legend = FALSE,
                    color = "black") + # Texto en negro para que se lea mejor
    
    labs(title = title, 
         subtitle = paste("Top", top_genes_label, "genes etiquetados por significancia"),
         x = "Log2 Fold Change", 
         y = "-log10(Adjusted P-value)",
         color = "Estado") +
    ylim(0,375) +
    theme_bw() + 
    theme(legend.position = "top")
  
  return(volcano_plot)
}



create_heatmap <- function(full_matrix, results_df, metadata, group_col = "TCGA_GTEX_main_category", comparison_name, top_n = 50, logFC_threshold = 1, p_threshold = 0.05, output_dir = "DEA_output") {
  # 1. Filtrar los genes de interés
  genes_sig <- results_df %>% dplyr::filter(abs(logFC) > logFC_threshold & adj.P.Val < p_threshold) %>% dplyr::arrange(adj.P.Val) %>% head(top_n)

  if (nrow(genes_sig) == 0) {
    message("No se encontraron genes significativos para el heatmap con los umbrales definidos.")
    return(NULL)
  }

  # 2. Preparar Matriz de Expresión
  # Subconjunto de la matriz completa solo con los genes seleccionados
  heatmap_matrix <- full_matrix[rownames(genes_sig), ]

  # Escalar por fila (Z-Score). Esto es necesario para que el color represente
  # la desviación del promedio de cada gen, no su valor absoluto.
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

  # 3. Preparar Anotación de Muestras (Columna)
  # La metadata debe estar ordenada exactamente como las columnas de la matriz (ya lo garantizamos en main_analysis.R)

  annotation_data <- metadata %>% dplyr::select(!!rlang::sym(group_col)) # Usamos rlang::sym para manejar nombres de columnas con espacios/puntos

  # Asegurar que los row/col names coincidan (aunque ya debería ser cierto)
  rownames(annotation_data) <- colnames(heatmap_matrix)

  # 4. Generar el Heatmap
  pheatmap(
       heatmap_matrix_scaled,
       color = colorRampPalette(c("blue", "white", "red"))(100), # Colores: Azul (Bajo) a Rojo (Alto)
       cluster_rows = TRUE,
       cluster_cols = TRUE,
       show_rownames = TRUE,
       show_colnames = FALSE, # Ocultamos IDs de muestra para que no sature
       annotation_col = annotation_data, # Añadimos la barra de Tumor/Normal
       main = paste("Heatmap Top", nrow(genes_sig), "DEGs:", comparison_name),
       fontsize = 6,
       filename = paste0("DEA_output/", gsub(" - ", "_vs_", comparison_name), "_heatmap_top", nrow(genes_sig), ".png"),
       width = 8,
       height = 10
  )

  cat("Heatmap guardado exitosamente para la comparación:", comparison_name, "\n")

  # Devolvemos el número de genes ploteados (al pedo)
  return(nrow(genes_sig))
}