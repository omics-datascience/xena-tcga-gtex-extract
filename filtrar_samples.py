import polars as pl
import sys

# Configuración de archivos
COUNTS_FILE = "cohort_TCGA_TARGET_GTEx/expected_counts_without_TARGET_samples.tsv"
META_FILE = "cohort_TCGA_TARGET_GTEx/TCGA_GTEX_category.txt"
OUTPUT_COUNTS = "filtered_datasets/filtered_counts.txt"
OUTPUT_META = "filtered_datasets/filtered_metadata.txt"

# Obtener keywords de argumentos
keywords = sys.argv[1:]
pattern = "|".join(keywords)

# Filtrar Metadatos
df_meta = pl.read_csv(META_FILE, separator='\t')

# Filtramos las filas que contienen alguna de las keywords en la columna 'TCGA_GTEX_main_category'
filtered_meta = df_meta.filter(
    pl.col("TCGA_GTEX_main_category").str.contains(pattern)
)
filtered_meta.write_csv(OUTPUT_META, separator='\t')

# Extraemos la lista de IDs a mantener
ids_to_keep = filtered_meta["sample"].to_list()

# Filtrar la Matriz de Conteos
lazy_counts = pl.scan_csv(COUNTS_FILE, separator='\t')

# Seleccionamos la columna de genes y las muestras filtradas
# collect_schema() para pedir los nombres en Lazy mode
all_columns = lazy_counts.collect_schema().names()
all_columns_set = set(all_columns)

# Buscamos la interseccion: columnas que están en el archivo Y en nuestros IDs
selected_cols = [all_columns[0]] + [c for c in ids_to_keep if c in all_columns_set]

# Ejecutar el filtrado y guardar el resultado
print(f"Escribiendo {len(selected_cols) - 1} muestras en el archivo de salida...")
lazy_counts.select(selected_cols).sink_csv(OUTPUT_COUNTS, separator='\t')

print("Proceso finalizado!")
