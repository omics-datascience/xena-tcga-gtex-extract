import polars as pl
import sys

# Configuración de archivos
COUNTS_FILE = "cohort_TCGA_TARGET_GTEx/expected_counts_without_TARGET_samples.tsv"
META_FILE = "cohort_TCGA_TARGET_GTEx/TCGA_GTEX_category.txt"
DB_FILE = "cohort_TCGA_TARGET_GTEx/probeMap_gencode.v23.annotation.gene.probemap"
OUTPUT_COUNTS = "filtered_datasets/filtered_counts.txt"
OUTPUT_META = "filtered_datasets/filtered_metadata.txt"


def aplicar_mapeo_genes(lazy_counts, db_path, id_col_name):
    """
    Realiza un Left Join con la base de datos de genes y reemplaza la columna de ID
    por el Gene Symbol,manteniendo el nombre original de la columna.
    """
    print(f"   -> Preparando mapeo de genes usando: {db_path}")

    # 1. Escanear la base de datos de mapeo (Lazy)
    # Asumimos que tiene cabeceras 'id' y 'gene'
    map_lf = pl.scan_csv(db_path, separator='\t').select(['id', 'gene'])

    # 2. Realizar el Join
    # Unimos la tabla de conteos con el mapa usando la columna de ID
    joined_lf = lazy_counts.join(
        map_lf,
        left_on=id_col_name,
        right_on='id',
        how='left'
    )

    # 3. Construir la selección de columnas para mantener el orden y nombre
    # Obtenemos todas las columnas del lazy frame original EXCEPTO la de ID
    cols_muestras = [c for c in lazy_counts.collect_schema().names() if c != id_col_name]

    # Definimos la nueva primera columna:
    # Usamos 'gene'. Si es nulo (no hubo match),rellenamos con el ID original ('id_col_name')
    # Finalmente,renombramos esta columna para que se llame igual que la original (ej. "")
    new_first_col = pl.col("gene").fill_null(pl.col(id_col_name)).alias(id_col_name)

    # 4. Retornar el LazyFrame con el orden correcto: [Gene,Muestra1,Muestra2...]
    return joined_lf.select([new_first_col] + cols_muestras)


# Obtener keywords de argumentos
keywords = sys.argv[1:]
if not keywords:
    print("Error: Faltan las categorías.")
    print(f"Uso: python {sys.argv[0]} categoria1 categoria2")
    print("Use comillas si la categoría tiene espacios.")
    print("Use bash listar_categorias.sh para ver las categorías disponibles.")
    sys.exit(1)

print(f"Filtrando coincidencias exactas en 'TCGA_GTEX_main_category' para: {keywords}")

# 1. Filtrar Metadatos
df_meta = pl.read_csv(META_FILE, separator='\t')

# Filtramos las filas que contienen alguna de las keywords
filtered_meta = df_meta.filter(
    pl.col("TCGA_GTEX_main_category").is_in(keywords)
)
filtered_meta.write_csv(OUTPUT_META, separator='\t')

num_muestras = filtered_meta.height
if num_muestras == 0:
    print(f"ERROR: No se encontraron muestras exactas para: {keywords}")
    sys.exit(1)

print(f"   -> Se encontraron {num_muestras} muestras.")

# Extraemos la lista de IDs a mantener
ids_to_keep = filtered_meta["sample"].to_list()

# 2. Configurar la Matriz de Conteos (Lazy)
lazy_counts = pl.scan_csv(COUNTS_FILE, separator='\t')

# collect_schema() para pedir los nombres en Lazy mode
all_columns = lazy_counts.collect_schema().names()
all_columns_set = set(all_columns)
id_col_name = all_columns[0]  # Detectamos el nombre de la primera columna (ej. "" o "id")

# Buscamos la interseccion: columnas que están en el archivo Y en nuestros IDs
# Mantenemos la primera columna siempre
selected_cols = [id_col_name] + [c for c in ids_to_keep if c in all_columns_set]

print(f"   -> Seleccionando {len(selected_cols) - 1} columnas de muestras...")

# Aplicamos la selección de columnas preliminar (para no cargar datos inútiles al join)
lazy_counts_filtered = lazy_counts.select(selected_cols)

# 3. Aplicar el cambio de ID a Gene Symbol
# Llamamos a nuestra función auxiliar
lazy_final = aplicar_mapeo_genes(lazy_counts_filtered, DB_FILE, id_col_name)

# 4. Ejecutar y guardar
print("Escribiendo archivo de conteos filtrado y mapeado...")
# Usamos sink_csv para procesar en streaming sin cargar todo en RAM
lazy_final.sink_csv(OUTPUT_COUNTS, separator='\t')

print(f"Proceso finalizado! Resultados en: {OUTPUT_COUNTS}")
