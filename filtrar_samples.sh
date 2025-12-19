#!/bin/bash

COUNTS_FILE="cohort_TCGA_TARGET_GTEx/expected_counts_without_TARGET_samples.tsv"
META_FILE="cohort_TCGA_TARGET_GTEx/TCGA_GTEX_category.txt"


mkdir -p filtered_datasets
OUTPUT_COUNTS="filtered_datasets/filtered_counts.txt"     # Matriz filtrada
OUTPUT_META="filtered_datasets/filtered_metadata.txt"     # Metadatos filtrados

IDS_FILE="ids_temp.tmp"                 # Archivo temporal

# --- VALIDAR ARGUMENTOS Y DEFINIR PATRÓN DE BÚSQUEDA ---
if [ $# -eq 0 ]; then
    echo "Error: Faltan las palabras clave."
    echo "Uso: $0 Keyword1 [Keyword2 ...]"
    exit 1
fi

# Unimos argumentos con OR (|) para grep
IFS='|'
PATTERN="$*"
unset IFS

echo "Filtro RNA-Seq: $PATTERN"
# --- GENERAR NUEVO ARCHIVO DE METADATOS ---
echo "1. Creando archivo de metadatos filtrado ($OUTPUT_META)..."

# Extraemos el encabezado del metadata original (primera línea)
head -n 1 "$META_FILE" > "$OUTPUT_META"
# Buscamos las filas que coincidan y las agregamos al nuevo archivo
grep -E "$PATTERN" "$META_FILE" >> "$OUTPUT_META"

NUM_MUESTRAS=$(($(wc -l < "$OUTPUT_META") - 1))
if [ "$NUM_MUESTRAS" -le 0 ]; then
    echo "ERROR: No se encontraron muestras!"
    rm "$OUTPUT_META"
    exit 1
fi

echo "   -> Se encontraron $NUM_MUESTRAS muestras."

# --- PREPARAR LISTA DE IDs ---
# Extraemos solo la columna 1 (Sample ID) del nuevo archivo de metadatos
# 'NR>1' salta el encabezado para no buscar "SampleID" como si fuera una muestra
awk 'NR>1 {print $1}' "$OUTPUT_META" > "$IDS_FILE"

# --- FILTRAR LA MATRIZ DE CONTEOS ---
echo "2. Filtrando la matriz de conteos ($OUTPUT_COUNTS)..."

awk -F'\t' '
    # Cargar IDs en memoria
    NR==FNR {
        keep[$1] = 1; 
        next
    }

    # Procesar encabezado de la matriz
    FNR==1 {
        printf "%s", "gen"; # Imprimir nombre columna gen
        
        for (i=2; i<=NF; i++) {
            if ($i in keep) {
                cols[i] = 1;      # Guardar indice de columna válida
                printf "\t%s", $i # Imprimir nombre de muestra
            }
        }
        printf "\n";
    }

    # Procesar filas de datos
    FNR>1 {
        printf "%s", $1; # Imprimir gen actual
        
        for (i=2; i<=NF; i++) {
            if (i in cols) {
                printf "\t%s", $i # Imprimir valor solo si la columna es valida
            }
        }
        printf "\n";
    }
' "$IDS_FILE" "$COUNTS_FILE" > "$OUTPUT_COUNTS"

rm "$IDS_FILE"

echo "Proceso finalizado!"
echo "   - Metadatos: $OUTPUT_META"
echo "   - Conteos:   $OUTPUT_COUNTS"