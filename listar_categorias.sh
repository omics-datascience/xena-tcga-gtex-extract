#!/bin/bash

# --- CONFIGURACIÓN ---
META_FILE="cohort_TCGA_TARGET_GTEx/TCGA_GTEX_category.txt"
TARGET_COL="TCGA_GTEX_main_category"

# --- VALIDACIÓN ---
if [ ! -f "$META_FILE" ]; then
    echo "❌ Error: No encuentro el archivo $META_FILE"
    exit 1
fi

echo "Conteo de muestras por categoría en: '$TARGET_COL'"
echo "--------------------------------------------------------"
echo -e " CANTIDAD\tCATEGORÍA"
echo "--------------------------------------------------------"

# --- PROCESAMIENTO ---
# 1. awk busca la columna y extrae los valores
# 2. sort agrupa (necesario para uniq)
# 3. uniq -c cuenta las ocurrencias
# 4. sort -nr ordena numéricamente (n) de mayor a menor (r)

awk -v col="$TARGET_COL" -F'\t' '
    NR==1 {
        for (i=1; i<=NF; i++) {
            if ($i == col) {
                target_idx = i;
                break;
            }
        }
        if (target_idx == 0) {
            # Si falla, imprimimos en error estándar
            print "ERROR: Columna no encontrada" > "/dev/stderr"
            exit 1
        }
    }
    NR>1 {
        if (target_idx > 0) {
            print $target_idx
        }
    }
' "$META_FILE" | sort | uniq -c 

echo "--------------------------------------------------------"
echo "Total de categorías únicas encontradas: $(awk -v col="$TARGET_COL" -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i==col)idx=i} NR>1{print $idx}' "$META_FILE" | sort | uniq | wc -l)"