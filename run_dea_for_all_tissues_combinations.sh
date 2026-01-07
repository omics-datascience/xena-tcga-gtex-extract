#!/bin/bash

ARCHIVO="all_tissues_combinations.tsv"

if [ ! -f "$ARCHIVO" ]; then
    echo "Error: El archivo $ARCHIVO no existe."
    exit 1
fi

# Leemos el archivo línea por línea
# tail -n +2 sirve para saltar la primera línea (encabezados)
# IFS=$'\t' define el delimitador como un tabulador
tail -n +2 "$ARCHIVO" | while IFS=$'\t' read -r sano cancer; do
    
    # Limpiar posibles espacios en blanco extra o saltos de línea
    sano=$(echo "$sano" | tr -d '\r' | xargs)
    cancer=$(echo "$cancer" | tr -d '\r' | xargs)

    # Si por alguna razón la línea está vacía, saltar
    [ -z "$sano" ] && continue

    echo "Procesando: $sano vs $cancer"

    # 1. Ejecutar el script de Python
    python3 filtrar_samples.py "$sano" "$cancer"

    # 2. Ejecutar el script de Bash con el formato "sano vs cancer"
    bash run_dea.sh "$sano vs $cancer"

    echo "Finalizado: $sano vs $cancer"
    echo "---------------------------------------"

done