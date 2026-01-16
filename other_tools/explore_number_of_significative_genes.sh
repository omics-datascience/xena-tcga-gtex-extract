#!/bin/bash

# Asignar los argumentos a variables
P_VAL=$1
LOGFC=$2

# Verificar argumentos
if [ -z "$P_VAL" ] || [ -z "$LOGFC" ]; then
    echo "Error: Faltan argumentos."
    echo "Uso: $0 <p-value> <logfc>"
    exit 1
fi

# Imprimir la cabecera con los umbrales
echo "Thresholds for significative genes:"
echo "p.value.adj threshold = $P_VAL"
echo "|LogFC| threshold = $LOGFC"
echo ""

# Imprimir columnas
printf "%-60s %s\n" "CARPETA" "GENES"
printf "%-60s %s\n" "-------" "-----"

for f in ../DEA_output/*/DEA_output/*_metrics.csv; do
    # Extraer el nombre de la carpeta
    folder=$(echo "$f" | cut -d'/' -f3)
    
    lines=$(wc -l < "$f")
    
    # Restar el header para obtener n genes
    genes=$((lines - 1))
    
    printf "%-60s %d\n" "$folder" "$genes"
done