#!/bin/bash

# Argumentos: se correspnden con los umbrales para decir que un gen fue DEA
# Deben ser los mismos que utilizamos en generate_expression_metrics.py
P_VAL=$1
LOGFC=$2
if [ -z "$P_VAL" ] || [ -z "$LOGFC" ]; then
    echo "Error: Faltan argumentos."
    echo "Uso: $0 <p-value> <logfc>"
    exit 1
fi

echo "Thresholds for significative genes:"
echo "p.value.adj threshold = $P_VAL"
echo "|LogFC| threshold = $LOGFC"
echo ""

# Columnas
printf "%-60s %s\n" "COMPARACION" "GENES"
printf "%-60s %s\n" "-------" "-----"

for f in ../DEA_output/*/DEA_output/*_metrics.csv; do
    # Extraer nombre carpeta (se corresponde con la comparacion sano enfermo)
    folder=$(echo "$f" | cut -d'/' -f3)
    lines=$(wc -l < "$f")
    # Resto linea de header
    genes=$((lines - 1))
    printf "%-60s %d\n" "$folder" "$genes"
done