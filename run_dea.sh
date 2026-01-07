#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Error: run-name parameter is mandatory"
    exit 1
fi

run_name=$1
cd DEA_limma
echo "Iniciando análisis de expresión diferencial con limma..."
Rscript main_analysis.R
cd ..
rm -rf DEA_output/"$run_name"
mkdir -p DEA_output/"$run_name"
mv DEA_limma/DEA_output/ DEA_output/"$run_name"/
echo "Análisis de expresión diferencial completado. Resultados en la carpeta DEA_output/$run_name/."