#!/bin/bash

cd DEA_limma
echo "Iniciando análisis de expresión diferencial con limma..."
Rscript main_analysis.R
cd ..
rm -rf DEA_output/
mv DEA_limma/DEA_output/ DEA_output/
echo "Análisis de expresión diferencial completado. Resultados en la carpeta 'DEA_output/'."