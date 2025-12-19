# TCGA & GTEx RNA-Seq Analysis Pipeline

https://github.com/omics-datascience/xena-tcga-gtex-extract

Este repositorio contiene un pipeline automatizado en Bash y R para realizar Análisis de Expresión Diferencial (DEA) comparando cohortes de **TCGA** (Cáncer) y **GTEx** (Tejido normal).

El flujo de trabajo permite limpiar los datos (eliminando muestras pediátricas de TARGET), explorar las categorías disponibles, filtrar subconjuntos de interés y ejecutar el análisis estadístico utilizando **Limma**.

## Estructura del Proyecto

```text
.
├── cohort_TCGA_TARGET_GTEx/
│   └── delete_target_samples.sh   # Script de limpieza de datasets
├── filtered_datasets/             # Directorio de salida para matrices filtradas
├── DEA_limma/                     # Directorio con scripts R para realizar el analisis de expresion diferencial
├── DEA_output/                    # Directorio de salida para resultados de Limma
├── listar_categorias.sh           # Script de exploración de metadatos
├── filtrar_samples.sh             # Script de generación de datasets
├── run_dea.sh                     # Script de ejecución del análisis (Limma)
└── README.md
```

## Uso del Pipeline

Sigue estos pasos secuenciales para ejecutar el análisis:

1. Preprocesamiento de Datos
Antes de comenzar, es necesario limpiar el dataset original eliminando las muestras pertenecientes a la base de datos TARGET (muestras pediátricas), para trabajar únicamente con TCGA y GTEx.

    ```Bash
    cd cohort_TCGA_TARGET_GTEx
    bash delete_target_samples.sh
    ```

2. Exploración de Categorías
Para definir qué grupos comparar, utiliza este script que lista todas las categorías disponibles en los metadatos. El output muestra el recuento de muestras por tejido/enfermedad para GTEx y TCGA.

    ```Bash
    bash listar_categorias.sh
    ```

3. Generación del Dataset
Una vez decididos los grupos a comparar (basado en el paso anterior), utiliza este script para crear la matriz de conteos filtrada. El script acepta múltiples categorías como argumentos.  
Sintaxis: bash filtrar_samples.sh "CATEGORIA_1" "CATEGORIA_2" ...  
Ejemplo de uso (Adenocarcinoma de Colon vs. Colon Normal):  

    ```Bash
    bash filtrar_samples.sh "TCGA Colon Adenocarcinoma" "GTEX Colon"
    ```

    Nota: Este paso generará dos archivos (filtered_metadata.txt y filtered_counts.txt) listos para el análisis de expresion diferencial en la carpeta filtered_datasets/.  

4. Análisis de Expresión Diferencial (DEA)
Finalmente, ejecuta el análisis estadístico. Este script toma el dataset generado en el paso anterior y utiliza Limma para encontrar genes diferencialmente expresados.

    ```Bash
    bash run_dea.sh
    ```

Los resultados finales (tablas de genes y volcanoplots) se guardarán automáticamente en la carpeta DEA_output/.

## Requisitos

- Entorno Unix/Linux (Bash) con Herramientas estándars awk, sed, grep
- R 4.5

Librerias R: Instalar usando:

```R
Rscript DEA_limma/requirements/check_and_install_packages.r
```
