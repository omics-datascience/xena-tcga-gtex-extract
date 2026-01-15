# TCGA & GTEx RNA-Seq Analysis Pipeline

Este repositorio contiene un pipeline automatizado en Bash, Python y R para realizar Análisis de Expresión Diferencial (DEA) comparando cohortes de **TCGA** (Cáncer) y **GTEx** (Tejido normal).

El flujo de trabajo permite limpiar los datos (eliminando muestras pediátricas de TARGET), explorar las categorías disponibles, filtrar subconjuntos de interés y ejecutar el análisis estadístico utilizando **Limma**.

## Estructura del Proyecto

```text
.
├── cohort_TCGA_TARGET_GTEx/
│   └── delete_target_samples.sh            # Script de limpieza de datasets
├── filtered_datasets/                      # Directorio de salida para matrices filtradas
├── DEA_limma/                              # Directorio con scripts R para realizar el analisis de expresion diferencial
├── DEA_output/                             # Directorio de salida para resultados de Limma
├── listar_categorias.sh                    # Script de exploración de metadatos
├── filtrar_samples.py                      # Script de generación de datasets y cambios de ID de genes de Ensembl a HUGO Symbol
├── run_dea.sh                              # Script de ejecución del análisis (Limma) para un dataset filtrado
├── run_dea.sh                              # Script de ejecución del análisis (Limma) para un dataset filtrado
├── map_genes.py                            # Convierte Identificadores de genes Ensembl a nomenclatura HUGO
├── run_dea_for_all_tissues_combinations.sh # Script automatizado para hacer DEA con todas las combinaciones descritas en el archivo all_tissues_combinations.tsv
├── all_tissues_combinations.tsv            # Archivo von todas las combinaciones logicas de tegidos sanos vs enfermos en el dataset
├── Interpretacion_de_resutlados.md         # Docuementacion: Explicacion de graficas de resultados
└── README.md                               # Documentacion
```

## Requisitos

- Entorno Unix/Linux (Bash) con Herramientas estándars awk, sed, grep
- R 4.5
- Python 3.12 o superior

Librerias R: Instalar usando:

```bash
Rscript DEA_limma/requirements/check_and_install_packages.r
```

Librerias Python: Instalar usando:

```bash
pip3 install polars pandas
```

## Descarga y Preprocesamiento del dataset

Antes de empezar a usar las herramientas de este repositorio, es necesario descargar el dataset que combina los datos de TCGA, GTEx y TARGET desde UCSC XENA Browser.

### Descargar

1. Descargar los siguientes archivos desde [A combined cohort of TCGA, TARGET and GTEx samples](https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443):
   1. [RSEM expected_count (n=19,109) UCSC Toil RNA-seq Recompute](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz)
   2. [ID Gene Mapping](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap)
   3. [TCGA GTEX main categories (n=17,221) UCSC Toil RNA-seq Recompute](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_GTEX_category.txt)
   4. [TCGA TARGET GTEX selected phenotypes (n=19.131) UCSC Toil RNA-seq Recompute](https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz)
2. Descomprimir en caso de que haga falta
3. Mover los 4 archivos dentro del directorio **cohort_TCGA_TARGET_GTEx**

### Procesar: Eliminar TARGET

El siguiente paso es correr por única vez el archivo delete_target_samples.sh para eliminar todas las muestras correspondientes a la base de datos TARGET (que viene incluida en este datasets). A fines de comparacion con analisis de expresion diferencial, no queremos incluir las muestras de origen pediatrico de esta base de datos.

```bash
cd cohort_TCGA_TARGET_GTEx
bash delete_target_samples.sh
```

Esto creara el archivo *expected_counts_without_TARGET_samples.tsv* el cual usaremos en todos los analisis de expresion diferencial para obtener las muestras correspondientes a la categoria deseada.

## Flujo del Pipeline

Sigue estos pasos secuenciales para ejecutar el análisis:

1. Exploración de Categorías
Para definir qué grupos comparar, utiliza este script que lista todas las categorías disponibles en los metadatos. El output muestra el recuento de muestras por tejido/enfermedad para GTEx y TCGA.

    ```Bash
    bash listar_categorias.sh
    ```

2. Generación del Dataset
Una vez decididos los grupos a comparar (basado en el paso anterior), utiliza este script para crear la matriz de conteos filtrada. El script acepta múltiples categorías como argumentos.  
Sintaxis: bash filtrar_samples.py "CATEGORIA_1" "CATEGORIA_2" ...  
Ejemplo de uso (Adenocarcinoma de Colon vs. Colon Normal):  

    ```Bash
    bash filtrar_samples.py "TCGA Colon Adenocarcinoma" "GTEX Colon"
    ```

    Notas:  

    - Este paso generará dos archivos (filtered_metadata.txt y filtered_counts.txt) listos para el análisis de expresion - diferencial en la carpeta **filtered_datasets/**.  
    - El script tambien cambia los identificadores de genes desde la nomenclatura ENSEMBL a la nomenclatura HUGO Symbol. Para ello usa al archivo [cohort_TCGA_TARGET_GTEx/probeMap_gencode.v23.annotation.gene.probemap](cohort_TCGA_TARGET_GTEx/probeMap_gencode.v23.annotation.gene.probemap) como referencias de mapeo.

3. Análisis de Expresión Diferencial (DEA)
Finalmente, ejecuta el análisis estadístico. Este script toma el dataset generado en el paso anterior y utiliza Limma para encontrar genes diferencialmente expresados.

    ```Bash
    bash run_dea.sh 'run-name'
    ```

'run-name' es un nombre para la corrida. Los resultados finales (PCA, tablas de genes, volcanoplots y heatmaps) se guardarán automáticamente en la carpeta DEA_output/run-name.

## Procesameinto de a pares de tejidos sanos vs enfermos

Alternativamente al proceso de 3 pasos descrito anteriormente, se creo el archivo *all_tissues_combinations.tsv* que contiene las 24 combinaciones logicas de comparaciones entre tegidos sanos y con cancer (para no comparar, por ejemplo, tegido sano de colon con tejido enfermo de mama). Usando el script *run_dea_for_all_tissues_combinations.sh* se corre el flujo completo de seleccion de categorias (obteniendolas del archivo antes mencionado), filtrado y analisis de expresion diferencial. Los resultados se guardaran en la carpeta DEA_output. 

```bash
bash run_dea_for_all_tissues_combinations.sh
```

## Resultados

Luego de elegir las categorias, filtrar las muestras y hacer el DEA, los resultados se almacenan en la carpeta DEA_output dentro de una nueva carpeta con el nombre usado como parametro al correr el script run_dea.sh.  

Dentro de este directorio encontraremos:  

- Archivo de control de calidad del dataset usado en el DEA: Contiene el grafico del analisis de componentes principales (PCA)
- Archivos de resultados:
  - Grafico VolcanoPlot
  - Grafico HeatMap
  - Tabla de valores estadisticos resultados del DEA (ordenados en primera instancia por p valor ajustado y en segunda instancia por valor absoluto de LogFC)
  - Tabla de valores estadisticos resultados del DEA con el top 50 de genes diferencialmente expresados.

En el archivo [Interpretacion de resultados](Interpretacion_de_resutlados.md) se explica como se genero y como interpretar cada resultado.
