import polars as pl
import os
import glob
import sys
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
COUNTS_FILE = os.path.join(PROJECT_ROOT, 'filtered_datasets', 'filtered_counts.txt')

# prametros
P_VAL_ADJ_THRESHOLD = 0.01
LOGFC_ABS_THRESHOLD = 3.0


def load_and_classify_counts():
    """Carga la matriz de counts y clasifica muestras una sola vez."""
    print("Cargando matriz de conteos filtrada...")
    if not os.path.exists(COUNTS_FILE):
        print(f"ERROR: No se encuentra el archivo de conteos en: {COUNTS_FILE}")
        sys.exit(1)

    try:
        counts_df = pl.read_csv(COUNTS_FILE, separator='\t', null_values=["NA", "nan"])

        # Normalizar nombre columna Gen (normalmente viene vacia)
        first_col = counts_df.columns[0]
        if first_col != "Gene":
            counts_df = counts_df.rename({first_col: "Gene"})

        # Clasificación de muestras sanas y enfermas
        # Uso prefijos GTEX y TCGA para determinar sano o enfermo
        all_samples = [c for c in counts_df.columns if c != "Gene"]
        gtex_samples = [s for s in all_samples if s.startswith("GTEX")]
        tcga_samples = [s for s in all_samples if s.startswith("TCGA")]

        print(f"Matriz cargada: {counts_df.height} genes. Muestras: {len(gtex_samples)} Sanas (GTEX), {len(tcga_samples)} Tumor (TCGA).")
        return counts_df, gtex_samples, tcga_samples

    except Exception as e:
        print(f"Error cargando matriz: {e}")
        sys.exit(1)


def process_dea_file(file_path, counts_df, gtex_samples, tcga_samples):
    filename = os.path.basename(file_path)
    # Defino nombre reemplazando el sufijo '_results.csv' por '_metrics.csv'
    output_path = file_path.replace('_results.csv', '_metrics.csv')

    print(f"\nProcesando: {filename}")

    # Cargo y filtro Genes DEA
    try:
        dea_df = pl.read_csv(
            file_path,
            infer_schema_length=10000,
            schema_overrides={
                "P.Value": pl.Float64,
                "adj.P.Val": pl.Float64,
                "logFC": pl.Float64,
                "t": pl.Float64,
                "B": pl.Float64,
                "AveExpr": pl.Float64
            }
        )
        # Normalización de columna Gene
        if dea_df.columns[0] == "":
            dea_df = dea_df.rename({"": "Gene"})
        elif "Gene" not in dea_df.columns:
            dea_df = dea_df.rename({dea_df.columns[0]: "Gene"})

        # Filtros
        dea_filtered = dea_df.filter(
            (pl.col("adj.P.Val") < P_VAL_ADJ_THRESHOLD) &
            (pl.col("logFC").abs() >= LOGFC_ABS_THRESHOLD)
        )

        genes_passed = dea_filtered["Gene"].to_list()
        print(f"  -> Genes significativos: {len(genes_passed)}")

        if len(genes_passed) == 0:
            print("  -> Ningun gen paso los filtros. Revisar Umbrales.")
            return

    except Exception as e:
        print(f"  Error leyendo CSV DEA: {e}")
        return

    # Calcular media y mediana
    # Me quedo solo con genes DEA
    subset_expr = counts_df.filter(pl.col("Gene").is_in(genes_passed))
    # Calculo Media y Mediana de sanos y de enfermos
    stats = subset_expr.select(
        pl.col("Gene"),
        pl.concat_list(gtex_samples).list.mean().round(6).alias("Mean_GTEX_Healthy"),
        pl.concat_list(gtex_samples).list.median().round(6).alias("Median_GTEX_Healthy"),
        pl.concat_list(tcga_samples).list.mean().round(6).alias("Mean_TCGA_Tumor"),
        pl.concat_list(tcga_samples).list.median().round(6).alias("Median_TCGA_Tumor")
    )

    # Guardo
    try:
        stats.write_csv(output_path)
        print(f"  -> Archivo generado: {os.path.basename(output_path)}")
    except Exception as e:
        print(f"  Error escribiendo archivo: {e}")


def main():
    # Configuración de argumentos de línea de comando
    parser = argparse.ArgumentParser(description="Genera metricas de expresión (media/mediana) para genes DEA.")
    parser.add_argument("input_dir", type=str, help="Carpeta raiz donde buscar los archivos _results.csv")
    args = parser.parse_args()
    input_root_dir = args.input_dir
    if not os.path.isdir(input_root_dir):
        print(f"ERROR: La carpeta '{input_root_dir}' no existe.")
        sys.exit(1)

    # Cargar Datos
    counts_df, gtex_samples, tcga_samples = load_and_classify_counts()

    # Busco archivo _results.csv
    results_file = None
    print(f"Buscando archivos en: {input_root_dir}")
    search_pattern = os.path.join(input_root_dir, "**", "*_results.csv")
    files = glob.glob(search_pattern, recursive=True)
    for f in files:
        fname = os.path.basename(f)
        if fname.endswith("top50_results.csv"):
            continue
        results_file = f

    # Proceso datos
    process_dea_file(results_file, counts_df, gtex_samples, tcga_samples)

    print("\nFinalizado!")


if __name__ == "__main__":
    main()