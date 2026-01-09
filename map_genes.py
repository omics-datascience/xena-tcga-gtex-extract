import pandas as pd
import argparse
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Reemplaza IDs por Gene Symbols usando un TSV y guarda el resultado.'
    )
    parser.add_argument(
        '-db', '--database',
        required=True,
        help='Ruta al archivo TSV de referencia (debe contener columnas "id" y "gene")'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Ruta al archivo CSV de resultados de expresión diferencial'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=False,
        default=None,
        help='Directorio donde guardar el archivo resultante (Opcional). Si no se indica, se guarda junto al input.'
    )

    args = parser.parse_args()

    if not os.path.exists(args.database):
        sys.exit(f"Error: No se encontró el archivo de base de datos: {args.database}")
    if not os.path.exists(args.input):
        sys.exit(f"Error: No se encontró el archivo de entrada: {args.input}")

    print(f"Procesando: {args.input} ...")
    db_df = pd.read_csv(args.database, sep='\t')

    # Diccionario de mapeo {id: gene}
    mapping_dict = dict(zip(db_df['id'], db_df['gene']))

    results_df = pd.read_csv(args.input)
    # Identificar columna ID (primera columna)
    id_col_name = results_df.columns[0]
    
    # Reemplazar IDs por Gene Symbols
    results_df['gene'] = results_df.iloc[:, 0].map(mapping_dict)
    
    # Eliminar columna vieja y reordenar 'gene' al principio
    results_df = results_df.drop(columns=[id_col_name])
    cols = ['gene'] + [c for c in results_df.columns if c != 'gene']
    results_df = results_df[cols]

    # Obtenemos el nombre del archivo sin la ruta y sin la extensión
    input_dir, input_filename = os.path.split(args.input)
    filename_no_ext, ext = os.path.splitext(input_filename)
    # Definimos el nuevo nombre del archivo con el sufijo gene_symbols
    new_filename = f"{filename_no_ext}_gene_symbols{ext}"

    if args.output_dir:
        # Si se especificó un directorio de salida
        out_dir = args.output_dir
        os.makedirs(out_dir, exist_ok=True)
        output_path = os.path.join(out_dir, new_filename)
    else:
        output_path = os.path.join(input_dir, new_filename)

    results_df.to_csv(output_path, index=False)

    print("-" * 30)
    print("Finalizado.")
    print(f"Archivo guardado en: {output_path}")


if __name__ == "__main__":
    main()