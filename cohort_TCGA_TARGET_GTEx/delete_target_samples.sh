
awk -F'\t' '$7 == "TARGET" {print $1}' TcgaTargetGTEX_phenotype.txt > target_sample_ids.txt
EXPECTED_COUNTS="TcgaTargetGtex_gene_expected_count"
echo "Número de muestras TARGET a eliminar:"
wc -l target_sample_ids.txt
awk -F'\t' '
    # 1. Leer el archivo de IDs a borrar (FNR==NR indica que es el primer archivo)
    FNR==NR {
        muestras_borrar[$1] = 1; 
        next
    }

    # 2. Procesar el archivo de conteos (segundo archivo)
    {
        linea_salida = ""
        
        # Recorremos todas las columnas de la fila actual
        for(i=1; i<=NF; i++) {
            
            # Si estamos en la cabecera (fila 1), chequeamos qué columnas eliminar
            if (FNR == 1) {
                # Si el nombre de la columna está en nuestra lista negra, marcamos el índice "i"
                if ($i in muestras_borrar) {
                    columnas_excluidas[i] = 1
                }
            }

            # Si la columna "i" NO está marcada como excluida, la imprimimos
            if (! (i in columnas_excluidas)) {
                # Añadimos tabulador solo si no es el primer elemento de la línea
                printf "%s%s", (linea_salida == "" ? "" : "\t"), $i
                linea_salida = "iniciado"
            }
        }
        # Salto de línea al terminar la fila
        printf "\n"
    }
' target_sample_ids.txt "$EXPECTED_COUNTS" > expected_counts_without_TARGET_samples.tsv
rm target_sample_ids.txt
echo "Archivo de conteos esperado sin muestras TARGET creado: expected_counts_without_TARGET_samples.tsv"
echo "Número de muestras en el nuevo archivo:"
head -1 expected_counts_without_TARGET_samples.tsv | tr '\t' '\n' | wc