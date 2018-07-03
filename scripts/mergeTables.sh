while read p; do grep -h $p salmonella_representatives* | paste - - >> salmonella_mergedTable; done < representatives_barcodes
