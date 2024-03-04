# scMultimodal

## Overview
To run scmulti:
```angular2html
python run.py -o results_dir \
        --rna scrna.anno.h5ad \
        --rds scrna.anno.rds \
        --atac atac.anno.h5ad \
        --name projectA \
        --fragment fragments.tsv.gz \
        --meta Metadata.tsv \
        --gtf species_gene_annotation.gtf \
        --fasta reference_genome.fa \
        --chromsize sizes.genome \
        --upstream 0 1500 \
        --downstream 100 15000
```
