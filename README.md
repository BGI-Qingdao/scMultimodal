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

To run label transfer:
```angular2html
Rscript label_transfer.R --rna scrna.anno.rds \
    --atac integrated.atac.rds \
    -m Metadata.tsv \
    --fragments fragments.tsv.gz \
    --gff species_gene_annotation.gtf \
    --name projectA \
    -o results_dir/
```

To integrate multiple ATAC datasets (by different libraries):
```angular2html
    Rscript integrate.R -n name --list datasets_ids.txt --data_path /path/to/samples/
```
Example of a `datasets_ids.txt`:<br>
```
sample_id_1
sample_id_2
sample_id_3
...
```
Example of a data path:<br>
```
 |--/path/to/samples/
    |--sample_id_1
       |--Peak_matrix
          |--barcodes.tsv
          |--matrix.mtx
          |--peak.bed
    |--sample_id_2
       |--Peak_matrix
          |--barcodes.tsv
          |--matrix.mtx
          |--peak.bed
```
