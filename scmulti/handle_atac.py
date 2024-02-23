#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 20 Feb 2024 17:02
# @Author: Yao LI
# @File: evo_fish/handle_atac.py
import subprocess


class RefGenome:
    def __init__(self, fasta_fn, chrom_size=None, ):
        self.fasta_fn = fasta_fn

        self.chromsize = chrom_size
        self.fasta = None

    def get_chromsize(self):
        pass


class scATAC:
    def __init__(self,
                 atac_fn=None,
                 scrna_fn=None,
                 gff=None,
                 fragment_dir=None,
                 fragment_metadata_fn=None,
                 rscript_path='/dellfsqd3/MGI_QINGDAO/USER/lishuangshuang/software/miniconda3/envs/scRNA_v2.4/lib/R/bin/Rscript',
                 project=None,
                 gene_name_col=None):
        self.atac_fn = atac_fn
        self.scrna_fn = scrna_fn
        self.gff = gff
        self.fragment_dir = fragment_dir
        self.metadata_fn = fragment_metadata_fn
        self.rscript_path = rscript_path
        self.name = project
        self.gene_name_col = gene_name_col

    def transfer_label(self):
        subprocess.call(f'{self.rscript_path} label_transfer.R '
                        f'--rna {self.scrna_fn} '
                        f'--peaks {self.atac_fn} '
                        f'--fragments {self.fragment_dir} '
                        f'--metadata_fn {self.metadata_fn} '
                        f'--gff {self.gff} '
                        f'--name {self.name} '
                        f'--gene_name_col {self.gene_name_col}',
                        shell=True)

    def get_atac_rds(self):
        subprocess.call(f'{self.rscript_path} 01.get_atac_rds.R '
                        f'-d {self.fragment_dir} '
                        f'-m {self.metadata_fn}',
                        shell=True)

    def rds2h5ad(self, fn):
        subprocess.call(f'{self.rscript_path} rds2h5ad.R {fn}', shell=True)
