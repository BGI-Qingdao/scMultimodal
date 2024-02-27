#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 20 Feb 2024 17:02
# @Author: Yao LI
# @File: evo_fish/handle_atac.py
import os
import subprocess
import pyranges as pr


class ComputeParams:
    def __init__(self,
                 n_topics=None,
                 n_cpu=5,
                 n_iter=500,
                 random_state=555,
                 alpha=50,
                 alpha_by_topic=True,
                 eta=0.1,
                 eta_by_topic=False,
                 select_model=16,
                 return_model=True,
                 metrics=None,
                 plot_metrics=False,
                 input_format='AUTO',
                 shift=73,
                 ext_size=146,
                 keep_dup='all',
                 q_value=0.05,
                 peak_half_width=250
                 ):
        if metrics is None:
            metrics = ['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood']
        if n_topics is None:
            n_topics = [2, 4, 10, 16, 32, 48]
        self.n_topics = n_topics,
        self.n_cpu = n_cpu,
        self.n_iter = n_iter,
        self.random_state = random_state,
        self.alpha = alpha,
        self.alpha_by_topic = alpha_by_topic,
        self.eta = eta,
        self.eta_by_topic = eta_by_topic,
        self.select_model = select_model,
        self.return_model = return_model,
        self.metrics = metrics,
        self.plot_metrics = plot_metrics,
        self.input_format = input_format,
        self.shift = shift,
        self.ext_size = ext_size,
        self.keep_dup = keep_dup,
        self.q_value = q_value,
        self.peak_half_width = peak_half_width


class ScMultiParams:
    def __init__(self, work_dir=None, output_dir=None):
        self.narrowPeaks_path = None
        self.consensus_bed_path = None
        self.bed_paths = None
        self.bw_paths = None
        self.output_dir = output_dir
        self.work_dir = work_dir
        self.path_to_regions = None


class RefGenome(ScMultiParams):
    def __init__(self, fasta_fn, chromsize_fn=None, work_dir=None, output_dir=None):
        super().__init__(work_dir, output_dir)
        self.fasta_fn = fasta_fn
        self.chromsize_fn = chromsize_fn

    def get_chromsize(self):
        subprocess.call(f'/dellfsqd2/ST_OCEAN/USER/liyao1/tools/anaconda3/envs/scenicplus/bin/faidx '
                        f'{self.fasta_fn} -i chromsizes > {self.work_dir}/sizes.genome', shell=True)
        # TODO: check if a process was finished
        if self.chromsize_fn is None:
            self.chromsize_fn = f'{self.work_dir}/sizes.genome'


class scATAC(ScMultiParams):
    def __init__(self,
                 atac_fn=None,
                 scrna_fn=None,
                 gff=None,
                 fragment_dir=None,
                 fragment_metadata_fn=None,
                 rscript_path='/dellfsqd3/MGI_QINGDAO/USER/lishuangshuang/software/miniconda3/envs/scRNA_v2.4/lib/R/bin/Rscript',
                 project=None,
                 gene_name_col=None,
                 work_dir=None,
                 output_dir=None):
        super().__init__(work_dir, output_dir)
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

    def rds2h5ad(self, fn: str):
        subprocess.call(f'{self.rscript_path} rds2h5ad.R {fn}', shell=True)
