#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 27 Feb 2024 10:08
# @Author: Yao LI
# @File: evo_fish/custom_databases.py
import os
import subprocess
import pandas as pd
from typing import Union, Optional
import argparse
import sys


class CustomDatabase:
    def __init__(self,
                 prefix,
                 output_dir=None,
                 fasta=None,
                 sub_fasta_fn: str = 'consensus_regions.fasta',
                 motifs_dir='/dellfsqd2/ST_OCEAN/USER/liyao1/tools/custom_cistarget/motif_cb',
                 motifs_list_filename='/dellfsqd2/ST_OCEAN/USER/liyao1/tools/custom_cistarget/motifs.lst',
                 nbr_threads=48,
                 bedtool_path='/dellfsqd2/ST_OCEAN/USER/liyao1/tools/bedtools',
                 create_cistarget_databases_dir='/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/00.custom/create_cisTarget_databases'):
        self.fasta = fasta
        self.motif_anno = None
        self.gene_anno = None
        self.ranking_db = None
        self.score_db = None
        self.tf_list = None
        self.sub_fasta = sub_fasta_fn
        self.motifs_dir = motifs_dir
        self.motifs_list_filename = motifs_list_filename
        self.db_prefix = prefix
        self.nbr_threads = nbr_threads
        self.create_cistarget_databases_dir = create_cistarget_databases_dir
        self.bedtool_path = bedtool_path

        self.human_tfs = '/dellfsqd2/ST_OCEAN/USER/liyao1/06.stereopy/resource/tfs/allTFs_hg38.txt'
        self.human_tbl = '/dellfsqd2/ST_OCEAN/USER/liyao1/06.stereopy/resource/motifs/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
        self.mouse_tfs = '/dellfsqd2/ST_OCEAN/USER/liyao1/06.stereopy/resource/tfs/allTFs_mm.txt'
        self.mouse_tbl = '/dellfsqd2/ST_OCEAN/USER/liyao1/06.stereopy/resource/motifs/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'

        if output_dir is None:
            output_dir = os.getcwd()
        self.saving_dir = os.path.join(output_dir, 'cistarget_database')
        if not os.path.exists(self.saving_dir):
            os.makedirs(self.saving_dir)

        self.species_options = ['hg', 'mm', 'custom']
        self.species_files = {'hg': {'tfs': self.human_tfs, 'tbl': self.human_tbl},
                              'mm': {'tfs': self.mouse_tfs, 'tbl': self.mouse_tbl}}

    def extact_fasta(self, consensus_regions_bed):
        """

        :param consensus_regions_bed:
        :return:
        """
        cmd = f'{self.bedtool_path}  getfasta  -fi  {self.fasta}  -fo {self.saving_dir}/{self.sub_fasta}  -bed  {consensus_regions_bed}'
        subprocess.call(cmd, shell=True)

    def extact_fasta_by_distance(self, chromsize_fn, upstream: int,
                                 downstream: int = 0, exon_fn=None):
        """

        :param chromsize_fn:
        :param upstream:
        :param downstream:
        :param exon_fn:
        :return:
        """
        if exon_fn is None:
            cmd1 = f'{self.bedtool_path} flank -l {upstream} -r {downstream} -g {chromsize_fn} > {self.saving_dir}/upstream-{upstream}.bed'
        else:
            cmd1 = f'{self.bedtool_path} flank -l {upstream} -r {downstream} -i {exon_fn} -g {chromsize_fn} > {self.saving_dir}/upstream-{upstream}.bed'
        cmd2 = f'{self.bedtool_path} getfasta  -fi {self.fasta}  -fo {self.sub_fasta}  -bed {self.saving_dir}/upstream-{upstream}.bed  -nameOnly'
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)

    def create_feathers(self):
        """

        :return:
        """
        # subprocess.run('source /dellfsqd2/ST_OCEAN/USER/liyao1/tools/anaconda3/bin/activate create_cistarget_databases', shell=True)
        cmd = f'/dellfsqd2/ST_OCEAN/USER/liyao1/tools/anaconda3/envs/create_cistarget_databases/bin/python ' \
              f'{self.create_cistarget_databases_dir}/create_cistarget_motif_databases.py ' \
              f'-f {self.saving_dir}/{self.sub_fasta} ' \
              f'-M {self.motifs_dir} ' \
              f'-m {self.motifs_list_filename} ' \
              f'-o {self.saving_dir}/{self.db_prefix} ' \
              f'-t {self.nbr_threads} '
        subprocess.call(cmd, shell=True)
        # subprocess.run('source /dellfsqd2/ST_OCEAN/USER/liyao1/tools/anaconda3/bin/deactivate', shell=True)

    def get_tfs(self, atac_fn, ortholog_groups_fn, ref_species='mm', ref_col: Union[str, int] = 'target_id',
                target_col: Union[str, int] = 'gene_name', header=None, delimiter='\t', ref_tfs=None):
        """

        :param atac_fn:
        :param ortholog_groups_fn:
        :param ref_species:
        :param ref_col:
        :param target_col:
        :param header:
        :param delimiter:
        :param ref_tfs:
        :return:
        """
        if ref_species not in self.species_options:
            print("Invalid species. Expected one of: %s" % self.species_options)
        if ref_species not in list(self.species_files.keys()) and ref_tfs is None:
            raise ValueError("A reference tbl file needs to be provided when species is set to custom.")
        elif ref_species in list(self.species_files.keys()):
            ref_tfs = self.species_files[ref_species]['tfs']

        import scanpy as sc
        adata = sc.read_h5ad(atac_fn)
        Agene_id_list = list(adata.var_names)

        hm = pd.read_csv(ortholog_groups_fn, delimiter=delimiter, header=header)
        sub_hm = hm[hm[ref_col].isin(list(Agene_id_list))]
        ah_anno = sub_hm[sub_hm[target_col] != '-']
        ah_anno[[ref_col, target_col]].to_csv(f'{self.saving_dir}/{self.db_prefix}_{ref_species}_gene_name.txt',
                                              index=False)  # TODO: is this necessary

        with open(ref_tfs, 'r') as f:
            reference_tfs = f.read().splitlines()
        homolog_tfs = list(
            set(ah_anno[ah_anno[target_col].isin(reference_tfs)][ref_col]))  # TODO: better way to get TFs
        self.tf_list = f'{self.db_prefix}_TFs.txt'
        with open(os.path.join(self.saving_dir, self.tf_list), 'w') as f:
            f.writelines('\n'.join(homolog_tfs))

    def make_motif_annotation(self, ref_species='mm', ref_col: Union[str, int] = 'target_id',
                              target_col: Union[str, int] = 'gene_name', ref_tbl_fn=None):
        """

        :param ref_species:
        :param ref_col:
        :param target_col:
        :param ref_tbl_fn:
        :return:
        """
        if ref_species not in self.species_options:
            print("Invalid species. Expected one of: %s" % self.species_options)
        target_homolog_fn = f'{self.saving_dir}/{self.db_prefix}_{ref_species}_gene_name.txt'

        if ref_species not in list(self.species_files.keys()) and ref_tbl_fn is None:
            raise ValueError("A reference tbl file needs to be provided when species is set to custom.")
        elif ref_species in list(self.species_files.keys()):
            ref_tbl_fn = self.species_files[ref_species]['tbl']

        ref_tbl_df = pd.read_csv(ref_tbl_fn, sep='\t', header=0, dtype=str, na_values='')
        columns = ref_tbl_df.columns.values

        homolog = pd.read_csv(target_homolog_fn, dtype=str, na_values='')  # TODO: rely on homology genes...
        homolog.columns = [ref_col, target_col]
        homolog = homolog.dropna(axis='index')

        ref_tbl_df = ref_tbl_df.merge(homolog, how='left', left_on='gene_name', right_on=target_col)
        ref_tbl_df = ref_tbl_df.dropna(axis='index', subset=ref_col)
        ref_tbl_df = ref_tbl_df.drop(axis='columns', labels=target_col)
        ref_tbl_df = ref_tbl_df.rename(columns={ref_col: target_col})
        ref_tbl_df = ref_tbl_df[columns]
        self.motif_anno = f'motifs-v10-nr.{self.db_prefix}-m0.001-o0.0.tbl'
        ref_tbl_df.to_csv(os.path.join(self.saving_dir, self.motif_anno), sep='\t', index=False)


def create_custom_database(prefix, fasta, output_dir, consensus_regions_path, atac_fn, ortholog_groups_fn,
                           ref_species='hg', bed='/dellfsqd2/ST_OCEAN/USER/liyao1/tools/bedtools'):
    cdb = CustomDatabase(prefix,
                         fasta=fasta,
                         output_dir=output_dir, bedtool_path=bed)
    cdb.extact_fasta(consensus_regions_path)
    cdb.get_tfs(atac_fn, ortholog_groups_fn, ref_species=ref_species, ref_col=0, target_col=1, header=None)
    cdb.make_motif_annotation(ref_species=ref_species, ref_col=0, target_col=1)
    cdb.create_feathers()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='scmulti',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='',
        add_help=True,
    )
    parser.add_argument('-p', '--prefix', help='species name')
    parser.add_argument('--fasta', help='reference genome fasta file')
    parser.add_argument('-o', '--output', help='output directory to save all the intermediate and final results.')
    parser.add_argument('-c', '--consensus_regions', help='consensus_regions.bed file')
    parser.add_argument('--ortholog',
                        help='files contains target species gene name/id and reference gene name/id (hg or mm)')
    parser.add_argument('-ref', '--ref_species', default='hg', help='human or mouse')
    parser.add_argument('--atac', help='scATAC-seq data h5ad file')
    parser.add_argument('--bedpath', default='/dellfsqd2/ST_OCEAN/USER/liyao1/tools/bedtools', help='path to bedtools')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    create_custom_database(args.prefix, args.fasta, args.output, args.consensus_regions, args.atac, args.ortholog,
                           ref_species=args.ref_species, bed=args.bedpath)
