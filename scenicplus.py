#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 18 Dec 2023 16:58
# @Author: Yao LI
# @File: evo_fish/scenicplus.py
import os
import scanpy as sc
import numpy as np
from pycisTopic.utils import (
    read_fragments_from_file,
    prepare_tag_cells
)
import pyranges as pr
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, export_pseudobulk_one_sample


def get_chromsizes(size_fn):
    """

    :param size_fn:
    :return: PyRanges object
    """
    import pyranges as pr
    import pandas as pd

    chromsizes = pd.read_csv(size_fn, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    # Exceptionally in this case, to agree with CellRangerARC annotations
    chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in
                                range(len(chromsizes['Chromosome']))]
    chromsizes['Chromosome'] = [
        chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else
        chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
    chromsizes = pr.PyRanges(chromsizes)
    return chromsizes


def main(fn, sample_id, bed_fn, size_fn, output_dir, tmp_dir, n_cpu: int = 12, anno_label='celltype'):
    adata = sc.read_h5ad(fn)
    cell_data = adata.obs
    cell_data['sample_id'] = sample_id
    cell_data[anno_label] = cell_data[anno_label].astype(str)

    fragments_dict = {sample_id: bed_fn}

    chromsizes = get_chromsizes(size_fn)

    # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
    bw_paths, bed_paths = export_pseudobulk(input_data=cell_data,
                                            variable=anno_label,
                                            # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                                            sample_id_col='sample_id',
                                            chromsizes=chromsizes,
                                            bed_path=os.path.join(output_dir, 'pseudobulk_bed_files/'),
                                            # specify where pseudobulk_bed_files should be stored
                                            bigwig_path=os.path.join(output_dir, 'pseudobulk_bw_files/'),
                                            # specify where pseudobulk_bw_files should be stored
                                            path_to_fragments=fragments_dict,  # location of fragment fiels
                                            n_cpu=n_cpu,
                                            # specify the number of cores to use, we use ray for multi processing
                                            normalize_bigwig=True,
                                            remove_duplicates=True,
                                            _temp_dir=os.path.join(tmp_dir, sample_id),
                                            split_pattern='-')


def get_pseudobulk(cell_data, anno_label, chromsizes, sample_id, output_dir, fragments_dict, n_cpu=12):
    # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
    from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
    bw_paths, bed_paths = export_pseudobulk(input_data=cell_data,
                                            variable=anno_label,
                                            # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                                            sample_id_col='sample_id',
                                            chromsizes=chromsizes,
                                            bed_path=os.path.join(output_dir, 'pseudobulk_bed_files/'),
                                            # specify where pseudobulk_bed_files should be stored
                                            bigwig_path=os.path.join(output_dir, 'pseudobulk_bw_files/'),
                                            # specify where pseudobulk_bw_files should be stored
                                            path_to_fragments=fragments_dict,  # location of fragment fiels
                                            n_cpu=n_cpu,
                                            # specify the number of cores to use, we use ray for multi processing
                                            normalize_bigwig=True,
                                            remove_duplicates=True,
                                            split_pattern='-')
    return bw_paths, bed_paths


def call_peaks(macs_path, bed_paths, work_dir):
    from pycisTopic.pseudobulk_peak_calling import peak_calling
    macs_path = 'macs2'
    # Run peak calling
    narrow_peaks_dict = peak_calling(macs_path,
                                     bed_paths,
                                     os.path.join(work_dir, 'scATAC/consensus_peak_calling/MACS/'),
                                     genome_size='hs',
                                     n_cpu=8,
                                     input_format='BEDPE',
                                     shift=73,
                                     ext_size=146,
                                     keep_dup='all',
                                     q_value=0.05)  # _temp_dir=os.path.join(tmp_dir, 'ray_spill')


def one_sample(group, chromsizes, cell_data, output_dir, split_pattern='-', sample_id_col='sample_id', variable='celltype'):
    fragments_df_dict = {}
    for sample_id in fragments_dict.keys():
        fragments_df = read_fragments_from_file(
            fragments_dict[sample_id],
            use_polars=True
        ).df
        fragments_df.Start = np.int32(fragments_df.Start)
        fragments_df.End = np.int32(fragments_df.End)
        if "Score" in fragments_df:
            fragments_df.Score = np.int32(fragments_df.Score)
        if "barcode" in cell_data:
            fragments_df = fragments_df.loc[
                fragments_df["Name"].isin(cell_data["barcode"].tolist())
            ]
        else:
            fragments_df = fragments_df.loc[
                fragments_df["Name"].isin(
                    prepare_tag_cells(cell_data.index.tolist(), split_pattern)
                )
            ]
        fragments_df_dict[sample_id] = fragments_df
    if "barcode" in cell_data:
        cell_data = cell_data.loc[:, [variable, sample_id_col, "barcode"]]
    else:
        cell_data = cell_data.loc[:, [variable, sample_id_col]]
    cell_data[variable] = cell_data[variable].replace(
        " ", "", regex=True)
    cell_data[variable] = cell_data[variable].replace(
        "[^A-Za-z0-9]+", "_", regex=True)
    groups = sorted(list(set(cell_data[variable])))
    if isinstance(chromsizes, pd.DataFrame):
        chromsizes = chromsizes.loc[:, ["Chromosome", "Start", "End"]]
        chromsizes = pr.PyRanges(chromsizes)

    try:
        export_pseudobulk_one_sample(
            cell_data=cell_data,
            group=group,
            fragments_df_dict=fragments_df_dict,
            chromsizes=chromsizes,
            bigwig_path=os.path.join(output_dir, 'pseudobulk_bw_files/'),
            bed_path=os.path.join(output_dir, 'pseudobulk_bed_files/'),
            sample_id_col=sample_id_col,
            normalize_bigwig=True,
            remove_duplicates=True,
            split_pattern=split_pattern
        )
    except ZeroDivisionError:
        print(f'{group}: ZeroDivisionError')


def preprocess(adata1):
    sc.pp.filter_cells(adata1, min_genes=200)
    sc.pp.filter_genes(adata1, min_cells=3)
    sc.external.pp.scrublet(adata1)  # estimates doublets
    adata = adata1[adata1.obs['predicted_doublet'] == False]  # do the actual filtering
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    mito_filter = 25
    n_counts_filter = 4300
    adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
    adata = adata[adata.obs.pct_counts_mt < mito_filter, :]
    return adata


if __name__ == '__main__':
    # scRNA-seq
    scrna_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/scrna'
    fn1 = 'a_liver.h5ad'
    adata_rna = sc.read_h5ad(os.path.join(scrna_dir, fn1))
    adata_atac = sc.read_h5ad('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/02.Derived_Lethenteron/merged.atac.h5ad')

    work_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/'
    output_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/atac'

    anno_label = 'predicted.id'
    cell_data = adata_atac.obs
    cell_data['sample_id'] = 'a_liver'
    cell_data[anno_label] = cell_data[anno_label].astype(str)
    # 2023-12-21
    cell_data['barcode'] = list(cell_data.index)
    cell_data['barcode'] = cell_data['barcode'].astype(str)
    # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

    # ATAC-seq fragments
    bed_fn = '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-zhangjin/07.project/5.liver_ATAC/0.qc/QC/lsqsm-youti-gan-2_1/output/lsqsm-youti-gan-2_1.fragments.tsv.gz'
    fragments_dict = {'a_liver': bed_fn}

    # Get chromosome sizes (for 七腮鳗 here)
    size_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/data/Lethenteron_reissneri_sizes.genome'
    chromsizes = pd.read_csv(size_fn, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    chromsizes = pr.PyRanges(chromsizes)

    # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
    celltypes = list(set(adata_atac.obs[anno_label]))
    for ct in celltypes:
        one_sample(ct, chromsizes, cell_data, output_dir, split_pattern=None, variable=anno_label)

    # bw_paths, bed_paths = export_pseudobulk(input_data=cell_data,
    #                                         variable=anno_label,
    #                                         # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
    #                                         sample_id_col='sample_id',
    #                                         chromsizes=chromsizes,
    #                                         bed_path=os.path.join(output_dir, 'pseudobulk_bed_files/'),
    #                                         # specify where pseudobulk_bed_files should be stored
    #                                         bigwig_path=os.path.join(output_dir, 'pseudobulk_bw_files/'),
    #                                         # specify where pseudobulk_bw_files should be stored
    #                                         path_to_fragments=fragments_dict,  # location of fragment fiels
    #                                         n_cpu=20,
    #                                         # specify the number of cores to use, we use ray for multi processing
    #                                         normalize_bigwig=True,
    #                                         remove_duplicates=False,
    #                                         split_pattern='-') #_temp_dir=os.path.join(tmp_dir, 'a_liver')
    #
    # import pickle
    # pickle.dump(bed_paths,
    #             open(os.path.join(output_dir, 'pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
    # pickle.dump(bw_paths,
    #             open(os.path.join(output_dir, 'pseudobulk_bw_files/bw_paths.pkl'), 'wb'))

    # get bed files

