#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 23 Nov 2023 16:07
# @Author: Yao LI
# @File: evo_fish/tutorial.py

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


def one_sample(group, chromsizes, cell_data):
    variable = 'celltype'
    sample_id_col = 'sample_id'
    chromsizes = chromsizes
    split_pattern = '-'
    output_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/02.tutorial/atac'

    print(cell_data)
    fragments_df_dict = {}
    for sample_id in fragments_dict.keys():
        fragments_df = read_fragments_from_file(
            fragments_dict[sample_id],
            use_polars=True
        ).df
        print(fragments_df)
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
        print(fragments_df)
    if "barcode" in cell_data:
        cell_data = cell_data.loc[:, [variable, sample_id_col, "barcode"]]
    else:
        cell_data = cell_data.loc[:, [variable, sample_id_col]]
    cell_data[variable] = cell_data[variable].replace(
        " ", "", regex=True)
    cell_data[variable] = cell_data[variable].replace(
        "[^A-Za-z0-9]+", "_", regex=True)
    print(cell_data)
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


if __name__ == '__main__':
    # scRNA-seq
    work_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/00.pbmc_tutorial'
    # fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/pbmc_tutorial/data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5'
    # # adata1 = sc.read_10x_h5(os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5'))
    # adata1 = sc.read_10x_h5(fn)
    #
    # sc.pp.filter_cells(adata1, min_genes=200)
    # sc.pp.filter_genes(adata1, min_cells=3)
    # sc.external.pp.scrublet(adata1)  # estimates doublets
    # adata = adata1[adata1.obs['predicted_doublet'] == False]  # do the actual filtering
    # adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    # sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # mito_filter = 25
    # n_counts_filter = 4300
    # adata = adata[adata.obs.n_genes_by_counts < n_counts_filter, :]
    # adata = adata[adata.obs.pct_counts_mt < mito_filter, :]
    # 
    # adata.raw = adata
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # adata = adata[:, adata.var.highly_variable]
    # sc.pp.scale(adata, max_value=10)
    # 
    # # adata_ref = sc.datasets.pbmc3k_processed()  # use the preprocessed data from the Scanpy tutorial as reference
    # adata_ref = sc.read_h5ad('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/pbmc_tutorial/data/adata_ref.h5ad')
    # var_names = adata_ref.var_names.intersection(adata.var_names)  # use genes which are present in both assays
    # adata_ref = adata_ref[:, var_names]
    # adata = adata[:, var_names]
    # sc.pp.pca(adata_ref)  # calculate PCA embedding
    # sc.pp.neighbors(adata_ref)  # calculate neighborhood graph
    # sc.tl.umap(adata_ref)  # calculate umap embedding
    # sc.tl.ingest(adata, adata_ref, obs='louvain')  # run label transfer
    # adata.obs.rename({'louvain': 'ingest_celltype_label'}, inplace=True, axis=1)
    # 
    # sc.tl.pca(adata, svd_solver='arpack')
    # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata, resolution=0.8, key_added='leiden_res_0.8')
    # tmp_df = adata.obs.groupby(['leiden_res_0.8', 'ingest_celltype_label']).size().unstack(fill_value=0)
    # tmp_df = (tmp_df / tmp_df.sum(0)).fillna(0)
    # leiden_to_annotation = tmp_df.idxmax(1).to_dict()
    # 
    # leiden_to_annotation['7'] = 'B cells 1'
    # leiden_to_annotation['11'] = 'B cells 2'
    # leiden_to_annotation = {cluster: leiden_to_annotation[cluster].replace(' ', '_') for cluster in
    #                         leiden_to_annotation.keys()}
    # 
    # adata.obs['celltype'] = [leiden_to_annotation[cluster_id] for cluster_id in adata.obs['leiden_res_0.8']]
    # del (leiden_to_annotation)
    # del (tmp_df)
    #
    # adata.write(os.path.join(work_dir, 'scRNA/adata.h5ad'), compression='gzip')

    adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/adata.h5ad'))
    # Annotation:
    anno_label = 'celltype'
    cell_data = adata.obs
    cell_data['sample_id'] = '10x_pbmc'
    cell_data[anno_label] = cell_data[anno_label].astype(str)
    # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

    # ATAC-seq
    fragments_dict = {'10x_pbmc': os.path.join(work_dir, 'data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz')}

    # Get chromosome sizes
    size_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/00.pbmc_tutorial/data/hg38.chrom.sizes'
    chromsizes = pd.read_csv(size_fn, sep='\t', header=None)
    print(chromsizes)
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
    print(chromsizes)

    # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
    # celltypes = ['CD14+_Monocytes', 'CD8_T_cells', 'B_cells_2', 'CD4_T_cells', 'FCGR3A+_Monocytes', 'B_cells_1']
    celltypes = ['NK_cells', 'Dendritic_cells']
    for ct in celltypes:
        one_sample(ct, chromsizes, cell_data)
        break
