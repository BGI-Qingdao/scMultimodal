#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 18 Dec 2023 16:58
# @Author: Yao LI
# @File: evo_fish/testrun_scenicplus.py
import os
import sys
import time
import scanpy as sc
import numpy as np
from pycisTopic.utils import (
    read_fragments_from_file,
    prepare_tag_cells
)
import pyranges as pr
import pandas as pd
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, export_pseudobulk_one_sample
import dill
import warnings

warnings.filterwarnings("ignore")


# Set stderr to null to avoid strange messages from ray


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


def one_sample(group, chromsizes, cell_data, output_dir, split_pattern='-', sample_id_col='sample_id',
               variable='celltype'):
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
    print('Starting...')
    # scRNA-seq
    scrna_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/scrna'
    fn1 = 'a_liver.h5ad'
    adata_rna = sc.read_h5ad(os.path.join(scrna_dir, fn1))
    adata_atac = sc.read_h5ad(
        '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/02.Derived_Lethenteron/merged.atac.h5ad')

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
    transfered_labels = pd.read_csv(
        '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/label_transfered_cell_data.csv', index_col=0)
    cell_data['celltype'] = transfered_labels['predicted.id']

    # ATAC-seq fragments
    fragment_fn = '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-zhangjin/07.project/5.liver_ATAC/0.qc/QC/lsqsm-youti-gan-2_1/output/lsqsm-youti-gan-2_1.fragments.tsv.gz'
    fragments_dict = {'a_liver': fragment_fn}

    # ------------------------------------------------
    # Get chromosome sizes (for 七腮鳗 here)
    # ------------------------------------------------
    size_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/data/Lethenteron_reissneri_sizes.genome'
    chromsizes = pd.read_csv(size_fn, sep='\t', header=None)
    chromsizes.columns = ['Chromosome', 'End']
    chromsizes['Start'] = [0] * chromsizes.shape[0]
    chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
    chromsizes = pr.PyRanges(chromsizes)

    # ------------------------------------------------
    # Create pseudo bulk data
    # ------------------------------------------------
    # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
    # celltypes = list(set(adata_atac.obs[anno_label]))
    # for ct in celltypes:
    #     one_sample(ct, chromsizes, cell_data, output_dir, split_pattern=None, variable=anno_label)

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
    import pickle

    bed_paths = pickle.load(
        open('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/atac/pseudobulk_bed_files/bed_paths.pkl',
             'rb'))
    bw_paths = pickle.load(
        open('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/atac/pseudobulk_bw_files/bw_paths.pkl',
             'rb'))

    # ------------------------------------------------
    from pycisTopic.pseudobulk_peak_calling import peak_calling

    macs_path = 'macs2'
    # ------------------------------------------------
    # Run peak calling
    # ------------------------------------------------
    # sys.stderr = open(os.devnull, "w")  # silence stderr
    # narrow_peaks_dict = peak_calling(macs_path,
    #                                  bed_paths,
    #                                  '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron',
    #                                  genome_size=1369200000,
    #                                  n_cpu=8,
    #                                  input_format='AUTO',
    #                                  shift=73,
    #                                  ext_size=146,
    #                                  keep_dup='all',
    #                                  q_value=0.05)
    # sys.stderr = sys.__stderr__  # unsilence stderr
    # import glob
    # files = glob.glob('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/narrow_peaks_dict/*.narrowPeak')
    # narrow_peaks_dict = {x: pr.read_bed(x) for x in files}
    # # Rename
    # for x in narrow_peaks_dict.keys():
    #     narrow_peaks_dict[x].columns = [
    #         'Chromosome',
    #         'Start',
    #         'End',
    #         'Name',
    #         'Score',
    #         'Strand',
    #         'FC_summit',
    #         '-log10_pval',
    #         '-log10_qval',
    #         'Summit']
    #
    # from pycisTopic.iterative_peak_calling import *
    # # Other param
    # peak_half_width = 250
    # # Get consensus peaks
    # sys.stderr = open(os.devnull, "w")  # silence stderr
    # consensus_peaks = get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes)
    # sys.stderr = sys.__stderr__  # unsilence stderr
    # consensus_peaks.to_bed(
    #     path='/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/consensus_regions.bed',
    #     keep=True,
    #     compression='infer',
    #     chain=False)

    # ------------------------------------------------
    # Creating a cisTopic object and topic modeling
    # ------------------------------------------------
    # generate a binary count matrix of ATAC-seq fragments over consensus peaks.
    # This matrix, along with metadata, will be stored in a cisTopic object and be used for topic modeling.
    path_to_regions = {
        'a_liver': '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/consensus_regions.bed'}

    # from pycisTopic.cistopic_class import *
    # cistopic_obj = create_cistopic_object_from_fragments(
    #     path_to_fragments=fragments_dict['a_liver'],
    #     path_to_regions=path_to_regions['a_liver'],
    #     n_cpu=1,
    #     project='a_liver')
    # pickle.dump(cistopic_obj,open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'wb'))
    # cistopic_obj = pickle.load(open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'rb'))

    # Add cell metadata
    # cell_data = pd.read_csv(os.path.join(work_dir, 'cell_data.csv'))
    # cistopic makes use of the sample_id to match the correct cell barcodes to the metadata, let's add the sample_id as a suffix to the cell barcodes
    # cell_data['barcode'] = cell_data['barcode'] + '___' + cell_data['sample_id']
    # print(cell_data['barcode'][0:5])
    # cell_data = cell_data.set_index('barcode')
    # print(cell_data)
    # print(cistopic_obj.cell_data.copy())
    # cistopic_obj.add_cell_data(cell_data)
    # # cistopic_obj.add_cell_data(cell_data[['a_liver']])
    # pickle.dump(cistopic_obj,open(os.path.join(work_dir, 'atac/cistopic_obj_models_500_iter_LDA.pkl'), 'wb'))
    # cistopic_obj = pickle.load(open(os.path.join(work_dir, 'atac/cistopic_obj_models_500_iter_LDA.pkl'), 'rb'))

    # Run topic modeling. The purpose of this is twofold:
    # To find sets of co-accessible regions (topics), this will be used downstream as candidate enhancers
    # (together with Differentially Accessible Regions (DARs)).
    # To impute dropouts.
    from pycisTopic.cistopic_class import *

    # models = run_cgs_models(cistopic_obj,
    #                         n_topics=[2, 4, 10, 16, 32, 48],
    #                         n_cpu=5,
    #                         n_iter=500,
    #                         random_state=555,
    #                         alpha=50,
    #                         alpha_by_topic=True,
    #                         eta=0.1,
    #                         eta_by_topic=False,
    #                         save_path=work_dir)
    # if not os.path.exists(os.path.join(work_dir, 'atac/models')):
    #     os.makedirs(os.path.join(work_dir, 'atac/models'))
    #
    # pickle.dump(models,
    #             open(os.path.join(work_dir, 'atac/models/a_liver_youti_2_1.pkl'), 'wb'))
    models = pickle.load(open(os.path.join(work_dir, 'atac/models/a_liver_youti_2_1.pkl'), 'rb'))

    # ------------------------------------------------
    # Analyze models
    # ------------------------------------------------
    from pycisTopic.lda_models import *
    # model = evaluate_models(models,
    #                         select_model=16,
    #                         return_model=True,
    #                         metrics=['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
    #                         plot_metrics=False)
    # cistopic_obj.add_LDA_model(model)
    # pickle.dump(cistopic_obj,
    #             open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'wb'))
    # cistopic_obj = pickle.load(open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'rb'))
    # cistopic_obj.add_cell_data(cell_data)
    # print(cistopic_obj.cell_data)
    # pickle.dump(cistopic_obj,
    #             open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'wb'))

    # Visualization
    from pycisTopic.clust_vis import *
    # run_umap(cistopic_obj, target='cell', scale=True)
    # plot_metadata(cistopic_obj, reduction_name='UMAP', variables=['celltype'], save=os.path.join(work_dir, 'figs/fig1'))
    # plot_topic(cistopic_obj, reduction_name='UMAP', num_columns=4, save=os.path.join(work_dir, 'figs/fig2'))

    # correct this kind of batch effects
    # Harmony
    # harmony(cistopic_obj, 'sample_id', random_state=555)
    # # UMAP
    # run_umap(cistopic_obj, reduction_name='harmony_UMAP',
    #          target='cell', harmony=True)
    # run_tsne(cistopic_obj, reduction_name='harmony_tSNE',
    #          target='cell', harmony=True)
    # Plot again
    # run_umap(cistopic_obj, target='cell', scale=True)
    # plot_metadata(cistopic_obj, reduction_name='UMAP', variables=['celltype'], save=os.path.join(work_dir, 'figs/fig1'))
    # plot_topic(cistopic_obj, reduction_name='UMAP', num_columns=4, save=os.path.join(work_dir, 'figs/fig2'))

    # cell_topic_heatmap(cistopic_obj,
    #                    variables=['celltype'],
    #                    scale=False,
    #                    legend_loc_x=1.05,
    #                    legend_loc_y=-1.2,
    #                    legend_dist_y=-1,
    #                    figsize=(10, 20),
    #                    save=work_dir + 'figs/heatmap_topic_contr.pdf')

    # ------------------------------------------------
    # Inferring candidate enhancer regions
    # ------------------------------------------------
    from pycisTopic.topic_binarization import *

    # region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu', save=work_dir + 'figs/otsu.pdf')
    # region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop=3000, save=work_dir + 'figs/top3k.pdf')
    # binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=True, num_columns=5,
    #                                        nbins=100, save=work_dir + 'figs/cell_topic.pdf')

    # from pycisTopic.diff_features import *
    # imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None,scale_factor=10 ** 6)
    # normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10 ** 4)
    # variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot=False)
    # markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='celltype', var_features=variable_regions,split_pattern='-')
    #
    # if not os.path.exists(os.path.join(work_dir, 'atac/candidate_enhancers')):
    #     os.makedirs(os.path.join(work_dir, 'atac/candidate_enhancers'))
    # pickle.dump(region_bin_topics_otsu,
    #             open(os.path.join(work_dir, 'atac/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
    # pickle.dump(region_bin_topics_top3k,
    #             open(os.path.join(work_dir, 'atac/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
    # pickle.dump(markers_dict, open(os.path.join(work_dir, 'atac/candidate_enhancers/markers_dict.pkl'), 'wb'))

    # ------------------------------------------------
    # Motif enrichment analysis using pycistarget
    # ------------------------------------------------
    # 1. Convert to dictionary of pyranges objects.
    region_bin_topics_otsu = pickle.load(
        open(os.path.join(work_dir, 'atac/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
    region_bin_topics_top3k = pickle.load(
        open(os.path.join(work_dir, 'atac/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
    markers_dict = pickle.load(
        open(os.path.join(work_dir, 'atac/candidate_enhancers/markers_dict.pkl'), 'rb'))

    from pycistarget.utils import region_names_to_coordinates
    # region_sets = {}
    # region_sets['topics_otsu'] = {}
    # region_sets['DARs_region_bin_topics_top3k'] = {}
    # region_sets['DARs_markers_dict'] = {}
    # for topic in region_bin_topics_otsu.keys():
    #     regions = region_bin_topics_otsu[topic].index[
    #         region_bin_topics_otsu[topic].index.str.startswith('Hic_chr')]  # only keep regions on known chromosomes
    #     region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
    # for DAR in region_bin_topics_top3k.keys():
    #     regions = region_bin_topics_top3k[DAR].index[
    #         region_bin_topics_top3k[DAR].index.str.startswith('Hic_chr')]  # only keep regions on known chromosomes
    #     region_sets['DARs_region_bin_topics_top3k'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    # for DAR in markers_dict.keys():  # TODO: why there are empty dataframes
    #     if not markers_dict[DAR].empty:
    #         regions = markers_dict[DAR].index[
    #             markers_dict[DAR].index.str.startswith('Hic_chr')]  # only keep regions on known chromosomes
    #         region_sets['DARs_markers_dict'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    # # for key in region_sets.keys():
    #     print(f'{key}: {region_sets[key].keys()}')

    # 2. CUSTOM DATABASE
    # db_fpath = "/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/00.custom"
    # # 2.1. make custom scoring database (feather)
    # # scores_db = os.path.join(db_fpath, 'Lethenteron.regions_vs_motifs.scores.feather') # peaks from ATAC data analysis by ? software
    # scores_db = os.path.join(db_fpath, 'Lethenteron.no_name.regions_vs_motifs.scores.feather')  # consensus peaks from previous step
    # # 2.2. make custom ranking database (feather)
    # # rankings_db = os.path.join(db_fpath, 'Lethenteron.regions_vs_motifs.rankings.feather')
    # rankings_db = os.path.join(db_fpath, 'Lethenteron.no_name.regions_vs_motifs.rankings.feather') # consensus peaks from previous step
    # # 2.3. make custom motif annotation metadata (tbl)
    motif_annot_fpath = "/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/00.custom"
    custom_annotation_fn = os.path.join(motif_annot_fpath, 'lethenteron_annotation.gtf')
    custom_annotation = pd.read_csv(custom_annotation_fn, sep='\t')
    # custom_annotation['Strand']
    # motif_annotation_fn = os.path.join(motif_annot_fpath, 'motifs-v10-nr.lethenteron-m0.001-o0.0.tbl')
    #
    # from scenicplus.wrappers.run_pycistarget import run_pycistarget
    # if not os.path.exists(os.path.join(work_dir, 'motifs')):
    #     os.makedirs(os.path.join(work_dir, 'motifs'))
    # run_pycistarget(
    #     region_sets=region_sets,
    #     species='custom',
    #     custom_annot=custom_annotation,
    #     save_path=os.path.join(work_dir, 'motifs'),
    #     ctx_db_path=rankings_db,
    #     dem_db_path=scores_db,
    #     path_to_motif_annotations=motif_annotation_fn,
    #     run_without_promoters=True,
    #     n_cpu=20,
    #     annotation_version="2024"
    # )

    # show the motifs found for topic 8 (specific to B-cells) using DEM
    # import dill
    # menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
    # menr['DEM_topics_otsu_All'].DEM_results('Topic8')


    # inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
    _stderr = sys.stderr
    null = open(os.devnull, 'wb')
    # atac_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/02.Derived_Lethenteron/merged.atac.h5ad'
    scrna_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/scrna/a_liver.h5ad'
    adata = sc.read_h5ad(scrna_fn)  # key_to_group_by, need to match with cistopic_obj
    cistopic_obj = pickle.load(open(os.path.join(work_dir, 'atac/cistopic_obj.pkl'), 'rb'))
    cistopic_obj.selected_model = models[2]
    print('cistopic_obj', type(cistopic_obj))
    menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

    # ------------------------------------------------
    #        Create the SCENIC+ object
    # ------------------------------------------------
    from scenicplus.scenicplus_class import create_SCENICPLUS_object
    # scplus_obj = create_SCENICPLUS_object(
    #     GEX_anndata=adata.raw.to_adata(),
    #     cisTopic_obj=cistopic_obj,
    #     menr=menr,
    #     bc_transform_func=lambda x: f'{x}___a_liver',
    #     nr_cells_per_metacells=5,
    #     meta_cell_split=' ',
    #     multi_ome_mode=False,
    #     key_to_group_by='celltype')
    #     # function to convert scATAC-seq barcodes to scRNA-seq ones
    # # scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
    # print(type(scplus_obj.X_EXP))
    # pickle.dump(scplus_obj, open(os.path.join(work_dir, 'scplus_obj.pkl'), 'wb'))
    # print(type(scplus_obj))
    # 2024-01-28 Done.

    biomart_host = "http://may2015.archive.ensembl.org"

    if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
        os.makedirs(os.path.join(work_dir, 'scenicplus'))

    scplus_obj = pickle.load(open(os.path.join(work_dir, 'scplus_obj.pkl'), 'rb'))

    import pyranges as pr
    custom_annot = custom_annotation.copy()
    custom_annot['Strand'] = custom_annot['Strand'].replace(1, '+')
    custom_annot['Strand'] = custom_annot['Strand'].replace(-1, '-')
    custom_annot['Transcription_Start_Site'] = custom_annot['Start']
    pr_annot = pr.PyRanges(custom_annot)
    # with custom annot a pandas DataFrame containing the extra Transcription_Start_Site column

    s_time = time.time()
    # for species different from human, mouse or fruit fly
    from scenicplus.enhancer_to_gene import get_search_space

    get_search_space(
        scplus_obj,
        pr_annot=pr_annot,  # see above
        pr_chromsizes=chromsizes,  # see above
        upstream=[1000, 150000], downstream=[1000, 150000],
        biomart_host='http://sep2019.archive.ensembl.org/')

    from scenicplus.wrappers.run_scenicplus import run_scenicplus
    # run the wrapper like usual
    try:
        sys.stderr = open(os.devnull, "w")  # silence stderr
        run_scenicplus(
            scplus_obj=scplus_obj,
            variable=['celltype'],
            species='custom',
            assembly='custom',
            tf_file='/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/00.custom/TFs_Lethenteron.txt',
            save_path=os.path.join(work_dir, 'scenicplus'),
            biomart_host=None,
            upstream=[1000, 150000],
            downstream=[1000, 150000],
            calculate_TF_eGRN_correlation=False,
            calculate_DEGs_DARs=False,
            export_to_loom_file=False,
            export_to_UCSC_file=False,
            n_cpu=24,
            _temp_dir=None)
        sys.stderr = sys.__stderr__  # unsilence stderr
    except Exception as e:
        # in case of failure, still save the object
        dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
        raise (e)

    e_time = time.time()
    print(f'finished in {e_time - s_time} seconds')
    print(scplus_obj.uns.keys())
