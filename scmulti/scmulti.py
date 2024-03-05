#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 18 Feb 2024 13:59
# @Author: Yao LI
# @File: evo_fish/scmulti.py

import os
import sys
import warnings

warnings.filterwarnings("ignore")

import glob
import dill
import pickle
import scanpy as sc
import numpy as np
import pyranges as pr
import pandas as pd
from typing import Optional, Union, List

from pycisTopic.utils import (
    read_fragments_from_file,
    prepare_tag_cells
)
from pycisTopic.pseudobulk_peak_calling import (
    export_pseudobulk,
    export_pseudobulk_one_sample
)

from pycisTopic.iterative_peak_calling import *
from pycisTopic.cistopic_class import *
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *

from custom_annotation import get_custom_annot


class ScMulti:
    def __init__(self,
                 refgenome_obj,
                 rna_fn,
                 fragment_fn,
                 scparams,
                 atac_fn=None,
                 atac_cell_data=None,
                 sample_id='sample_one',
                 fragments_dict=None,
                 sample_id_col='sample_id',
                 split_pattern=None,  # '-',
                 variable='celltype',
                 atac_anno_label='predicted.id',
                 gff_fn=None,
                 macs_path='macs2',
                 n_cpu=24,
                 c_value=None):  # params_obj=None
        # input exp data
        self.rna_data = sc.read_h5ad(rna_fn)
        self.atac_data = sc.read_h5ad(atac_fn)

        # species reference genome
        self.fasta = refgenome_obj.fasta_fn
        self.chromsizes_fn = refgenome_obj.chromsize_fn
        if self.chromsizes_fn is None:
            raise FileNotFoundError('Chromosome size file not found')
        self.chromsizes = None
        self.genome_size = None
        self.gff_fn = gff_fn
        self.c_value = c_value

        # ATAC meta data
        self.cell_data = atac_cell_data
        self.sample_id = sample_id
        self.sample_id_col = sample_id_col
        self.variable = variable
        self.atac_anno_label = atac_anno_label

        # ATAC data attributes
        self.fragment_fn = fragment_fn
        if fragments_dict is None:
            self.fragments_dict = {self.sample_id: self.fragment_fn}
        else:
            self.fragments_dict = fragments_dict
        self.split_pattern = split_pattern

        # Other attributes
        self.macs_path = macs_path
        self.n_cpu = n_cpu

        # saving results
        self.work_dir = scparams.work_dir
        self.output_dir = scparams.output_dir
        self.bed_paths = None
        self.bw_paths = None
        self.narrow_peaks_paths = None
        self.path_to_regions = None

        # inter results
        self.cistopic_obj = None
        self.scplus_obj = None
        self.models = None

        # enrich
        self.scores_db_fn = ''
        self.rankings_db_fn = ''
        self.custom_annotation_fn = ''
        self.motif_annotation_fn = ''
        self.tf_fn = ''

        # methods
        # self.create_atac_cell_data()
        # self.create_chromesize()
        # self.get_genome_size(c_value=c_value)

    def create_atac_cell_data(self, celltype_name_seps: Optional[Union[str, List[str]]]):
        """
        Extract needed annotation information (cell types, barcodes, sample ids...) from ATAC data.
        :param inplace:
        :param celltype_name_seps: delimiter in cell type names. could be a string '-' or a list of str ['-','.']
        :return:
        """
        cell_data = self.atac_data.obs.copy()
        cell_data[self.sample_id_col] = self.sample_id
        cell_data[self.atac_anno_label] = cell_data[self.atac_anno_label].astype(str)
        cell_data['celltype_label'] = cell_data[self.atac_anno_label]
        # Do not allow symbols other than _ exist in cell type name
        cell_data[self.atac_anno_label].str.replace("[^A-Za-z0-9]+", "_", regex=True, inplace=True)
        # if isinstance(celltype_name_seps, str):
        #     celltype_name_seps = list(celltype_name_seps)
        # if celltype_name_seps is None:
        #     celltype_name_seps = ['\s+', '-']
        # # special characters
        # pattern = r'[!@#$%^&*.]'
        # cell_data[self.atac_anno_label] = cell_data[self.atac_anno_label].str.replace(pattern, '_', regex=True)
        # for sep in celltype_name_seps:
        #     if sep not in pattern:
        #         cell_data[self.atac_anno_label].replace(sep, '_', regex=True, inplace=True)
        # # remove + in cell type name
        # cell_data[self.atac_anno_label].replace('[+]', '', regex=True, inplace=True)

        cell_data['barcode'] = list(cell_data.index)
        cell_data['barcode'] = cell_data['barcode'].astype(str)
        # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
        cell_data[self.variable] = cell_data[self.atac_anno_label]
        print(cell_data)
        self.cell_data = cell_data
        return cell_data

    def create_chromesize(self):
        """
        Generate chrome size object (PyRange)
        :return:
        """
        chromsizes = pd.read_csv(self.chromsizes_fn, sep='\t', header=None)
        chromsizes.columns = ['Chromosome', 'End']
        chromsizes['Start'] = [0] * chromsizes.shape[0]
        chromsizes = chromsizes.loc[:, ['Chromosome', 'Start', 'End']]
        chromsizes = pr.PyRanges(chromsizes)
        self.chromsizes = chromsizes
        print('Chrome Size:')
        print(self.chromsizes)
        return chromsizes

    def get_genome_size(self, c_value=None):
        """
        Total genome size for the species. (not effective genome size)
        :param c_value: C-value is the amount, in picograms (pg), of DNA contained within a haploid nucleus
        :return:
        """
        if c_value is None:
            chromsizes = pd.read_csv(self.chromsizes_fn, sep='\t', header=None)
            self.genome_size = int(sum(chromsizes[1]))
        else:
            self.genome_size = c_value_to_bp(c_value)
        return self.genome_size

    def one_sample(self, cell_type):
        """
        Generate pseudobulk peaks data for given cell type
        :param cell_type:
        :return:
        """
        fragments_df_dict = {}
        for sample_id in self.fragments_dict.keys():
            fragments_df = read_fragments_from_file(
                self.fragments_dict[sample_id],
                use_polars=True
            ).df
            fragments_df.Start = np.int32(fragments_df.Start)
            fragments_df.End = np.int32(fragments_df.End)
            if "Score" in fragments_df:
                fragments_df.Score = np.int32(fragments_df.Score)
            if "barcode" in self.cell_data:
                fragments_df = fragments_df.loc[
                    fragments_df["Name"].isin(self.cell_data["barcode"].tolist())
                ]
            else:
                fragments_df = fragments_df.loc[
                    fragments_df["Name"].isin(
                        prepare_tag_cells(self.cell_data.index.tolist(), self.split_pattern)
                    )
                ]
            fragments_df_dict[sample_id] = fragments_df
        if "barcode" in self.cell_data:
            self.cell_data = self.cell_data.loc[:, [self.variable, self.sample_id_col, "barcode"]]
        else:
            self.cell_data = self.cell_data.loc[:, [self.variable, self.sample_id_col]]
        self.cell_data[self.variable] = self.cell_data[self.variable].replace(
            " ", "", regex=True)
        self.cell_data[self.variable] = self.cell_data[self.variable].replace(
            "[^A-Za-z0-9]+", "_", regex=True)
        # groups = sorted(list(set(cell_data[self.variable])))
        if isinstance(self.chromsizes, pd.DataFrame):
            self.chromsizes = self.chromsizes.loc[:, ["Chromosome", "Start", "End"]]
            self.chromsizes = pr.PyRanges(self.chromsizes)

        # set up saving directories
        if not os.path.exists(os.path.join(self.output_dir, 'pseudobulk_bw_files')):
            os.makedirs(os.path.join(self.output_dir, 'pseudobulk_bw_files'))
        if not os.path.exists(os.path.join(self.output_dir, 'pseudobulk_bed_files')):
            os.makedirs(os.path.join(self.output_dir, 'pseudobulk_bed_files'))
        try:
            export_pseudobulk_one_sample(
                cell_data=self.cell_data,
                group=cell_type,
                fragments_df_dict=fragments_df_dict,
                chromsizes=self.chromsizes,
                bigwig_path=os.path.join(self.output_dir, 'pseudobulk_bw_files/'),
                bed_path=os.path.join(self.output_dir, 'pseudobulk_bed_files/'),
                sample_id_col=self.sample_id_col,
                normalize_bigwig=True,
                remove_duplicates=True,
                split_pattern=self.split_pattern
            )
        except ZeroDivisionError:
            print(f'{cell_type}: ZeroDivisionError')

    def generate_pseudo_peaks(self, celltypes: Optional[list] = None):
        """
        Generate pseudobulk peaks for all cell types in data
        :return:
        """
        # if not os.path.exists(os.path.join(self.output_dir, 'pseudobulk_bed_files')):
        #     os.makedirs(os.path.join(self.output_dir, 'pseudobulk_bed_files'))
        # if not os.path.exists(os.path.join(self.output_dir, 'pseudobulk_bw_files')):
        #     os.makedirs(os.path.join(self.output_dir, 'pseudobulk_bw_files'))
        # to create pseudo bulk data, needs cell type annotation from scRNA-seq data
        if celltypes is None:
            celltypes = list(set(self.cell_data[self.atac_anno_label]))
        for celltype in celltypes:
            self.one_sample(celltype)
        self.bed_paths = os.path.join(self.output_dir, 'pseudobulk_bed_files/bed_paths.pkl')
        self.bw_paths = os.path.join(self.output_dir, 'pseudobulk_bed_files/bw_paths.pkl')

    def custom_call_peaks(self,
                          convert_strand=False,
                          genome_size=None,
                          input_format='AUTO',
                          shift=73,
                          ext_size=146,
                          keep_dup='all',
                          q_value=0.05):
        """

        :param convert_strand:
        :param :
        :param input_format:
        :param shift:
        :param ext_size:
        :param keep_dup:
        :param q_value:
        :return:
        """
        import glob
        import subprocess
        if genome_size is None:
            genome_size = self.genome_size
        if self.genome_size is None:
            self.get_genome_size(c_value=self.c_value)
        paths = os.path.join(self.output_dir, 'pseudobulk_bed_files')
        out = os.path.join(self.output_dir, 'consensus_peak_calling/narrow_peaks')
        self.narrow_peaks_paths = out
        if not os.path.exists(out):
            os.makedirs(out)
        for bed_file in glob.glob(f'{paths}/*.bed.gz'):
            if convert_strand:
                add_strand(bed_file)
            cmd = f'{self.macs_path} callpeak --treatment {bed_file} --name {get_celltype_name(bed_file)} --outdir {out}  --format {input_format} --gsize {genome_size} --qvalue {q_value} --nomodel --shift {shift} --extsize {ext_size} --keep-dup {keep_dup} --call-summits --nolambda'
            print(cmd)
            try:
                subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    "command '{}' return with error (code {}): {}".format(
                        e.cmd, e.returncode, e.output
                    )
                )

    def custom_consensus_peaks(self, narrow_peaks_paths=None, peak_half_width=250):
        """

        :param narrow_peaks_paths:
        :param peak_half_width:
        :return:
        """
        if narrow_peaks_paths is None:
            narrow_peaks_paths = self.narrow_peaks_paths
        narrow_peaks_dict = read_narrowPeaks(narrow_peaks_paths)
        # Rename
        for x in narrow_peaks_dict.keys():
            narrow_peaks_dict[x].columns = [
                'Chromosome',
                'Start',
                'End',
                'Name',
                'Score',
                'Strand',
                'FC_summit',
                '-log10_pval',
                '-log10_qval',
                'Summit']

        # Get consensus peaks
        print('self.chromsizes')
        print(self.chromsizes)
        sys.stderr = open(os.devnull, "w")  # silence stderr
        consensus_peaks = get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=self.chromsizes)
        sys.stderr = sys.__stderr__  # unsilence stderr
        consensus_bed_path = os.path.join(self.output_dir, 'consensus_peak_calling/consensus_regions.bed')
        consensus_peaks.to_bed(
            path=consensus_bed_path,
            keep=True,
            compression='infer',
            chain=False)
        self.path_to_regions = {self.sample_id: consensus_bed_path}

    def create_cistopic_obj(self, n_cpu=None):
        """

        :param n_cpu:
        :return:
        """
        if n_cpu is None:
            n_cpu = self.n_cpu
        cistopic_obj = create_cistopic_object_from_fragments(
            path_to_fragments=self.fragments_dict[self.sample_id],
            path_to_regions=self.path_to_regions[self.sample_id],
            n_cpu=n_cpu,
            project=self.sample_id)
        # Add cell metadata
        cell_data = self.cell_data.copy()
        if cell_data.index.name is None:
            # cistopic makes use of the sample_id to match the correct cell barcodes to the metadata,
            # let's add the sample_id as a suffix to the cell barcodes
            cell_data['barcode'] = cell_data['barcode'] + '___' + cell_data['sample_id']
            cell_data = cell_data.set_index('barcode')
        cistopic_obj.add_cell_data(cell_data)
        self.cell_data = cell_data
        # cistopic_obj.add_cell_data(cell_data[['a_liver']])

        if not os.path.exists(os.path.join(self.output_dir, 'atac')):
            os.makedirs(os.path.join(self.output_dir, 'atac'))
        pickle.dump(cistopic_obj, open(os.path.join(self.output_dir, 'atac/cistopic_obj.pkl'), 'wb'))

        self.cistopic_obj = cistopic_obj
        return cistopic_obj

    def get_topics(self,
                   n_cpu=None,
                   n_topics=None,
                   n_iter=500,
                   random_state=555,
                   alpha=50,
                   alpha_by_topic=True,
                   eta=0.1,
                   eta_by_topic=False,
                   select_model=16,
                   return_model=True,
                   metrics=None,
                   plot_metrics=False
                   ):
        """

        :param n_cpu:
        :param n_topics:
        :param n_iter:
        :param random_state:
        :param alpha:
        :param alpha_by_topic:
        :param eta:
        :param eta_by_topic:
        :return:
        """
        if metrics is None:
            metrics = ['Arun_2010', 'Cao_Juan_2009', 'Minmo_2011', 'loglikelihood']
        if n_topics is None:
            n_topics = [2, 4, 10, 16, 32, 48]
        if n_cpu is None:
            n_cpu = self.n_cpu

        # Run topic modeling. The purpose of this is twofold:
        # To find sets of co-accessible regions (topics), this will be used downstream as candidate enhancers
        # (together with Differentially Accessible Regions (DARs)).
        # To impute dropouts.
        self.models = run_cgs_models(self.cistopic_obj,
                                     n_topics=n_topics,
                                     n_cpu=n_cpu,
                                     n_iter=n_iter,
                                     random_state=random_state,
                                     alpha=alpha,
                                     alpha_by_topic=alpha_by_topic,
                                     eta=eta,
                                     eta_by_topic=eta_by_topic,
                                     save_path=self.output_dir)
        if not os.path.exists(os.path.join(self.output_dir, 'atac/models')):
            os.makedirs(os.path.join(self.output_dir, 'atac/models'))
        pickle.dump(self.models,
                    open(os.path.join(self.output_dir, f'atac/models/{self.sample_id}_models_500_iter_LDA.pkl'), 'wb'))

        # ------------------------------------------------
        # Analyze models
        # ------------------------------------------------
        model = evaluate_models(self.models,
                                select_model=select_model,
                                return_model=return_model,
                                metrics=metrics,
                                plot_metrics=plot_metrics)
        print('evaluate models...')
        print(model)
        self.cistopic_obj.add_LDA_model(model)
        self.cistopic_obj.add_cell_data(self.cell_data)
        pickle.dump(self.cistopic_obj, open(os.path.join(self.output_dir, 'atac/cistopic_obj_models.pkl'), 'wb'))

    def auto_best_topic(self,
                        metrics=None,
                        min_topics_coh: Optional[int] = 5):
        """
        Automatically select the best model based on model quality metrics
        :param metrics: model quality metrics
        :param min_topics_coh: Minimum number of topics on a topic to use its coherence for model selection.
        :return:
        """
        if metrics is None:
            metrics = ["Minmo_2011", "loglikelihood", "Cao_Juan_2009", "Arun_2010"]
        metrics_dict = {}
        all_topics = sorted([self.models[x].n_topic for x in range(0, len(self.models))])
        in_index = [i for i in range(len(all_topics)) if all_topics[i] >= min_topics_coh]
        for metric in metrics:
            try:
                metric_var = [
                    self.models[index].metrics.loc["Metric", metric]
                    for index in range(0, len(all_topics))]
                metric_negative = [-x for x in metric_var]
                metric_rescale = (metric_negative - min(metric_negative)) / (
                        max(metric_negative) - min(metric_negative)
                )
                metrics_dict[metric] = np.array(
                    subset_list(metric_rescale, in_index)
                )
            except KeyError:
                print(f'{metric} not found')

        combined_metric = sum(metrics_dict.values())
        best_model_index = all_topics[combined_metric.tolist().index(max(combined_metric))]
        best_model = self.models[all_topics.index(best_model_index)]
        return best_model

    def visualize_topics(self):
        """

        :return:
        """
        #
        if not os.path.exists(os.path.join(self.output_dir, 'figs')):
            os.makedirs(os.path.join(self.output_dir, 'figs'))
        out_path = os.path.join(self.output_dir, 'figs')
        # Visualization
        run_umap(self.cistopic_obj, target='cell', scale=True)
        plot_metadata(self.cistopic_obj, reduction_name='UMAP', variables=[self.variable],
                      save=os.path.join(out_path, 'metadata1'))
        plot_topic(self.cistopic_obj, reduction_name='UMAP', num_columns=4, save=os.path.join(out_path, 'topic_plot1'))

        # correct this kind of batch effects
        # Harmony
        harmony(self.cistopic_obj, 'sample_id', random_state=555)
        # # UMAP
        run_umap(self.cistopic_obj, reduction_name='harmony_UMAP',
                 target='cell', harmony=True)
        run_tsne(self.cistopic_obj, reduction_name='harmony_tSNE',
                 target='cell', harmony=True)
        # Plot again
        run_umap(self.cistopic_obj, target='cell', scale=True)
        plot_metadata(self.cistopic_obj, reduction_name='UMAP', variables=[self.variable],
                      save=os.path.join(out_path, 'metadata'))
        plot_topic(self.cistopic_obj, reduction_name='UMAP', num_columns=4, save=os.path.join(out_path, 'topic_plot'))

        cell_topic_heatmap(self.cistopic_obj,
                           variables=[self.variable],
                           scale=False,
                           legend_loc_x=1.05,
                           legend_loc_y=-1.2,
                           legend_dist_y=-1,
                           figsize=(10, 20),
                           save=out_path + '/heatmap_topic_contr.pdf')

    def candidate_enhancer_regions(self, split_pattern='-'):
        """

        :param split_pattern:
        :return:
        """
        if not os.path.exists(os.path.join(self.output_dir, 'figs')):
            os.makedirs(os.path.join(self.output_dir, 'figs'))
        out_path = os.path.join(self.output_dir, 'figs')

        region_bin_topics_otsu = binarize_topics(self.cistopic_obj, method='otsu', save=out_path + '/otsu.pdf')
        region_bin_topics_top3k = binarize_topics(self.cistopic_obj, method='ntop', ntop=3000,
                                                  save=out_path + '/top3k.pdf')
        binarized_cell_topic = binarize_topics(self.cistopic_obj, target='cell', method='li', plot=True, num_columns=5,
                                               nbins=100, save=out_path + '/cell_topic.pdf')

        imputed_acc_obj = impute_accessibility(self.cistopic_obj, selected_cells=None, selected_regions=None,
                                               scale_factor=10 ** 6)
        normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10 ** 4)
        variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot=False)
        markers_dict = find_diff_features(self.cistopic_obj, imputed_acc_obj, variable=self.variable,
                                          var_features=variable_regions, split_pattern=split_pattern)

        if not os.path.exists(os.path.join(self.output_dir, 'atac/candidate_enhancers')):
            os.makedirs(os.path.join(self.output_dir, 'atac/candidate_enhancers'))
        pickle.dump(region_bin_topics_otsu,
                    open(os.path.join(self.output_dir, 'atac/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
        pickle.dump(region_bin_topics_top3k,
                    open(os.path.join(self.output_dir, 'atac/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
        pickle.dump(markers_dict,
                    open(os.path.join(self.output_dir, 'atac/candidate_enhancers/markers_dict.pkl'), 'wb'))

    def create_region_sets(self, chrom_header: tuple):
        """

        :param chrom_header:
        :return:
        """
        # 1. Convert to dictionary of pyranges objects.
        region_bin_topics_otsu = pickle.load(
            open(os.path.join(self.output_dir, 'atac/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
        region_bin_topics_top3k = pickle.load(
            open(os.path.join(self.output_dir, 'atac/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
        markers_dict = pickle.load(
            open(os.path.join(self.output_dir, 'atac/candidate_enhancers/markers_dict.pkl'), 'rb'))

        from pycistarget.utils import region_names_to_coordinates
        region_sets = {}
        region_sets['topics_otsu'] = {}
        region_sets['DARs_region_bin_topics_top3k'] = {}
        region_sets['DARs_markers_dict'] = {}
        for topic in region_bin_topics_otsu.keys():
            regions = region_bin_topics_otsu[topic].index[
                region_bin_topics_otsu[topic].index.str.startswith(
                    chrom_header)]  # only keep regions on known chromosomes
            region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
        for DAR in region_bin_topics_top3k.keys():
            regions = region_bin_topics_top3k[DAR].index[
                region_bin_topics_top3k[DAR].index.str.startswith(
                    chrom_header)]  # only keep regions on known chromosomes
            region_sets['DARs_region_bin_topics_top3k'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
        for DAR in markers_dict.keys():  # TODO: why there are empty dataframes
            if not markers_dict[DAR].empty:
                regions = markers_dict[DAR].index[
                    markers_dict[DAR].index.str.startswith(chrom_header)]  # only keep regions on known chromosomes
                region_sets['DARs_markers_dict'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
        for key in region_sets.keys():
            print(f'{key}: {region_sets[key].keys()}')

        return region_sets

    def make_custom_databases(self, db_path, prefix):
        """

        :param db_path:
        :param prefix:
        :return:
        """
        self.scores_db_fn = os.path.join(db_path, f'{prefix}.regions_vs_motifs.scores.feather')
        self.rankings_db_fn = os.path.join(db_path, f'{prefix}.regions_vs_motifs.rankings.feather')
        self.custom_annotation_fn = os.path.join(db_path, f'{prefix}_custom_annotation.gtf')
        self.motif_annotation_fn = os.path.join(db_path, f'motifs-v10-nr.{prefix}-m0.001-o0.0.tbl')
        self.tf_fn = os.path.join(db_path, f'{prefix}_TFs.txt')

    def get_pr_annot(self):
        """
        with custom annot a pandas DataFrame containing the extra Transcription_Start_Site column
        :return:
        """
        custom_annotation = get_custom_annot(self.gff_fn, self.custom_annotation_fn)
        custom_annot = custom_annotation.copy()
        custom_annot['Strand'] = custom_annot['Strand'].replace('1', '+')  # TODO: is there a more robust way
        custom_annot['Strand'] = custom_annot['Strand'].replace('-1', '-')
        custom_annot['Transcription_Start_Site'] = custom_annot['Start']
        pr_annot = pr.PyRanges(custom_annot)
        return pr_annot

    def enrichment(self, chrom_header: tuple = ('chr'), species='custom', n_cpu=None):
        """

        :param chrom_header:
        :param species:
        :return:
        """
        if n_cpu is None:
            n_cpu = self.n_cpu
        if self.gff_fn is None:
            raise ValueError("GTF/GFF file not provided")
        custom_annotation = get_custom_annot(self.gff_fn, self.custom_annotation_fn)

        from scenicplus.wrappers.run_pycistarget import run_pycistarget
        if not os.path.exists(os.path.join(self.work_dir, 'motifs')):
            os.makedirs(os.path.join(self.work_dir, 'motifs'))
        region_sets = self.create_region_sets(chrom_header)
        run_pycistarget(
            region_sets=region_sets,
            species=species,
            custom_annot=custom_annotation,
            save_path=os.path.join(self.output_dir, 'motifs'),
            ctx_db_path=self.rankings_db_fn,
            dem_db_path=self.scores_db_fn,
            path_to_motif_annotations=self.motif_annotation_fn,
            run_without_promoters=True,
            n_cpu=n_cpu,
            annotation_version="2024"
        )

    def network(self,
                meta_cell_split=' ',
                multi_ome_mode=False,
                tf_file=None,
                upstream=[1000, 150000],
                downstream=[1000, 150000],
                use_gene_boundaries=True,
                n_cpu=None):
        """

        :param meta_cell_split:
        :param multi_ome_mode:
        :param tf_file:
        :param upstream:
        :param downstream:
        :param use_gene_boundaries:
        :return:
        """
        if n_cpu is None:
            n_cpu = self.n_cpu
        if tf_file is None:
            tf_file = self.tf_fn
        _stderr = sys.stderr
        null = open(os.devnull, 'wb')
        # ensure selected model exists
        if not self.cistopic_obj.selected_model:  # TODO: 为什么后面的步骤会改变这个数值
            self.cistopic_obj.selected_model = self.auto_best_topic()

        menr = dill.load(open(os.path.join(self.output_dir, 'motifs/menr.pkl'), 'rb'))

        # ------------------------------------------------
        #        Create the SCENIC+ object
        # ------------------------------------------------
        from scenicplus.scenicplus_class import create_SCENICPLUS_object
        if self.rna_data.raw is None:  # TODO: why use raw?
            gex_data = self.rna_data
        else:
            gex_data = self.rna_data.raw.to_adata()
        self.scplus_obj = create_SCENICPLUS_object(
            GEX_anndata=gex_data,
            cisTopic_obj=self.cistopic_obj,
            menr=menr,
            bc_transform_func=lambda x: f'{x}___{self.sample_id}',
            nr_cells_per_metacells=5,
            meta_cell_split=meta_cell_split,
            multi_ome_mode=multi_ome_mode,
            key_to_group_by=self.variable)
        # function to convert scATAC-seq barcodes to scRNA-seq ones
        # scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())

        from scenicplus.enhancer_to_gene import get_search_space
        pr_annot = self.get_pr_annot()

        get_search_space(
            self.scplus_obj,
            pr_annot=pr_annot,  # see above
            pr_chromsizes=self.chromsizes,  # see above
            upstream=upstream, downstream=downstream,
            use_gene_boundaries=use_gene_boundaries,
            biomart_host='http://sep2019.archive.ensembl.org/')

        from scenicplus.wrappers.run_scenicplus import run_scenicplus
        # run the wrapper like usual
        try:
            # sys.stderr = open(os.devnull, "w")  # silence stderr
            run_scenicplus(
                scplus_obj=self.scplus_obj,
                variable=[self.variable],
                species='custom',
                assembly='custom',
                tf_file=tf_file,
                save_path=os.path.join(self.output_dir, 'scenicplus'),
                biomart_host=None,
                upstream=upstream,
                downstream=downstream,
                calculate_TF_eGRN_correlation=False,
                calculate_DEGs_DARs=False,
                export_to_loom_file=False,
                export_to_UCSC_file=False,
                n_cpu=n_cpu,
                _temp_dir=None)
            # sys.stderr = sys.__stderr__  # unsilence stderr
        except Exception as e:
            # in case of failure, still save the object
            raise (e)
        finally:
            dill.dump(self.scplus_obj, open(os.path.join(self.output_dir, 'scenicplus/scplus_obj.pkl'), 'wb'),
                      protocol=-1)


def get_celltype_name(fn):
    """ Only suitable for bed files in this package"""
    bn = os.path.basename(fn)
    ct = bn.replace('.bed.gz', '')
    return ct


def read_narrowPeaks(narrowPeaks_path):
    """
    Create narrow peaks dictionary
    :param narrowPeaks_path:
    :return:
    """
    files = glob.glob(f'{narrowPeaks_path}/*.narrowPeak')
    if files:
        narrow_peaks_dict = {x: pr.read_bed(x) for x in files}
    else:
        print('Did not find narrow peak files')
        narrow_peaks_dict = {}
    return narrow_peaks_dict


def add_strand(fn):  # TODO:
    """

    :param fn:
    :return:
    """
    df = pd.read_csv(fn, compression='gzip', sep='\t', header=None)
    df[5] = '+'
    df.to_csv(fn, index=False, compression='gzip', header=None, sep='\t')


def c_value_to_bp(c_value: float) -> int:
    """
    Convert C-value (in picograms/pg) into base pairs (bp)
    Animal genome size data are available at https://www.genomesize.com/index.php
    :param c_value:
    :return:
    """
    genome_size = c_value * 0.978 * 1e9
    return int(genome_size)


def subset_list(target_list, index_list):
    X = list(map(target_list.__getitem__, index_list))
    return X
