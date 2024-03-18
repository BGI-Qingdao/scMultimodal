#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 29 Jan 2024 17:03
# @Author: Yao LI
# @File: evo_fish/plotting.py
import os
import sys
import dill
import pickle
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# from scmulti import ScMulti
from scenicplus.RSS import *
from scenicplus.plotting.correlation_plot import *
from scmulti import match_key
from scenicplus.cistromes import *

COLORS = ["#065143", "#70B77E", "#E0A890", "#F56476", "#CE1483",
          '#d60000', '#e2afaf', '#018700', '#a17569', '#e6a500', '#004b00',
          '#6b004f', '#573b00', '#005659', '#5e7b87', '#0000dd', '#00acc6',
          '#bcb6ff', '#bf03b8', '#645472', '#790000', '#0774d8', '#729a7c',
          '#8287ff', '#ff7ed1', '#8e7b01', '#9e4b00', '#8eba00', '#a57bb8',
          '#5901a3', '#8c3bff', '#a03a52', '#a1c8c8', '#f2007b', '#ff7752',
          '#bac389', '#15e18c', '#60383b', '#546744', '#380000', '#e252ff',
          ]


class DownScMulti:
    def __init__(self, scplus_obj, work_dir, variable='celltype', celltypes=None):
        if celltypes is None:
            celltypes = []
        self.scplus_obj = scplus_obj
        self.work_dir = work_dir
        self.saving_dir = os.path.join(self.work_dir, 'downstream')
        if not os.path.exists(self.saving_dir):
            os.makedirs(self.saving_dir)
        self.variable = variable
        self.celltypes = celltypes

        self.selected_eRegulons = None

    # ---------------------------------------------------------------
    # Calculate correlation coefficient
    # ---------------------------------------------------------------
    def calculate(self):
        from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
        apply_std_filtering_to_eRegulons(self.scplus_obj)
        from scenicplus.eregulon_enrichment import score_eRegulons
        region_ranking = dill.load(open(os.path.join(self.work_dir, 'scenicplus/region_ranking.pkl'),
                                        'rb'))  # load ranking calculated using the wrapper function
        gene_ranking = dill.load(open(os.path.join(self.work_dir, 'scenicplus/gene_ranking.pkl'),
                                      'rb'))  # load ranking calculated using the wrapper function
        score_eRegulons(self.scplus_obj,
                        ranking=region_ranking,
                        eRegulon_signatures_key='eRegulon_signatures_filtered',
                        key_added='eRegulon_AUC_filtered',
                        enrichment_type='region',
                        auc_threshold=0.05,
                        normalize=False,
                        n_cpu=5)
        score_eRegulons(self.scplus_obj,
                        gene_ranking,
                        eRegulon_signatures_key='eRegulon_signatures_filtered',
                        key_added='eRegulon_AUC_filtered',
                        enrichment_type='gene',
                        auc_threshold=0.05,
                        normalize=False,
                        n_cpu=5)

        from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

        generate_pseudobulks(
            scplus_obj=self.scplus_obj,
            variable='celltype',
            auc_key='eRegulon_AUC_filtered',
            signature_key='Gene_based',
            nr_cells=1)
        generate_pseudobulks(
            scplus_obj=self.scplus_obj,
            variable='celltype',
            auc_key='eRegulon_AUC_filtered',
            signature_key='Region_based',
            nr_cells=1)

        TF_cistrome_correlation(
            self.scplus_obj,
            use_pseudobulk=False,  # 2024-03-18: use_pseudobulk=True,
            variable='celltype',
            auc_key='eRegulon_AUC_filtered',
            signature_key='Gene_based',
            out_key='filtered_gene_based')
        TF_cistrome_correlation(
            self.scplus_obj,
            use_pseudobulk=False,  # 2024-03-18: use_pseudobulk=True,
            variable='celltype',
            auc_key='eRegulon_AUC_filtered',
            signature_key='Region_based',
            out_key='filtered_region_based')

        n_targets = [int(x.split('(')[1].replace('r)', '')) for x in
                     self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
        rho = self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
        adj_pval = self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

        # thresholds = {
        #     'rho': [-0.75, 0.70],
        #     'n_targets': 0
        # }
        thresholds = {
            'rho': [-0.10, 0.10],
            'n_targets': 0
        }
        return n_targets, rho, adj_pval, thresholds

    def calculate2(self):
        from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
        apply_std_filtering_to_eRegulons(self.scplus_obj)
        from scenicplus.eregulon_enrichment import score_eRegulons
        region_ranking = dill.load(open(os.path.join(self.work_dir, 'scenicplus/region_ranking.pkl'),
                                        'rb'))  # load ranking calculated using the wrapper function
        gene_ranking = dill.load(open(os.path.join(self.work_dir, 'scenicplus/gene_ranking.pkl'),
                                      'rb'))  # load ranking calculated using the wrapper function
        score_eRegulons(self.scplus_obj,
                        ranking=region_ranking,
                        eRegulon_signatures_key='eRegulon_signatures_filtered',
                        key_added='eRegulon_AUC_filtered',
                        enrichment_type='region',
                        auc_threshold=0.05,
                        normalize=False,
                        n_cpu=5)
        score_eRegulons(self.scplus_obj,
                        gene_ranking,
                        eRegulon_signatures_key='eRegulon_signatures_filtered',
                        key_added='eRegulon_AUC_filtered',
                        enrichment_type='gene',
                        auc_threshold=0.05,
                        normalize=False,
                        n_cpu=5)

        eregulon_correlation(self.scplus_obj,
                             auc_key='eRegulon_AUC',
                             signature_key1='Gene_based',
                             signature_key2='Region_based',
                             nSignif=3,
                             out_key='Unfiltered')

        n_targets = [int(x.split('(')[1].replace('r)', '')) for x in
                     self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
        rho = self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
        adj_pval = self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

        thresholds = {
            'rho': [-0.30, 0.30],
            'n_targets': 0
        }
        return n_targets, rho, adj_pval, thresholds

    # ---------------------------------------------------------------
    # 1. Plot correlation coefficient
    # ---------------------------------------------------------------
    def plot_correlation(self, rho, n_targets, adj_pval, thresholds):
        fig, ax = plt.subplots(figsize=(10, 5))
        sc = ax.scatter(rho, n_targets, c=-np.log10(adj_pval), s=5)
        ax.set_xlabel('Correlation coefficient')
        ax.set_ylabel('nr. target regions')
        # ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
        ax.vlines(x=thresholds['rho'], ymin=0, ymax=max(n_targets), color='black', ls='dashed', lw=1)
        ax.text(x=thresholds['rho'][0], y=max(n_targets), s=str(thresholds['rho'][0]))
        ax.text(x=thresholds['rho'][1], y=max(n_targets), s=str(thresholds['rho'][1]))
        sns.despine(ax=ax)
        fig.colorbar(sc, label='-log10(adjusted_pvalue)', ax=ax)
        plt.savefig(os.path.join(self.saving_dir, 'cor_coe.png'))

        selected_cistromes = self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
            np.logical_or(
                self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
                self.scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
            )]['Cistrome'].to_list()
        selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
        self.selected_eRegulons = selected_eRegulons
        print('selected_eRegulons')
        print(selected_eRegulons)
        print(type(selected_eRegulons))
        selected_eRegulons_gene_sig = [
            x for x in self.scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
            if x.split('_(')[0] in selected_eRegulons]
        selected_eRegulons_region_sig = [
            x for x in self.scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
            if x.split('_(')[0] in selected_eRegulons]
        # save the results in the scenicplus object
        self.scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig,
                                                    'Region_based': selected_eRegulons_region_sig}
        print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')

    # ---------------------------------------------------------------
    # 2. heatmap
    # ---------------------------------------------------------------
    def plot_heatmap(self):
        from scenicplus.plotting.dotplot import heatmap_dotplot
        heatmap_dotplot(
            scplus_obj=self.scplus_obj,
            size_matrix=self.scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'],
            # specify what to plot as dot sizes, target region enrichment in this case
            color_matrix=self.scplus_obj.to_df('EXP'),  # specify  what to plot as colors, TF expression in this case
            scale_size_matrix=True,
            scale_color_matrix=True,
            group_variable='celltype',
            subset_eRegulons=self.scplus_obj.uns['selected_eRegulon']['Gene_based'],
            figsize=(50, 50),
            orientation='vertical',
            save=os.path.join(self.saving_dir, 'heatmap.pdf'),
            index_order=self.celltypes)

    # ---------------------------------------------------------------
    # 3. rss plot
    # ---------------------------------------------------------------
    def plot_erss(self, **kwargs):
        regulon_specificity_scores(
            self.scplus_obj,
            variable=self.variable,
            auc_key='eRegulon_AUC_filtered',
            signature_keys=['Region_based'],
            selected_regulons=[x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
            out_key_suffix='_filtered')

        plot_rss(
            self.scplus_obj,
            f'{self.variable}_filtered',
            save=os.path.join(self.saving_dir, 'rss.pdf'),
            **kwargs)

    # ---------------------------------------------------------------
    # 4. Jaccard Heatmap
    # ---------------------------------------------------------------
    def plot_jaccard_heatmap(self, **kwargs):
        flat_list = lambda t: [item for sublist in t for item in sublist]
        selected_markers = list(set(flat_list(
            [self.scplus_obj.uns['RSS'][f'{self.variable}_filtered'].loc[celltype].sort_values(ascending=False).head(
                30).index.to_list()
             for celltype in self.scplus_obj.uns['RSS'][f'{self.variable}_filtered'].index])))
        jaccard_heatmap(
            self.scplus_obj,
            method='intersect',
            gene_or_region_based='Region_based',
            use_plotly=False,
            selected_regulons=selected_markers,
            signature_key='eRegulon_signatures_filtered',
            figsize=(10, 10),
            return_data=False,
            vmax=0.5,
            cmap='plasma',
            save=os.path.join(self.saving_dir, 'jaccard_heatmap.pdf'),
            **kwargs)

    # ---------------------------------------------------------------
    # 5. Plot dimensionality reductions
    # ---------------------------------------------------------------
    def plot_dim_reduction(self, color_dict=None):
        from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
        run_eRegulons_umap(
            scplus_obj=self.scplus_obj,
            auc_key='eRegulon_AUC_filtered',
            reduction_name='eRegulons_UMAP',  # overwrite previously calculated UMAP
        )
        run_eRegulons_tsne(
            scplus_obj=self.scplus_obj,
            auc_key='eRegulon_AUC_filtered',
            reduction_name='eRegulons_tSNE',  # overwrite previously calculated tSNE
        )
        from scenicplus.dimensionality_reduction import plot_metadata_given_ax
        # specify color_dictionary
        if color_dict is None:
            color_dict = dict(zip(self.celltypes, COLORS))
        fig, axs = plt.subplots(ncols=2, figsize=(16, 8))
        plot_metadata_given_ax(
            scplus_obj=self.scplus_obj,
            ax=axs[0],
            reduction_name='eRegulons_UMAP',
            variable=self.variable,
            color_dictionary={self.variable: color_dict}
        )
        plot_metadata_given_ax(
            scplus_obj=self.scplus_obj,
            ax=axs[1],
            reduction_name='eRegulons_tSNE',
            variable=self.variable,
            color_dictionary={self.variable: color_dict}
        )
        fig.tight_layout()
        sns.despine(ax=axs[0])  # remove top and right edge of axis border
        sns.despine(ax=axs[1])  # remove top and right edge of axis border
        plt.savefig(os.path.join(self.saving_dir, 'dim_reduc.pdf'), format='pdf')

    # ---------------------------------------------------------------
    # 7. Plot activity / expression of an eRegulon
    # ---------------------------------------------------------------
    def plot_eregulon(self, selected_regulons=None, **kwargs):
        from scenicplus.dimensionality_reduction import plot_eRegulon
        if selected_regulons is None:
            selected_regulons = self.selected_eRegulons
        plot_eRegulon(
            scplus_obj=self.scplus_obj,
            reduction_name='eRegulons_tSNE',
            selected_regulons=selected_regulons,
            # ['MSTRG.2388_extended_+', 'MSTRG.8533_extended_+', 'nbisL1-mrna-14072_extended_+','MSTRG.18385_extended_+', 'MSTRG.14903_extended_+'],
            scale=True,
            auc_key='eRegulon_AUC_filtered',
            save=os.path.join(self.saving_dir, 'eRegulon_activity.pdf'), **kwargs)

    # ---------------------------------------------------------------
    # 8. Plot Network
    # ---------------------------------------------------------------
    def plot_network(self, target_celltype, target_regulons: list, n_top_features=1000, figsize=(20, 20)):
        # plot hvg and hvr
        from pycisTopic.diff_features import find_highly_variable_features
        hvr = find_highly_variable_features(
            self.scplus_obj.to_df('ACC').loc[list(set(self.scplus_obj.uns['eRegulon_metadata_filtered']['Region']))],
            n_top_features=n_top_features,
            plot=False, save=os.path.join(self.saving_dir, 'hvr.pdf'))
        hvg = find_highly_variable_features(
            self.scplus_obj.to_df('EXP')[list(set(self.scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T,
            n_top_features=n_top_features,
            plot=False, save=os.path.join(self.saving_dir, 'hvg.pdf'))


        from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape
        nx_tables = create_nx_tables(
            scplus_obj=self.scplus_obj,
            eRegulon_metadata_key='eRegulon_metadata_filtered',
            subset_regions=hvr,  # subset_eRegulons=None, 使用target_regulons
            subset_genes=hvg,
            add_differential_gene_expression=True,
            add_differential_region_accessibility=True,
            differential_variable=[self.variable])
        # print(nx_tables)
        # pickle.dump(nx_tables, open(os.path.join(self.saving_dir, 'nx_tables.pkl'), 'wb'))

        # assign regulons color
        # network_color = ['Orange', 'Purple', 'Red', 'Yellow', 'Blue', 'Green', 'Pink', 'Gray', 'Brown']
        # if len(network_color) < len(target_regulons):
        #     print('Warning! not enough color for each regulon')
        #     target_regulons = target_regulons[:len(network_color)]
        #     print(target_regulons)
        # category_color = dict(zip(target_regulons, network_color))
        category_color = None
        # make the big dictionary
        G, pos, edge_tables, node_tables = create_nx_graph(nx_tables,
                                                           use_edge_tables=['TF2R', 'R2G'],
                                                           color_edge_by={'TF2R': {'variable': 'TF',
                                                                                   'category_color': category_color},
                                                                          'R2G': {'variable': 'R2G_rho',
                                                                                  'continuous_color': 'viridis',
                                                                                  'v_min': -1,
                                                                                  'v_max': 1}},
                                                           transparency_edge_by={
                                                               'R2G': {'variable': 'R2G_importance', 'min_alpha': 0.1,
                                                                       'v_min': 0}},
                                                           width_edge_by={
                                                               'R2G': {'variable': 'R2G_importance', 'max_size': 1.5,
                                                                       'min_size': 1}},
                                                           color_node_by={'TF': {'variable': 'TF',
                                                                                 'category_color': category_color},
                                                                          'Gene': {
                                                                              'variable': f'celltype_Log2FC_{target_celltype}',
                                                                              'continuous_color': 'bwr'},
                                                                          'Region': {
                                                                              'variable': f'celltype_Log2FC_{target_celltype}',
                                                                              'continuous_color': 'viridis'}},
                                                           transparency_node_by={
                                                               'Region': {
                                                                   'variable': f'celltype_Log2FC_{target_celltype}',
                                                                   'min_alpha': 0.1},
                                                               'Gene': {
                                                                   'variable': f'celltype_Log2FC_{target_celltype}',
                                                                   'min_alpha': 0.1}},
                                                           size_node_by={
                                                               'TF': {'variable': 'fixed_size', 'fixed_size': 30},
                                                               'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
                                                               'Region': {'variable': 'fixed_size',
                                                                          'fixed_size': 10}},
                                                           shape_node_by={
                                                               'TF': {'variable': 'fixed_shape',
                                                                      'fixed_shape': 'ellipse'},
                                                               'Gene': {'variable': 'fixed_shape',
                                                                        'fixed_shape': 'ellipse'},
                                                               'Region': {'variable': 'fixed_shape',
                                                                          'fixed_shape': 'diamond'}},
                                                           label_size_by={
                                                               'TF': {'variable': 'fixed_label_size',
                                                                      'fixed_label_size': 20.0},
                                                               'Gene': {'variable': 'fixed_label_size',
                                                                        'fixed_label_size': 10.0},
                                                               'Region': {'variable': 'fixed_label_size',
                                                                          'fixed_label_size': 0.0}},
                                                           layout='kamada_kawai_layout',
                                                           scale_position_by=250)

        plt.figure(figsize=figsize)
        plot_networkxx(G, pos, save=os.path.join(self.saving_dir, f'{target_celltype}_network.pdf'))
        export_to_cytoscape(G, pos, out_file=os.path.join(self.saving_dir, f'{target_celltype}_network.cys'))


def select_regulons(scplus_obj, thresholds=None):
    if thresholds is None:
        thresholds = {'rho': [-0.10, 0.10], 'n_targets': 0}
    selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
        np.logical_or(
            scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
            scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
        )]['Cistrome'].to_list()
    selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
    return selected_eRegulons


def plot_networkxx(G, pos, save='network.pdf'):
    """
    A function to plot networks with networkx

    Parameters
    ---------
    G: Graph
        A networkx graph
    pos: Dict
        Position values
    """
    import networkx as nx
    nx.draw_networkx_nodes(G, pos, node_color=nx.get_node_attributes(G, 'color').values(),
                           node_size=list(nx.get_node_attributes(G, 'size').values()),
                           node_shape='D')
    nx.draw_networkx_edges(G, pos, edge_color=nx.get_edge_attributes(G, 'color').values(),
                           width=list(nx.get_edge_attributes(G, 'width').values()))
    fontsize_d = {y: x['size'] for x, y in zip(list(nx.get_node_attributes(G, 'font').values()),
                                               list(nx.get_node_attributes(G, 'label').values())) if x['size'] != 0.0}
    fontcolor_d = {y: x['color'] for x, y in zip(list(nx.get_node_attributes(G, 'font').values()),
                                                 list(nx.get_node_attributes(G, 'label').values())) if x['size'] != 0.0}
    for node, (x, y) in pos.items():
        if node in fontsize_d.keys():
            plt.text(x, y, node, fontsize=fontsize_d[node], color=fontcolor_d[node], ha='center', va='center')
    ax = plt.gca()
    ax.margins(0.11)
    plt.tight_layout()
    plt.axis("off")
    plt.savefig(save, format='pdf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='scmulti',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='',
        add_help=True,
    )
    parser.add_argument('-o', '--output', help='output directory to save all the intermediate and final results.')
    parser.add_argument('--scenic', help='scplus object')
    parser.add_argument('-v', '--variable', default='celltype', help='Meta data label storing cell type information')
    parser.add_argument('-ct', '--celltypes', nargs='*', type=str, default=None, help='Interested cell types')
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    # scplus_obj = dill.load(open(args.scenic, 'rb'))
    scplus_obj = dill.load(open(os.path.join(args.output, 'scplus_obj_downstream.pkl'), 'rb'))

    if args.celltypes is None:
        celltypes = sorted(list(set(scplus_obj.metadata_cell['celltype'])))
    else:
        celltypes = sorted(args.celltypes)
    print(f'Using Cell Types: {celltypes}')

    dmulti = DownScMulti(scplus_obj, args.output, celltypes=celltypes)
    n_targets, rho, adj_pval, thresholds = dmulti.calculate()
    dmulti.plot_correlation(rho, n_targets, adj_pval, thresholds)
    dmulti.plot_heatmap()
    dmulti.plot_erss(num_columns=2, top_n=10, figsize=(20, 30))
    dmulti.plot_jaccard_heatmap()
    dmulti.plot_dim_reduction()
    dmulti.plot_eregulon()

    if dmulti.selected_eRegulons is None:
        selected_regulons = select_regulons(dmulti.scplus_obj)
        target_regulons = [i.split('_')[0] for i in selected_regulons]
    else:
        target_regulons = dmulti.selected_eRegulons
    for celltype in celltypes:
        dmulti.plot_network(celltype, target_regulons, n_top_features=1000, figsize=(20, 20))
    # except Exception as e:
    #     raise (e)
    # finally:
    #     dill.dump(dmulti.scplus_obj, open(os.path.join(args.output, 'scplus_obj_downstream.pkl'), 'wb'), protocol=-1)
