#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 29 Jan 2024 17:03
# @Author: Yao LI
# @File: evo_fish/downstream.py
import os
import sys
import dill
import pickle
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

work_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron'
# scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))
# scrna_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/scrna/a_liver.h5ad'
# adata = sc.read_h5ad(scrna_fn)  # key_to_group_by, need to match with cistopic_obj
# celltypes = sorted(list(set(adata.obs['celltype'])))
celltypes = ['Adipocyte', 'Cholangiocytes', 'Endothelial', 'Hepatocyte', 'Red']
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj_downstream.pkl'), 'rb'))

# ---------------------------------------------------------------
# Calculate correlation coefficient
# ---------------------------------------------------------------
# from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
# apply_std_filtering_to_eRegulons(scplus_obj)
# from scenicplus.eregulon_enrichment import score_eRegulons
# region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/region_ranking.pkl'),
#                                 'rb'))  # load ranking calculated using the wrapper function
# gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/gene_ranking.pkl'),
#                               'rb'))  # load ranking calculated using the wrapper function
# score_eRegulons(scplus_obj,
#                 ranking=region_ranking,
#                 eRegulon_signatures_key='eRegulon_signatures_filtered',
#                 key_added='eRegulon_AUC_filtered',
#                 enrichment_type='region',
#                 auc_threshold=0.05,
#                 normalize=False,
#                 n_cpu=5)
# score_eRegulons(scplus_obj,
#                 gene_ranking,
#                 eRegulon_signatures_key='eRegulon_signatures_filtered',
#                 key_added='eRegulon_AUC_filtered',
#                 enrichment_type='gene',
#                 auc_threshold=0.05,
#                 normalize=False,
#                 n_cpu=5)
# 
# from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks
# 
# generate_pseudobulks(
#     scplus_obj=scplus_obj,
#     variable='celltype',
#     auc_key='eRegulon_AUC_filtered',
#     signature_key='Gene_based')
# generate_pseudobulks(
#     scplus_obj=scplus_obj,
#     variable='celltype',
#     auc_key='eRegulon_AUC_filtered',
#     signature_key='Region_based')
# 
# TF_cistrome_correlation(
#     scplus_obj,
#     use_pseudobulk=True,
#     variable='celltype',
#     auc_key='eRegulon_AUC_filtered',
#     signature_key='Gene_based',
#     out_key='filtered_gene_based')
# TF_cistrome_correlation(
#     scplus_obj,
#     use_pseudobulk=True,
#     variable='celltype',
#     auc_key='eRegulon_AUC_filtered',
#     signature_key='Region_based',
#     out_key='filtered_region_based')
# 
# n_targets = [int(x.split('(')[1].replace('r)', '')) for x in
#              scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Cistrome']]
# rho = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'].to_list()
# adj_pval = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Adjusted_p-value'].to_list()

# thresholds = {
#     'rho': [-0.75, 0.70],
#     'n_targets': 0
# }
# thresholds = {
#     'rho': [-0.30, 0.30],
#     'n_targets': 0
# }
# ---------------------------------------------------------------
# 1. Plot correlation coefficient
# ---------------------------------------------------------------
# fig, ax = plt.subplots(figsize=(10, 5))
# sc = ax.scatter(rho, n_targets, c=-np.log10(adj_pval), s=5)
# ax.set_xlabel('Correlation coefficient')
# ax.set_ylabel('nr. target regions')
# # ax.hlines(y = thresholds['n_targets'], xmin = min(rho), xmax = max(rho), color = 'black', ls = 'dashed', lw = 1)
# ax.vlines(x=thresholds['rho'], ymin=0, ymax=max(n_targets), color='black', ls='dashed', lw=1)
# ax.text(x=thresholds['rho'][0], y=max(n_targets), s=str(thresholds['rho'][0]))
# ax.text(x=thresholds['rho'][1], y=max(n_targets), s=str(thresholds['rho'][1]))
# sns.despine(ax=ax)
# fig.colorbar(sc, label='-log10(adjusted_pvalue)', ax=ax)
# plt.savefig('cor_coe.png')

# selected_cistromes = scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].loc[
#     np.logical_or(
#         scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] > thresholds['rho'][1],
#         scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based']['Rho'] < thresholds['rho'][0]
#     )]['Cistrome'].to_list()
# selected_eRegulons = [x.split('_(')[0] for x in selected_cistromes]
# selected_eRegulons_gene_sig = [
#     x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Gene_based'].keys()
#     if x.split('_(')[0] in selected_eRegulons]
# selected_eRegulons_region_sig = [
#     x for x in scplus_obj.uns['eRegulon_signatures_filtered']['Region_based'].keys()
#     if x.split('_(')[0] in selected_eRegulons]
# # save the results in the scenicplus object
# scplus_obj.uns['selected_eRegulon'] = {'Gene_based': selected_eRegulons_gene_sig,
#                                        'Region_based': selected_eRegulons_region_sig}
# print(f'selected: {len(selected_eRegulons_gene_sig)} eRegulons')
# dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj_downstream.pkl'), 'wb'), protocol=-1)

# ---------------------------------------------------------------
# 2. heatmap
# ---------------------------------------------------------------
# from scenicplus.plotting.dotplot import heatmap_dotplot
# heatmap_dotplot(
#     scplus_obj=scplus_obj,
#     size_matrix=scplus_obj.uns['eRegulon_AUC_filtered']['Region_based'],
#     # specify what to plot as dot sizes, target region enrichment in this case
#     color_matrix=scplus_obj.to_df('EXP'),  # specify  what to plot as colors, TF expression in this case
#     scale_size_matrix=True,
#     scale_color_matrix=True,
#     group_variable='celltype',
#     subset_eRegulons=scplus_obj.uns['selected_eRegulon']['Gene_based'],
#     figsize=(50, 50),
#     orientation='vertical',
#     save='heatmap.pdf',
#     index_order=celltypes)


# ---------------------------------------------------------------
# 3. rss plot
# ---------------------------------------------------------------
# from scenicplus.RSS import *
#
# regulon_specificity_scores(
#     scplus_obj,
#     variable='celltype',
#     auc_key='eRegulon_AUC_filtered',
#     signature_keys=['Region_based'],
#     selected_regulons=[x for x in scplus_obj.uns['selected_eRegulon']['Region_based'] if '-' not in x],
#     out_key_suffix='_filtered')
# plot_rss(scplus_obj, 'celltype_filtered', num_columns=2, top_n=10, figsize=(20, 30), save='rss.pdf')

# ---------------------------------------------------------------
# 4. Jaccard Heatmap
# ---------------------------------------------------------------
# flat_list = lambda t: [item for sublist in t for item in sublist]
# selected_markers = list(set(flat_list(
#     [scplus_obj.uns['RSS']['celltype_filtered'].loc[celltype].sort_values(ascending=False).head(30).index.to_list()
#      for celltype in scplus_obj.uns['RSS']['celltype_filtered'].index])))
# from scenicplus.plotting.correlation_plot import *
# region_intersetc_data, Z = jaccard_heatmap(
#     scplus_obj,
#     method='intersect',
#     gene_or_region_based='Region_based',
#     use_plotly=False,
#     selected_regulons=selected_markers,
#     signature_key='eRegulon_signatures_filtered',
#     figsize=(10, 10), return_data=True, vmax=0.5, cmap='plasma', save='jaccard_heatmap.pdf')

# ---------------------------------------------------------------
# 5. Plot dimensionality reductions
# ---------------------------------------------------------------
# from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
# run_eRegulons_umap(
#     scplus_obj=scplus_obj,
#     auc_key='eRegulon_AUC_filtered',
#     reduction_name='eRegulons_UMAP',  # overwrite previously calculated UMAP
# )
# run_eRegulons_tsne(
#     scplus_obj=scplus_obj,
#     auc_key='eRegulon_AUC_filtered',
#     reduction_name='eRegulons_tSNE',  # overwrite previously calculated tSNE
# )
# from scenicplus.dimensionality_reduction import plot_metadata_given_ax
# # specify color_dictionary
# color_dict = {
#     'Adipocyte': "#065143",
#     'Cholangiocytes': "#70B77E",
#     'Endothelial': "#E0A890",
#     'Hepatocyte': "#F56476",
#     'Red': "#CE1483",
# }
#
# fig, axs = plt.subplots(ncols=2, figsize=(16, 8))
# plot_metadata_given_ax(
#     scplus_obj=scplus_obj,
#     ax=axs[0],
#     reduction_name='eRegulons_UMAP',
#     variable='celltype',
#     color_dictionary={'celltype': color_dict}
# )
# plot_metadata_given_ax(
#     scplus_obj=scplus_obj,
#     ax=axs[1],
#     reduction_name='eRegulons_tSNE',
#     variable='celltype',
#     color_dictionary={'celltype': color_dict}
# )
# fig.tight_layout()
# sns.despine(ax=axs[0])  # remove top and right edge of axis border
# sns.despine(ax=axs[1])  # remove top and right edge of axis border
# plt.savefig('dim_reduc.pdf', format='pdf')

# ---------------------------------------------------------------
# 7. Plot activity / expression of an eRegulon
# ---------------------------------------------------------------
# from scenicplus.dimensionality_reduction import plot_eRegulon
# plot_eRegulon(
#     scplus_obj=scplus_obj,
#     reduction_name='eRegulons_tSNE',
#     selected_regulons=['MSTRG.2388_extended_+', 'MSTRG.8533_extended_+', 'nbisL1-mrna-14072_extended_+',
#                        'MSTRG.18385_extended_+', 'MSTRG.14903_extended_+'],
#     scale=True,
#     auc_key='eRegulon_AUC_filtered',
#     save='eregulon_activity.pdf')
#
# dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj_downstream.pkl'), 'wb'), protocol=-1)
# ---------------------------------------------------------------
# 8. Plot Network
# ---------------------------------------------------------------
# from pycisTopic.diff_features import find_highly_variable_features
# hvr = find_highly_variable_features(
#     scplus_obj.to_df('ACC').loc[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Region']))], n_top_features=1000,
#     plot=False)
# hvg = find_highly_variable_features(
#     scplus_obj.to_df('EXP')[list(set(scplus_obj.uns['eRegulon_metadata_filtered']['Gene']))].T, n_top_features=1000,
#     plot=False)
# from scenicplus.networks import create_nx_tables, create_nx_graph, plot_networkx, export_to_cytoscape

# nx_tables = create_nx_tables(
#     scplus_obj=scplus_obj,
#     eRegulon_metadata_key='eRegulon_metadata_filtered',
#     subset_regions=hvr,
#     subset_genes=hvg,
#     add_differential_gene_expression=True,
#     add_differential_region_accessibility=True,
#     differential_variable=['celltype'])
# pickle.dump(nx_tables, open(os.path.join(work_dir, 'nx_tables.pkl'), 'wb'))
# nx_tables = pickle.load(open(os.path.join(work_dir, 'nx_tables.pkl'), 'rb'))
# G, pos, edge_tables, node_tables = create_nx_graph(nx_tables,
#                                                    use_edge_tables=['TF2R', 'R2G'],
#                                                    color_edge_by={'TF2R': {'variable': 'TF',
#                                                                            'category_color': {'MSTRG.13254': 'Orange',
#                                                                                               'MSTRG.1399': 'Purple',
#                                                                                               'nbisL1-mrna-3223': 'Red',
#                                                                                               'MSTRG.14903': 'Yellow',
#                                                                                               'MSTRG.2388': 'Blue',
#                                                                                               'MSTRG.15165': 'Green',
#                                                                                               'MSTRG.9306': 'Pink',
#                                                                                               'nbisL1-mrna-7336': 'Gray',
#                                                                                               'MSTRG.12536': 'Brown'}},
#                                                                   'R2G': {'variable': 'R2G_rho',
#                                                                           'continuous_color': 'viridis', 'v_min': -1,
#                                                                           'v_max': 1}},
#                                                    transparency_edge_by={
#                                                        'R2G': {'variable': 'R2G_importance', 'min_alpha': 0.1,
#                                                                'v_min': 0}},
#                                                    width_edge_by={'R2G': {'variable': 'R2G_importance', 'max_size': 1.5,
#                                                                           'min_size': 1}},
#                                                    color_node_by={'TF': {'variable': 'TF',
#                                                                          'category_color': {'MSTRG.13254': 'Orange',
#                                                                                             'MSTRG.1399': 'Purple',
#                                                                                             'nbisL1-mrna-3223': 'Red',
#                                                                                             'MSTRG.14903': 'Yellow',
#                                                                                             'MSTRG.2388': 'Blue',
#                                                                                             'MSTRG.15165': 'Green',
#                                                                                             'MSTRG.9306': 'Pink',
#                                                                                             'nbisL1-mrna-7336': 'Gray',
#                                                                                             'MSTRG.12536': 'Brown'}},
#                                                                   'Gene': {'variable': 'celltype_Log2FC_Hepatocyte',
#                                                                            'continuous_color': 'bwr'},
#                                                                   'Region': {
#                                                                       'variable': 'celltype_Log2FC_Hepatocyte',
#                                                                       'continuous_color': 'viridis'}},
#                                                    transparency_node_by={
#                                                        'Region': {'variable': 'celltype_Log2FC_Hepatocyte',
#                                                                   'min_alpha': 0.1},
#                                                        'Gene': {'variable': 'celltype_Log2FC_Hepatocyte',
#                                                                 'min_alpha': 0.1}},
#                                                    size_node_by={'TF': {'variable': 'fixed_size', 'fixed_size': 30},
#                                                                  'Gene': {'variable': 'fixed_size', 'fixed_size': 15},
#                                                                  'Region': {'variable': 'fixed_size',
#                                                                             'fixed_size': 10}},
#                                                    shape_node_by={
#                                                        'TF': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                                        'Gene': {'variable': 'fixed_shape', 'fixed_shape': 'ellipse'},
#                                                        'Region': {'variable': 'fixed_shape', 'fixed_shape': 'diamond'}},
#                                                    label_size_by={
#                                                        'TF': {'variable': 'fixed_label_size', 'fixed_label_size': 20.0},
#                                                        'Gene': {'variable': 'fixed_label_size',
#                                                                 'fixed_label_size': 10.0},
#                                                        'Region': {'variable': 'fixed_label_size',
#                                                                   'fixed_label_size': 0.0}},
#                                                    layout='kamada_kawai_layout',
#                                                    scale_position_by=250)
#
# plt.figure(figsize=(20, 20))
# plot_networkx(G, pos)
# export_to_cytoscape(G, pos, out_file=os.path.join(work_dir, 'scenicplus/network_Hepatocyte.cys'))
