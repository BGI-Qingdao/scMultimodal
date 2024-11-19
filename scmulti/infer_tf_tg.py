#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 16 Oct 2024 15:20
# @Author: Yao LI
# @File: scmulti/infer_tf_tg.py
# Use HMM to infer TFs and Target Genes relationship
import os
import sys
import time

import gzip
import pickle
import argparse
import textwrap
import numpy as np
import pandas as pd
import scanpy as sc
from hmmlearn import hmm


class fTRM:
    """Fish Transcriptional regulatory module (TRM)
    A TRM contains TF, motifs and target genes."""

    def __init__(self):
        pass


def read_fragments(fn):
    df = pd.read_csv(fn, compression='gzip', header=None, sep='\t', quotechar='"')
    df.columns = ['chrom', 'start', 'end', 'cell', 'count']
    return df


def get_region_reads(fragments: pd.DataFrame, chrom: str, start, end):
    sub_df = fragments[fragments.chrom == chrom]
    regions_reads = sub_df[sub_df.start >= start]
    regions_reads = regions_reads[regions_reads.end <= end]
    return regions_reads


def train_sequences(regions_reads, ref_genome, chrom):
    seqs = []
    for index, row in regions_reads.iterrows():
        seq = ref_genome[chrom][row.start - 1:row.end]
        seqs.append(seq)
    return seqs


def get_regions_meta(scplus_obj: sc.AnnData):
    """Extract TF and its corresponding regions (containing motif sequence) info"""
    meta = scplus_obj.uns['eRegulon_metadata'][['TF', 'Region', 'Gene', 'R2G_importance']].drop_duplicates(inplace=True)
    return meta


def get_seq(region_id, ref_genome):
    chrom, pos = region_id.split(':')
    start, end = pos.split('-')
    seq = ref_genome[chrom][int(start) - 1:int(end)]  # TODO: starts from 1?
    return chrom, seq


def train4tf(tf_meta, ref_genome):
    """Train a HMM for one TF on each chromosome"""
    regions = list(set(tf_meta.Region))
    known_regions = [get_seq(region, ref_genome) for region in regions]
    train_sequences = [encode_sequence(seq) for seq in known_regions]
    # Train HMM
    n_states = 5  # Choose a suitable number of states
    model = hmm.MultinomialHMM(n_components=n_states, n_iter=100, random_state=42)

    chrome_seq = ref_genome[chrom]
    chrome_seq = genome[chrom][int(start) - 1000:int(end) + 1000]
    encoded_chr¢®ƒ√om = encode_sequence(chrome_seq)
    window_size = 500  # Length of binding region
    step = 1

    print(f'Scaning chromosome {chrom}...')
    start_time = time.time()
    scores = []
    for i in range(0, len(encoded_chrom) - window_size + 1, step):
        window = encoded_chrom[i:i + window_size]
        window_counts = sequence_to_counts(window)
        score = model.score(window_counts)
        scores.append((i, i + window_size, score))
    end_time = time.time()
    print(f'Done, used {round(end_time - start_time, 2)} seconds.')

    # Filter and output potential binding sites
    threshold = -800  # Adjust threshold based on your data
    predicted_sites = [(start, end) for start, end, score in scores if score > threshold]
    print("Predicted binding sites:", predicted_sites)

    # Or top 1%
    scores_array = np.array([score for _, _, score in scores])
    dynamic_threshold = np.percentile(scores_array, 99)  # Top 1% scores
    predicted_sites = [(start, end) for start, end, score in scores if score > dynamic_threshold]
    # predicted_sites: not right?


def create_parser():
    parser = argparse.ArgumentParser(
        prog='scmulti_hmm',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True,
    )
    parser.add_argument('-o', '--output', help='output directory to save all the intermediate and final results.')
    parser.add_argument('-i', '--input', help='scenic+ results')
    return parser


def read_fasta(fasta_file):
    """
    Read fa/fasta file
    :param fasta_file: protome fa file name
    :param min_len: amino acid sequence length cutoff (only seq longer than the len_cutoff will be kept)
    :return:
    """
    sequences = {}
    with open(fasta_file, "r") as f:
        current_sequence_id = None
        current_sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):  # header line
                # check if max_len > length > than cutoff
                if current_sequence_id:
                    sequences[current_sequence_id] = current_sequence
                current_sequence_id = line[1:].split(' ')[0]
                current_sequence = ""
            else:  # sequence line
                current_sequence += line
        # save last sequence (if meet length requirement)
        if current_sequence_id:
            sequences[current_sequence_id] = current_sequence
    return sequences


def encode_sequence(seq):
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3, 'N': 5, 'n': 5}
    return np.array([encoding[base] for base in seq])


def encode_sequence(seq, vocab="ACGTatcgNn"):
    encoding = {char: idx for idx, char in enumerate(vocab)}
    return [encoding[base] for base in seq]


def sequence_to_counts(sequence, vocab_size=10):
    counts = np.zeros((len(sequence), vocab_size), dtype=int)
    for i, symbol in enumerate(sequence):
        counts[i, symbol] = 1
    return counts


def main():
    # Parse arguments.
    theparser = create_parser()
    args = theparser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    scplus_obj = pickle.load(open(args.input, 'rb'))
    meta_data = get_regions_meta(scplus_obj)
    meta_tf = meta_data[meta_data.TF == 'LOC112142958']
    meta_tf['Chromosome'] = meta_tf.Region.str.split(':', expand=True)[0]
    print(meta_tf[meta_tf.chrom == 'NC_050513.1'])
    tfs = list(set(meta_data.TF))

    # states = tfs + ['null']  # n TFs + 1 None
    # print(f'Using {len(states)} hidden states.')

    species_genome_fn = ''
    genome = read_fasta(species_genome_fn)


if __name__ == '__main__':
    observations = ['Genes']  # m target genes in scRNA-seq data
    start_probability = {'TF1': 0.5, 'TF2': 0.5}  # based on eRegulons of OM fish
    transition_probability = {'TF1': {'TF1': 1, 'TF2': 0},
                              'TF2': {'TF1': 0,
                                      'TF2': 1}}  # state aka TF won't transfer to another TF, state never changes
    emission_probability = {'TF1': {'Gene1': 0.7, 'Gene2': 0.3},
                            'TF2': {'Gene1': 0.7, 'Gene2': 0.3}}  # adjacency value / importance

    # 读取scRNA-seq数据
    scrna_data_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.multi_modal/DATA/03.Oryzias_melastigma/04.annoted_scRNA'
    scrna_data_fn = 'OM_liver_anno.h5ad'
    rna_data = sc.read_h5ad(os.path.join(scrna_data_dir, scrna_data_fn))
    # 假设数据集中有基因表达矩阵
    gene_expression = rna_data.X
    gene_names = rna_data.var_names

    # 读取scATAC-seq数据
    atac_data = sc.read_h5ad(
        '/dellfsqd2/ST_OCEAN/USER/liyao1/11.multi_modal/exp/03.Oryzias_melastigma/02.label_transfer/liver/liver.merged.atac.h5ad')  # 请根据你的数据文件路径修改
    # 假设数据集中包含列: cell_id, peak_id, accessibility
    # 处理ATAC-seq数据，计算每个基因的平均可访问性
    accessibility_df = atac_data.groupby('gene_id')['accessibility'].mean().reset_index()

    # 结合表达和可访问性数据
    # 创建一个DataFrame，包含每个基因的表达和可访问性
    combined_df = pd.DataFrame({
        'gene_id': gene_names,
        'expression': gene_expression.mean(axis=0),  # 取平均值
        'accessibility': accessibility_df['accessibility']
    })

    # 对数据进行标准化处理
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    X = scaler.fit_transform(combined_df[['expression', 'accessibility']])

    # 定义HMM模型
    # 假设我们有两个隐藏状态，例如 'Active' 和 'Inactive'
    model = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=100)

    # 训练HMM模型
    model.fit(X)

    # 预测隐藏状态
    hidden_states = model.predict(X)

    # 将预测的状态添加到原始DataFrame中
    combined_df['predicted_state'] = hidden_states

    # 打印预测结果
    print(combined_df[['gene_id', 'expression', 'accessibility', 'predicted_state']])
