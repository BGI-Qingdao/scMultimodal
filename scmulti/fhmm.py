#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 21 Nov 2024 17:26
# @Author: Yao LI
# @File: scmulti/fhmm.py
import os
import sys
import time
from dataclasses import dataclass
from typing import List

import pickle
import argparse
import numpy as np
import pandas as pd
from hmmlearn import hmm

from infer_tf_tg import fRegion, fRegulons, fTRM
from infer_tf_tg import read_fasta, read_fragments, read_regulon_meta, handle_eRegulons_meta, get_region_reads


@dataclass
class ModelContainer:
    """
    For each TF, models = {'Chrom': [HMM]}
    """
    tf: str
    models: dict


class fHMM:
    """A HMM for each chromosome a TF has ..."""

    def __init__(self, TF: fTRM, chrom: str, states: list = None, n_iter: int = 100, random_state: int = 42):
        # Model Parameters
        if states is None:
            states = ['A', 'a', 'T', 't', 'C', 'c', 'G', 'g', 'N', 'n']
        self.states = states
        self.n_states = len(states)
        self.n_iter = n_iter
        self.random_state = random_state
        self.init_em = None  # the initial emission matrix aka start probability matrix
        self.emission_matrix = None  # the emission matrix b/w hidden states and the observed sequence
        self.trans_matrix = None  # the transitional probability matrix b/w hidden states

        # Data Attributes
        # 1. reference data
        self.tf_name = TF.tf  # tf name
        self.enhancer = TF.region  # designated enhancer! of the TF
        self.region_length = len(TF.region)
        self.chrom = chrom
        # 2. Training data
        self.aligned_sequences = []

    def init_emission_matrix(self, trm: fTRM, base_prob: float = 0.1, state_prob: float = 0.7) -> np.array:
        """
        Create a start probability matrix based on known regulons results
        a matrix where
        rows corresponding to positions in the region and columns corresponding to the states A, C, G, T, N
        :param trm:
        :param base_prob:
        :param state_prob:
        :return:
        """
        region_list = trm[self.chrom]
        if len(region_list) == 1:
            region = region_list[0]
            assert len(region.seq) > 0
            # Fill with 0.1 as the base probability
            start_prob_matrix = np.full((len(region), self.n_states), base_prob)
            # Set 0.7 probability for the most probable state at each position
            for i, state in enumerate(region.seq):
                start_prob_matrix[i, state] = state_prob

            # normalize? matrix
            start_prob_matrix = start_prob_matrix / start_prob_matrix.sum(axis=1, keepdims=True)

            self.init_em = start_prob_matrix
            return start_prob_matrix
        # for region in region_list:
        else:
            print(f'{self.tf_name} has {len(region_list)} regions on chromosome {self.chrom}!')

    def cal_emission_matrix(self):
        from collections import Counter
        # 初始化参考区域碱基频率矩阵
        base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        n_states = 4  # A, C, G, T
        base_counts = np.zeros((self.region_length, n_states))  # 初始化计数矩阵
        # 填充参考区域的碱基频率
        for aln_seq in self.aligned_sequences:
            seq = aln_seq["sequence"]
            start = aln_seq["start"]

            for i, base in enumerate(seq):
                if base in base_map:
                    position = start + i  # 参考区域的位置
                    base_counts[position, base_map[base]] += 1
        # 归一化为概率
        emission_probs = base_counts / base_counts.sum(axis=1, keepdims=True)
        return emission_probs

    def cal_transitional_matrix(self, stay_prob=0.1, trans_prob=0.9, last_prob=1.0):
        # init
        transition_probs = np.zeros((self.region_length, self.region_length))
        # TODO: 定义转移规则
        for i in range(self.region_length - 1):
            transition_probs[i, i] = stay_prob  # 停留在当前状态
            transition_probs[i, i + 1] = trans_prob  # 转移到下一个状态
        transition_probs[-1, -1] = last_prob  # 最后一个状态只能停留
        print("Transition Probability Matrix (sample):")
        print(transition_probs[:5, :5])

    def match_positions(self, regions: List[fRegion], normalize=False):
        new_nuleotides = np.zeros((len(regions), len(self.enhancer)), dtype=int)
        for region in regions:
            # # check if region is in
            # if self.enhancer.start <= region.start <= self.enhancer.end:
            #     pass
            region.relative_start = region.start - self.enhancer.start
            region.relative_end = region.end - self.enhancer.start

    def train(self, trm, emission_probs, chrom, n_iter=500, random_state=42):
        """
        Train a HMM for one TF on each chromosome
        :param emission_probs:
        :param chrom:
        :param n_iter:
        :param random_state:
        :return:
        """
        ref_length = 500
        # encoded_aligned_sequences = tf_training_dataset[chrom].encoded_sequences
        # Train HMM
        model = hmm.MultinomialHMM(n_components=ref_length, n_iter=n_iter, random_state=random_state)
        # initialize HMM parameters
        model.startprob_ = self.init_emission_matrix(trm)  # init start probability matrix based on eRegulons
        model.transmat_ = np.eye(ref_length, k=1) + np.eye(ref_length, k=0) * 0.1  # TODO: 简单线性模型
        model.transmat_ /= model.transmat_.sum(axis=1, keepdims=True)  # normalize
        model.emissionprob_ = emission_probs  # calculate emission probability matrix

        print("HMM successfully built with custom probabilities.")
        return model


def encode_sequence(seq) -> np.array:
    """S"""
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3, 'N': 5, 'n': 5}
    return np.array([encoding[base] for base in seq])


def get_sequence(genome: dict, chrom: str, start: int, end: int) -> np.array:
    seq = genome[chrom][start - 1:end]
    encoded_seq = encode_sequence(seq)
    return encoded_seq


def sequence_to_counts(sequence, vocab_size=10) -> np.array:
    """S"""
    counts = np.zeros((len(sequence), vocab_size), dtype=int)
    for i, symbol in enumerate(sequence):
        counts[i, symbol] = 1
    return counts


def predict(model, chrom, genome, window_size=500, step=1, threshold=-800) -> list:
    """

    :param model:
    :param chrom:
    :param genome:
    :param window_size: Length of binding region
    :param step:
    :param threshold: Adjust threshold based on your data
    :return:
    """
    chrome_seq = genome[chrom]  # TODO
    encoded_chrom = encode_sequence(chrome_seq)

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
    predicted_sites = [(start, end) for start, end, score in scores if score > threshold]
    print("Predicted binding sites:", predicted_sites)

    # Or top 1%
    # scores_array = np.array([score for _, _, score in scores])
    # dynamic_threshold = np.percentile(scores_array, 99)  # Top 1% scores
    # predicted_sites = [(start, end) for start, end, score in scores if score > dynamic_threshold]
    return predicted_sites


# -------------------------------------------------------
# OTHERS:
# -------------------------------------------------------
def create_parser():
    parser = argparse.ArgumentParser(
        prog='scmulti_hmm',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=True,
    )
    parser.add_argument('-o', '--output', help='output directory to save all the intermediate and final results.')
    parser.add_argument('-s', '--scenic', help='scenic+ results pickle file name')
    parser.add_argument('-f', '--fragments', help='scATAC-seq fragments.tsv.gz file name')
    parser.add_argument('-g', '--genome', help='reference genome file name')
    return parser


def main():
    # Parse arguments.
    theparser = create_parser()
    args = theparser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # reference genome
    genome = read_fasta(args.genome)
    # associated sequencing data: scRNA-seq, scATAC-seq, etc.
    fragments = read_fragments(args.fragments)
    # Scenic+ results
    # meta_data = get_regulon_meta(args.scenic)
    meta_data = read_regulon_meta(args.scenic)

    # Parse fish Transcriptional regulatory modules (TRMs)
    itemlist = handle_eRegulons_meta(meta_data)
    # with open(f'{args.output}/regulons.pkl', 'wb') as pickle_file:
    #     pickle.dump(itemlist, pickle_file)

    # Create training sequences
    total_training_dataset = {}
    for tf, fregulons in itemlist.items():
        training_dataset = get_training_data(fregulons, fragments, genome)
        print(training_dataset)
        total_training_dataset[tf] = training_dataset
        break


# -------------------------------------------------------
# MODEL:
# -------------------------------------------------------
@dataclass(frozen=True)
class TrainingSequences:
    """T"""
    tf: str
    chrom: str
    encoded_sequences: list


def get_training_data(regulons: fRegulons, fragments: pd.DataFrame, genome: dict) -> dict:
    """

    :param tf:
    :param regulons:
    :return:
    """
    training_dataset = {}
    # for one TF only
    start_time = time.time()
    for trm in regulons[:3]:
        regions_reads = get_region_reads(fragments, trm.region.chrom, trm.region.start, trm.region.end)
        training_seqs = list(
            regions_reads.apply(lambda row: get_sequence(genome, row.chrom, row.start, row.end), axis=1)
        )
        # save data in a dictionary
        # training_dataset[trm.chrom] = training_seqs
        ts = TrainingSequences(chrom=trm.chrom, tf=trm.tf, encoded_sequences=training_seqs)
        training_dataset[trm.chrom] = ts
    end_time = time.time()
    print(f'Generating training dataset done in {round((end_time-start_time)/60, 2)} minutes.')
    return training_dataset


if __name__ == '__main__':
    main()
