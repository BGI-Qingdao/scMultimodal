#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 16 Oct 2024 15:20
# @Author: Yao LI
# @File: scmulti/infer_tf_tg.py
# Use HMM to infer TFs and Target Genes relationship
import os
import sys
import time
from dataclasses import dataclass

import pickle
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from hmmlearn import hmm


# Load & extracting sequencing data
def read_fragments(fn):
    df = pd.read_csv(fn, compression='gzip', header=None, sep='\t', quotechar='"')
    df.columns = ['chrom', 'start', 'end', 'cell', 'count']
    return df


def get_region_reads(fragments: pd.DataFrame, chrom: str, start: int, end: int):
    sub_df = fragments[fragments.chrom == chrom]
    regions_reads = sub_df[sub_df.start >= start]
    regions_reads = regions_reads[regions_reads.end <= end]
    return regions_reads


# -------------------------------------------------------
# DATA:
# -------------------------------------------------------
class fRegion:
    def __init__(self, region_id: str):
        self.region_id = region_id
        self.chrom, self.start, self.end = get_region_info(region_id)
        self.seq = ''

    def get_seq(self, genome: dict):
        """Extract region sequences when needed."""
        self.seq = genome[self.chrom][self.start - 1:self.end]  # TODO: starts from 1?
        return self.seq

    def __repr__(self):
        return f"fRegion: {self.region_id}"


class fTRM:
    """Fish Transcriptional regulatory module (TRM)
        A TRM contains TF, regions and target genes on a single chromosome."""

    def __init__(self, tf, region_id, targets=None, celltype=None, species=None):
        if targets is None:
            targets = []
        self.tf = tf  # tf name
        self.region = fRegion(region_id)  # a region id, only one
        self.chrom = self.region.chrom
        self.targets = targets
        # optional
        self.celltype = celltype
        self.species = species

    def __repr__(self):
        return f"fTRM(TF: {self.tf}, {self.region}, Targets: {self.targets})"


class fRegulons:
    """
    A list of fTRMs
    """

    def __init__(self):
        self.regulons = {}  # Dictionary to store lists of items

    def add_trm(self, trm):
        """
        Add an Item to the dictionary.
        If the name already exists, append the item to the list for that name.
        """
        if not isinstance(trm, fTRM):
            raise ValueError("Only instances of the Item class can be added.")
        if trm.tf not in self.regulons:
            self.regulons[trm.tf] = []
        self.regulons[trm.tf].append(trm)

    def remove_trm(self, tf_name, region: fRegion = None):
        """
        Remove items by name.
        If a value is provided, remove only the items with that specific value.
        """
        if tf_name in self.regulons:
            if region is None:
                del self.regulons[tf_name]  # Remove all items with this name
            else:
                self.regulons[tf_name] = [
                    trm for trm in self.regulons[tf_name] if trm.region != region
                ]
                if not self.regulons[tf_name]:  # Remove the key if the list is empty
                    del self.regulons[tf_name]

    def __repr__(self):
        return f"ItemList({self.regulons})"

    def __getitem__(self, tf_name):
        return self.regulons.get(tf_name, [])  # Return the list of items or an empty list

    def __setitem__(self, tf_name, items):
        if not all(isinstance(trm, fTRM) for trm in items):
            raise ValueError("Value must be a list of Item instances.")
        self.regulons[tf_name] = items

    def __delitem__(self, tf_name):
        del self.regulons[tf_name]

    def __len__(self):
        return len(self.regulons)  # Number of unique keys

    def __contains__(self, tf_name):
        return tf_name in self.regulons

    def keys(self):
        """Return all keys (names of items)."""
        return self.regulons.keys()

    def values(self):
        """Return all values (lists of items)."""
        return self.regulons.values()

    def items(self):
        """Return all key-value pairs as (name, list of items)."""
        return self.regulons.items()


# Load in reference genome
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


# Load eRegulon results
def get_regulon_meta(scplus_fn: str) -> pd.DataFrame:
    """
    Extract eRegulon metadata from Scenic+ output (pickle file)
    (TF and its corresponding regions (containing motif sequence) info)
    :param scplus_fn:
    :return:
    """
    scplus_obj = pickle.load(open(scplus_fn, 'rb'))
    meta_data = scplus_obj.uns['eRegulon_metadata'][['TF', 'Region', 'Gene', 'R2G_importance']].drop_duplicates(
        inplace=True)
    meta_data['Chromosome'] = meta_data.Region.str.split(':', expand=True)[0]
    return meta_data


def read_regulon_meta(fn):
    """Read eRegulon metadata csv file"""
    meta_data = pd.read_csv(fn, sep='\t')
    meta_data = meta_data[['TF', 'Region', 'Gene', 'R2G_importance']].drop_duplicates(inplace=False)
    meta_data['Chromosome'] = meta_data.Region.str.split(':', expand=True)[0]
    return meta_data


def get_region_info(region_id: str):
    """
    Extract region information from region id in Scenic+ eRegulon metadata
    :param region_id:
    :return: chromosome id, region start position, region end position
    """
    chrom, pos = region_id.split(':')
    start, end = pos.split('-')
    return chrom, int(start), int(end)


def handle_eRegulons_meta(eRegulons_meta: pd.DataFrame) -> fRegulons:
    """
    For each TF-region-target, create a fTRM object
    :param eRegulons_meta:
    :return: rRegulons object, aka a dictionary of fTRMs, key is TF name, value is a list of fTRMs with the same TF name
    """
    regulons = fRegulons()
    TFs = list(set(eRegulons_meta.TF))
    for tf in TFs:
        tf_df = eRegulons_meta[eRegulons_meta.TF == tf]
        # for each chromosome:
        chroms = list(set(tf_df.Chromosome))
        for chrom in chroms:
            sub_df = tf_df[tf_df.Chromosome == chrom]
            # on the same chromosome, for each region (TF binding region, containing motifs)
            for index, row in sub_df.iterrows():
                # for each TF, get its target gene names on one chromosome
                target_genes = list(set(sub_df.Gene))
                ftrm = fTRM(tf, row.Region, targets=target_genes)
                regulons.add_trm(ftrm)
    return regulons


# -------------------------------------------------------
# MODEL:
# -------------------------------------------------------
@dataclass(frozen=True)
class TrainingSequences:
    """T"""
    chrom: str
    tf: str
    encoded_sequences: list


def get_training_data(tf, regulons: fRegulons, fragments: pd.DataFrame, genome: dict):
    """

    :param tf:
    :param regulons:
    :return:
    """
    training_dataset = {}
    # for one TF only
    tf_regulons = regulons[tf]  # a list
    for trm in tf_regulons:
        regions_reads = get_region_reads(fragments, trm.region.chrom, trm.region.start, trm.region.end)
        training_seqs = list(
            regions_reads.apply(lambda row: get_sequence(genome, row.chrom, row.start, row.end), axis=1))
        # save data in a dictionary
        # training_dataset[trm.chrom] = training_seqs
        ts = TrainingSequences(chrom=trm.chrom, tf=trm.tf, encoded_sequences=training_seqs)
        training_dataset[trm.chrom] = ts
    return training_dataset


def encode_sequence(seq) -> np.array:
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3, 'N': 5, 'n': 5}
    return np.array([encoding[base] for base in seq])


def get_sequence(genome: dict, chrom: str, start: int, end: int) -> np.array:
    seq = genome[chrom][start - 1:end]
    encoded_seq = encode_sequence(seq)
    return encoded_seq


# def encode_sequence(seq, vocab="ACGTatcgNn"):
#     encoding = {char: idx for idx, char in enumerate(vocab)}
#     return [encoding[base] for base in seq]


def sequence_to_counts(sequence, vocab_size=10):
    counts = np.zeros((len(sequence), vocab_size), dtype=int)
    for i, symbol in enumerate(sequence):
        counts[i, symbol] = 1
    return counts


def train4tf(tf_training_dataset, chrom, n_states=5, n_iter=100, random_state=42):
    """Train a HMM for one TF on each chromosome"""
    encoded_training_sequences = tf_training_dataset[chrom].encoded_sequences
    # Train HMM
    model = hmm.MultinomialHMM(n_components=n_states, n_iter=n_iter, random_state=random_state)

    # WRONG! TODO:
    lengths = [len(seq) for seq in encoded_training_sequences]
    X = np.concatenate(encoded_training_sequences).reshape(-1, 1)
    model.fit(X, lengths)
    return model


def predict_region(model, chrom, genome, fragments, window_size=500, step=1, threshold=-800):
    # ts = tf_training_dataset[chrom]
    chrome_seq = genome[chrom]  # TODO
    encoded_chrom = encode_sequence(chrome_seq)
    # window_size = 500  # Length of binding region

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
    # threshold = -800  # Adjust threshold based on your data
    predicted_sites = [(start, end) for start, end, score in scores if score > threshold]
    print("Predicted binding sites:", predicted_sites)

    # Or top 1%
    scores_array = np.array([score for _, _, score in scores])
    dynamic_threshold = np.percentile(scores_array, 99)  # Top 1% scores
    predicted_sites = [(start, end) for start, end, score in scores if score > dynamic_threshold]
    # predicted_sites: not right?
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
    TFS = list(set(meta_data.TF))

    # Parse fish Transcriptional regulatory modules (TRMs)
    itemlist = handle_eRegulons_meta(meta_data)
    with open(f'{args.output}/regulons.pkl', 'wb') as pickle_file:
        pickle.dump(itemlist, pickle_file)

    # Create training sequences


if __name__ == '__main__':
    main()
