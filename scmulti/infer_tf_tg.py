#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 16 Oct 2024 15:20
# @Author: Yao LI
# @File: scmulti/infer_tf_tg.py
# Use HMM to infer TFs and Target Genes relationship
import os
import sys
from dataclasses import dataclass

import pickle
import argparse
import pandas as pd


# Load & extracting sequencing data
def read_fragments(fn) -> pd.DataFrame:
    """Read fragments.tsv.gz file generated by scATAC-seq"""
    df = pd.read_csv(fn, compression='gzip', header=None, sep='\t', quotechar='"')
    df.columns = ['chrom', 'start', 'end', 'cell', 'count']
    return df


def get_region_reads(fragments: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    """Given an interested region, extract scATAC fragments mapped to the region"""
    sub_df = fragments[fragments.chrom == chrom]
    regions_reads = sub_df[sub_df.start >= start]
    regions_reads = regions_reads[regions_reads.end <= end]
    return regions_reads


# -------------------------------------------------------
# DATA:
# -------------------------------------------------------
class fRegion:
    """
    A region on genome, mostly likely to be an Enhancer, 500bp long
    """

    def __init__(self, region_id: str):
        self.region_id = region_id
        self.chrom, self.start, self.end = get_region_info(region_id)
        self.seq = ''
        self.relative_start = None
        self.relative_end = None

    def get_seq(self, genome: dict) -> str:
        """Extract region sequences when needed."""
        self.seq = genome[self.chrom][self.start - 1:self.end]  # TODO: starts from 1?
        return self.seq

    def __repr__(self):
        return f"fRegion: {self.region_id}"

    def __len__(self):
        return len(self.seq)  # Number of bps


class fTRM:
    """Fish Transcriptional regulatory module (TRM)
        A TRM contains TF, regions and target genes on a single chromosome."""

    def __init__(self, tf, region_id, targets=None, celltype=None, species=None):
        self.tf = tf  # tf name

        self.dic = {}

        self.region = fRegion(region_id)  # a region id, only one
        self.chrom = self.region.chrom
        if targets is None:
            targets = []
        self.targets = targets
        # optional
        self.celltype = celltype
        self.species = species

    def __repr__(self):
        return f"fTRM(TF: {self.tf}, {self.region}, Targets: {self.targets})"


class fRegulons:
    """
    A dictionary of Fish Transcriptional Regulatory Modules
    """

    def __init__(self):
        self.regulons = {}  # Dictionary to store lists of trms

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
def read_fasta(fasta_file) -> dict:
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

