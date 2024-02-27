#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 24 Jan 2024 14:11
# @Author: Yao LI
# @File: evo_fish/custom_annotation.py
import pandas as pd


def get_custom_annot(gtf_fn: str, out_fn: str) -> pd.DataFrame:
    """
    :param gtf_fn:
    :param out_fn:
    :return:
    Object: to create a dataframe like this:
    pandas DataFrame with genome annotation for custom species (i.e. for a species other than homo_sapiens, mus_musculus,
    drosophila_melanogaster or gallus_gallus).
    This DataFrame should (minimally) look like the example below, and only contains protein coding genes:
    >>> custom_annot
            Chromosome      Start  Strand     Gene Transcript_type
        8053         chrY   22490397       1      PRY  protein_coding
        8153         chrY   12662368       1    USP9Y  protein_coding
        8155         chrY   12701231       1    USP9Y  protein_coding
        8158         chrY   12847045       1    USP9Y  protein_coding
        8328         chrY   22096007      -1     PRY2  protein_coding
        ...           ...        ...     ...      ...             ...
        246958       chr1  181483738       1  CACNA1E  protein_coding
        246960       chr1  181732466       1  CACNA1E  protein_coding
        246962       chr1  181776101       1  CACNA1E  protein_coding
        246963       chr1  181793668       1  CACNA1E  protein_coding
        246965       chr1  203305519       1     BTG2  protein_coding
        [78812 rows x 5 columns]
    """
    df = pd.read_csv(gtf_fn, header=None, sep='\t')
    # extract gene names
    df[9] = df[8].str.split(';', expand=True)[1]
    df[9] = df[9].str.split(' ', expand=True)[2]
    df[9] = df[9].str.replace('"', '')
    # subset only essential columns
    new = df[[0, 3, 4, 6, 9, 2]]
    # select only protein_coding genes/entries
    new = new[new[2] == 'exon']
    new.columns = ['Chromosome', 'Start', 'End', 'Strand', 'Gene', 'Transcript_type']
    new['Transcript_type'] = new['Transcript_type'].str.replace('exon', 'protein_coding')
    # convert Strand anno from 1/-1 to +/-
    new['Strand'] = new['Strand'].str.replace('+', '1')
    new['Strand'] = new['Strand'].str.replace('-', '-1')
    # save custom annotation
    print('Custom Annotation:')
    print(new)
    new.to_csv(out_fn, sep='\t', index=False)
    return new


if __name__ == '__main__':
    gtf_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/01.Lethenteron/data/LRfinal.gtf'
    out_fn = 'lethenteron_annotation.gtf'
    get_custom_annot(gtf_fn, out_fn)
