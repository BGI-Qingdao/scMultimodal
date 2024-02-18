#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 30 Jan 2024 08:56
# @Author: Yao LI
# @File: evo_fish/make_TF_list.py
import pandas as pd

gene_ortho = pd.read_csv('Lethenteron_mm_gene_name.txt')
with open('/dellfsqd2/ST_OCEAN/USER/liyao1/06.stereopy/resource/tfs/allTFs_mm.txt', 'r') as f:
    mm_tfs = f.read().splitlines()

target_tfs = list(gene_ortho[gene_ortho.gene_name.isin(mm_tfs)].lamprey)
with open('/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/01.Lethenteron/00.custom/TFs_Lethenteron.txt', 'w') as f:
    f.writelines('\n'.join(target_tfs))
