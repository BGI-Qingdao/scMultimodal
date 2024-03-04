#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 23 Feb 2024 10:06
# @Author: Yao LI
# @File: evo_fish/run.py
import os
import pickle
import scanpy as sc
import os
import argparse
import textwrap

from handle_atac import ScMultiParams, RefGenome, scATAC
from scmulti import ScMulti

descrip = '''
   _____ ________  _____  ____  __________
  / ___// ____/  |/  / / / / / /_  __/  _/
  \__ \/ /   / /|_/ / / / / /   / /  / /  
 ___/ / /___/ /  / / /_/ / /___/ / _/ /   
/____/\____/_/  /_/\____/_____/_/ /___/   

A tool                    
'''

epilog = '''
--------------------------------BGI Qingdao
'''


def create_parser():
    parser = argparse.ArgumentParser(
        prog='scmulti',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(descrip),
        add_help=True,
        epilog=epilog,
    )
    parser.add_argument('-o', '--output', help='output directory to save all the intermediate and final results.')
    # parser.add_argument('-d', '--data', help='directory where BGI ATAC QC output stored')
    parser.add_argument('--rna', help='scRNA-seq data h5ad file')
    parser.add_argument('--rds', help='scRNA-seq data RDS file')
    parser.add_argument('--atac', help='scATAC-seq data h5ad file')
    parser.add_argument('-n', '--name', help='species name/project name')
    parser.add_argument('--fragment', help='')
    parser.add_argument('--meta', help='')
    parser.add_argument('--gtf', help='')
    parser.add_argument('-f', '--fasta', help='')
    parser.add_argument('--chromsize', help='')
    parser.add_argument('--macs_path', help='')
    parser.add_argument('-t', '--tf', help='')

    parser.add_argument('--convert_strand', default=True, help='')
    parser.add_argument('--chrom_header', default=('NC', 'NW'), help='')
    parser.add_argument('--c_value', default=None, help='')
    parser.add_argument('--genome_size', default=None, help='')
    parser.add_argument('--upstream', default=[1000, 150000], nargs=2, type=int, action='store',
                        metavar=('up_bp', 'down_bp'), help='')
    parser.add_argument('--downstream', default=[1000, 150000], nargs=2, type=int, action='store',
                        metavar=('up_bp', 'down_bp'), help='')
    parser.add_argument('--use_gene_boundaries', default=True, help='')

    parser.add_argument('--celltype_name_seps', default=None, help='')
    parser.add_argument('--celltypes', default=None, help='')
    parser.add_argument('--input_format', default='AUTO', help='')
    parser.add_argument('--shift', default=73, help='')
    parser.add_argument('--ext_size', default=146, help='')
    parser.add_argument('--keep_dup', default='all', help='')
    parser.add_argument('--q_value', default=0.05, help='')
    parser.add_argument('--narrow_peaks_paths', default=None, help='')
    parser.add_argument('--peak_half_width', default=250, help='')
    parser.add_argument('--n_cpu', default=24, help='')
    parser.add_argument('--n_topics', default=None, help='')
    parser.add_argument('--n_iter', default=500, help='')
    parser.add_argument('--random_state', default=555, help='')
    parser.add_argument('--alpha', default=50, help='')
    parser.add_argument('--alpha_by_topic', default=True, help='')
    parser.add_argument('--eta', default=0.1, help='')
    parser.add_argument('--eta_by_topic', default=False, help='')
    parser.add_argument('--select_model', default=16, help='')
    parser.add_argument('--return_model', default=True, help='')
    parser.add_argument('--metrics', default=None, help='')
    parser.add_argument('--plot_metrics', default=False, help='')
    parser.add_argument('--min_topics_coh', default=5, help='')
    parser.add_argument('--split_pattern', default=None, help='')
    parser.add_argument('--db_path', default=None, help='')
    parser.add_argument('--species', default='custom', choices=['custom'], help='')
    parser.add_argument('--meta_cell_split', default=' ', help='')
    parser.add_argument('--multi_ome_mode', default=False, help='')

    # subparsers = parser.add_subparsers(help="sub-command help")
    # parser_atac = subparsers.add_parser(
    #     "atac", help="Handle scATAC-seq data."
    # )
    # parser_atac.add_argument(
    #     "expression_mtx_fname",
    #     type=str,
    #     help="The name of the file that contains the expression matrix for the single cell experiment."
    #          " Two file formats are supported: csv (rows=cells x columns=genes) or loom (rows=genes x columns=cells).",
    # )
    # add_receptor_parameters(parser_atac)

    return parser


def main(argv=None):
    # Parse arguments.
    parser = create_parser()
    args = parser.parse_args(args=argv)
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        args.func(args)

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    # 0.
    scparams = ScMultiParams(work_dir=args.output, output_dir=args.output)
    # 1. Reference Genome
    refgenome = RefGenome(work_dir=args.output, output_dir=args.output, fasta_fn=args.fasta,
                          chromsize_fn=args.chromsize)
    refgenome.get_chromsize()
    # 2. Integration
    proj = ScMulti(refgenome,
                   args.rna,
                   args.fragment,
                   scparams,
                   atac_fn=args.atac,
                   sample_id=args.name,
                   macs_path=args.macs_path,
                   gff_fn=args.gtf,
                   n_cpu=args.n_cpu)
    if args.tf == 'all':
        all_genes = list(proj.atac_data.var_names)
        with open(os.path.join(args.output, 'cistarget_database/all_genes.txt'), 'w') as f:
            f.writelines('\n'.join(all_genes))

    try:
        proj.create_atac_cell_data()
        proj.create_chromesize()
        proj.generate_pseudo_peaks()
        proj.custom_call_peaks(convert_strand=args.convert_strand)
        proj.custom_consensus_peaks(
            narrow_peaks_paths=os.path.join(args.output, 'consensus_peak_calling/narrow_peaks'))
        proj.create_cistopic_obj()
        proj.get_topics()
        proj.visualize_topics()
        proj.candidate_enhancer_regions()
        proj.make_custom_databases(db_path=os.path.join(args.output, 'cistarget_database'), prefix=args.name)
        proj.enrichment(chrom_header=args.chrom_header)
        proj = pickle.load(open(os.path.join(args.output, 'proj.pkl'), 'rb'))
        proj.get_genome_size(c_value=args.c_value)
        proj.network(tf_file=args.tf, upstream=args.upstream,
                     downstream=args.downstream, use_gene_boundaries=args.use_gene_boundaries)
    except Exception as e:
        raise (e)
    finally:
        pickle.dump(proj, open(os.path.join(args.output, 'proj.pkl'), 'wb'))


if __name__ == '__main__':
    # # basic params
    # work_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma'
    # output_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/derived'
    # # if not os.path.exists(work_dir):
    # #     os.makedirs(work_dir)
    # scparams = ScMultiParams(work_dir=work_dir, output_dir=output_dir)
    #
    # # exp1
    # data_dir = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/00.QC/DP8480004851TR_L01_2/output/'
    # frags = 'DP8480004851TR_L01_2.fragments.tsv.gz'
    # meta = 'DP8480004851TR_L01_2.Metadata.tsv'
    # gtf_fn = '/dellfsqd1/ST_OCEAN/ST_OCEAN/USRS/liuqun/project/SingleCell/06.Atlas/01.Oryzias_melastigma/01.genome/GCF_002922805.2_ASM292280v2_genomic.gtf'
    # scrna_rds = '/dellfsqd1/ST_OCEAN/C_OCEAN/USERS/c-tongbaichuan/fish_stomach/04.annotation_v2/01.Oryzias_melastigma/02.intestine/01.Oryzias_melastigma_intestine.anno.rds'
    # project_name = 'DP8480004851TR_L01_2'
    #
    # # Create a scATAC object
    # atacobj = scATAC(atac_fn=os.path.join(data_dir, 'Peak_matrix'),
    #                  scrna_fn=scrna_rds,
    #                  gff=gtf_fn,
    #                  fragment_dir=os.path.join(data_dir, frags),
    #                  fragment_metadata_fn=os.path.join(data_dir, meta),
    #                  project=project_name)
    # # Annotate scATAC data via label transferring
    # atacobj.transfer_label()
    # atacobj.rds2h5ad(os.path.join(output_dir, f'{project_name}.merged.atac.rds'))
    #
    # # Reference Genome
    # fasta_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/DATA/03.Oryzias_melastigma/GCF_002922805.2_ASM292280v2_genomic.fna'
    # chrom_size_fn = '/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/sizes.genome'
    # refgenome = RefGenome(work_dir=work_dir, output_dir=output_dir, fasta_fn=fasta_fn, chromsize_fn=chrom_size_fn)
    # refgenome.get_chromsize()
    #
    # # Multimodal data integration
    # scdata_h5ad = os.path.join(output_dir, f'{project_name}.merged.atac.h5ad')
    # scrna_h5ad = '/dellfsqd2/ST_OCEAN/USER/guolidong/FishGut_V1/10.AnnotationH5ad/01.Oryzias_melastigma/01.Oryzias_melastigma_intestine.anno.h5ad'
    #
    # proj = ScMulti(refgenome,
    #                scrna_h5ad,
    #                os.path.join(data_dir, frags),
    #                scparams,
    #                atac_fn=scdata_h5ad,
    #                sample_id=project_name,
    #                macs_path='macs2',
    #                gff_fn=gtf_fn,
    #                n_cpu = 24)
    # # all_genes = list(proj.atac_data.var_names)
    # # with open(os.path.join(output_dir, 'cistarget_database/all_genes.txt'), 'w') as f:
    # #     f.writelines('\n'.join(all_genes))
    # # proj.create_atac_cell_data()
    # # proj.create_chromesize()
    # # proj.generate_pseudo_peaks()
    # # print('call pesudo peaks')
    # # proj.custom_call_peaks(convert_strand=True)
    # # proj.custom_consensus_peaks(
    # #     narrow_peaks_paths='/dellfsqd2/ST_OCEAN/USER/liyao1/11.evo_fish/exp/03.Oryzias_melastigma/derived/consensus_peak_calling/narrow_peaks')
    # # print('get topics')
    # # proj.create_cistopic_obj()
    # # proj.get_topics()
    # # # proj.visualize_topics()
    # # proj.candidate_enhancer_regions()
    # # proj.make_custom_databases(db_path=os.path.join(output_dir, 'cistarget_database'), prefix='OM')
    # # proj.enrichment(chrom_header=('NC', 'NW'))
    # try:
    #     proj = pickle.load(open(os.path.join(work_dir, 'proj.pkl'), 'rb'))
    #     proj.get_genome_size(c_value=0.87)
    #     proj.network(tf_file=os.path.join(output_dir, 'cistarget_database/all_genes.txt'), upstream=[0, 15000], downstream=[10, 15000], use_gene_boundaries=True)
    # except Exception as e:
    #     raise (e)
    # finally:
    #     pickle.dump(proj, open(os.path.join(work_dir, 'proj.pkl'), 'wb'))

    main()
