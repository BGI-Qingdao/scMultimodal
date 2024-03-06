#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 23 Feb 2024 10:06
# @Author: Yao LI
# @File: evo_fish/run.py
import sys
import os
import pickle
import scanpy as sc
import os
import argparse
import textwrap

from handle_atac import ScMultiParams, RefGenome, scATAC
from scmulti import ScMulti
from custom_databases import create_custom_database


class TupleType:
    def __init__(self, sep=','):
        self.sep = sep

    def __call__(self, value):
        return tuple(value.split(self.sep))


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
    parser.add_argument('--rna', help='scRNA-seq data h5ad file')
    parser.add_argument('--rds', help='scRNA-seq data RDS file')
    parser.add_argument('--atac', help='scATAC-seq data h5ad file')
    parser.add_argument('-n', '--name', help='species name/project name')
    parser.add_argument('--fragment', help='')
    parser.add_argument('--meta', help='')
    parser.add_argument('--gtf', help='')
    parser.add_argument('-f', '--fasta', help='')
    parser.add_argument('--chromsize', help='')
    parser.add_argument('--macs_path', default='macs2', help='')
    parser.add_argument('-t', '--tf', default='all', help='')
    parser.add_argument('-v', '--variable', default='celltype', help='')

    parser.add_argument('--make_custom_database', default=True, action='store_true')

    parser.add_argument('--convert_strand', default=True, action='store_true', help='')
    parser.add_argument('--chrom_header', type=TupleType(), default=('NC', 'NW'), help='')
    parser.add_argument('--c_value', default=None, type=float, help='')
    parser.add_argument('--genome_size', default=None, type=int, help='')
    parser.add_argument('--upstream', default=[1000, 150000], nargs=2, type=int, action='store',
                        metavar=('up_bp', 'down_bp'), help='')
    parser.add_argument('--downstream', default=[1000, 150000], nargs=2, type=int, action='store',
                        metavar=('up_bp', 'down_bp'), help='')
    parser.add_argument('--use_gene_boundaries', default=True, action='store_true', help='')

    # parser.add_argument('--celltype_name_seps', default=None, nargs='*', type=str,
    #                     help='Delimiter in cell type names. e.g. Blood-1, sep is -')
    parser.add_argument('--celltypes', default=None, nargs='*', type=str, action='store',
                        help='Interested cell types. When not provided, will use all the cell types in data to compute')
    parser.add_argument('--input_format', default='AUTO', help='')
    parser.add_argument('--shift', default=73, type=int, help='')
    parser.add_argument('--ext_size', default=146, type=int, help='')
    parser.add_argument('--keep_dup', default='all', help='')
    parser.add_argument('--q_value', default=0.05, type=float, help='')
    parser.add_argument('--narrow_peaks_paths', default=None, help='')
    parser.add_argument('--peak_half_width', default=250, type=int, help='')
    parser.add_argument('--n_cpu', default=24, type=int, help='')
    parser.add_argument('--n_topics', default=None, help='')
    parser.add_argument('--n_iter', default=500, type=int, help='')
    parser.add_argument('--random_state', default=555, type=int, help='')
    parser.add_argument('--alpha', default=50, type=int, help='')
    parser.add_argument('--alpha_by_topic', default=True, action='store_true', help='')
    parser.add_argument('--eta', default=0.1, type=float, help='')
    parser.add_argument('--eta_by_topic', default=False, action='store_false', help='')
    parser.add_argument('--select_model', default=16, type=int, help='')
    parser.add_argument('--return_model', default=True, action='store_true', help='')
    parser.add_argument('--metrics', default=None, help='')
    parser.add_argument('--plot_metrics', default=False, action='store_false', help='')
    parser.add_argument('--min_topics_coh', default=5, type=int, help='')
    parser.add_argument('--split_pattern', default=None, help='')
    parser.add_argument('--db_path', default=None, help='')
    parser.add_argument('--species', default='custom', choices=['custom'], help='')
    parser.add_argument('--meta_cell_split', default=' ', help='')
    parser.add_argument('--multi_ome_mode', default=False, action='store_false', help='')

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
    # subcommand_parser.set_defaults(func=subcommand_function)

    return parser


def main():
    # Parse arguments.
    parser = create_parser()
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if not os.path.exists(args.output):
        os.makedirs(args.output)
    # 0.
    # scparams = ScMultiParams(work_dir=args.output, output_dir=args.output)
    # # 1. Reference Genome
    # refgenome = RefGenome(work_dir=args.output, output_dir=args.output, fasta_fn=args.fasta,
    #                       chromsize_fn=args.chromsize)
    # refgenome.get_chromsize()
    # 2. Integration
    # proj = ScMulti(refgenome,
    #                args.rna,
    #                args.fragment,
    #                scparams,
    #                atac_fn=args.atac,
    #                sample_id=args.name,
    #                macs_path=args.macs_path,
    #                gff_fn=args.gtf,
    #                n_cpu=args.n_cpu,
    #                c_value=args.c_value,
    #                variable=args.variable)
    # if args.tf == 'all':
    #     all_genes = list(proj.atac_data.var_names)
    #     with open(os.path.join(args.output, 'all_genes.txt'), 'w') as f:
    #         f.writelines('\n'.join(all_genes))
    #     tf_file = os.path.join(args.output, 'all_genes.txt')
    # else:
    #     tf_file = args.tf
    tf_file = os.path.join(args.output, 'all_genes.txt')
    proj = pickle.load(open(os.path.join(args.output, 'proj.pkl'), 'rb'))
    try:
        # proj.create_atac_cell_data()  # celltype_name_seps=args.celltype_name_seps
        # proj.create_chromesize()
        # proj.get_genome_size()
        #
        # proj.generate_pseudo_peaks(celltypes=args.celltypes)
        # proj.custom_call_peaks(convert_strand=args.convert_strand,
        #                        input_format=args.input_format,
        #                        shift=args.shift,
        #                        ext_size=args.ext_size,
        #                        keep_dup=args.keep_dup,
        #                        q_value=args.q_value)
        # if args.narrow_peaks_paths is None:
        #     narrow_peaks_paths = os.path.join(args.output, 'consensus_peak_calling/narrow_peaks')
        # else:
        #     narrow_peaks_paths = args.narrow_peaks_paths
        # proj.custom_consensus_peaks(
        #     narrow_peaks_paths=narrow_peaks_paths,
        #     peak_half_width=args.peak_half_width)
        # proj.create_cistopic_obj()
        # proj.get_topics(n_topics=args.n_topics,
        #                 n_iter=args.n_iter,
        #                 random_state=args.random_state,
        #                 alpha=args.alpha,
        #                 alpha_by_topic=args.alpha_by_topic,
        #                 eta=args.eta,
        #                 eta_by_topic=args.eta_by_topic,
        #                 select_model=args.select_model,
        #                 return_model=args.return_model,
        #                 metrics=args.metrics,
        #                 plot_metrics=args.plot_metrics)
        # proj.visualize_topics()
        proj.candidate_enhancer_regions()
        # if args.make_custom_database:
        #     create_custom_database(args.name, args.fasta, args.output, args.consensus_regions, args.atac, args.ortholog,
        #                    ref_species=args.ref_species, bed=args.bedpath)
        proj.parse_custom_databases_name(db_path=os.path.join(args.output, 'cistarget_database'), prefix=args.name)
        proj.enrichment(chrom_header=args.chrom_header)
        proj.network(tf_file=tf_file, upstream=args.upstream,
                     downstream=args.downstream, use_gene_boundaries=args.use_gene_boundaries)
    except Exception as e:
        raise (e)
    finally:
        pickle.dump(proj, open(os.path.join(args.output, 'proj.pkl'), 'wb'))


if __name__ == '__main__':
    main()
