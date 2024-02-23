#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date: Created on 19 Feb 2024 16:33
# @Author: Yao LI
# @File: evo_fish/pipeline.py
import argparse
import textwrap

descrip = '''
   _____ ________  _____  ____  __________
  / ___// ____/  |/  / / / / / /_  __/  _/
  \__ \/ /   / /|_/ / / / / /   / /  / /  
 ___/ / /___/ /  / / /_/ / /___/ / _/ /   
/____/\____/_/  /_/\____/_____/_/ /___/   
                    
Hello!                      
'''

epilog = '''
--------------------------------BGI Qingdao
'''


def create_argument_parser():
    parser = argparse.ArgumentParser(
        prog='SCMULTI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(descrip),
        epilog=epilog,
    )
    parser.add_argument('--outPath', type=str, help='''input the outpath''',)
    return parser


def main(argv=None):
    # Parse arguments.
    parser = create_argument_parser()
    args = parser.parse_args(args=argv)
    if not hasattr(args, "func"):
        parser.print_help()
    else:
        args.func(args)


if __name__ == "__main__":
    main()
