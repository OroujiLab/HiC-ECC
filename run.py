"""
run.py — HiC-ECC CLI Entry Point

Usage:
    python run.py --config config/config.yaml --module all
    python run.py --config config/config.yaml --module enhance
    python run.py --config config/config.yaml --module compare
    python run.py --config config/config.yaml --module compare --tissue-idx 0
    python run.py --config config/config.yaml --module call
"""

import argparse
import yaml
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from enhance import run_enhancement
from compare import run_comparison
from call    import run_calling

MODULES = ['all', 'enhance', 'compare', 'call']


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(
        description='HiC-ECC: Hi-C Enhancement, Comparison, and Calling Pipeline'
    )
    parser.add_argument(
        '--config', required=True,
        help='Path to config YAML file (e.g. config/config.yaml)'
    )
    parser.add_argument(
        '--module', required=True, choices=MODULES,
        help='Module to run. Choices: all, enhance, compare, call'
    )
    parser.add_argument(
        '--sample-idx', type=int, default=None,
        help='For enhance module: index of sample in samples list (0-based). '
             'If set, only processes this sample. '
             'Use this to parallelize across SLURM jobs.'
    )
    parser.add_argument(
        '--tissue-idx', type=int, default=None,
        help='For compare module: index of reference tissue in samples list (0-based). '
             'If set, only compares this tissue against all others. '
             'Use this to parallelize across SLURM jobs.'
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    print(f'Loaded config: {args.config}')
    print(f'Samples: {cfg["samples"]}')
    print(f'Module:  {args.module}')
    if args.sample_idx is not None:
        print(f'Sample: {cfg["samples"][args.sample_idx]} (idx={args.sample_idx})')
    if args.tissue_idx is not None:
        print(f'Reference tissue: {cfg["samples"][args.tissue_idx]} (idx={args.tissue_idx})')
    print()

    if args.module in ('all', 'enhance'):
        print('=' * 50)
        print('STEP 1: Enhancement')
        print('=' * 50)
        run_enhancement(cfg, sample_idx=args.sample_idx)

    if args.module in ('all', 'compare'):
        print('=' * 50)
        print('STEP 2: Comparison')
        print('=' * 50)
        run_comparison(cfg, tissue_idx=args.tissue_idx)

    if args.module in ('all', 'call'):
        print('=' * 50)
        print('STEP 3: Calling')
        print('=' * 50)
        run_calling(cfg)

    print('\nHiC-ECC pipeline complete.')


if __name__ == '__main__':
    main()
