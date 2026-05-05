"""
run.py — HiC-ECC CLI Entry Point

Usage:
    python run.py --config config/config.yaml --steps enhance compare call
    python run.py --config config/config.yaml --steps enhance
    python run.py --config config/config.yaml  # runs all steps
"""

import argparse
import yaml
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from enhance import run_enhancement
from compare import run_comparison
from call    import run_calling

STEPS = ['enhance', 'compare', 'call']


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
        '--steps', nargs='+', choices=STEPS, default=STEPS,
        help=f'Steps to run (default: all). Choices: {STEPS}'
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    print(f'Loaded config: {args.config}')
    print(f'Samples: {cfg["samples"]}')
    print(f'Steps:   {args.steps}\n')

    if 'enhance' in args.steps:
        print('=' * 50)
        print('STEP 1: Enhancement')
        print('=' * 50)
        run_enhancement(cfg)

    if 'compare' in args.steps:
        print('=' * 50)
        print('STEP 2: Comparison')
        print('=' * 50)
        run_comparison(cfg)

    if 'call' in args.steps:
        print('=' * 50)
        print('STEP 3: Calling')
        print('=' * 50)
        run_calling(cfg)

    print('\nHiC-ECC pipeline complete.')


if __name__ == '__main__':
    main()
