"""
call.py — HiC-ECC Call Module
Loop and TAD calling on enhanced .cool files.
"""

import os
import subprocess


def run_loop_calling(tissue, cool, out_dir, threads):
    """Run Chromosight loop detection."""
    out = os.path.join(out_dir, f'{tissue}_loops')
    print(f'[loops] {tissue}')
    subprocess.run(
        ['chromosight', 'detect', '-t', str(threads), cool, out],
        check=True
    )


def run_tad_calling(tissue, cool, out_dir):
    """Run hicFindTADs TAD calling."""
    prefix = os.path.join(out_dir, f'{tissue}_tads')
    print(f'[TADs] {tissue}')
    subprocess.run(
        ['hicFindTADs', '-m', cool,
         '--outPrefix', prefix,
         '--correctForMultipleTesting', 'fdr'],
        check=True
    )


def run_calling(cfg):
    """Main entry point. Runs loop and TAD calling for all samples."""
    tissues    = cfg['samples']
    resolution = cfg['resolution']
    output_dir = cfg['output_dir']

    cool_dir    = os.path.join(output_dir, 'cool_files', str(resolution))
    cool_suffix = f"_deephic.{resolution//1000}kb.cool"
    out_dir     = os.path.join(output_dir, 'call_output')
    threads     = cfg['calling']['loops']['threads']

    os.makedirs(out_dir, exist_ok=True)

    for tissue in tissues:
        cool = os.path.join(cool_dir, tissue, tissue + cool_suffix)

        if not os.path.isfile(cool):
            print(f'[SKIP] Missing: {cool}')
            continue

        run_loop_calling(tissue, cool, out_dir, threads)
        run_tad_calling(tissue, cool, out_dir)
        print(f'[{tissue}] done')

    print('\nCalling complete.')
