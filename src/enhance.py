"""
enhance.py — HiC-ECC Enhancement Module
Converts HiC-Pro matrices to DeepHiC-enhanced .cool files.
"""

import os
import re
import subprocess
import numpy as np


def patch_deephic_root(deephic_py, deephic_root):
    """Patch all_parser.py root_dir to point to the correct output directory."""
    parser_path = os.path.join(deephic_py, 'all_parser.py')
    with open(parser_path, 'r') as f:
        content = f.read()
    content = re.sub(r"root_dir = '.*'", f"root_dir = '{deephic_root}/'", content)
    with open(parser_path, 'w') as f:
        f.write(content)
    print(f'Patched all_parser.py: root_dir = {deephic_root}/')


def run_hicpro2deephic(tissue, sample_dir, resolution, deephic_py, deephic_root):
    """Convert HiC-Pro matrix to DeepHiC input format."""
    res_str = str(resolution)
    tissue_dir = os.path.join(sample_dir, tissue, 'hic_results', 'matrix', tissue, 'raw', res_str)
    subprocess.run([
        'python', os.path.join(deephic_py, 'hicpro2deephic.py'),
        '--bed', os.path.join(tissue_dir, f'{tissue}_{res_str}_abs.bed'),
        '--mat', os.path.join(tissue_dir, f'{tissue}_{res_str}.matrix'),
        '-r', res_str,
        '-o', os.path.join(deephic_root, 'mat', f'{tissue}_{res_str}'),
    ], check=True)


def run_data_generate(tissue, resolution, species, deephic_py, deephic_root, chunk, stride, bound, lrc):
    """Generate DeepHiC input data chunks."""
    res_str = str(resolution)
    subprocess.run([
        'python', os.path.join(deephic_py, 'data_generate.py'),
        '-hr', res_str, '-lr', res_str,
        '-lrc', str(lrc), '-s', species,
        '-chunk', str(chunk), '-stride', str(stride),
        '-bound', str(bound), '-scale', '1',
        '-c', f'{tissue}_{res_str}',
    ], check=True, cwd=deephic_root)


def run_data_predict(tissue, resolution, deephic_py, deephic_root, checkpoint):
    """Run DeepHiC prediction."""
    res_str = str(resolution)
    subprocess.run([
        'python', os.path.join(deephic_py, 'data_predict.py'),
        '-lr', res_str,
        '-ckpt', checkpoint,
        '-c', f'{tissue}_{res_str}',
    ], check=True, cwd=deephic_root)


def convert_to_cool(tissue, resolution, deephic_root, cool_dir, chrom_sizes):
    """Convert DeepHiC predictions (.npz) to .cool format."""
    res_kb = resolution // 1000
    tissue_cool_dir = os.path.join(cool_dir, tissue)
    os.makedirs(tissue_cool_dir, exist_ok=True)

    bg2    = os.path.join(tissue_cool_dir, f'{tissue}_deephic_{res_kb}kb.bg2')
    cool   = os.path.join(tissue_cool_dir, f'{tissue}_deephic.{res_kb}kb.cool')
    sr_dir = os.path.join(deephic_root, 'predict', f'{tissue}_{resolution}', 'sr16')

    # npz → bg2
    with open(bg2, 'w') as out:
        for fname in sorted(os.listdir(sr_dir)):
            if not fname.endswith(f'_{resolution}.npz'):
                continue
            chrom  = fname.replace('predict_', '').replace(f'_{resolution}.npz', '')
            matrix = np.load(os.path.join(sr_dir, fname))['deephic']
            for i in range(matrix.shape[0]):
                for j in range(i, matrix.shape[0]):
                    if matrix[i, j] > 0:
                        out.write(
                            f'{chrom}\t{i*resolution}\t{(i+1)*resolution}\t'
                            f'{chrom}\t{j*resolution}\t{(j+1)*resolution}\t{matrix[i,j]}\n'
                        )

    # bg2 → cool
    subprocess.run([
        'cooler', 'load', '-f', 'bg2',
        f'{chrom_sizes}:{resolution}',
        bg2, cool,
        '--count-as-float', '--input-copy-status', 'duplex',
    ], check=True)

    # balance
    subprocess.run(['cooler', 'balance', cool, '--force', '--max-iters', '500'], check=True)
    print(f'[{tissue}] done → {cool}')


def run_enhancement(cfg):
    """Main entry point. Runs full enhancement pipeline for all samples."""
    tissues     = cfg['samples']
    genome      = cfg['genome']
    species     = cfg['species']
    resolution  = cfg['resolution']
    chrom_sizes = cfg['chrom_sizes']
    sample_dir  = cfg['hicpro_dir']
    output_dir  = cfg['output_dir']

    enh         = cfg['enhancement']
    deephic_py  = enh['deephic_path']
    checkpoint  = enh['checkpoint']
    chunk       = enh['chunk']
    stride      = enh['stride']
    bound       = enh['bound']
    lrc         = enh['lrc']

    deephic_root = os.path.join(output_dir, 'DeepHiC', genome)
    cool_dir     = os.path.join(output_dir, 'cool_files', str(resolution))

    os.makedirs(os.path.join(deephic_root, 'mat'), exist_ok=True)
    patch_deephic_root(deephic_py, deephic_root)

    # Section 1: enhance
    for tissue in tissues:
        predict_dir = os.path.join(deephic_root, 'predict', f'{tissue}_{resolution}', 'sr16')
        if os.path.isdir(predict_dir) and len(os.listdir(predict_dir)) > 0:
            print(f'[{tissue}] prediction exists, skipping')
        else:
            print(f'\n[{tissue}] hicpro2deephic')
            run_hicpro2deephic(tissue, sample_dir, resolution, deephic_py, deephic_root)
            print(f'[{tissue}] data_generate')
            run_data_generate(tissue, resolution, species, deephic_py, deephic_root, chunk, stride, bound, lrc)
            print(f'[{tissue}] data_predict')
            run_data_predict(tissue, resolution, deephic_py, deephic_root, checkpoint)

    # Section 2: convert to cool
    for tissue in tissues:
        res_kb = resolution // 1000
        cool   = os.path.join(cool_dir, tissue, f'{tissue}_deephic.{res_kb}kb.cool')
        if os.path.isfile(cool) and os.path.getsize(cool) > 0:
            print(f'[{tissue}] cool exists, skipping')
        else:
            print(f'\n[{tissue}] converting to cool')
            convert_to_cool(tissue, resolution, deephic_root, cool_dir, chrom_sizes)

    print('\nEnhancement complete.')
