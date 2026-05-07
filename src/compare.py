"""
compare.py — HiC-ECC Comparison Module
Pairwise Hi-C comparison with CHESS + unique region extraction.
"""

import os
import subprocess
import itertools
import pandas as pd


def generate_bed_files(genome, window, resolution, chromosomes, bed_dir):
    """Generate CHESS region pair BED files for each chromosome."""
    os.makedirs(bed_dir, exist_ok=True)
    for chrom in chromosomes:
        bed = os.path.join(bed_dir, f'{genome}_{chrom}_{window}_win_{resolution}_step.bed')
        if os.path.isfile(bed):
            print(f'[skip] BED exists: {chrom}')
            continue
        print(f'[pairs] {chrom}')
        subprocess.run(
            ['chess', 'pairs', genome, str(window), str(resolution), bed, '--chromosome', chrom],
            check=True
        )


def run_chess_sim(tissues, cool_dir, cool_suffix, chromosomes, genome, window,
                  resolution, out_dir, bed_dir, threads):
    """Run pairwise CHESS sim for all unique tissue pairs."""
    pairs = list(itertools.combinations(tissues, 2))
    print(f'{len(pairs)} pairs × {len(chromosomes)} chromosomes = {len(pairs)*len(chromosomes)} comparisons')

    for ref, comp in pairs:
        ref_cool  = os.path.join(cool_dir, ref,  ref  + cool_suffix)
        comp_cool = os.path.join(cool_dir, comp, comp + cool_suffix)

        for cool, label in [(ref_cool, ref), (comp_cool, comp)]:
            if not os.path.isfile(cool):
                print(f'[SKIP] Missing: {cool}')
                continue

        pair_dir = os.path.join(out_dir, ref)
        os.makedirs(pair_dir, exist_ok=True)

        for chrom in chromosomes:
            bed = os.path.join(bed_dir, f'{genome}_{chrom}_{window}_win_{resolution}_step.bed')
            out = os.path.join(pair_dir, f'{chrom}_{ref}_vs_{comp}{cool_suffix.replace(".cool", "")}.tsv')
            print(f'[sim] {ref} vs {comp} | {chrom}')
            subprocess.run(
                ['chess', 'sim', ref_cool, comp_cool, bed, out, '-p', str(threads)],
                check=True
            )


def extract_unique_regions(tissues, chromosomes, genome, window, resolution,
                           out_dir, bed_dir, uniq_dir, sn_thr, zsim_thr):
    """Extract tissue-unique regions from CHESS output for each chromosome."""
    os.makedirs(uniq_dir, exist_ok=True)
    win_label = f'{window // 1000}kb'

    for chrom in chromosomes:
        bed_file = os.path.join(bed_dir, f'{genome}_{chrom}_{window}_win_{resolution}_step.bed')
        if not os.path.isfile(bed_file):
            print(f'[skip] no BED for {chrom}')
            continue
        regions = pd.read_csv(bed_file, sep='\t', header=None)

        for tissue in tissues:
            samp_sims = []
            for tissue_dir in tissues:
                tsv_dir = os.path.join(out_dir, tissue_dir)
                if not os.path.isdir(tsv_dir):
                    continue
                for f in os.listdir(tsv_dir):
                    if f.endswith('.tsv') and chrom in f:
                        if f'_{tissue}_vs_' in f or f'_vs_{tissue}_' in f:
                            sim  = pd.read_csv(os.path.join(tsv_dir, f), sep='\t', index_col=0)
                            filt = sim[(sim['SN'] >= sn_thr) & (sim['z_ssim'] <= zsim_thr)]
                            samp_sims.append(filt)

            if not samp_sims:
                print(f'[{tissue}] {chrom}: no data')
                continue

            uniq_idx = set(samp_sims[0].index)
            for s in samp_sims[1:]:
                uniq_idx.intersection_update(s.index)

            out_bed = os.path.join(uniq_dir, f'{genome}_{tissue}_{chrom}_{win_label}_{zsim_thr}.bed')
            tmp = regions.iloc[list(uniq_idx), :3].merge(
                samp_sims[0]['ssim'], how='inner', left_index=True, right_index=True
            ).rename(columns={0: 'chr', 1: 'start', 2: 'end', 'ssim': 'ssim_score'})
            tmp.to_csv(out_bed, sep='\t', header=False, index=False)
            print(f'[{tissue}] {chrom}: {len(tmp)} unique regions → {out_bed}')


def run_comparison(cfg, tissue_idx=None):
    """Main entry point. Runs CHESS comparison + unique region extraction.
    
    Args:
        cfg: parsed config dict
        tissue_idx: if set, only run comparisons where this tissue is the reference.
                    Use for SLURM parallelization (one job per tissue).
    """
    tissues     = cfg['samples']
    genome      = cfg['genome']
    resolution  = cfg['resolution']
    chromosomes = cfg['chromosomes']
    output_dir  = cfg['output_dir']

    cmp         = cfg['comparison']
    window      = cmp['window']
    threads     = cmp['threads']
    sn_thr      = cmp['sn_threshold']
    zsim_thr    = cmp['zsim_threshold']

    cool_dir    = os.path.join(output_dir, 'cool_files', str(resolution))
    cool_suffix = f"_deephic.{resolution//1000}kb.cool"
    out_dir     = os.path.join(output_dir, 'chess_output')
    bed_dir     = os.path.join(out_dir, 'bed_files')
    uniq_dir    = os.path.join(out_dir, 'unique_regions')

    # if tissue_idx set, only run for that reference tissue
    if tissue_idx is not None:
        ref_tissue = tissues[tissue_idx]
        print(f'Running in parallel mode: reference tissue = {ref_tissue}')
        tissues_to_run = [ref_tissue]
    else:
        tissues_to_run = tissues

    generate_bed_files(genome, window, resolution, chromosomes, bed_dir)
    all_tissues = tissues if tissue_idx is not None else None
    run_chess_sim(tissues_to_run, cool_dir, cool_suffix, chromosomes, genome,
                  window, resolution, out_dir, bed_dir, threads, all_tissues=all_tissues)
    extract_unique_regions(tissues, chromosomes, genome, window, resolution,
                           out_dir, bed_dir, uniq_dir, sn_thr, zsim_thr)

    print('\nComparison complete.')
