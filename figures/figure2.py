import cooler
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import fanc
import fanc.plotting as fancplot
import cooltools
from cooltools.api.insulation import insulation
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import os

def load_and_process_cool(cool_path, chrom, start_pos, end_pos, resolution):
    """Load region from .cool file using cooler"""
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=False).fetch(f'{chrom}:{start_pos}-{end_pos}')
    return np.array(matrix)

def calculate_insulation(cool_path, chrom, start_pos, end_pos, window=250000):
#def calculate_insulation(cool_path, chrom, start_pos, end_pos, window=300000):
    """Calculate insulation score for a region"""
    clr = cooler.Cooler(cool_path)
    insulation_table = insulation(clr, [window], ignore_diags=2, append_raw_scores=True)    
    
    mask = (insulation_table['chrom'] == chrom) & \
           (insulation_table['start'] >= start_pos) & \
           (insulation_table['end'] <= end_pos)
    
    return insulation_table[mask]

def plot_comparison(ax, matrix_reg, matrix_pred, start_pos, end_pos, resolution, tad_path_enh, tad_path_reg, chrom):
    # Normalize matrices
    matrix_reg = matrix_reg / np.max(matrix_reg)
    matrix_pred = matrix_pred / np.max(matrix_pred)

    vmax = max(np.max(matrix_reg), np.max(matrix_pred))

    # Lower triangle: DeepHiC (blue)
    lower_triangle = np.tril(np.ones_like(matrix_pred), -1)
    matrix_pred_masked = np.ma.masked_where(lower_triangle == 0, matrix_pred)
    im1 = ax.imshow(matrix_pred_masked, cmap='Blues', origin='upper', vmax=vmax)

    # Upper triangle: HiC-Pro (red)
    upper_triangle = np.triu(np.ones_like(matrix_reg), 1)
    matrix_reg_masked = np.ma.masked_where(upper_triangle == 0, matrix_reg)
    im2 = ax.imshow(matrix_reg_masked, cmap='Reds', origin='upper', vmax=vmax)

    # Add TAD visualization for enhanced (lower triangle - blue)
    if os.path.exists(tad_path_enh):
        with open(tad_path_enh, 'r') as tad_file:
            for line in tad_file:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                chrom_tad, start, end = parts[0], int(parts[1]), int(parts[2])

                if chrom_tad == chrom and start < end_pos and end > start_pos:
                    plot_start = max(0, (start - start_pos) / resolution)
                    plot_end = min((end - start_pos) / resolution, (end_pos - start_pos) / resolution)

                    # Lower triangle TAD
                    verts = [
                        (plot_start, plot_start),
                        (plot_start, plot_end),
                        (plot_end, plot_end),
                        (plot_start, plot_start)
                    ]
                    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                    path = Path(verts, codes)
                    patch = PathPatch(path, facecolor='none', edgecolor='blue', 
                                    alpha=0.7, lw=1.5, ls='-')
                    ax.add_patch(patch)

    # Add TAD visualization for regular (upper triangle - red)
    if os.path.exists(tad_path_reg):
        with open(tad_path_reg, 'r') as tad_file:
            for line in tad_file:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                chrom_tad, start, end = parts[0], int(parts[1]), int(parts[2])

                if chrom_tad == chrom and start < end_pos and end > start_pos:
                    plot_start = max(0, (start - start_pos) / resolution)
                    plot_end = min((end - start_pos) / resolution, (end_pos - start_pos) / resolution)

                    # Upper triangle TAD
#                    verts = [
#                        (plot_end, plot_start),
#                        (plot_end, plot_end),
#                        (plot_start, plot_end),
#                        (plot_end, plot_start)
#                    ]
                    verts = [
                        (plot_start, plot_start),
                        (plot_end, plot_end),
                        (plot_end, plot_end),
                        (plot_start, plot_start)
                    ]
                    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                    path = Path(verts, codes)
                    patch = PathPatch(path, facecolor='none', edgecolor='red',
                                    alpha=0.7, lw=1.5, ls='-')
                    ax.add_patch(patch)

    # Add colorbars
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes('left', size='5%', pad=0.05)
    cax2 = divider.append_axes('right', size='5%', pad=0.05)

    cbar1 = plt.colorbar(im1, cax=cax1)
    cbar2 = plt.colorbar(im2, cax=cax2)

    cax1.yaxis.set_label_position('left')
    cbar1.ax.yaxis.set_ticks_position('left')

    ax.set_xticks([])
    ax.set_yticks([])

    return im1, im2
def plot_insulation(ax, insulation_enh, start_pos, end_pos,window=250000):
#def plot_insulation(ax, insulation_enh, start_pos, end_pos,window=300000):
    """Plot insulation score for DeepHiC only"""
    pos_enh = (insulation_enh['start'] + insulation_enh['end']) / 2
#    score_enh = insulation_enh['log2_insulation_score_100000']

    score_col = f'log2_insulation_score_{window}'
    score_enh = insulation_enh[score_col]

    
    ax.plot(pos_enh, score_enh, color='blue', linewidth=1.5, label='DeepHiC')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax.set_xlim(start_pos, end_pos)
    ax.set_ylabel('Insulation', fontsize=8)
    ax.legend(loc='upper right', fontsize=6)
    ax.tick_params(labelsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# Parameters
SAMP = 'Kidney'
CHR = 'chr1'
Chr = '1'
RESOLUTION = 100000
COOL_BASE = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Hi-C/Mouse/cool_files_all/100kb'

SPANS = [
    (10000000, 10500000, "5kb"),
    (14500000, 15500000, "1Mb")
]

SPANS = [
    (10000000, 10500000, "5kb"),
    (10000000, 20000000, "1Mb")
]

#SPANS = [
#    (124500000, 125500000, "5kb"),
#    (14500000, 15500000, "1Mb")
#]

# Paths to .cool files
cool_reg = f'{COOL_BASE}/{SAMP}/{SAMP}_hicpro.100kb.cool'
cool_enh = f'{COOL_BASE}/{SAMP}/{SAMP}_deephic.100kb.cool'

# TAD paths
tad_reg = f'{COOL_BASE}/{SAMP}/findTADs/{SAMP}_hicpro_tads_domains.bed'
tad_enh = f'{COOL_BASE}/{SAMP}/findTADs/{SAMP}_deephic_tads_domains.bed'

# Create figure
fig = plt.figure(figsize=(20, 25))
gs = fig.add_gridspec(2, 1, hspace=0.4)

for i, (start_pos, end_pos, span_label) in enumerate(SPANS):
    # Load matrices
    matrix_reg = load_and_process_cool(cool_reg, CHR, start_pos, end_pos, RESOLUTION)
    matrix_pred = load_and_process_cool(cool_enh, CHR, start_pos, end_pos, RESOLUTION)
    
    # Calculate insulation
#    insulation_reg = calculate_insulation(cool_reg, CHR, start_pos, end_pos)
    insulation_enh = calculate_insulation(cool_enh, CHR, start_pos, end_pos, window=300000)

    # Create subplot
    ax = fig.add_subplot(gs[i, 0])
    im1, im2 = plot_comparison(ax, matrix_reg, matrix_pred, start_pos, end_pos, 
                               RESOLUTION, tad_enh, tad_reg, CHR)

    # Add insulation plot
    axins_insulation = ax.inset_axes((0, -0.15, 1, .08))
    plot_insulation(axins_insulation, insulation_enh, start_pos, end_pos, window=300000)

    # Add gene plot
#    bed_path = '/cluster/home/t111631uhn/HiC_ECC/genes.gtf'
    bed_path = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Genomes/gencode.vM25.annotation.gtf'
    bed = fanc.load(bed_path)
    axins_bed = ax.inset_axes((0, -0.25, 1, .08))
    bedplot = fancplot.GenePlot(bed, ax=axins_bed, n_ticks=5)
    bedplot.plot(f'Chr{Chr}:{start_pos}-{end_pos}')

    # Add title
    ax.set_title(f'{span_label} at 100kb resolution', pad=20)

# Add overall title
fig.suptitle(f'DeepHiC Comparison for {SAMP} {CHR}', fontsize=16, y=0.95)

# Save plot
plt.savefig(f'{SAMP}_deephic_comparison_with_TADs_100kb_fixed.pdf', format='pdf', bbox_inches='tight', dpi=500)
plt.show()
