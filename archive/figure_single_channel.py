import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from sklearn.preprocessing import normalize
import os
import fanc
import fanc.plotting as fancplot
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.colors import ListedColormap
import pandas as pd
import cooler

def load_cool_region(cool_path, chrom, start_pos, end_pos):
    """Load region from .cool file"""
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=False).fetch(f'{chrom}:{start_pos}-{end_pos}')
    return np.array(matrix)

def get_unique_regions(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep='\t', header=None)
    df_filt = df[(df[0]==f'chr{CHR}') & (df[1]<END) & (df[2]>START)]
    df_filt.loc[df_filt[1]<START, 1] = START
    df_filt.loc[df_filt[2]>END, 2] = END
    return df_filt[[1,2]].values.tolist()

def get_uniq_matrix(region_start, region_end):
    mat_size = int((END - START) / RESOLUTION)
    matrix = np.zeros((mat_size, mat_size))
    start_idx = int((region_start - START) / RESOLUTION)
    end_idx = int((region_end - START) / RESOLUTION)
    matrix[start_idx:end_idx, start_idx:end_idx] = 1
    return matrix

def mergeIntervals(intervals):
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = [last[0], max(last[1], current[1])]
        else:
            merged.append(current)
    return merged

def get_matrix_sum(*matrices, ones=False):
    result = matrices[0].copy()
    for mat in matrices[1:]:
        result += mat
    if ones:
        result = (result > 0).astype(int)
    return result

def combine_uniq_regions(uniq_reg_list):
    if not uniq_reg_list:
        mat_size = int((END - START) / RESOLUTION)
        return np.zeros((mat_size, mat_size))
    
    if MERGE_INTERVALS:
        intervals = mergeIntervals(uniq_reg_list)
    else:
        intervals = uniq_reg_list
    
    matrices = []
    for region in intervals:
        uniq_mat = get_uniq_matrix(region[0], region[1])
        matrices.append(uniq_mat)
    
    if not matrices:
        mat_size = int((END - START) / RESOLUTION)
        return np.zeros((mat_size, mat_size))
    
    return get_matrix_sum(*matrices, ones=True)

def combine_uniq_regions(uniq_reg_list):
    # Get the actual matrix size from cooler
    clr = cooler.Cooler(f'{COOL_BASE}/{TISSUES[0]}/{TISSUES[0]}_deephic.10kb.cool')
    test_matrix = clr.matrix(balance=False).fetch(f'chr{CHR}:{START}-{END+RESOLUTION}')
    mat_size = test_matrix.shape[0]
    
    if not uniq_reg_list:
        return np.zeros((mat_size, mat_size))
    
    if MERGE_INTERVALS:
        intervals = mergeIntervals(uniq_reg_list)
    else:
        intervals = uniq_reg_list
    
    matrices = []
    for region in intervals:
        uniq_mat = get_uniq_matrix(region[0], region[1])
        matrices.append(uniq_mat)
    
    if not matrices:
        return np.zeros((mat_size, mat_size))
    
    return get_matrix_sum(*matrices, ones=True)






def get_uniq_interaction(cool_path, uniq_mat, chrom, start, end, resolution):
    """Load from cool and apply unique region mask"""
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=False).fetch(f'{chrom}:{start}-{end+resolution}')
    return uniq_mat * np.array(matrix)


def get_uniq_interaction(cool_path, uniq_mat):
    """Load from cool and apply unique region mask"""
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=False).fetch(f'chr{CHR}:{START}-{END+RESOLUTION}')
    matrix = np.array(matrix)
    
    # Resize uniq_mat to match matrix if needed
    if uniq_mat.shape != matrix.shape:
        print(f"Warning: Resizing uniq_mat from {uniq_mat.shape} to {matrix.shape}")
        if uniq_mat.shape[0] > matrix.shape[0]:
            uniq_mat = uniq_mat[:matrix.shape[0], :matrix.shape[1]]
        else:
            new_mat = np.zeros(matrix.shape, dtype=uniq_mat.dtype)
            new_mat[:uniq_mat.shape[0], :uniq_mat.shape[1]] = uniq_mat
            uniq_mat = new_mat
    
    return uniq_mat * matrix


def get_nth_percentile(matrix, percentile):
    non_zero = matrix[matrix > 0]
    if len(non_zero) == 0:
        return matrix
    threshold = np.percentile(non_zero, percentile)
    return np.where(matrix >= threshold, matrix, 0)

def calculate_contact_density(matrix, window_size=5):
    size = matrix.shape[0]
    density = np.zeros(size)
    half_window = window_size // 2

    for i in range(size):
        start = max(0, i - half_window)
        end = min(size, i + half_window + 1)
        density[i] = np.sum(matrix[start:end, start:end])

    if np.max(density) > 0:
        density = density / np.max(density)
    return density

# Configuration
RESOLUTION = 10000
CHR = 17
START = 14500000
END = 15500000

GENE = 'Mtor5'
PERCENTILE = 70
WIN_SZ = '310kb'
MERGE_INTERVALS = True
TISSUES = ['Kidney', 'Liver', 'Large_Intestine', 'Pancreas', 'Small_Intestine', 'Brain', 'Lung', 'Spleen']
size = int((END - START) / RESOLUTION)
extent = [0, size, size, 0]
alpha = 1

# Base paths
COOL_BASE = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Hi-C/Mouse/cool_files_all'
UNIQ_BASE = '/cluster/home/t111631uhn/HiC_ECC/new_scripts/chess_output/unique_regions'

# Colors
colors_left = ["#89288F", "#F47D2B", "#272E6A", "#8A9FD1", "#C06CAB", "#D51F26", "#FEE500", "#F9B712"]
cmaps = [ListedColormap([color]) for color in colors_left]

# Create figure
fig, axes = plt.subplots(4, 2, figsize=(10, 20), gridspec_kw={'wspace': 0.1, 'hspace': 0.5})
fig.suptitle(f'Channel-wise Interactions - Chr{CHR}:{START}-{END}', fontsize=16, y=0.95)

# Plot each tissue
for i, tis in enumerate(TISSUES):
    ax = axes[i // 2, i % 2]
    ax.set_title(tis, fontsize=12)
    ax.plot([0, size], [0, size], color="black", linewidth=1)

    # Get matrix
    uniq_path = f'{UNIQ_BASE}/mm10_{tis}_chr{CHR}_{WIN_SZ}_0.bed'
    uniq = get_unique_regions(uniq_path)
    
    if uniq is not None:
        uniq_mat = combine_uniq_regions(uniq)
        cool_enh = f'{COOL_BASE}/{tis}/{tis}_deephic.10kb.cool'
        
        #enh_uniq_data = get_uniq_interaction(cool_enh, uniq_mat, f'chr{CHR}', START, END, RESOLUTION)
        enh_uniq_data = get_uniq_interaction(cool_enh, uniq_mat)
        enh_uniq_nth_data = get_nth_percentile(enh_uniq_data, PERCENTILE)
        enh_mat = normalize(enh_uniq_nth_data, axis=1, norm='l1')
        enh_mat = np.maximum(enh_mat, enh_mat.T)
    else:
        mat_size = int((END - START) / RESOLUTION)
        enh_mat = np.zeros((mat_size, mat_size))

    # Plot matrix
    ax.imshow(enh_mat, alpha=alpha, cmap=cmaps[i], extent=extent,
              norm=colors.LogNorm(vmin=0.01, vmax=0.1))

    # Add TADs
    tad_path = f'{COOL_BASE}/{tis}/findTADs/{tis}_deephic_tads_domain.bed'
    if os.path.exists(tad_path):
        with open(tad_path, 'r') as tad_file:
            for line in tad_file:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                if chrom == f'chr{CHR}' and start < END and end > START:
                    plot_start = max(0, (start - START) / RESOLUTION)
                    plot_end = min((end - START) / RESOLUTION, (END - START) / RESOLUTION)
                    verts = [
                        (plot_start, plot_start),
                        (plot_end, plot_start),
                        (plot_end, plot_end),
                        (plot_start, plot_start),
                    ]
                    codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                    path = Path(verts, codes)
                    patch = PathPatch(path, facecolor='none', edgecolor=cmaps[i](256),
                                    alpha=0.7, lw=1.5, ls='-')
                    ax.add_patch(patch)

    ax.tick_params(left=False, right=False, labelleft=False, 
                   labelbottom=False, bottom=False)

    # Contact density heatmap
    axins_bottom = ax.inset_axes((0, -0.15, 1, .05))
    density = calculate_contact_density(enh_mat)
    density_2d = density.reshape(1, -1)
    im = axins_bottom.imshow(density_2d, aspect='auto', cmap='coolwarm', 
                            interpolation='nearest')
    axins_bottom.set_xticks([])
    axins_bottom.set_yticks([])
    axins_bottom.set_ylabel(' ', fontsize=5)

    # Gene track
#    bed_path = '/cluster/home/t111631uhn/HiC_ECC/genes.gtf'
    bed_path = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Genomes/gencode.vM25.annotation.gtf'
    bed = fanc.load(bed_path)
    axins_bed = ax.inset_axes((0, -0.25, 1, .05))
    bedplot = fancplot.GenePlot(bed, ax=axins_bed, n_ticks=5)
    bedplot.plot(f'Chr{CHR}:{START}-{END}')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])

save_dir = f'{COOL_BASE}/'
os.makedirs(save_dir, exist_ok=True)
image_name = f"HiC_channel_{CHR}_{GENE}_density_percentile{PERCENTILE}.pdf"
plt.savefig(image_name, format='pdf', bbox_inches='tight', dpi=500)
print(f"Saved: {image_name}")
plt.close()
