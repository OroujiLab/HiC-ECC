import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from sklearn.preprocessing import normalize
from matplotlib.patches import Arc, Rectangle
from matplotlib.colors import ListedColormap
import fanc
import fanc.plotting as fancplot
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from scipy.ndimage import label
from scipy import stats
from scipy.stats import ranksums
import sys, numpy as np, pandas as pd, cooler
import os

# Parameters
RESOLUTION = 10000
CHR = 4
START = 148000000
END = 149200000

GENE = 'Test_Region'
PERCENTILE = 95
WIN_SZ = '310kb'
MERGE_INTERVALS = True

TISSUES = ['Liver','Small_Intestine','Lung','Kidney','Brain','Spleen','Pancreas','Large_Intestine']
COLORS = [
    "#F47D2B",  # Liver
    "#C06CAB",  # Small_Intestine
    "#FEE500",  # Lung
    "#89288F",  # Kidney
    "#D51F26",  # Brain
    "#F9B712",  # Spleen
    "#8A9FD1",  # Pancreas
    "#272E6A"   # Large_Intestine
]

COOL_BASE = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Hi-C/Mouse/cool_files_all'
UNIQ_BASE = '/cluster/home/t111631uhn/HiC_ECC/new_scripts/chess_output/unique_regions'

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
    clr = cooler.Cooler(f'{COOL_BASE}/{TISSUES[0]}/{TISSUES[0]}_deephic.10kb.cool')
    test_matrix = clr.matrix(balance=False).fetch(f'chr{CHR}:{START}-{END+RESOLUTION}')
    mat_size = test_matrix.shape[0]
    
    start_idx = int((region_start - START) / RESOLUTION)
    end_idx = int((region_end - START) / RESOLUTION)
    
    m = np.zeros((mat_size, mat_size), dtype=np.int8)
    bins = list(range(start_idx, min(end_idx+1, mat_size)))
    for i in bins:
        for j in bins:
            if i < mat_size and j < mat_size:
                m[i, j] = 1
    return m

def calculate_tissue_specificity_track(ref_tissue, tissue_raw_data, tissues, start, end, resolution):
    """
    Calculate tissue-specificity score for each genomic bin
    Combined score = log2FC * -log10(p-value)
    Higher score = more specific to ref_tissue
    """
    # Get matrix size
    mat_size = tissue_raw_data[ref_tissue].shape[0]
    specificity = np.zeros(mat_size)
    
    print(f"    Calculating tissue-specificity track for {ref_tissue}...")
    
    # For each genomic bin
    for bin_idx in range(mat_size):
        # Get interaction frequencies for this bin (average across its row/column)
        ref_row = tissue_raw_data[ref_tissue][bin_idx, :]
        ref_col = tissue_raw_data[ref_tissue][:, bin_idx]
        ref_signal = np.mean(np.concatenate([ref_row, ref_col]))
        
        # Get signals from all other tissues
        other_signals = []
        for tissue in tissues:
            if tissue != ref_tissue and tissue in tissue_raw_data:
                other_row = tissue_raw_data[tissue][bin_idx, :]
                other_col = tissue_raw_data[tissue][:, bin_idx]
                other_signal = np.mean(np.concatenate([other_row, other_col]))
                other_signals.append(other_signal)
        
        if len(other_signals) > 0:
            other_mean = np.mean(other_signals)
            
            # Calculate log2 fold-change
            log2fc = np.log2((ref_signal + 0.01) / (other_mean + 0.01))
            
            # Statistical test if we have enough other tissues
            if len(other_signals) >= 3:
                try:
                    # Test if ref_signal is significantly different from other tissues
                    stat, p_val = ranksums(other_signals, [ref_signal] * len(other_signals))
                    
                    # Weight by significance
                    if p_val > 0 and p_val < 1:
                        sig_weight = -np.log10(p_val)
                    elif p_val == 0:
                        sig_weight = 10
                    else:
                        sig_weight = 0
                    
                    # Combined score (capped at reasonable range)
                    combined_score = log2fc * sig_weight
                    specificity[bin_idx] = np.clip(combined_score, -10, 10)
                    
                except:
                    # If statistical test fails, just use fold-change
                    specificity[bin_idx] = np.clip(log2fc, -5, 5)
            else:
                # Just use fold-change if not enough tissues
                specificity[bin_idx] = np.clip(log2fc, -5, 5)
        else:
            specificity[bin_idx] = 0
    
    print(f"      Score range: {np.min(specificity):.2f} to {np.max(specificity):.2f}")
    print(f"      Mean score: {np.mean(specificity):.2f}")
    
    return specificity

def calculate_insulation_score(matrix, window_size=10):
    n = matrix.shape[0]
    insulation_score = np.zeros(n)

    for i in range(n):
        upstream_start = max(0, i - window_size)
        upstream_end = i
        downstream_start = i
        downstream_end = min(n, i + window_size)

        if upstream_end > upstream_start and downstream_end > downstream_start:
            insulation_score[i] = np.nanmean(
                matrix[upstream_start:upstream_end, downstream_start:downstream_end]
            )
        else:
            insulation_score[i] = np.nan

    mean_score = np.nanmean(insulation_score)
    std_score = np.nanstd(insulation_score)
    if std_score > 0:
        insulation_score = (insulation_score - mean_score) / std_score
    return insulation_score

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
    if not matrices:
        return None
    result = np.zeros(matrices[0].shape)
    for mat in matrices:
        m = np.matrix(mat)
        result += m
    if ones:
        result[result > 0] = 1
    return result

def combine_uniq_regions(uniq_reg_list):
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

def get_uniq_interaction(cool_path, uniq_mat):
    """Load from cool and apply unique region mask"""
    clr = cooler.Cooler(cool_path)
    matrix = clr.matrix(balance=False).fetch(f'chr{CHR}:{START}-{END+RESOLUTION}')
    matrix = np.array(matrix)
    
    if uniq_mat.shape != matrix.shape:
        if uniq_mat.shape[0] > matrix.shape[0]:
            uniq_mat = uniq_mat[:matrix.shape[0], :matrix.shape[1]]
        else:
            new_mat = np.zeros(matrix.shape, dtype=uniq_mat.dtype)
            new_mat[:uniq_mat.shape[0], :uniq_mat.shape[1]] = uniq_mat
            uniq_mat = new_mat
    
    return uniq_mat * matrix

def get_nth_percentile(matrix, n):
    result = matrix.copy()
    nth_val = np.percentile(result, n)
    result[result < nth_val] = 0
    return result

def calculate_region_purity(region_coords, tissue_matrices, tissues):
    start, end = region_coords
    total_positions = (end - start) * (end - start)

    tissue_values = {}
    for tissue in tissues:
        if tissue in tissue_matrices:
            region_slice = tissue_matrices[tissue][start:end, start:end]
            non_zero = np.sum(region_slice != 0)
            total = np.sum(region_slice)

            tissue_values[tissue] = {
                'total': total,
                'coverage': non_zero/total_positions * 100
            }

    total_signal = sum(t['total'] for t in tissue_values.values())

    purity_scores = {}
    for tissue, values in tissue_values.items():
        purity = (values['total']/total_signal * 100) if total_signal > 0 else 0
        coverage = values['coverage']
        purity_scores[tissue] = {
            'purity': purity,
            'coverage': coverage
        }

    if purity_scores:
        dominant_tissue = max(purity_scores.items(), key=lambda x: x[1]['purity'])
        return {
            'purity_scores': purity_scores,
            'dominant_tissue': dominant_tissue[0],
            'dominant_purity': dominant_tissue[1]['purity'],
            'dominant_coverage': dominant_tissue[1]['coverage']
        }
    else:
        return {
            'purity_scores': {},
            'dominant_tissue': None,
            'dominant_purity': 0,
            'dominant_coverage': 0
        }

def merge_overlapping_regions(regions, max_gap=3):
    """Merge regions that overlap or are very close together"""
    if not regions:
        return []
    
    regions = sorted(regions, key=lambda x: x['start'])
    merged = [regions[0]]
    
    for current in regions[1:]:
        last = merged[-1]
        
        if current['start'] <= last['end'] + max_gap:
            merged[-1] = {
                'start': min(last['start'], current['start']),
                'end': max(last['end'], current['end']),
                'size_kb': (max(last['end'], current['end']) - min(last['start'], current['start'])) * RESOLUTION / 1000,
            }
        else:
            merged.append(current)
    
    return merged

def load_loops(bedpe_path, chrom, start, end):
    """Load significant loops from BEDPE file"""
    if not os.path.exists(bedpe_path):
        return None
    
    loops = []
    with open(bedpe_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            
            chr1, start1, end1 = parts[0], int(parts[1]), int(parts[2])
            chr2, start2, end2 = parts[3], int(parts[4]), int(parts[5])
            
            if chr1 == f'chr{chrom}' and chr2 == f'chr{chrom}':
                anchor1_mid = (start1 + end1) // 2
                anchor2_mid = (start2 + end2) // 2
                
                if start <= anchor1_mid <= end and start <= anchor2_mid <= end:
                    loops.append({
                        'anchor1': anchor1_mid,
                        'anchor2': anchor2_mid,
                        'start1': start1,
                        'end1': end1,
                        'start2': start2,
                        'end2': end2
                    })
    
    return loops

def load_loops(bedpe_path, chrom, start, end):
    """Load loops from BEDPE file"""
    if not os.path.exists(bedpe_path):
        return []
    
    loops = []
    with open(bedpe_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:  # Need 8 columns for score and qvalue
                continue
            
            chr1, start1, end1 = parts[0], int(parts[1]), int(parts[2])
            chr2, start2, end2 = parts[3], int(parts[4]), int(parts[5])
            score = float(parts[6])      # ADD THIS
            qvalue = float(parts[7])     # ADD THIS
            
            if chr1 == f'chr{chrom}' and chr2 == f'chr{chrom}':
                anchor1_mid = (start1 + end1) // 2
                anchor2_mid = (start2 + end2) // 2
                
                if start <= anchor1_mid <= end and start <= anchor2_mid <= end:
                    loops.append({
                        'anchor1': anchor1_mid,
                        'anchor2': anchor2_mid,
                        'score': score,      # ADD THIS
                        'qvalue': qvalue     # ADD THIS
                    })
    
    return loops



def filter_loops_by_matrix_signal(loops, matrix, start, resolution, min_signal_threshold=0.01):
    """Only keep loops where BOTH anchors have visible signal in the filtered matrix"""
    filtered_loops = []
    
    for loop in loops:
        bin1 = (loop['anchor1'] - start) // resolution
        bin2 = (loop['anchor2'] - start) // resolution
        
        if bin1 < matrix.shape[0] and bin2 < matrix.shape[1]:
            if matrix[bin1, bin2] > min_signal_threshold:
                filtered_loops.append(loop)
    
    return filtered_loops

def n_cmaps(tissues, num_colormaps=8):
    cmaps = []
    for color in COLORS[:num_colormaps]:
        cmap = ListedColormap([color])
        cmaps.append(cmap)
    return cmaps

# ===== LOAD TISSUE MATRICES =====
tissue_matrices = {}
tissue_raw_data = {}

print("Loading tissue data...")
for tis in TISSUES:
    uniq_path = f'{UNIQ_BASE}/mm10_{tis}_chr{CHR}_{WIN_SZ}_0.bed'
    uniq = get_unique_regions(uniq_path)
    
    if uniq is not None:
        uniq_mat = combine_uniq_regions(uniq)
        cool_path = f'{COOL_BASE}/{tis}/{tis}_deephic.10kb.cool'
        
        enh_uniq_data = get_uniq_interaction(cool_path, uniq_mat)
        tissue_raw_data[tis] = enh_uniq_data.copy()
        
        enh_uniq_nth_data = get_nth_percentile(enh_uniq_data.copy(), PERCENTILE)
        enh_mat = normalize(enh_uniq_nth_data, axis=1, norm='l1')
        enh_mat = np.maximum(enh_mat, enh_mat.T)
        tissue_matrices[tis] = enh_mat
        print(f"  {tis} loaded: {enh_mat.shape}, nonzero: {np.count_nonzero(enh_uniq_data)}")
def plot_hic_with_density(density_tissue):
    max_arc_height = 0

    fig, ax = plt.subplots(figsize=(15, 20))
    size = tissue_matrices[list(tissue_matrices.keys())[0]].shape[0]
    ax.plot([0, size], [0, size], color="white", linewidth=3)
    extent = [0, size, size, 0]
    alpha = 0.5
    cmaps = n_cmaps(TISSUES)

    axins_arc = ax.inset_axes((0, -0.6, 1, .4))
    axins_arc.plot([0])

    axins_insulation = ax.inset_axes((0, -0.8, 1, .15))

    density_plotted = False

    for i, tis in enumerate(TISSUES):
        if tis in tissue_matrices:
            mat = tissue_matrices[tis]
            
            # Plot CHESS regions
            ax.imshow(mat, alpha=alpha, cmap=cmaps[i], extent=extent,
                     norm=colors.LogNorm(vmin=0.01, vmax=0.1))

            # TADs on UPPER triangle
            tad_path = f'{COOL_BASE}/{tis}/findTADs/{tis}_deephic_tads_domains.bed'
            if os.path.exists(tad_path):
                with open(tad_path, 'r') as tad_file:
                    for line in tad_file:
                        parts = line.strip().split()
                        if len(parts) < 3:
                            continue
                        chrom, start, end = parts[0], int(parts[1]), int(parts[2])

                        if chrom == f'chr{CHR}' and start < END and end > START:
                            plot_start = max(0, (start - START) // RESOLUTION)
                            plot_end = min((end - START) // RESOLUTION, size)

                            verts = [
                                (plot_start, plot_start),
                                (plot_start, plot_end),
                                (plot_end, plot_end),
                                (plot_end, plot_start),
                                (plot_start, plot_start)
                            ]
                            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]

                            path = Path(verts, codes)
                            patch = PathPatch(path, facecolor='none', edgecolor=COLORS[i],
                                            alpha=0.7, lw=1.5, ls='--')
                            
                            clip_verts = [(0, 0), (size, 0), (size, size), (0, 0)]
                            clip_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                            clip_path = Path(clip_verts, clip_codes)
                            patch.set_clip_path(PathPatch(clip_path, transform=ax.transData))
                            
                            ax.add_patch(patch)

            # Loops on LOWER triangle
            loop_path = f'{COOL_BASE}/chromosight/Spleen_loops_highconf.bedpe'
            loops = load_loops(loop_path, CHR, START, END)
            loops = filter_loops_by_matrix_signal(loops, mat, START, RESOLUTION, min_signal_threshold=0.01)
        
            if loops:
                for loop in loops:
                    bin1 = (loop['anchor1'] - START) // RESOLUTION
                    bin2 = (loop['anchor2'] - START) // RESOLUTION
                
                    if bin1 < bin2:
                        bin1, bin2 = bin2, bin1
                    
                    square_size = 2
                    
                    rect = Rectangle((bin2 - square_size/2, bin1 - square_size/2), 
                                   square_size, square_size,
                                   fill=False,
                                   edgecolor='black',
                                   linewidth=1.5,
                                   linestyle='--',
                                   alpha=0.7)
                    
                    clip_verts = [(0, 0), (size, size), (0, size), (0, 0)]
                    clip_codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
                    clip_path = Path(clip_verts, clip_codes)
                    rect.set_clip_path(PathPatch(clip_path, transform=ax.transData))
                    
                    ax.add_patch(rect)
                    if bin1 > bin2:
                        ax.text(bin2 + square_size, bin1, '*', 
                                fontsize=8, color='black', fontweight='bold')            
#                        ax.text(bin2 + square_size, bin1, f"*\nscore={loop['score']:.2f}\nq={loop['qvalue']:.2e}", 
#                                fontsize=5, color='black', verticalalignment='center',fontweight='bold',bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9))
                print(f'\t {tis}: Plotted {len(loops)} loops')

            # Purity regions
            threshold = np.mean(mat) + 0.3*np.std(mat)
            binary_mat = mat > threshold
            labeled_regions, num_regions = label(binary_mat)

            detected_regions = []
            min_region_size = 3

            for region_idx in range(1, num_regions + 1):
                region_mask = labeled_regions == region_idx
                rows, cols = np.where(region_mask)
                
                if len(rows) < min_region_size:
                    continue

                start_r, end_r = min(rows), max(rows)
                region_size_bp = (end_r - start_r) * RESOLUTION
                
                if region_size_bp >= 30000:
                    purity_info = calculate_region_purity((start_r, end_r), tissue_matrices, TISSUES)
                    
                    if purity_info['dominant_tissue'] == tis and purity_info['dominant_purity'] > 60:
                        detected_regions.append({
                            'start': start_r,
                            'end': end_r,
                            'size_kb': region_size_bp/1000,
                        })

            merged_regions = merge_overlapping_regions(detected_regions, max_gap=3)

            for region in merged_regions:
                purity_info = calculate_region_purity((region['start'], region['end']), tissue_matrices, TISSUES)
                
                rect = Rectangle((region['start']-0.5, region['start']-0.5), 
                                region['end']-region['start']+1, region['end']-region['start']+1,
                                fill=False,
                                edgecolor='black',
                                linewidth=1,
                                linestyle='-')
                ax.add_patch(rect)

#                ax.text(region['end']+1, region['start'],
                ax.text(region['end']-16, region['start']-1,
                       f"{purity_info['dominant_tissue']}\n"
                       f"Purity: {purity_info['dominant_purity']:.1f}%\n"
                       f"Coverage: {purity_info['dominant_coverage']:.1f}%\n"
                       f"Span: {region['size_kb']:.0f}kb",
                       color='black',
                       fontsize=9,
                       fontweight='bold',
                       verticalalignment='bottom',
                       bbox=dict(boxstyle='round,pad=0.5',
                                facecolor='white',
                                edgecolor='black',
                                alpha=0.8))

            ax.tick_params(left=False, right=False, labelleft=False,
                          labelbottom=False, bottom=False)
            print(f'\t {tis} Enhanced Matrix Plotted')
# Tissue-specificity track
            if tis == density_tissue and not density_plotted:
                specificity_score = calculate_tissue_specificity_track(
                    tis, tissue_raw_data, TISSUES, START, END, RESOLUTION
                )
                
                axins_specificity = ax.inset_axes((0, -0.05, 1, .02))
                specificity_2d = specificity_score.reshape(1, -1)
                
                # Calculate actual data range
                actual_min = np.min(specificity_score)
                actual_max = np.max(specificity_score)
                max_abs = max(abs(actual_min), abs(actual_max))
                
                # Use capped range for better color visibility
                if max_abs < 2:
                    vmin, vmax = -2, 2
                elif max_abs > 3:
                    vmin, vmax = -3, 3
                else:
                    vmin, vmax = -max_abs, max_abs
                
                print(f"      Data range: {actual_min:.2f} to {actual_max:.2f}")
                print(f"      Display range: {vmin:.2f} to {vmax:.2f}")
                
                im = axins_specificity.imshow(specificity_2d, aspect='auto',
                                              cmap='RdBu_r', interpolation='nearest',
                                              vmin=vmin, vmax=vmax)
                axins_specificity.set_xticks([])
                axins_specificity.set_yticks([])
                axins_specificity.set_ylabel('Specificity', fontsize=8, rotation=90, 
                                            ha='center', va='bottom')
                
                # Create a second inset axes specifically for the colorbar on the right
                axins_cbar = ax.inset_axes((1.01, -0.05, 0.01, .02))
                cbar = plt.colorbar(im, cax=axins_cbar, orientation='vertical')
                cbar.set_label('log2FC × -log10(p)', rotation=270, labelpad=15, fontsize=8)
                cbar.set_ticks([vmin, 0, vmax])
                cbar.set_ticklabels([f'{vmin:.1f}', '0', f'{vmax:.1f}'])
                cbar.ax.tick_params(labelsize=8)
                
                density_plotted = True
            # Insulation score
            insulation_score = calculate_insulation_score(mat, window_size=10)
            axins_insulation.plot(insulation_score, color=COLORS[i], label=tis, linewidth=1.5)

            # Arc plotting
            if tis in tissue_raw_data:
                arc_data = get_nth_percentile(tissue_raw_data[tis].copy(), 99)
                arc_data_tril = np.tril(arc_data, k=0)
                arc_data_tril_norm = normalize(arc_data_tril, axis=1, norm='l1')

                intervals = np.nonzero(arc_data_tril_norm)
                for j in range(len(intervals[0])):
                    interval_begin = intervals[0][j]
                    interval_end = intervals[1][j]

                    width = (interval_end - interval_begin)
                    center = interval_begin + width / 2
                    half_height = width
                    interaction_freq = arc_data_tril_norm[interval_begin][interval_end]
                    linewidth = 2*np.sqrt(interaction_freq)

                    if abs(half_height) > abs(max_arc_height):
                        max_arc_height = half_height
                    if abs(width) > 0:
                        axins_arc.add_patch(Arc((center, 0),
                                            width,
                                            2*half_height,
                                            theta1=0,
                                            theta2=180,
                                            color=cmaps[i](256),
                                            linewidth=linewidth))
                print(f'\t {tis} Arcs Plotted')

    # Simplified legend - only tissue colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS[i], label=TISSUES[i]) 
        for i in range(len(TISSUES))
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

    # Insulation
    axins_insulation.set_xlabel('')
    axins_insulation.set_ylabel('Insulation score', fontsize=10)
#    axins_insulation.legend(fontsize='small', ncol=8, loc="lower center")
    axins_insulation.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    axins_insulation.set_xlim(0, size)
    axins_insulation.set_xticks([])

    all_scores = []
    for tis in TISSUES:
        if tis in tissue_matrices:
            mat = tissue_matrices[tis]
            insulation_score = calculate_insulation_score(mat, window_size=10)
            valid_scores = insulation_score[~np.isnan(insulation_score)]
            all_scores.extend(valid_scores)

    if len(all_scores) > 0:
        y_min = np.percentile(all_scores, 5)
        y_max = np.percentile(all_scores, 95)
        y_range = max(abs(y_min), abs(y_max))
        axins_insulation.set_ylim(-y_range, y_range)

    if max_arc_height > 0:
        axins_arc.set_ylim(max_arc_height, 0)
        axins_arc.tick_params(left=False, right=False, labelleft=False,
                            labelbottom=False, bottom=False)
        axins_arc.set_ylabel('Enhanced - Top 1%', fontsize='medium')
    
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])

    bed_path = '/cluster/projects/epigenomics/BACKUP_31032025/EpigenomeLab/Aminnn/Genomes/gencode.vM25.annotation.gtf'
    if os.path.exists(bed_path):
        bed = fanc.load(bed_path)
        axins_bed = ax.inset_axes((0, -0.14, 1, .05))
        bedplot = fancplot.GenePlot(bed, ax=axins_bed, n_ticks=5, group_by="gene_name", squash=True, show_labels=True)
        bedplot.plot(f'chr{CHR}:{START}-{END}')
    
    fig.subplots_adjust(hspace=0.3)

    image_name = f"HiC_CHESS_Mouse_{GENE}_chr{CHR}_{START}_{END}_with_specificity_{PERCENTILE}_percentile.pdf"
    fig.savefig(image_name, format='pdf', bbox_inches='tight', dpi=500)
    plt.close(fig)
    print(f"\nSaved: {image_name}")
plot_hic_with_density('Spleen')
