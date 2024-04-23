import argparse
import shutil
import h5py
import numpy as np
import pandas as pd
import cooler
import cooltools
import cooltools.lib.plotting
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import json
import pickle
from tqdm import tqdm
import os, subprocess
import bioframe
import pyBigWig
import warnings
from sklearn.decomposition import PCA
import sys
from collections import defaultdict

## from concurrent.futures import ProcessPoolExecutor, as_completed
## import multiprocessing

###### utils functions for creating pairs from BAM files ######

def is_valid_record(chrom1, pos1, chrom2, pos2, chrom_sizes):
    """
    Checks if the positions are within the valid range for their chromosomes
    """

    return (chrom1 in chrom_sizes and pos1 <= chrom_sizes[chrom1]) and \
           (chrom2 in chrom_sizes and pos2 <= chrom_sizes[chrom2])

def filter_pairs(chrom_sizes_file, input_file, output_file):
    """
    Filters the pairs file to remove invalid records
    """

    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom_sizes[parts[0]] = int(parts[1])

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        line_number = 0  # Initialize line number counter
        for line in f_in:
            # Skip lines starting with #
            if line.startswith('#'):
                continue
            line_number += 1  # Increment line number
            try:
                parts = line.strip().split()
                chrom1, pos1, chrom2, pos2 = parts[1], int(parts[2]), parts[3], int(parts[4])            
                if is_valid_record(chrom1, pos1, chrom2, pos2, chrom_sizes):
                    f_out.write(line)
            except ValueError as e:
                print(f"Error processing line {line_number}: {line.strip()} - {e}")
                continue
                
def bin_reference_genome(chromsize_file, res=[10000], tmp_dir='tmp'):
    """
    Bin the reference genome hg19/hg38.fa into bed files at specific resolutions; can pass multires list
    """
    # Ensure the temporary directory exists
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    
    chromsize = pd.read_csv(chromsize_file, sep='\t', header=None, index_col=0).to_dict()[1]
    #can bin the reference genome at different resolutions from res-list
    for resolution in res:
        res_dir = os.path.join(tmp_dir, f'{resolution}_resolution')
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)
        for c in chromsize:
            ngene = int(chromsize[c]//resolution) + 1
            bed = [[c, i*resolution, (i+1)*resolution] for i in range(ngene)]
            bed[-1][-1] = chromsize[c]
            np.savetxt(os.path.join(res_dir, f'{c}.bed'), bed, fmt='%s', delimiter='\t')
            print(resolution, c)
    return 0
    
def assign_bins_to_bed(bed, bins):
    """
    Assign cooler bins at a specific res to bed files by overlapping the start positions
    """
    #ensure that the starts are ints and not strs
    bins['start'] = pd.to_numeric(bins['start'], errors='coerce').fillna(0).astype(int)
    bed['start'] = pd.to_numeric(bed['start'], errors='coerce').fillna(0).astype(int)
    #overlap the starts in bins and the bed df
    bed['bin_id'] = np.digitize(bed['start'], bins['start'])
    return bed

def assign_bins_to_loops(loops, bins):
    """
    Assign cooler bins at a specific res to loop-anchors by overlapping the start positions
    """
     #convert strs to ints
    bins['start'] = pd.to_numeric(bins['start'], errors='coerce').fillna(0).astype(int)
    loops['start1'] = pd.to_numeric(loops['start1'], errors='coerce').fillna(0).astype(int)
    loops['start2'] = pd.to_numeric(loops['start2'], errors='coerce').fillna(0).astype(int)
    #assign bins to loop anchors
    loops['bin1'] = np.digitize(loops['start1'], bins['start'])
    loops['bin2'] = np.digitize(loops['start2'], bins['start'])
    return loops

def parsebed(chiafile, lower=50000, upper=5000000):
    """
    Parse bedfiles of different resolutions and non-standard formats to the standard (chrA: regions) format
    """
    coords = defaultdict(set)
    with open(chiafile, 'r') as o:
        for line in o:
            # Skip header or comment lines
            if line.startswith('#'):
                continue

            p = line.strip().split()

            # Skip lines with mitochondrial DNA or any unassembled contigs or scaffolds
            if 'M' in p[0] or '_' in p[0]:
                continue
            
            # Convert start and end positions to integers
            s1, e1, s2, e2 = int(p[1]), int(p[2]), int(p[4]), int(p[5])
            # Ensure that the region is ordered from smallest to largest
            if s1 > s2:
                s1, s2 = s2, s1
                e1, e2 = e2, e1
            
            # Skip regions outside the specified bounds
            if s2 - s1 > upper or s2 - s1 < lower:
                continue
            
            # Format chromosome to include 'chr' prefix if it is not present
            chrom1 = f'chr{p[0].lstrip("chr")}'
            chrom2 = f'chr{p[3].lstrip("chr")}'

            # chrom1 and chrom2 are expected to be the same:
            if chrom1 != chrom2:
                continue  # Or handle the case if needed
            
            # Add to coordinates, assuming symmetric interactions (only one chromosomal identifier needed)
            coords[chrom1].add((s1, e1, s2, e2))

    # Sort regions for each chromosome
    for c in coords:
        coords[c] = sorted(coords[c])
    return coords
    
###### functions for feature extraction from HiC matrices ######

def sqrt_norm(matrix):
	coverage = (np.sqrt(np.sum(matrix, axis=-1)))
	with np.errstate(divide='ignore', invalid='ignore'):
		matrix = matrix / coverage.reshape((-1, 1)); matrix = matrix / coverage.reshape((1, -1))
	matrix[np.isnan(matrix)] = 0.0; matrix[np.isinf(matrix)] = 0.0
	return matrix
    
def pearson(matrix):
	return np.corrcoef(matrix)
    
def get_oe_logtrans(M,threshold=1):
    """ 
    The O/E matrix is calculated as the log2 ratio of the raw contact matrix to the expected contact matrix.
    The expected contact matrix is calculated by filling in the average value of the diagonals of the raw contact matrix.
    Remove the NaN bins before calculating O/E so that interpolated edges aren't used
    """
    #process M to mask out the centeromeric bins (NANs currently)
    M = np.nan_to_num(M)
    
    #construct expected matrix
    E = np.zeros_like(M).astype(float)
    l = len(M)
    sums = []
    for i in range(M.shape[0]):
        contacts = np.diag(M,i)
		#using on non zero diagonal elements to get the denominator
        non_zero_indices = np.nonzero(contacts)[0]
        if len(non_zero_indices) > 0:
              #chr wide expected, not factorized by number of chrs for genome-wide comparison
              expected = contacts.sum() / len(non_zero_indices) 
            
        else:
              expected = 0
        sums.append(expected)
        #uniform distribution of contacts across diagonals assumed
        x_diag,y_diag = np.diag_indices(M.shape[0]-i)
        x,y = x_diag,y_diag+i
        E[x,y] = expected
    E += E.T
    eps=1e-5
    E = np.nan_to_num(E) + eps
    OE = M / E 
    OE[OE == 0] = 1 #to avoid neg inf in log
    OE = np.log(OE) #log transform the OE to get equal-sized bins, as the expected values need to be log binned
	# threshold elements based on M (TODO: decide on an adaptive threshold)
    OE_filtered = np.where(OE > threshold, OE, 0)
    return OE , OE_filtered, E, sums

def get_A_B(matrix):
    """
    Raw -> normalised -> O/E -> Pearson -> PCA gives A/B
    """
    
    #ensure diagonals are 1
    np.fill_diagonal(matrix, 1)
    #get pearson matrix
    matrix = np.corrcoef(matrix)
    np.fill_diagonal(matrix, 1)
    matrix[np.isnan(matrix)] = 0.0
    pca = PCA(n_components=1)
    y = pca.fit_transform(matrix)
    return y, pca


def e1_bigwig(bins, e1_values, chromsizes, output_file='e1.bw'):
    """
    Save e1 as a bigwig track for one chr
    """
    
    e1_values = e1_values.flatten()
    chroms = bins['chrom'].values
    chroms = np.array([chroms[i].encode() for i in range(len(chroms))])
    starts = bins['start'].values.astype(int)
    ends = bins['end'].values.astype(int)
    #adding chromsize header to bigwig file
    bw = pyBigWig.open(output_file, "w")
    bw.addHeader(list(chromsizes.items())) #dict of 'chr' and 'size'
    #adding entries (bulk addition as each chroms, starts, ends, values can be numpy arrays)
    bw.addEntries(chroms, starts, ends=ends, values=e1_values)  
    bw.close()
    print(f'BigWig file saved to {output_file}')

def insulation_score(m, windowsize=500000, res=10000):
    """
    needs <10kb bin size to detect tads using diamond score method and >100M total filtered reads mapped to the genome
    """
    windowsize_bin = int(windowsize / res)
    m = np.nan_to_num(m)
    score = np.ones((m.shape[0]))
    for i in range(0, m.shape[0]):
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.sum(m[max(0, i - windowsize_bin): i, i + 1: min(m.shape[0] - 1, i + windowsize_bin + 1)]) / (np.sum(
                 m[max(0, i - windowsize_bin):min(m.shape[0], i + windowsize_bin + 1),
                   max(0, i - windowsize_bin):min(m.shape[0], i + windowsize_bin + 1)]))
            if np.isnan(v):
                v = 1.0
        score[i] = v
		#get log2 of the score
		#score[score == 0] = 1
		#score = np.log2(score)
    return score

def ins_bigwig(bins, ins_values, chromsizes, output_file='ins.bw'):
    ins_values = ins_values.flatten()
    chroms = bins['chrom'].values
    chroms = np.array([chroms[i].encode() for i in range(len(chroms))])
    starts = bins['start'].values.astype(int)
    ends = bins['end'].values.astype(int)
    #adding chromsize header to bigwig file
    bw = pyBigWig.open(output_file, "w")
    bw.addHeader(list(chromsizes.items())) #dict of 'chr' and 'size'
    #adding entries (bulk addition as each chroms, starts, ends, values can be numpy arrays)
    bw.addEntries(chroms, starts, ends=ends, values=ins_values)  
    bw.close()
    print(f'BigWig file saved to {output_file}')

######## Sub-compartments (ground truth:gm12878, hg19)
    
def get_subcompartments_intra(subcompartment_file, bins, resolution=100000):
    """
    INTRA-CHR ::
    Load subcompartment annotations of varying bin length from bed file, then transfer the labels to current resolution HiC bins
    Assuming that the chromosome assembly for subcompartment_file matches the bins (hg19/hg38)
    """
    #add header to subcompartment file
    subcompartments = pd.read_csv(subcompartment_file, sep='\t', header=None)
    #only select the first 4 columns
    subcompartments = subcompartments.iloc[:, 0:4]
    subcompartments.columns = ['chrom', 'start', 'end', 'subcompartment']
    subcompartments['chrom'] = subcompartments['chrom'].astype(str)
    subcompartments['start'] = subcompartments['start'].astype(int) 
    subcompartments['end'] = subcompartments['end'].astype(int)
    subcompartments['subcompartment'] = subcompartments['subcompartment'].astype(str)
    #transfer the subcompartment labels to current resolution HiC bins
    subC_labels = [] 
    #loop through each bin's end position to get the subcompartment label
    bins = bins.reset_index(drop=True)
    for i in range(len(bins)):
        bin_end = bins['end'][i]
        chrom = bins['chrom'][i]
        subC = subcompartments[subcompartments['chrom'] == chrom] #get all labels of the current chromosome
        subC = subC.reset_index(drop=True)
        bin_label = ''
        for j in range(len(subC)):
            subC_end = subC['end'][j]
            if bin_end < subC_end:
                bin_label = subC['subcompartment'][j]
                break
            else:
                x = subC_end - bin_end
                y = resolution - x
                if x > y:
                    bin_label = subC['subcompartment'][j]
                else:
                    bin_label = subC['subcompartment'][j+1]
            if j == len(subC)-1:
                bin_label = subC['subcompartment'][j]
        subC_labels.append(bin_label)
    #append the subcompartment labels to the bins df
    subC_labels = pd.Series(subC_labels).replace('nan', 'N').values
    bins['subcompartment'] = subC_labels
    return bins

######## Loop related functions #######

def loop1D_distance(loops, res):
    """
    Get the 1-D genomic distance between loop anchors
    """
    dis = []
    for c in loops:
        for s1, e1, s2, e2 in loops[c]:
            a = (s1 + e1) // (2 * res)
            b = (s2 + e2) // (2 * res)
            #convert dist from bins to bp by multiplying with res
            dis.append((b-a)*res)
    return dis

def compare_loops(gt_loops, pred_loops):
    """
    Compare the loops called by different tools. Input loops have bin_ids in the last two columns.
    A loop is overlapping if either bin_id (bin1 or bin2) coincides in the two loops.
    Calculate the % overlap of pred_loops to gt_loops and return both the percentage and the DataFrame of overlapping loops.
    """
    if len(pred_loops) == 0:
        return 0, pd.DataFrame()

    # Perform an inner merge on 'bin1'
    coinciding_loops_bin1 = pd.merge(gt_loops, pred_loops, on='bin1')
    # Perform an inner merge on 'bin2'
    coinciding_loops_bin2 = pd.merge(gt_loops, pred_loops, on='bin2')
    # Concatenate results from both merges and remove duplicates
    coinciding_loops = pd.concat([coinciding_loops_bin1, coinciding_loops_bin2]).drop_duplicates().reset_index(drop=True)
    # Calculate percentage overlap
    percent_overlap = len(coinciding_loops) / len(pred_loops)
    return percent_overlap, coinciding_loops

