#code to assign bins to anchors etc, other things to be done to all feature calls etc

#need a function to make all the temp directories per cell line

#function for creating bins on the reference genome using bedtools python wrapper

import subprocess
import os
import pandas as pd
import numpy as np
from collections import defaultdict

# relative and absolute imports for running module and script respectively
try:
    from ..configs.config1 import Config
except ImportError:
    from configs.config1 import Config

class ValidContact:
    """
    Base class for validation of 3d genome feature calls 
    """

    def __init__(self, config):
        self.config = config
        self.temp_dir = config.paths.temp_dir #should exist by this step
        self.res_list = config.genomic_params.resolutions

class ValidLoop(ValidContact):
    """
    Assessing the reliability of loop calls by using reference databases and size distribution plots
    """

    def __init__(self, config, res):
        super().__init__(config)
        #self.loop_file = #set post configuring hiccups output params for file name
        #self.loop_list = #set from reading file a df containing chr 
        self.resolution = res

    def parse_bedpe(self, bedpe_file, lower=50000, upper=5000000):
        """
        Parse bedfiles of different resolutions and non-standard formats to the standard (chrN: start end) format
        """
        coords = defaultdict(set)
        with open(bedpe_file, 'r') as o:
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

    # def loop1D_distance(self):
    #     """
    #     Get the 1-D genomic distance between loop anchor centroids a and b
    #     """
    #     dis = []
    #     for c in loops:
    #         for s1, e1, s2, e2 in loops[c]:
    #             a = (s1 + e1) // (2 * res)
    #             b = (s2 + e2) // (2 * res)
    #             #convert dist from bins to bp by multiplying with res
    #             dis.append((b-a)*res)
    #     return dis
        

def bedtools_makewindows(chromsizes_file, resolution, tmp_dir):
    """runs this command: $ bedtools makewindows -g hg19.txt -w 1000000
        can add additional if statements to be able to run -b input.bed and -n 10 command
    """

    #create sub-directory for resoltuion/window sizes 
    resolution_dir = os.path.join(tmp_dir, str(resolution))
    if not os.path.exists(resolution_dir):
        os.makedirs(resolution_dir)

    #extract name of output file from chromsizes_file (example: bins_hg19.bed)
    output_file = f'bins_{os.path.basename(chromsizes_file).split(".")[0]}.bed'
    output_path = os.path.join(resolution_dir, output_file)

    # Construct the bedtools command
    cmd = f'bedtools makewindows -g "{chromsizes_file}" -w {resolution}'

    # Run the command and redirect output to a file
    with open(output_path, 'w') as file:
        process = subprocess.run(cmd, shell=True, stdout=file, stderr=subprocess.PIPE, text=True)

    # Handle command output
    if process.returncode == 0:
        print(f"Successfully created bins file: {output_path}")
    else:
        print(f"Error creating bins file: {output_file}. Error: {process.stderr}")

def gen_coords_to_bin_index(gen_coords, res):
    """
    Convert genomic coordinates (base pairs) to bin positions (indices used in matrix representation) of a given resolution
    Assumption: genomic coordinates are from the same chromosome
    pos1 = np.array(data['pos1'])
	pos2 = np.array(data['pos2'])
	bin1 = np.floor(pos1 / res).astype('int')
	bin2 = np.floor(pos2 / res).astype('int')
    """
    
    region_indices = [0]*len(gen_coords)
    for index, coord in enumerate(gen_coords):
        region_indices[index] = coord // res
    return region_indices

def standardize_chromosome(input_chrom, chrom_keys):
    """
    Normalize input chromosome name to match keys in the chrom_indices.
    Handles common chromosome naming conventions.
    """
    # Check direct match
    if input_chrom in chrom_keys:
        return input_chrom
    # Check with 'chr' prefix
    prefixed_chrom = f"chr{input_chrom}"
    if prefixed_chrom in chrom_keys:
        return prefixed_chrom
    # Remove 'chr' prefix if present and check
    if input_chrom.startswith('chr'):
        stripped_chrom = input_chrom[3:]
        if stripped_chrom in chrom_keys:
            return stripped_chrom
    # If no match found, return None or raise an error
    return None

#test function
if __name__ == "__main__":
    config = Config("GM12878", 10000)
    chromsizes_file = config.paths.chrom_sizes_file
    tmp_dir = config.paths.temp_dir
    resolution = config.resolution
    print(resolution)
    #bedtools_makewindows(chromsizes_file, resolution, tmp_dir)