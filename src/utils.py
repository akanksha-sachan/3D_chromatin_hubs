"""Validation of interaction calls, and other common utils"""

# pylint: disable=all
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

# relative and absolute imports for running module and script respectively
# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from ..configs.config_local import Config
except ImportError:
    from configs.config_local import Config


def bedtools_makewindows(chromsizes_file, resolution, tmp_dir):
    """runs this command: $ bedtools makewindows -g hg19.txt -w 1000000
    can add additional if statements to be able to run -b input.bed and -n 10 command
    """

    # create sub-directory for resoltuion/window sizes
    resolution_dir = os.path.join(tmp_dir, str(resolution))
    if not os.path.exists(resolution_dir):
        os.makedirs(resolution_dir)

    # extract name of output file from chromsizes_file (example: bins_hg19.bed)
    output_file = f'bins_{os.path.basename(chromsizes_file).split(".")[0]}.bed'
    output_path = os.path.join(resolution_dir, output_file)

    # Construct the bedtools command
    cmd = f'bedtools makewindows -g "{chromsizes_file}" -w {resolution}'

    # Run the command and redirect output to a file
    with open(output_path, "w") as file:
        process = subprocess.run(
            cmd, shell=True, stdout=file, stderr=subprocess.PIPE, text=True
        )

    # Handle command output
    if process.returncode == 0:
        print(f"Successfully created bins file: {output_path}")
    else:
        print(f"Error creating bins file: {output_file}. Error: {process.stderr}")

def bedpe_to_bed():
        """
        Convert bedpe to bed format for liftover utility
        """
        pass

def gen_coords_to_bin_index(gen_coords, res):
    """
    Convert genomic coordinates (base pairs) to bin positions (indices used in matrix representation) of a given resolution
    Assumption: genomic coordinates are from the same chromosome
    pos1 = np.array(data['pos1'])
        pos2 = np.array(data['pos2'])
        bin1 = np.floor(pos1 / res).astype('int')
        bin2 = np.floor(pos2 / res).astype('int')
    """

    region_indices = [0] * len(gen_coords)
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
    if input_chrom.startswith("chr"):
        stripped_chrom = input_chrom[3:]
        if stripped_chrom in chrom_keys:
            return stripped_chrom
    # If no match found, return None or raise an error
    return None


# test function
if __name__ == "__main__":
    config = Config("GM12878", 10000)
    chromsizes_file = config.paths.chrom_sizes_file
    tmp_dir = config.paths.temp_dir
    res = config.paths.resolution
    
