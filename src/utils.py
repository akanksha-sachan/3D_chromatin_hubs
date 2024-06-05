"""Validation of interaction calls, and other common utils"""

# pylint: disable=all
import os
import subprocess
import sys
from collections import defaultdict
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import h5py

# relative and absolute imports for running module and script respectively
# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from .configs.config_local import Config
except ImportError:
    from configs.config_local import Config


def oe_from_cooler(obs_numpy_matrix, threshold=1):
    """
    The O/E matrix is calculated as the log2 ratio of the raw contact matrix to the expected contact matrix.
    The expected contact matrix is calculated by filling in the average value of the diagonals of the raw contact matrix.
    Remove the NaN bins before calculating O/E so that interpolated edges aren't used

    Input: normalised observed counts as numpy matrix from cooler file
    Output: O/E matrix, O/E matrix with thresholded values, expected matrix, sums of expected values
    """

    # process matrix to mask out the centeromeric bins (NANs currently)
    matrix = np.nan_to_num(obs_numpy_matrix)

    # construct expected matrix
    expected_matrix = np.zeros_like(matrix).astype(float)
    sums = []
    for i in range(matrix.shape[0]):
        contacts = np.diag(matrix, i)
        # using on non zero diagonal elements to get the denominator
        non_zero_indices = np.nonzero(contacts)[0]
        if len(non_zero_indices) > 0:
            # chr wide expected, not factorized by number of chrs for genome-wide comparison
            expected_strength = contacts.sum() / len(non_zero_indices)
        else:
            expected_strength = 0
        sums.append(expected_strength)
        # uniform distribution of contacts across diagonals assumed
        x_diag, y_diag = np.diag_indices(matrix.shape[0] - i)
        x, y = x_diag, y_diag + i
        expected_matrix[x, y] = expected_strength
    expected_matrix += expected_matrix.T
    eps = 1e-5
    expected_matrix = np.nan_to_num(expected_matrix) + eps
    obs_over_expected = matrix / expected_matrix
    obs_over_expected[obs_over_expected == 0] = 1  # to avoid neg inf in log
    obs_over_expected = np.log(
        obs_over_expected
    )  # log transform the OE to get equal-sized bins, as the expected values need to be log binned
    # threshold elements based on M (TODO: decide on an adaptive threshold)
    obs_over_expected_filtered = np.where(
        obs_over_expected > threshold, obs_over_expected, 0
    )
    return obs_over_expected, obs_over_expected_filtered, expected_matrix, sums


def bedtools_makewindows(chromsizes_file, res, outfile):
    """runs this command: $ bedtools makewindows -g hg19.txt -w 1000000
    can add additional if statements to be able to run -b input.bed and -n 10 command
    """

    # Construct the bedtools command
    cmd = f'bedtools makewindows -g "{chromsizes_file}" -w {res}'

    # Run the command and redirect output to a file
    with open(outfile, "w") as file:
        process = subprocess.run(
            cmd, shell=True, stdout=file, stderr=subprocess.PIPE, text=True
        )

    # Handle command output
    if process.returncode == 0:
        print(f"Successfully created bins file: {outfile}")
    else:
        print(f"Error creating bins file: {outfile}. Error: {process.stderr}")


def read_null_terminated_string(binary_file) -> str:
    """
    Read null terminated string from a binary file
    """
    string_buffer = ""
    while True:
        byte = binary_file.read(1)
        if not byte:
            # EOF (end of file) or read error
            return string_buffer

        decoded_byte = byte.decode("utf-8", "backslashreplace")
        if decoded_byte == "\0":
            # Null terminator found, end the string
            return string_buffer

        string_buffer += decoded_byte


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


def format_loci_string(start, end, resolution_str):
    if "Mb" in resolution_str:
        start_unit = start // 1_000_000
        end_unit = end // 1_000_000
        unit = "Mb"
    elif "kb" in resolution_str:
        start_unit = start // 1_000
        end_unit = end // 1_000
        unit = "kb"
    else:
        start_unit = start
        end_unit = end
        unit = "bp"
    return f"{start_unit}-{end_unit}{unit}"


# Helper function for plotting
def plot_hic_map(dense_matrix, cmap, vmin=0, vmax=30, filename=None, title=""):
    d2 = dense_matrix
    d2[np.isnan(d2)] = 0
    d2[np.isinf(d2)] = 0
    plt.matshow(dense_matrix, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.title(title)
    if filename:
        plt.savefig(filename)
    plt.close()


def inspect_h5_file(file_path):
    with h5py.File(file_path, "r") as f:

        def print_attrs(name, obj):
            print(f"{name}: {dict(obj.attrs)}")

        print("HDF5 file structure:")
        f.visititems(print_attrs)

        print("\nReading all datasets:")
        for chrom in f.keys():
            for res_str in f[chrom].keys():
                for dataset in f[chrom][res_str].keys():
                    key = f"{chrom}/{res_str}/{dataset}"
                    data = pd.read_hdf(file_path, key=key)
                    print(f"Data for {key}:")
                    print(data.head())


if __name__ == "__main__":
    config = Config()
    bedtools_makewindows(
        config.paths.chrom_sizes_infile, 100000, config.paths.ref_genome_bins
    )
