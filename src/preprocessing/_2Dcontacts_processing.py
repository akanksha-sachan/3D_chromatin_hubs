######### CONTACT MATRIX (.hic/.mcool) QUERYING for Observed and OE counts #########
######### Calling loops, compartments and TADs #########
######### Creating edge lists of local and global interactions for making graphs #########

# pylint: disable=all
# import io
import json
import os
import struct
import subprocess
import sys

import cooler
import hicstraw
import numpy as np
from memory_profiler import profile
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA

# parallel processing
# import multiprocessing
# import numba as njit

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# relative and absolute imports for running module and script respectively
# to test functions of this script inside this script itself, without NB
try:
    from ..configs.config_local import Config
except ImportError:
    from configs.config_local import Config


class Query:
    """
    Base class for querying .mcool/.hic files for observed and OE counts
    Store queried contacts as sparse csr matrices fro scalability of adding edges
    """

    def __init__(self, config):
        self.config = config
        self.chromosomes = config.genomic_params.chromosomes
        self.res_list = config.genomic_params.resolutions
        self.res_strings = config.genomic_params.res_strs
        self.temp_dir = config.paths.temp_dir
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)


class HiCQuery(Query):
    """
    Querying .hic files for observed and OE counts

    Raises:
        ValueError: If queried resolution not found
    """

    def __init__(self, config):
        super().__init__(config)  # instantiate parent class
        self.hic_file = config.paths.hic_file  # path to .hic file
        self.hic_norm = config.genomic_params.hic_norm  # normalization method to query
        self.hic = hicstraw.HiCFile(self.hic_file)  # hic object from straw
        print("HiC file loaded")

        # checking for data availability
        if not self.resolutions_present():
            raise ValueError("Queried resolutions are not present in .hic file")

        # TODO: add check for normalisation presence in .hic file

    def resolutions_present(self) -> bool:
        """
        Check if the resolution is present in the .hic file
        """
        available_resolutions = self.hic.getResolutions()
        return all(res in available_resolutions for res in self.res_list)

    @profile
    def observed_intra(self, chrom, res):
        """
        returns csr sparse matrix of contact records for one chromosome
        straw object : .binX [0] .binY [1] .counts [2] as attributes
        """
        chrom = chrom[3:]
        res = int(res)
        observed_list = hicstraw.straw(
            "observed", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        )
        return observed_list

    @profile
    def oe_intra(self, chrom, res):
        """
        returns csr sparse matrix of contact records for one chromosome
        straw object : .binX [0] .binY [1] .counts [2] as attributes
        """
        chrom = chrom[3:]
        res = int(res)
        oe_list = hicstraw.straw(
            "oe", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        )
        return oe_list

    @profile
    def straw_to_csr(self, straw_obj, res):
        """
        Convert straw object to csr matrix
        """
        # convert to numpy
        straw_array = np.array(
            [(i.binX, i.binY, i.counts) for i in straw_obj],
            dtype=[("binX", np.int32), ("binY", np.int32), ("counts", np.float32)],
        )

        # use vectorized ops
        row = straw_array["binX"] // res
        col = straw_array["binY"] // res
        value = straw_array["counts"]
        dimension = max(row.max(), col.max()) + 1
        csr_mat = csr_matrix(
            (value, (row, col)), shape=(dimension, dimension), dtype=float
        )
        return csr_mat

    def read_null_terminated_string(self, binary_file) -> str:
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

    def read_hic_header(self, hic_file):
        """
        Read header of .hic file to get info
        """
        hic_header = {}
        with open(hic_file, "rb") as f:
            magic_string = struct.unpack("<3s", f.read(3))[0]
            f.read(1)
            if magic_string != b"HIC":
                return None  # this is not a valid .hic file
            version = struct.unpack("<i", f.read(4))[0]
            master_index = struct.unpack("<q", f.read(8))[0]
            hic_header["version"] = str(version)
            hic_header["master_index"] = str(master_index)
            genome = ""
            c = f.read(1).decode("utf-8")
            while c != "\0":
                genome += c
                c = f.read(1).decode("utf-8")
            hic_header["genome_id"] = str(genome)
            num_attributes = struct.unpack("<i", f.read(4))[0]
            attrs = {}
            for _ in range(num_attributes):
                key = struct.unpack("<s", f.read(1))[0]
                value = struct.unpack("<s", f.read(1))[0]
                attrs[key] = value
            hic_header["attributes"] = attrs
            # num_chrs = struct.unpack("<i", f.read(4))[0]
            # chroms = []
            # for _ in range(num_chrs):
            #     name = struct.unpack("<s", f.read(1))[0]
            #     length = struct.unpack("<i", f.read(4))[0]
            #     chroms.append((name, length))
            # hic_header["chromosomes"] = chroms

        return hic_header


# class McoolQuery(Query):
#     """
#     Querying .mcool files for observed and OE counts

#     Raises:
#         ValueError: If queried resolution not found
#     """

#     def __init__(self, config):
#         super().__init__(config)  # instantiate parent class
#         self.mcool_file = config.paths.mcool_file  # path to .mcool file
#         self.mcool = cooler.Cooler(self.mcool_file)  # cooler object from cooler
#         print("Mcool file loaded")


##### 3D A/B compartment calls #####


def get_oe_logtrans(matrix, threshold=1):
    """
    The O/E matrix is calculated as the log2 ratio of the raw contact matrix to the expected contact matrix.
    The expected contact matrix is calculated by filling in the average value of the diagonals of the raw contact matrix.
    Remove the NaN bins before calculating O/E so that interpolated edges aren't used
    """
    # process matrix to mask out the centeromeric bins (NANs currently)
    matrix = np.nan_to_num(matrix)

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


def get_ab_compartment(matrix):
    """
    Raw -> normalised -> O/E -> Pearson -> PCA gives A/B
    """
    # ensure diagonals are 1
    np.fill_diagonal(matrix, 1)
    # get pearson matrix
    matrix = np.corrcoef(matrix)
    np.fill_diagonal(matrix, 1)
    matrix[np.isnan(matrix)] = 0.0
    pca = PCA(n_components=1)
    projected_matrix = pca.fit_transform(matrix)
    return projected_matrix, pca


##### 3D loop calls #####


def run_hiccups(
    hic_file,
    output_dir,
    juicer_path,
    matrix_size=None,
    chromosomes=None,
    resolutions=None,
    verbose=False,
):
    """
    Runs the hiccups command from juicer_tools.

    Parameters:
    - hic_file: Path to the HiC file.
    - output_dir: Directory where the output will be saved.
    - juicer_path: Path to the juicer_tools jar file.
    - matrix_size: Optional. Size of the matrix.
    - chromosomes: Optional. List of chromosomes to include.
    - resolutions: Optional. List of resolutions to use.
    - verbose: Optional. If True, print additional output.
    """
    # Ensure the output directory exists
    try:
        os.stat(output_dir)
    except:
        os.mkdir(output_dir)

    # strip chromosome list of the 'chr' prefix
    if chromosomes:
        chromosomes = [str(chrom[3:]) for chrom in chromosomes]

    # Build the command
    cmd = ["java", "-jar", juicer_path, "hiccups"]

    if matrix_size:
        cmd.extend(["-m", str(matrix_size)])

    if chromosomes:
        cmd.extend(["-c", ",".join(map(str, chromosomes))])

    if resolutions:
        cmd.extend(["-r", ",".join(map(str, resolutions))])

    cmd.extend([hic_file, output_dir])

    # Execute the command
    if verbose:
        print("Running command: " + " ".join(cmd))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            if verbose:
                print("Error in command execution:")
                print(result.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd)

        if verbose:
            print("Command output:")
            print(result.stdout)

        return result.stdout, result.stderr
    except Exception as e:
        if verbose:
            print("Exception occurred while running the command:")
            print(e)
        raise


# class PeakachuLoops:
#     """
#     Peakachu loop calls
#     """
#     def __init__(self, config):
#         self.config = config
#         self.temp_dir = config["temp_dir"]
#         self.peakachu_dir = config["peakachu_dir"]
#         if not os.path.exists(self.temp_dir):
#             os.mkdir(self.temp_dir)


#     def call_loops(self, chrom):
#         """
#         Call loops using Peakachu wrapper
#         """
#         pass

##### 3D insulation score #####


def insulation_score(m, windowsize=500000, res=10000):
    """
    needs <10kb bin size to detect tads using diamond score method and >100M total filtered reads mapped to the genome
    """
    windowsize_bin = int(windowsize / res)
    m = np.nan_to_num(m)
    score = np.ones((m.shape[0]))
    for i in range(0, m.shape[0]):
        with np.errstate(divide="ignore", invalid="ignore"):
            v = np.sum(
                m[
                    max(0, i - windowsize_bin) : i,
                    i + 1 : min(m.shape[0] - 1, i + windowsize_bin + 1),
                ]
            ) / (
                np.sum(
                    m[
                        max(0, i - windowsize_bin) : min(
                            m.shape[0], i + windowsize_bin + 1
                        ),
                        max(0, i - windowsize_bin) : min(
                            m.shape[0], i + windowsize_bin + 1
                        ),
                    ]
                )
            )
            if np.isnan(v):
                v = 1.0
        score[i] = v
    # get log2 of the score
    # score[score == 0] = 1
    # score = np.log2(score)
    return score
