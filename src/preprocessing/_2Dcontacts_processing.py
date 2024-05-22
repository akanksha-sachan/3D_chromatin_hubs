######### CONTACT MATRIX (.hic/.mcool) QUERYING for Observed and OE counts #########
# pylint: disable=all
# import io
import json
import os
import struct
import subprocess
import sys

import cooler
import hicstraw

# import numba as njit
import numpy as np
from memory_profiler import profile
from scipy.sparse import csr_matrix

# parallel processing
# import multiprocessing

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# relative and absolute imports for running module and script respectively
try:
    from ..configs.config1 import Config
except ImportError:
    from configs.config1 import Config


###### .mcool/.hic file querying, we know which format is the input  ######


class Query:
    """
    Base class for querying .mcool/.hic files for observed and OE counts
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
        Returns a list of intrachr contacts
        """
        chrom = chrom[3:]
        res = int(res)
        contacts_observed = hicstraw.straw(
            "observed", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        ) #hicstraw object with bin.X bin.Y and counts as attrs of one row 
        return contacts_observed

    @profile
    def oe_intra(self, chrom, res):
        """
        Returns a list of intrachr contacts
        """
        chrom = chrom[3:]
        res = int(res)
        contacts_oe = hicstraw.straw(
            "oe", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        )
        return contacts_oe

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


class Process:
    """
    Post querying processing of observed/oe counts
    """

    def __init__(self, config):
        self.config = config
        self.temp_dir = config["temp_dir"]
        self.resolutions = config["resolutions"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)

    def sqrt_norm(self, matrix):
        """
        Square root normalization of a matrix
        """
        coverage = np.sqrt(np.sum(matrix, axis=-1))
        with np.errstate(divide="ignore", invalid="ignore"):
            matrix = matrix / coverage.reshape((-1, 1))
            matrix = matrix / coverage.reshape((1, -1))
        matrix[np.isnan(matrix)] = 0.0
        matrix[np.isinf(matrix)] = 0.0
        return matrix

    def pearson(self, matrix):
        """
        Pearson correlation coefficient between two matrices
        """
        return np.corrcoef(matrix)


# write main

if __name__ == "__main__":
    config = Config("GM12878", 1000)
    # want to check hic_geader function on the hic file mentioned in paths
    query = HiCQuery(config)
    # check if hic file is readable
    inform = query.read_hic_header(config.paths.hic_file); print(inform)
    contacts = query.oe_intra(
        config.genomic_params.chromosomes[0], config.paths.resolution
    ); print(contacts[0:10])
