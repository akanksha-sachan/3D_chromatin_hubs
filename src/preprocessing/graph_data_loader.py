######### CONTACT MATRIX (.hic/.mcool) -> EDGELIST (HDF5) of graph #########
######### Calling loops, compartments and TADs #########
######### Creating edge lists of local and global interactions; plotting edge_reliabilities #########

##TODO: fix variables from chrom 2 currents

# pylint: disable=all

import os
import struct
import sys
from functools import partial
from multiprocessing import Pool

import hicstraw
import numpy as np
import pandas as pd
import pyBigWig
from memory_profiler import profile
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# relative and absolute imports for running module and script respectively
# to test functions of this script inside this script itself, without NB
try:
    from ..configs.config_local import Config
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from utils import *

############ Query Hic data and call 3D genomic features to create HDF5 edgelist object ############


class Query:
    """
    Base class for querying both .mcool/.hic files for observed, OE counts, edges (not needed if only using .hic as input)
    Inherit: to return egdelist as hdf5 with local and global interactions
    """

    def __init__(self, config, chrom, res, res_str):
        self.config = config
        self.temp_dir = config.paths.temp_dir
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)
        self.chrom = chrom
        self.res = res
        self.res_str = res_str


class HiCQuery(Query):
    """
    Querying .hic files for observed and OE counts

    Raises:
        ValueError: If queried resolution not found
    """

    def __init__(self, config, chrom, res, res_str):
        super().__init__(config, chrom, res, res_str)  # instantiate parent class
        self.hic_file = config.paths.hic_infile  # path to .hic file
        self.hic_norm = config.genomic_params.hic_norm  # normalization method to query
        self.hic = hicstraw.HiCFile(self.hic_file)  # hic object from straw
        print("HiC file loaded")

        # instantiate nested classes
        self.ab_comp = self.ab_comp(
            config
        )  # passing config if they need specifc attrs from there
        self.loop = self.loop(config)
        self.tad = self.tad(config)

        # checking for data availability
        self.res_list = config.genomic_params.resolutions
        if not self.resolutions_present():
            raise ValueError("Queried resolutions are not present in .hic file")

        # TODO: add check for normalisation presence in .hic file

    def resolutions_present(self) -> bool:
        """
        Check if the resolution is present in the .hic file
        """
        available_resolutions = self.hic.getResolutions()

        return all(res in available_resolutions for res in self.res_list)

    def read_hic_header(self):
        """
        Read header of .hic file to get info
        """
        hic_file = self.hic_file
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

    @profile
    def observed_intra(self):
        """
        returns csr sparse matrix of contact records for one chromosome
        straw object : .binX [0] .binY [1] .counts [2] as attributes
        """
        chrom = self.chrom[3:]
        res = int(self.res)
        observed_list = hicstraw.straw(
            "observed", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        )

        return observed_list

    @profile
    def oe_intra(self, threshold=0):
        """
        returns contact records for one chromosome as straw object, thresholded if needed
        straw object : .binX [0] .binY [1] .counts [2] as attributes of list
        """
        chrom = self.chrom[3:]
        res = int(self.res)
        oe_list = hicstraw.straw(
            "oe", self.hic_norm, self.hic_file, chrom, chrom, "BP", res
        )
        if threshold != 0:
            oe_list = [record for record in oe_list if record.counts >= threshold]

        return oe_list

    def oe_intra_numpy(self, start, end, threshold=0):
        """
        Purpose: visualization of a slice of chromosomal OE matrix
        query OE as zoom data object from HiCFile class in straw
        mzd object : get records of form .binX .binY counts
        """
        chrom = self.chrom[3:]
        oe_mzd = self.hic.getMatrixZoomData(
            chrom, chrom, "oe", self.hic_norm, "BP", self.res
        )
        oe_numpy = oe_mzd.getRecordsAsMatrix(start, end, start, end)
        if threshold != 0:
            oe_numpy_thresh = np.where(oe_numpy > threshold, oe_numpy, 0)
            return oe_numpy_thresh
        return oe_numpy

    @profile
    def oe_intra_df(self, threshold=0):
        """
        returns DataFrame of contact records for one chromosome
        straw object : .binX [0] .binY [1] .counts [2] as attributes
        """
        chrom = self.chrom
        res = self.res
        oe_list = self.oe_intra(threshold)

        # Preallocate lists using list comprehensions
        x1 = [record.binX for record in oe_list]
        x2 = [record.binX + res for record in oe_list]
        y1 = [record.binY for record in oe_list]
        y2 = [record.binY + res for record in oe_list]
        counts = [record.counts for record in oe_list]

        # Use single value for chromosome columns
        chr_column = [chrom] * len(oe_list)

        # Create df from lists
        df = pd.DataFrame(
            {
                "chr1": chr_column,
                "x1": x1,
                "x2": x2,
                "chr2": chr_column,
                "y1": y1,
                "y2": y2,
                "counts": counts,
            }
        )

        # remove self-edges where x1 == y1 and x2 == y2
        df = df[(df["x1"] != df["y1"]) | (df["x2"] != df["y2"])]
        # sort the df by x1 (node1 starts) then y1 (node2 starts)
        df.sort_values(by=["x1", "y1"], inplace=True)
        df.reset_index(drop=True, inplace=True)

        return df

    @profile
    def oe_straw_to_csr(self, threshold=0):
        """
        Convert straw object to csr matrix
        """
        res = self.res
        straw_obj = self.oe_intra(threshold)
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

    class ab_comp:
        """Nested class for calling A/B compartments from .hic data and reliability checks"""

        def __init__(self, config):
            self.ab_bigwig_file = config.paths.compartments_infile

        def ab_score(self, matrix):
            """
            Raw -> normalised -> O/E -> Pearson -> PCA gives A/B

            Input: OE matrix
            Output: A/B compartment calls
            """

            matrix = matrix.toarray()
            # ensure diagonals are 1
            np.fill_diagonal(matrix, 1)
            # get pearson matrix
            matrix = np.corrcoef(matrix)
            np.fill_diagonal(matrix, 1)
            matrix[np.isnan(matrix)] = 0.0
            pca = PCA(n_components=1)
            projected_matrix = pca.fit_transform(matrix)
            return projected_matrix, pca

        def e1_to_bigwig(bins, e1_values, chromsizes, output_file="e1.bw"):
            """
            Save e1 as a bigwig track for one chr
            """

            e1_values = e1_values.flatten()
            chroms = bins["chrom"].values
            chroms = np.array([chroms[i].encode() for i in range(len(chroms))])
            starts = bins["start"].values.astype(int)
            ends = bins["end"].values.astype(int)
            # adding chromsize header to bigwig file
            bw = pyBigWig.open(output_file, "w")
            bw.addHeader(list(chromsizes.items()))  # dict of 'chr' and 'size'
            # adding entries (bulk addition as each chroms, starts, ends, values can be numpy arrays)
            bw.addEntries(chroms, starts, ends=ends, values=e1_values)
            bw.close()
            print(f"BigWig file saved to {output_file}")

        def get_ab_bigwig(self):
            #TODO (present in jupyter NB)
            """
            Get the PCA score for the compartment calls
            """
            pass

        def ab_score_correlation(self):
            #TODO
            """
            Reliability: Compare the compartment calls with reference databases using R2
            """
            pass

    class loop:
        """Nested class for analysis of loops from .hic data"""

        def __init__(self, config):
            self.hiccups_merged_infile = config.paths.hiccups_merged_infile
            self.loops_txt_infile = config.paths.loops_txt_infile
            self.loops_bedpe_outfile = config.paths.loops_bedpe_outfile
            self.geo_loops = config.paths.geo_loops_infile

        def looplist_to_bedpe(self):
            """
            Convert looplist.txt to bedpe format for visualization and comparison
            """
            import csv

            rows = []
            with open(self.loops_txt_infile, "r") as infile:
                reader = csv.reader(infile, delimiter="\t")

                # Skip the first line as it is a header
                next(reader)

                for row in reader:
                    # Convert chromosome names to have 'chr' prefix
                    chr1 = "chr" + row[0]
                    chr2 = "chr" + row[3]

                    # Create the new row for BEDPE format
                    new_row = [
                        chr1,
                        row[1],
                        row[2],
                        chr2,
                        row[4],
                        row[5],
                        row[6],
                        row[7],
                        row[8],
                        row[9],
                        row[10],
                        row[11],
                        row[12],
                        row[13],
                        row[14],
                        row[15],
                        row[16],
                        row[17],
                        row[18],
                        row[19],
                    ]
                    rows.append(new_row)

            # Sort rows lexicographically by the first column (chr1)
            rows.sort(key=lambda x: (x[0], int(x[1])))

            # Write the sorted rows to the output file
            with open(self.loops_bedpe_outfile, "w", newline="") as outfile:
                writer = csv.writer(outfile, delimiter="\t")

                # Write the BEDPE header
                writer.writerow(
                    [
                        "#chr1",
                        "x1",
                        "x2",
                        "chr2",
                        "y1",
                        "y2",
                        "color",
                        "o",
                        "e_bl",
                        "e_donut",
                        "e_h",
                        "e_v",
                        "fdr_bl",
                        "fdr_donut",
                        "fdr_h",
                        "fdr_v",
                        "num_collapsed",
                        "centroid1",
                        "centroid2",
                        "radius",
                    ]
                )

                # Write the sorted rows
                writer.writerows(rows)
                pass

        def parse_bedpe(self, bedpe_file, lower=50000, upper=5000000):
            """
            Parse bedfiles of different resolutions and non-standard formats to the standard (chrN: start end) format
            """
            from collections import defaultdict

            genomic_coords = defaultdict(set)
            with open(bedpe_file, "r") as o:
                for line in o:
                    # Skip header or comment lines
                    if line.startswith("#"):
                        continue

                    p = line.strip().split()

                    # Skip lines with mitochondrial DNA or any unassembled contigs or scaffolds
                    if "M" in p[0] or "_" in p[0]:
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
                    genomic_coords[chrom1].add((s1, e1, s2, e2))

            # Sort regions for each chromosome
            for chr in genomic_coords:
                genomic_coords[chr] = sorted(genomic_coords[chr])
            return genomic_coords

        def get_loop_size(self, loop_genomic_coords, res):
            """
            Get the 1-D genomic distance between loop anchor centroids a and b
            """
            dis = []
            for chr in loop_genomic_coords:
                for s1, e1, s2, e2 in loop_genomic_coords[chr]:
                    a = (s1 + e1) // (2 * res)
                    b = (s2 + e2) // (2 * res)
                    # convert dist from bins to bp by multiplying with res
                    dis.append((b - a) * res)
            return dis

        def loop_anchor_overlap_percent(self):
            #TODO
            """
            Reliability: Check for overlap between loop anchors within a 1-d dist threshold
            """
            ground_truth_loops = self.geo_loops
            hiccups_loops = self.hiccups_merged_infile

    class tad:
        """Nested class for insulation score calculation from .hic data"""

        def __init__(self, config):
            pass

        def get_insulation_score(self, m, windowsize=500000, res=10000):
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

        def ins_to_bigwig(bins, ins_values, chromsizes, output_file="ins.bw"):
            ins_values = ins_values.flatten()
            chroms = bins["chrom"].values
            chroms = np.array([chroms[i].encode() for i in range(len(chroms))])
            starts = bins["start"].values.astype(int)
            ends = bins["end"].values.astype(int)
            # adding chromsize header to bigwig file
            bw = pyBigWig.open(output_file, "w")
            bw.addHeader(list(chromsizes.items()))  # dict of 'chr' and 'size'
            # adding entries (bulk addition as each chroms, starts, ends, values can be numpy arrays)
            bw.addEntries(chroms, starts, ends=ends, values=ins_values)
            bw.close()
            print(f"BigWig file saved to {output_file}")


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

class DataLoader(HiCQuery):
    """
    Class for creating outputs for the data_loader script
    Returns:
    edgelists in .h5 format
    plots for oe matrices (thresholded matrices, ab_scores, correlations)
    """

    def __init__(self, config, chrom, res, res_str):
        super().__init__(
            config, chrom, res, res_str
        )  # instantiate parent class attributes
        self.edgelist_outfile = config.paths.edgelist_outfile

    def oe_intra_edgelist_single_chr(self, threshold=0, mode="a"):
        """Save pandas df oe edges in .h5 format"""
        oe_intra_df = self.oe_intra_df(threshold)
        thresh_str = str(threshold).replace(".", "_")
        
        with pd.HDFStore(self.edgelist_outfile, mode=mode) as store:
            store.put(f"{self.chrom}/_{self.res_str}/oe_intra_{thresh_str}", oe_intra_df, format="table")

    def oe_plot_single_chr(
        self, output_dir_oe_plot, start, end, threshold=0, cmap="bwr", vmin=0, vmax=1
    ):
        """save whole genome OE plots of a selected region and threshold"""
        oe_numpy_thresholded = self.oe_intra_numpy(start, end, threshold)
        filename = os.path.join(output_dir_oe_plot, f"{self.chrom}_oe_{threshold}.png")
        region_str = format_loci_string(start, end, self.res_str)
        plot_hic_map(
            oe_numpy_thresholded,
            cmap,
            vmin,
            vmax,
            filename,
            title=f"{self.chrom}:{region_str} OE",
        )


def whole_genome_edgelist(config, chromosomes, res, res_str, threshold=0):
    """
    multiprocess methods to run on the whole genome
    writing to .h5 using multiprocessing requires file locking
    """
    for chrom in chromosomes:  #use for loop to write sequentially to .h5
        loader = DataLoader(config, chrom, res, res_str)
        loader.oe_intra_edgelist_single_chr(threshold)


def oe_plots(
    chrom, config, res, res_str, output_dir_oe_plot, start, end, threshold=0
):
    """
    multiprocess to run on the whole genome
    writing to .h5 using multiprocessing requires file locking
    """
    loader = DataLoader(config, chrom, res, res_str)
    loader.oe_plot_single_chr(output_dir_oe_plot, start, end, threshold) #can pass cmap, vmin, vmax

def run_parallel_oe_plots(config, chromosomes, current_res, current_res_str, output_dir_oe_plot, start, end, threshold):
    # multiprocessing on whole genome
    with Pool() as pool:
        pool.map(
            partial(
                oe_plots,
                config=config,
                res=current_res,
                res_str=current_res_str,
                output_dir_oe_plot=output_dir_oe_plot,
                start=start,
                end=end,
                threshold=threshold,
            ),
            chromosomes,
        )

    print(f"All plots saved to {output_dir_oe_plot}")

if __name__ == "__main__":

    # whole genome run test
    config = Config()
    chromosomes = config.genomic_params.chromosomes
    current_res = config.genomic_params.resolutions[0]  # 1Mb for OE part
    current_res_str = config.genomic_params.res_strs[0]  # 1Mb for OE part

    # params for OE matrix visualisation
    threshold = 0  # tweak based on single chr viz
    start = 0
    end = 72000000

    # directory to save plots
    output_dir_oe_plot = os.path.join(config.paths.temp_dir, f"oe_plots_{threshold}")
    os.makedirs(output_dir_oe_plot, exist_ok=True)
    # custom colormap
    REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1, 1, 1), (1, 0, 0)])

    #write edgelist file for whole genome
    whole_genome_edgelist(config, chromosomes, current_res, current_res_str, threshold)
    
    
