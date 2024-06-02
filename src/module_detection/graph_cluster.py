######### EGDELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from ..configs.config_local import Config
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from utils import *


class Graph:
    """graph constructor class: building csr matrix graph data object from edglelists"""

    def __init__(
        self, config, current_chrom, current_res_str, nodeset_key="oe_intra_0"
    ):
        self.edgelist_h5_infile = (
            config.paths.edgelist_outfile
        )  # outfile from dataloader
        self.query_group_chrom = current_chrom  # current chrom to build graph for
        self.query_subgroup_res = f"_{current_res_str}"
        self.query_key_edge = nodeset_key
        self.bins_bed = config.paths.ref_genome_bins

        self.edge_df = None  # extracted df for quick access
        self.nodeset_attrs = None  # {start:idx} : start indicates the genomic loci of the node, and idx is nodeset order
        self.affinity_matrix = None  # graph matrix

    def load_edges(self):
        """extract pandas dataframe dataset for building graph from h5"""
        dataset_path = (
            f"{self.query_group_chrom}/{self.query_subgroup_res}/{self.query_key_edge}"
        )
        with pd.HDFStore(self.edgelist_h5_infile, mode="r") as store:
            self.edge_df = store[dataset_path]

    def edgelist_to_csr_affinity_matrix(self):
        """convert edgelist pandas df graph to a CSR matrix graph (affinity matrix) for clustering"""
        nodeset = pd.unique(
            self.edge_df[["x1", "y1"]].values.ravel("K")
        )  # node set: stores unique starts
        self.nodeset_attrs = {
            start: set_idx for set_idx, start in enumerate(nodeset)
        }  # save in format (loci start: idx of nodeset)
        rows_node_id = self.edge_df["x1"].map(
            self.nodeset_attrs
        )  # assigning nodeset id to edgelist nodes using start
        cols_node_id = self.edge_df["y1"].map(self.nodeset_attrs)
        weights = self.edge_df["counts"].values
        num_nodes = len(
            nodeset
        )  # create affinity matrix of only current edges and not whole chr.
        self.affinity_matrix = csr_matrix(
            (weights, (rows_node_id, cols_node_id)), shape=(num_nodes, num_nodes)
        )
        # make upper triangular matrix symmetric (not needed for spectral clustering but for visualization of graph)
        #self.affinity_matrix = affinity_matrix + affinity_matrix.T - csr_matrix((affinity_matrix.diagonal(), (range(num_nodes), range(num_nodes))), shape=(num_nodes, num_nodes))

    def construct_hub_edgelist():
        """
        1. find local overlapping nodes in the global nodeset, remove all other global nodes (and non-overlapping local nodes)
        2. borrow gobal edge between 2 oe nodes and distribute it among the local nodes in a fully connected manner
        """
        pass


# module detection class
class Cluster(Graph):
    """
    spectral clustering on intra chr graph built from
    a) single scale OE edges
    b) single scale loop edges
    c) multi scale OE (global) + loop (local) edges
    """

    def __init__(self, config, chrom, current_res_str):
        super().__init__(config, chrom, current_res_str)
        self.load_edges()  # call function from parent to load edges
        self.edgelist_to_csr_affinity_matrix()
        self.cluster_labels = None
        self.graphs_outdir = config.paths.gexf_dir

    def perform_spectral_clustering(self, n_clusters):
        """perform spectral clustering on the affinity matrix"""
        affinity_matrix = self.affinity_matrix
        spectral = SpectralClustering(
            n_clusters=n_clusters, affinity="precomputed", assign_labels="discretize"
        )  #cluster_labels coming from image segmentation algo
        self.cluster_labels = spectral.fit_predict(affinity_matrix)

    def append_labels_to_nodeset(self):
        """ add lables as attrs to the nodeset dict """
        nodeset_dict = self.nodeset_attrs #{start:idx}
        # append cluster labels as {start: (idx, cluster_label)}
        self.nodeset_attrs = {start: (set_idx, self.cluster_labels[set_idx]) for start, set_idx in nodeset_dict.items()}

    def create_gexf(self):
        """store the clusters in a gexf format for viz"""
        pass


if __name__ == "__main__":

    config = Config()
    # inspect_h5_file(config.paths.edgelist_outfile) #utils func

    modules = Cluster(
        config, config.genomic_params.chromosomes[0], config.genomic_params.res_strs[0]
    )
    modules.perform_spectral_clustering(2)
    modules.append_labels_to_nodeset()
