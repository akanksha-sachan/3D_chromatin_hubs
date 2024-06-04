######### EGDELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import sys
from functools import partial
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
    from ..preprocessing import *
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from preprocessing import *
    from utils import *


class Graph:
    """graph constructor class: building csr matrix graph data object from edglelists"""

    def __init__(
        self, config, current_chrom, current_res_str, nodeset_key="oe_intra_0"
    ):
        #basic params init
        self.config = config
        self.chrom = current_chrom
        self.res_str = current_res_str

        #query data_loader object init
        self.edgelist_h5_infile = (
            config.paths.edgelist_outfile
        ) 
        self.query_group_chrom = current_chrom  # current chrom to build graph for
        self.query_subgroup_res = f"_{current_res_str}"
        self.query_key_edge = nodeset_key
        self.bins_bed = config.paths.ref_genome_bins

        #graph data attrs init
        self.edge_df = None  # extracted df for quick access
        self.nodeset_attrs = None  # {start:(idx, attrs)} : start indicates the genomic loci of the node, and idx is nodeset order
        self.affinity_matrix = None  # affinity matrix to represent constructed graph

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
        weights = self.edge_df["counts"].values  # edge weights
        num_nodes = len(
            nodeset
        )  # create affinity matrix of only current edges and not whole chr.
        affinity_matrix = csr_matrix(
            (weights, (rows_node_id, cols_node_id)), shape=(num_nodes, num_nodes)
        )
        # make upper triangular matrix symmetric (not needed for spectral clustering but for visualization of graph)
        self.affinity_matrix = (
            affinity_matrix
            + affinity_matrix.T
            - csr_matrix(
                (affinity_matrix.diagonal(), (range(num_nodes), range(num_nodes))),
                shape=(num_nodes, num_nodes),
            )
        )

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

    def __init__(self, config, chrom, current_res, current_res_str, n_clusters=2):
        super().__init__(config, chrom, current_res_str)
        self.load_edges()  # call function from parent to load edges
        self.edgelist_to_csr_affinity_matrix()
        self.res = current_res
        self.number_of_clusters = n_clusters
        self.cluster_labels = None
        self.graphs_outdir = config.paths.gexf_dir

    def spectral_clustering(self):
        """perform spectral clustering on the affinity matrix"""
        affinity_matrix = self.affinity_matrix
        spectral = SpectralClustering(
            n_clusters=self.number_of_clusters,
            affinity="precomputed",
            assign_labels="discretize",
        )  # cluster_labels coming from image segmentation algo
        self.cluster_labels = spectral.fit_predict(affinity_matrix)

    def append_cluster_labels_to_nodeset(self):
        """add lables as attrs to the nodeset dict"""
        nodeset_dict = self.nodeset_attrs  # {start:set_idx}
        # append cluster labels as {start: (set_idx, cluster_label)}
        self.nodeset_attrs = {
            start: (set_idx, self.cluster_labels[set_idx])
            for start, set_idx in nodeset_dict.items()
        }
    
    def create_gexf(self):
        """store the clusters in a gexf format for viz"""
        G = nx.Graph()
        # add nodes from nodeset
        for start, (set_idx, cluster_label) in self.nodeset_attrs.items():
            G.add_node(
                set_idx,
                start=str(start),
                cluster_label=cluster_label,
            )
        # add edges from edgelist, map nodes in edgelist to nodes in the nx graph using start
        for _, row in self.edge_df.iterrows():
            start_x, start_y = row["x1"], row["y1"]
            start_x_set_idx, start_y_set_idx = (
                self.nodeset_attrs[start_x][0],
                self.nodeset_attrs[start_y][0],
            )
            G.add_edge(start_x_set_idx, start_y_set_idx, weight=row["counts"])
        # write to gexf
        outfile = (
            f"{self.graphs_outdir}/{self.query_group_chrom}_{self.query_key_edge}.gexf"
        )
        nx.write_gexf(G, outfile)
    
    def oe_confusion_matrix(self):
        """calculate confusion matrix for clustering OE edges using AB compartments as ground truth"""
        query = HiCQuery(self.config, self.chrom, self.res, self.res_str)
        ab_bed_dict = query.ab_comp.load_bigwig_chromosomal_ab()
        # append the A/B labels to the nodeset_attrs dict by mapping 'start' key to get {start: (set_idx, cluster_label, ab_label)}
        self.nodeset_attrs = {
            start: (set_idx, cluster_label, ab_bed_dict.get(start, (None, None))[1]) #get [1] from (signal, a/b label)
            for start, (set_idx, cluster_label) in self.nodeset_attrs.items()
        }
        #get the confusion matrix between cluster_labels as predicted and ab_labels as ground truth
        cluster_labels = np.array([cluster_label for _, (_, cluster_label, _) in self.nodeset_attrs.items()])
        ab_labels = np.array([ab_label for _, (_, _, ab_label) in self.nodeset_attrs.items()])
        mapped_ab_labels = np.where(ab_labels == 'A', 0, 1) #map A to 0 and B to 1
        confusion_matrix = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):
                confusion_matrix[i, j] = np.sum((cluster_labels == i) & (mapped_ab_labels == j))
        return confusion_matrix 

def single_chrom_clustering(chrom, config, res, res_str):
    """perform spectral clustering on single intra-chromosomal graph"""
    modules = Cluster(config, chrom, res, res_str, n_clusters=2)
    modules.spectral_clustering()
    modules.append_cluster_labels_to_nodeset()
    #modules.create_gexf()
    conf_mtx = modules.oe_confusion_matrix()
    return conf_mtx

def run_parallel_clustering(config):
    """run spectral clustering on all chromosomes in parallel"""
    with Pool() as pool:
        pool.map(
            partial(
                single_chrom_clustering,
                config=config,
                res = config.genomic_params.resolutions[0],
                res_str=config.genomic_params.res_strs[0],
            ),
            config.genomic_params.chromosomes,
        )


if __name__ == "__main__":

    config = Config()
    #single_chrom_clustering("chr1", config, 1000000, "1Mb")
    run_parallel_clustering(config)
