######### EDGELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import pickle
import sys
from functools import partial
from multiprocessing import Pool

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, dia_matrix
from sklearn.cluster import SpectralClustering
from sklearn.decomposition import PCA
from sklearn.metrics import f1_score, precision_score, recall_score

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


class HiCGraph:
    """graph constructor class: building csr matrix graph data object from edglelists"""

    def __init__(
        self, config, current_chrom, current_res, current_res_str, nodeset_key
    ):
        # basic params init
        self.config = config
        self.chrom = current_chrom
        self.current_res = current_res
        self.current_res_str = current_res_str

        # query data_loader object init
        self.edgelist_h5_infile = os.path.join(
            self.config.paths.edgelist_outdir, f"{self.chrom}.h5"
        )
        self.query_group_res = f"_{current_res_str}"
        self.query_key_edge = nodeset_key  # oe_intra_0, etc, other edge types
        self.bins_bed = config.paths.ref_genome_bins

        # graph data attrs init
        self.edge_df = None  # extracted df for quick access
        self.nodeset_attrs = None  # {start:(idx, attrs)} : start indicates the genomic loci of the node, and idx is nodeset order
        self.affinity_matrix = None  # affinity matrix to represent constructed graph
        self.affinity_key = config.affinity_key
        self.affinity_plot_dir = os.path.join(
            config.paths.temp_dir, f"affinity_matrices/{self.affinity_key}"
        )
        os.makedirs(self.affinity_plot_dir, exist_ok=True)  # create if not exists

        # long-range window filter parameters
        self.long_range_window_width = config.long_range_window_width
        self.long_range_min_distance = config.long_range_min_distance

        # set chromosome length
        self.chromosome_length = self.get_current_chromosome_length()

    def get_current_chromosome_length(self):
        """Read chromosome length from the ref_genome.chrom.sizes file."""
        chrom_sizes_infile = self.config.paths.chrom_sizes_infile
        chrom_sizes = pd.read_csv(
            chrom_sizes_infile, sep="\t", header=None, names=["chrom", "size"]
        )  # file has no header
        chrom_length = chrom_sizes.loc[
            chrom_sizes["chrom"] == self.chrom, "size"
        ].values[0]
        return chrom_length

    def load_edges(self):
        """extract pandas dataframe dataset for building graph from h5 per chr

        Returns: "chr1": chr_column,
                "x1": x1,
                "x2": x2,
                "chr2": chr_column,
                "y1": y1,
                "y2": y2,
                "counts": counts
        """
        dataset_path = f"{self.query_group_res}/{self.query_key_edge}"
        with pd.HDFStore(self.edgelist_h5_infile, mode="r") as store:
            self.edge_df = store[dataset_path]

    def apply_genomic_distance_filter(self, current_filter_start):
        """apply genomic distance filter to the edgelist df to select only long range interactions"""
        window_width = self.long_range_window_width
        if current_filter_start + window_width > self.chromosome_length:
            raise ValueError(
                "The start position and window width exceed the chromosome length."
            )
        end_distance = current_filter_start + window_width
        filtered_df = self.edge_df[
            (abs(self.edge_df["x1"] - self.edge_df["y1"]) >= current_filter_start)
            & (abs(self.edge_df["x1"] - self.edge_df["y1"]) < end_distance)
        ]
        self.edge_df = filtered_df

    def edgelist_to_csr_affinity_matrix(self):
        """convert edgelist pandas df graph to a CSR matrix graph (affinity matrix) for clustering
        Returns: affinity_matrix: csr_matrix (W); nodeset_attrs: {start: idx} dict
        """
        # (loci start: idx of nodeset) dict for adding nodes
        nodeset = pd.unique(self.edge_df[["x1", "y1"]].values.ravel("K"))
        self.nodeset_attrs = {start: set_idx for set_idx, start in enumerate(nodeset)}
        # assigning nodeset id to edgelist nodes using start
        rows_node_id = self.edge_df["x1"].map(self.nodeset_attrs)
        cols_node_id = self.edge_df["y1"].map(self.nodeset_attrs)
        weights = self.edge_df["counts"].values  # edge weights
        # create affinity matrix of only current edges and not whole chr
        num_nodes = len(nodeset)
        affinity_matrix = csr_matrix(
            (weights, (rows_node_id, cols_node_id)), shape=(num_nodes, num_nodes)
        )
        # make the affinity matrix symmetric by adding the transpose of it
        self.affinity_matrix = (
            affinity_matrix
            + affinity_matrix.T
            - csr_matrix(
                (affinity_matrix.diagonal(), (range(num_nodes), range(num_nodes))),
                shape=(num_nodes, num_nodes),
            )
        )
        return self.affinity_matrix, self.nodeset_attrs

    def compute_diagonal_matrix(self):
        """
        Compute the diagonal matrix D for the CSR affinity matrix.
        Returns
        D : dia_matrix; The diagonal matrix of the graph. D[i, i] is the sum of weights of all edges incident on i
        """
        if self.affinity_matrix is None:
            raise ValueError("Affinity matrix has not been computed yet.")
        # sum of weights of all edges (all columns) incident on each node (row)
        degrees = np.array(self.affinity_matrix.sum(axis=1)).flatten()
        D = dia_matrix((degrees, [0]), shape=self.affinity_matrix.shape)
        return D

    def save_thresholded_graph_to_neo4j_csv(self):
        """save the thresholded edgelist to disk for plotting in neo4j"""
        # create nodes dataframe
        nodes = [{"id": idx, "label": loci} for loci, idx in self.nodeset_attrs.items()]
        nodes_df = pd.DataFrame(nodes)
        # create edges dataframe
        row, col = self.affinity_matrix.nonzero()
        edges = [
            {
                "source": row[i],
                "target": col[i],
                "weight": self.affinity_matrix[row[i], col[i]],
            }
            for i in range(len(row))
        ]
        edges_df = pd.DataFrame(edges)
        # save nodes and edges to CSV
        nodes_csv_path = f"{self.affinity_plot_dir}/{self.chrom}_nodes.csv"
        edges_csv_path = f"{self.affinity_plot_dir}/{self.chrom}_edges.csv"
        nodes_df.to_csv(nodes_csv_path, index=False)
        edges_df.to_csv(edges_csv_path, index=False)

        print(f"Nodes CSV saved to {nodes_csv_path}")
        print(f"Edges CSV saved to {edges_csv_path}")

    def save_graph_as_gexf(self):
        """save the graph as gexf for visualization
        Params: nodeset_attrs: {start: idx} dict for adding nodes
                edge_df: edgelist df for adding edge weights
        """
        G = nx.Graph()
        # add nodes from nodeset
        for start, set_idx in self.nodeset_attrs.items():
            node_label = f"{self.chrom}:{start}"
            G.add_node(set_idx, label=node_label)
        # add edges from edgelist, map nodes in edgelist to nodes in the nx graph using start
        for _, row in self.edge_df.iterrows():
            start_x, start_y = row["x1"], row["y1"]
            start_x_set_idx, start_y_set_idx = (
                self.nodeset_attrs[start_x],
                self.nodeset_attrs[start_y],
            )
            G.add_edge(start_x_set_idx, start_y_set_idx, weight=row["counts"])
        # write to gexf
        dist_thresh = format_loci_string(
            config.min_distance_threshold, config.max_distance_threshold, "1Mb"
        )  # get string for dist thresh
        outfile = f"{self.affinity_plot_dir}/{self.chrom}_{dist_thresh}_affinity.gexf"
        nx.write_gexf(G, outfile)


class SpectralCluster(HiCGraph):
    """
    spectral clustering class using scikit-learn's SpectralClustering implementation
    Usage: two partitioning HiC graph to see if compartments are clustered
    """

    def __init__(
        self, config, chrom, current_res, current_res_str, nodeset_key, n_clusters=2
    ):
        super().__init__(config, chrom, current_res, current_res_str, nodeset_key)
        self.load_edges()  # call function from parent to load edges
        self.apply_genomic_distance_filter()  # call function from parent to filter edges
        self.edgelist_to_csr_affinity_matrix()
        self.number_of_clusters = n_clusters
        self.cluster_labels = None
        self.graphs_outdir = config.paths.gexf_dir

        # instantiate nested classes
        self.ab_evaluation = self.ab_evaluation(self)

    def spectral_clustering(self):
        """perform spectral clustering on the affinity matrix, add cluster labels to nodeset_attrs dict"""
        affinity_matrix = self.affinity_matrix
        spectral = SpectralClustering(
            n_clusters=self.number_of_clusters,
            affinity="precomputed",
            assign_labels="discretize",
        )  # cluster_labels coming from image segmentation algo
        self.cluster_labels = spectral.fit_predict(affinity_matrix)
        nodeset_dict = self.nodeset_attrs  # {start:set_idx}
        # append cluster labels as {start: (set_idx, cluster_label)}
        self.nodeset_attrs = {
            start: (set_idx, self.cluster_labels[set_idx])
            for start, set_idx in nodeset_dict.items()
        }

    def save_nodeset_attrs(self):
        # save the nodeset_attrs to disk as pkl in the same dir as affinity matrix pkl
        outfile = f"{self.affinity_plot_dir}/{self.chrom}_nodeset_attrs.pkl"
        with open(outfile, "wb") as file:
            pickle.dump(self.nodeset_attrs, file)

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
        outfile = f"{self.graphs_outdir}/{self.chrom}_{self.query_key_edge}.gexf"
        nx.write_gexf(G, outfile)

    class ab_evaluation:
        """class to calculate accuracy metrics for A/B clustering"""

        def __init__(self, parent):
            self.parent = parent
            self.metrics_csv_file = os.path.join(
                parent.config.paths.temp_dir,
                f"accuracy_metrics_{self.parent.affinity_key}.csv",
            )

        def oe_confusion_matrix(self, nodeset_attrs=None):
            ## Not applicable for hubs clustering ##
            """calculate confusion matrix for clustering OE edges using AB compartments as ground truth"""
            query = HiCQuery(
                self.parent.config,
                self.parent.chrom,
                self.parent.current_res,
                self.parent.current_res_str,
            )
            ab_bed_dict = query.ab_comp.load_bigwig_chromosomal_ab()
            # append the A/B labels to the nodeset_attrs dict by mapping 'start' key to get {start: (set_idx, cluster_label, ab_label)}

            if self.parent.nodeset_attrs is None:
                nodeset_attrs = nodeset_attrs
            else:
                nodeset_attrs = self.parent.nodeset_attrs

            nodeset_attrs = {
                start: (
                    set_idx,
                    cluster_label,
                    ab_bed_dict.get(start, (None, None))[1],
                )  # get [1] from (signal, a/b label)
                for start, (set_idx, cluster_label) in nodeset_attrs.items()
            }
            # get the confusion matrix between cluster_labels as predicted and ab_labels as ground truth
            cluster_labels = np.array(
                [cluster_label for _, (_, cluster_label, _) in nodeset_attrs.items()]
            )
            ab_labels = np.array(
                [ab_label for _, (_, _, ab_label) in nodeset_attrs.items()]
            )
            mapped_ab_labels = np.where(ab_labels == "A", 0, 1)  # map A to 0 and B to 1
            self.confusion_matrix = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    self.confusion_matrix[i, j] = np.sum(
                        (cluster_labels == i) & (mapped_ab_labels == j)
                    )
            return self.confusion_matrix

        def accuracy_metrics_single_chr(self):
            ## Not applicable for hubs clustering ##
            """calculate accuracy metrics (F1 score from the confusion matrix) for single chromosomal clustering"""
            confusion_matrix = self.confusion_matrix.astype(int)
            true_labels = np.concatenate(
                [
                    [0] * confusion_matrix[0, 0],
                    [1] * confusion_matrix[0, 1],
                    [0] * confusion_matrix[1, 0],
                    [1] * confusion_matrix[1, 1],
                ]
            )
            pred_labels = np.concatenate(
                [
                    [0] * confusion_matrix[0, 0],
                    [0] * confusion_matrix[0, 1],
                    [1] * confusion_matrix[1, 0],
                    [1] * confusion_matrix[1, 1],
                ]
            )
            precision = precision_score(true_labels, pred_labels)
            recall = recall_score(true_labels, pred_labels)
            f1 = f1_score(true_labels, pred_labels)
            accuracy_metrics_dict = {
                "chrom": self.parent.chrom,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            }
            #save to csv
            df = pd.DataFrame([accuracy_metrics_dict])
            df.to_csv(
                self.metrics_csv_file,
                mode="a",
                header=not os.path.exists(self.metrics_csv_file),
                index=False,
            )
            return (precision, recall, f1)


class RecursiveNormCut(HiCGraph):
    """
    Implementation of recursive bipartitioning (two-way) using Normalized Cut
    Usage: recursively partitioning the HiC graph to find hubs

    Reference
    ---------
    Jianbo Shi and Jitendra Malik. Normalized Cuts and Image Segmentation.
    IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
    VOL. 22, NO. 8, AUGUST 2000.
    """

    def __init__(
        self, config, chrom, current_res, current_res_str, nodeset_key, n_clusters=2
    ):
        ## loading HiC data into graph object ##
        super().__init__(config, chrom, current_res, current_res_str, nodeset_key)
        self.load_edges()  # call function from parent to load edges
        self.edgelist_to_csr_affinity_matrix()  # call function from parent to create W
        self.compute_diagonal_matrix()  # call function from parent to create D
        self.subgraph_labels = None
        self.graphs_outdir = (
            config.paths.gexf_dir
        )  

        # attributes for recursive norm cut 
        self.subgraph_list = []  # list to store subgraphs

        # instantiate nested classes
        self.evaluation = self.evaluation(self)

    def bipartition_norm_cut(self):
        """perform two way normalized cut/bipartitioning on the affinity matrix
        usage: base function for recursive bipartitioning
        """
        pass

    class evaluation:
        """class to calculate accuracy metrics for NormCut using stability and scores"""

        def __init__(self, parent):
            self.parent = parent


def run_single_chrom_spectral(chrom, config, res, res_str, nodeset_key):
    """perform spectral clustering on single intra-chromosomal graph"""
    modules = SpectralCluster(
        config, chrom, res, res_str, nodeset_key, n_clusters=2
    )  # initializes graph object as well
    modules.spectral_clustering()
    modules.save_nodeset_attrs()  # gives nodeset_attrs as pkl in the affinity matrix dir
    modules.create_gexf()  # gives graph objects


def run_parallel(config):
    """run spectral clustering on all chromosomes in parallel"""
    with Pool() as pool:
        pool.map(
            partial(
                run_single_chrom_spectral,
                config=config,
                res=config.current_res,
                res_str=config.current_res_str,
                nodeset_key=config.nodeset_key,
            ),
            config.param_lists.chromosomes,
        )


def run_single_chrom_ab_eval(chrom, config, res, res_str, nodeset_key):
    """use for loop for evaluation, multiprocessing this gives errors"""
    modules = SpectralCluster(
        config, chrom, res, res_str, nodeset_key, n_clusters=2
    )  # initializes graph object as well
    modules.spectral_clustering()
    modules.save_nodeset_attrs()  # gives nodeset_attrs as pkl in the affinity matrix dir
    conf_mtx = modules.ab_evaluation.oe_confusion_matrix()
    accuracy_metrics_tuple = modules.ab_evaluation.accuracy_metrics_single_chr()
    return conf_mtx, accuracy_metrics_tuple


def whole_genome_ab_eval(config):
    """calculate evaluation metrics for whole genome together"""
    for chrom in config.param_lists.chromosomes:
        _, accuracy_metrics_tuple = run_single_chrom_ab_eval(
            chrom,
            config,
            config.current_res,
            config.current_res_str,
            config.nodeset_key,
        )
        print(f"Chrom: {chrom}, Accuracy Metrics: {accuracy_metrics_tuple}")


if __name__ == "__main__":
    config = Config()
    chrom = config.param_lists.chromosomes[0]
    current_res = config.current_res
    current_res_str = config.current_res_str
    
