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
from sklearn.metrics import f1_score, precision_score, recall_score
import pickle

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
        self, config, current_chrom, current_res, current_res_str, nodeset_key
    ):
        #basic params init
        self.config = config
        self.chrom = current_chrom
        self.current_res = current_res
        self.current_res_str = current_res_str

        #query data_loader object init
        self.edgelist_h5_infile = os.path.join(
            self.config.paths.edgelist_outdir, f"{self.chrom}.h5"
        )
        self.query_group_res = f"_{current_res_str}"
        self.query_key_edge = nodeset_key #oe_intra_0, etc, other edge types
        self.bins_bed = config.paths.ref_genome_bins

        #graph data attrs init
        self.edge_df = None  # extracted df for quick access
        self.nodeset_attrs = None  # {start:(idx, attrs)} : start indicates the genomic loci of the node, and idx is nodeset order
        self.affinity_matrix = None  # affinity matrix to represent constructed graph
        self.affinity_key = config.genomic_params.affinity_key
        self.affinity_plot_dir = os.path.join(config.paths.temp_dir, f"affinity_matrices/{self.affinity_key}")
        os.makedirs(self.affinity_plot_dir, exist_ok=True) #create if not exists

    def load_edges(self):
        """extract pandas dataframe dataset for building graph from h5 per chr"""
        dataset_path = (
            f"{self.query_group_res}/{self.query_key_edge}"
        )
        with pd.HDFStore(self.edgelist_h5_infile, mode="r") as store:
            self.edge_df = store[dataset_path]

    def apply_genomic_distance_filter(self, distance_threshold):
        """apply genomic distance filter to the edgelist df to select only long range interactions"""
        pass

    def edgelist_to_csr_affinity_matrix(self):
        """convert edgelist pandas df graph to a CSR matrix graph (affinity matrix) for clustering
                "chr1": chr_column,
                "x1": x1,
                "x2": x2,
                "chr2": chr_column,
                "y1": y1,
                "y2": y2,
                "counts": counts
        """
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
        # save affinity matrix to disk as pkl for plotting input 
        outfile = f"{self.affinity_plot_dir}/{self.chrom}_affinity.pkl"
        with open(outfile, 'wb') as file:
            pickle.dump(self.affinity_matrix, file) 


class SpectralCluster(Graph):
    """
    spectral clustering class using scikit-learn's SpectralClustering implementation
    """

    def __init__(self, config, chrom, current_res, current_res_str, nodeset_key, n_clusters=2):
        super().__init__(config, chrom, current_res, current_res_str, nodeset_key)
        self.load_edges()  # call function from parent to load edges
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
        #save the nodeset_attrs to disk as pkl in the same dir as affinity matrix pkl
        outfile = f"{self.affinity_plot_dir}/{self.chrom}_nodeset_attrs.pkl"
        with open(outfile, 'wb') as file:
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
        outfile = (
            f"{self.graphs_outdir}/{self.chrom}_{self.query_key_edge}.gexf"
        )
        nx.write_gexf(G, outfile)

    class ab_evaluation:
        """class to calculate accuracy metrics for A/B clustering"""

        def __init__(self, parent):
            self.parent = parent
            self.metrics_csv_file = os.path.join(parent.config.paths.temp_dir, f"accuracy_metrics_{self.parent.affinity_key}.csv")
    
        def oe_confusion_matrix(self, nodeset_attrs=None):
            ## Not applicable for hubs clustering ##
            """calculate confusion matrix for clustering OE edges using AB compartments as ground truth"""
            query = HiCQuery(self.parent.config, self.parent.chrom, self.parent.current_res, self.parent.current_res_str)
            ab_bed_dict = query.ab_comp.load_bigwig_chromosomal_ab()
            # append the A/B labels to the nodeset_attrs dict by mapping 'start' key to get {start: (set_idx, cluster_label, ab_label)}
            
            if self.parent.nodeset_attrs is None:
                nodeset_attrs = nodeset_attrs
            else:
                nodeset_attrs = self.parent.nodeset_attrs

            nodeset_attrs = {
                start: (set_idx, cluster_label, ab_bed_dict.get(start, (None, None))[1]) #get [1] from (signal, a/b label)
                for start, (set_idx, cluster_label) in nodeset_attrs.items()
            }
            #get the confusion matrix between cluster_labels as predicted and ab_labels as ground truth
            cluster_labels = np.array([cluster_label for _, (_, cluster_label, _) in nodeset_attrs.items()])
            ab_labels = np.array([ab_label for _, (_, _, ab_label) in nodeset_attrs.items()])
            mapped_ab_labels = np.where(ab_labels == 'A', 0, 1) #map A to 0 and B to 1
            self.confusion_matrix = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    self.confusion_matrix[i, j] = np.sum((cluster_labels == i) & (mapped_ab_labels == j))
            return self.confusion_matrix 
        
        def accuracy_metrics_single_chr(self):
            ## Not applicable for hubs clustering ##
            """calculate accuracy metrics (F1 score from the confusion matrix) for single chromosomal clustering"""
            confusion_matrix = self.confusion_matrix.astype(int)
            true_labels = np.concatenate([[0] * confusion_matrix[0, 0], [1] * confusion_matrix[0, 1],
                                        [0] * confusion_matrix[1, 0], [1] * confusion_matrix[1, 1]])
            pred_labels = np.concatenate([[0] * confusion_matrix[0, 0], [0] * confusion_matrix[0, 1],
                                        [1] * confusion_matrix[1, 0], [1] * confusion_matrix[1, 1]])

            precision = precision_score(true_labels, pred_labels)
            recall = recall_score(true_labels, pred_labels)
            f1 = f1_score(true_labels, pred_labels)

            accuracy_metrics_dict = {
            'chrom': self.parent.chrom,
            'precision': precision,
            'recall': recall,
            'f1': f1
        }

            df = pd.DataFrame([accuracy_metrics_dict])
            df.to_csv(self.metrics_csv_file, mode='a', header=not os.path.exists(self.metrics_csv_file), index=False)

            return (precision, recall, f1)


class NormCut(Graph):
    """
    Implementation of recursive bipartitioning using Normalized Cut for finding modules
    
    Reference
    ---------
    Jianbo Shi and Jitendra Malik. Normalized Cuts and Image Segmentation.
    IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
    VOL. 22, NO. 8, AUGUST 2000.
    """

    def __init__(self, config, chrom, current_res, current_res_str, nodeset_key, n_clusters=2):
        super().__init__(config, chrom, current_res, current_res_str, nodeset_key)
        self.load_edges()  # call function from parent to load edges
        self.edgelist_to_csr_affinity_matrix() #call function from parent to load affinity matrix from edges
        self.cluster_labels = None
        self.graphs_outdir = config.paths.gexf_dir

        # instantiate nested classes
        self.evaluation = self.evaluation(self)  

    class evaluation:
        """class to calculate accuracy metrics for NormCut using stability and scores"""

        def __init__(self, parent):
            self.parent = parent

def run_single_chrom(chrom, config, res, res_str, nodeset_key):
    """perform spectral clustering on single intra-chromosomal graph"""
    modules = SpectralCluster(config, chrom, res, res_str, nodeset_key, n_clusters=2) #initializes graph object as well
    modules.spectral_clustering() 
    modules.save_nodeset_attrs() #gives nodeset_attrs as pkl in the affinity matrix dir
    modules.create_gexf() #gives graph objects

def run_parallel(config):
    """run spectral clustering on all chromosomes in parallel"""
    with Pool() as pool:
        pool.map(
            partial(
                run_single_chrom,
                config=config,
                res = config.current_res,
                res_str=config.current_res_str,
                nodeset_key=config.genomic_params.nodeset_key
            ),
            config.genomic_params.chromosomes,
        )

def run_single_chrom_eval(chrom, config, res, res_str, nodeset_key):
    """use for loop for evaluation, multiprocessing this gives errors"""
    modules = SpectralCluster(config, chrom, res, res_str, nodeset_key, n_clusters=2) #initializes graph object as well
    modules.spectral_clustering() 
    modules.save_nodeset_attrs() #gives nodeset_attrs as pkl in the affinity matrix dir
    conf_mtx = modules.ab_evaluation.oe_confusion_matrix()
    accuracy_metrics_tuple = modules.ab_evaluation.accuracy_metrics_single_chr()
    return conf_mtx, accuracy_metrics_tuple

def whole_genome_evaluation(config):
    """calculate evaluation metrics for whole genome together"""
    for chrom in config.genomic_params.chromosomes:
        _ , accuracy_metrics_tuple = run_single_chrom_eval(chrom, config, config.current_res, config.current_res_str, config.genomic_params.nodeset_key)
        print(f"Chrom: {chrom}, Accuracy Metrics: {accuracy_metrics_tuple}")

if __name__ == "__main__":
    config = Config()
   
