######### EGDELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import sys
from functools import partial
from multiprocessing import Pool

import h5py
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
from memory_profiler import profile
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from ..configs.config_local import Config
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from utils import *

#graph constructor class
class Graph:
    """intra chr graph built from a) single scale OE edges b) single scale loop edges c) multi scale OE (global) + loop (local) edges"""
    def __init__(self, config, chrom, current_res_str):
        self.edgelist_h5 = config.paths.edgelist_outfile #infile from dataloader
        self.chrom_group = chrom #current chrom to build graph of
        self.res_str_key = current_res_str #key for getting resolution specific edgelist

        self.intra_chrom_oe_csr = None #global oe csr matrix graph
        self.intra_chrom_loop_csr = None #local loop csr matrix graph
        self.intra_chrom_multi_scale_csr = None #multiscale oe + loop csr matrix graph
        
        self.load_h5() #new attrs added in this function
        self.node_indices = None

    def load_h5(self):
        """extract the pandas dataframe from the h5 file"""
        with h5py.File(self.edgelist_h5, 'r') as f:
            self.edgelist_df = f[self.chrom_group][self.res_str_key][:]
    
    def convert_to_csr(self):
        nodes = pd.unique(self.edgelist_df[['start', 'end']].values.ravel('K'))
        self.node_indices = {node: idx for idx, node in enumerate(nodes)}
        rows = self.edgelist_df['start'].map(self.node_indices)
        cols = self.edgelist_df['end'].map(self.node_indices)
        weights = self.edgelist_df['counts'].values
        num_nodes = len(nodes)
        self.csr_matrix = csr_matrix((weights, (rows, cols)), shape=(num_nodes, num_nodes))

#module detection class
class Cluster(Graph):
    def __init__(self, config, chrom, current_res_str):
        super().__init__(config, chrom, current_res_str)
        self.clusters_csr = None
        self.cluster_labels = None

    def perform_spectral_clustering(self, csr_matrix, n_clusters):
        spectral = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', assign_labels='discretize')
        self.cluster_labels = spectral.fit_predict(csr_matrix)
    
    def store_clusters(self, node_indices, cluster_labels, output_path):
        G = nx.Graph()
        reverse_node_indices = {v: k for k, v in node_indices.items()}
        
        for node, cluster in zip(reverse_node_indices.keys(), cluster_labels):
            G.add_node(reverse_node_indices[node], cluster=cluster)
        
        for i, j, weight in zip(*csr_matrix.nonzero(), csr_matrix.data):
            G.add_edge(reverse_node_indices[i], reverse_node_indices[j], weight=weight)
        
        nx.write_gexf(G, output_path)

# Example usage:
# graph = Graph()
# graph.load_h5('path_to_file.h5', 'chromosome_group')
# graph.convert_to_csr()

# cluster = Cluster()
# cluster.perform_spectral_clustering(graph.csr_matrix, n_clusters=5)
# cluster.store_clusters(graph.node_indices, cluster.cluster_labels, 'output_file.gexf')

if __name__ == "__main__":

    config = Config()
    inspect_h5_file(config.paths.edgelist_outfile)