######### EGDELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import sys
from functools import partial
from multiprocessing import Pool

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
    """building csr matrix graph data object from h5 edglelist file"""
    def __init__(self, config, current_chrom, current_res_str, nodeset_key='oe_intra_0'):
        self.edgelist_h5_infile = config.paths.edgelist_outfile #outfile from dataloader
        self.query_group_chrom = current_chrom #current chrom to build graph for
        self.query_subgroup_res = f'_{current_res_str}' 
        self.query_key_edge = nodeset_key
 
        self.edge_df = None #extracted df for quick access
        self.node_idx_loci = None #idx mapped from df rows to csr rows, dict: {idx: chrN:start} for each node
        self.csr_matrix = None #graph matrix

    def load_edges(self):
        """extract pandas dataframe dataset for building graph from h5"""
        dataset_path = f'{self.query_group_chrom}/{self.query_subgroup_res}/{self.query_key_edge}'
        with pd.HDFStore(self.edgelist_h5_infile, mode='r') as store:
            self.edge_df = store[dataset_path]
    
    def df_to_affinity_csr(self):
        """convert edgelist pandas df graph to a CSR matrix graph (affinity matrix) for clustering"""
        nodes = pd.unique(self.edge_df[['x1', 'y1']].values.ravel('K')) 
        #TODO: match with number of bins in the chrom

        self.node_idx_loci = {idx: start for idx, start in enumerate(nodes)} #save in format idx: chrN:start 
        rows = self.edge_df['x1'].map(self.node_idx_loci) #mapping csr row idx to df row idx of bin starts
        cols = self.edge_df['y1'].map(self.node_idx_loci) #mapping csr row idx to df row idx of bin starts
        weights = self.edge_df['counts'].values
        num_nodes = len(nodes)
        self.csr_matrix = csr_matrix((weights, (rows, cols)), shape=(num_nodes, num_nodes))


#module detection class
class Cluster(Graph):
    """
    spectral clustering on intra chr graph built from 
    a) single scale OE edges 
    b) single scale loop edges 
    c) multi scale OE (global) + loop (local) edges
    """

    def __init__(self, config, chrom, current_res_str):
        super().__init__(config, chrom, current_res_str)
        self.clusters_csr = None
        self.cluster_labels = None
        self.intra_chrom_oe_csr = None #global oe csr matrix graph
        self.intra_chrom_loop_csr = None #local loop csr matrix graph
        self.intra_chrom_multi_scale_csr = None #multiscale oe + loop csr matrix graph

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
    #inspect_h5_file(config.paths.edgelist_outfile) #utils func
    graph = Graph(config, config.genomic_params.chromosomes[0], config.genomic_params.res_strs[0])
    graph.load_edges()
    graph.df_to_csr()