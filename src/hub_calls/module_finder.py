### Creates and clusters a graph object using called edges and nodes ###

#create a class to make a graphg object and another spectral clustering class to find modules from it

import networkx as nx
import numpy as np
from sklearn.cluster import KMeans

class Graph:
    def __init__(self, edges):
        self.edges = edges
        self.graph = self.make_graph()
        
    def make_graph(self):
        G = nx.Graph()
        G.add_edges_from(self.edges)
        return G
    
    def get_graph(self):
        return self.graph
    

class SpectralClustering:
    def __init__(self, graph, n_clusters):
        self.graph = graph
        self.n_clusters = n_clusters
        
    def kmeans_cluster(self):
        A = nx.adjacency_matrix(self.graph)
        L = nx.laplacian_matrix(self.graph)
        eigvals, eigvecs = np.linalg.eig(L.toarray())
        eigvecs = eigvecs[:, eigvals.argsort()]
        kmeans = KMeans(n_clusters=self.n_clusters)
        kmeans.fit(eigvecs[:, 1:self.n_clusters+1])
        return kmeans.labels_
    
   