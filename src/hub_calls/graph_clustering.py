### Creates and clusters a graph object using called edges and nodes ###
# save the graph states as npy files (adj matrix/list) and node ID to bin name mapping as a csv

# create a class to make a graph object and another spectral clustering class to find modules from it

import networkx as nx
import numpy as np
from tqdm import tqdm
from sklearn.cluster import KMeans


class DataGenerator:
    def __init__(
        self,
        edges,
        edge_weight,
        batch_size,
        num_batch_per_iter,
        min_size=2,
        max_size=2,
        flag=False,
    ):
        self.edges = [[] for i in range(max_size + 1)]
        self.edge_weight = [[] for i in range(max_size + 1)]

        for e, ew in zip(edges, edge_weight):
            self.edges[len(e)].append(e)
            self.edge_weight[len(e)].append(ew)

        self.batch_size = batch_size
        self.num_batch_per_iter = num_batch_per_iter
        self.min_size = min_size
        self.max_size = max_size
        self.flag = flag
        for i in range(min_size, max_size + 1):
            self.edges[i] = np.array(self.edges[i])
            self.edge_weight[i] = np.array(self.edge_weight[i])

            while len(self.edges[i]) <= self.num_batch_per_iter * self.batch_size:
                self.edges[i] = np.concatenate([self.edges[i], self.edges[i]])
                self.edge_weight[i] = np.concatenate(
                    [self.edge_weight[i], self.edge_weight[i]]
                )
            self.shuffle(i)
        self.pointer = np.zeros((len(self.edges)), dtype="int")

    def shuffle(self, i):
        if self.flag:
            print("reach end, shuffling")
        index = np.random.permutation(len(self.edges[i]))
        self.edges[i] = (self.edges[i])[index]
        self.edge_weight[i] = (self.edge_weight[i])[index]

    def next_iter(self):

        return_edges = []
        return_edge_weight = []

        for i in range(self.min_size, self.max_size + 1):
            self.pointer[i] += self.num_batch_per_iter * self.batch_size

            if self.pointer[i] <= len(self.edges[i]):
                index = range(
                    self.pointer[i] - self.num_batch_per_iter * self.batch_size,
                    min(self.pointer[i], len(self.edges[i])),
                )
                edges = (self.edges[i])[index]
                edge_weight = (self.edge_weight[i])[index]
                return_edges += list(edges)
                return_edge_weight += list(edge_weight)
            else:
                index = range(
                    self.pointer[i] - self.num_batch_per_iter * self.batch_size,
                    min(self.pointer[i], len(self.edges[i])),
                )
                edges = (self.edges[i])[index]
                edge_weight = (self.edge_weight[i])[index]

                self.shuffle(i)
                left = self.num_batch_per_iter * self.batch_size - len(index)
                self.pointer[i] = 0
                self.pointer[i] += left
                index = range(0, self.pointer[i])
                edges, edge_weight = np.concatenate(
                    [edges, self.edges[i][index]]
                ), np.concatenate([edge_weight, self.edge_weight[i][index]])
                return_edges += list(edges)
                return_edge_weight += list(edge_weight)
        return np.asarray(return_edges), np.asarray(return_edge_weight)

    def balance_num(self, edges):
        cell = edges[:, 0]
        final = []
        choice, counts_ = np.unique(cell, return_counts=True)
        # num = int(np.mean(counts_))
        num = 50
        for c in tqdm(choice):
            final.append(np.random.choice(np.where(cell == c)[0], num, replace=True))
        final = np.concatenate(final, axis=-1)
        return final


class GraphGenerator:
    def __init__(self, config):
        self.local_edges = config.paths.edgelist[0] #pass edgelist of a specific nodeset
        self.global_edges = config.paths.edgelist[1] #pass edgelist of a specific nodeset

    def assign_local_nodes_to_global(self):
        """
        Assign local nodes to global nodes based on overlapping loop anchors to O/E bins 
        """
        pass

    def distribute_global_edges(self):
        """
        Distribute global edges to small scale local nodes
        """
        pass

    def normalize_edges(self):
        """
        Normalize edge weights between local and global edges
        """
        pass

    def csr_to_nx_graph(self):
        """
        Convert a csr matrix to a networkx graph object for getting gexf for viz
        """
        pass



class SpectralClustering:
    def __init__(self, graph):
        self.graph = graph
    
    def find_optimal_number_of_clusters(self):
        A = nx.adjacency_matrix(self.graph)
        L = nx.laplacian_matrix(self.graph)
        eigvals, eigvecs = np.linalg.eig(L.toarray())
        eigvals = eigvals[eigvals.argsort()]
        return eigvals

    def kmeans_cluster(self, num_clusters):
        A = nx.adjacency_matrix(self.graph)
        L = nx.laplacian_matrix(self.graph)
        eigvals, eigvecs = np.linalg.eig(L.toarray())
        eigvecs = eigvecs[:, eigvals.argsort()]
        kmeans = KMeans(num_clusters)
        kmeans.fit(eigvecs[:, 1 : num_clusters + 1])
        return kmeans.labels_


