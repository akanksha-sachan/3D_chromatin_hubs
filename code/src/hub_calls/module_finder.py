import community as community_louvain
import networkx as nx
import numpy as np
import pandas as pd

# implement spectral clustering


class Graph:
    def __init__(
        self, edge_list=None, matrix=None, projection=None, sample="", parent=None
    ):
        """
        Initializes a Graph object with an edge list or a connectivity matrix and a projection
        """
        if edge_list is not None:
            self.edge2matrix(edge_list)  # Converts edge list to connectivity matrix
            self.projection = np.arange(self.connect.shape[0])
            print("Connection Initialized")
        elif matrix is not None and projection is not None:
            self.connect = matrix
            self.projection = projection
        else:
            raise ValueError(
                "Either provide an edge list or a connectivity matrix with projection."
            )

        self.parent = parent
        self.sample = sample
        self.initialize_inheritance()

    def initialize_inheritance(self):
        """
        Initializes the inheritance of the parent graph
        """
        # Try inheriting attribute (bin_indices {interval_table}, etc) from parent
        try:
            self.interval_table = self.parent.interval_table
        except AttributeError:
            # If no parent or parent lacks interval_table, initialize for root graph
            print("Root Graph")
            if self.cell:
                # Load interval table specific to the sample from a file
                self.interval_table = pd.read_table(
                    f"../Temp/{self.sample}/chrom_bin_index.txt", sep="\t"
                )
            else:
                # If no sample specified, use an empty DataFrame
                self.interval_table = pd.DataFrame()
