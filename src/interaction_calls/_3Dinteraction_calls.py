import os
import sys
import numpy as np
import pandas as pd
import cooler
from sklearn.decomposition import PCA

##### 3D loop calls #####

class HICCUPSLoops:
    """
    HICCUPS loop calls
    """
    def __init__(self, config):
        self.config = config
        self.temp_dir = config["temp_dir"]
        self.hiccups_dir = config["hiccups_dir"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)
        

    def juicer_hiccups(self, chrom):
        """
        Call loops using HICCUPS wrapper
        """
        pass

class PeakachuLoops:
    """
    Peakachu loop calls
    """
    def __init__(self, config):
        self.config = config
        self.temp_dir = config["temp_dir"]
        self.peakachu_dir = config["peakachu_dir"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)
        

    def call_loops(self, chrom):
        """
        Call loops using Peakachu wrapper
        """
        pass

class LoadLoops:
    """
    Load pre-existing loops 
    """
    def __init__(self, config):
        self.config = config
        self.temp_dir = config["temp_dir"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)
        

    def load_loops(self, chrom):
        """
        Load pre-existing loops
        """
        pass

##### 3D A/B compartment calls #####

def get_ab_compartment(matrix):
    """
    Raw -> normalised -> O/E -> Pearson -> PCA gives A/B
    """
    
    #ensure diagonals are 1
    np.fill_diagonal(matrix, 1)
    #get pearson matrix
    matrix = np.corrcoef(matrix)
    np.fill_diagonal(matrix, 1)
    matrix[np.isnan(matrix)] = 0.0
    pca = PCA(n_components=1)
    y = pca.fit_transform(matrix)
    return y, pca

##### 3D insulation score #####

def insulation_score(m, windowsize=500000, res=10000):
    """
    needs <10kb bin size to detect tads using diamond score method and >100M total filtered reads mapped to the genome
    """
    windowsize_bin = int(windowsize / res)
    m = np.nan_to_num(m)
    score = np.ones((m.shape[0]))
    for i in range(0, m.shape[0]):
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.sum(m[max(0, i - windowsize_bin): i, i + 1: min(m.shape[0] - 1, i + windowsize_bin + 1)]) / (np.sum(
                 m[max(0, i - windowsize_bin):min(m.shape[0], i + windowsize_bin + 1),
                   max(0, i - windowsize_bin):min(m.shape[0], i + windowsize_bin + 1)]))
            if np.isnan(v):
                v = 1.0
        score[i] = v
		#get log2 of the score
		#score[score == 0] = 1
		#score = np.log2(score)
    return score