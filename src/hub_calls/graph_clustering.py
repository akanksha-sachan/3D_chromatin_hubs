######### EGDELIST (HDF5) -> clusters of graph #########

# pylint: disable=all

import os
import sys
from functools import partial
from multiprocessing import Pool

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from memory_profiler import profile
from scipy.sparse import csr_matrix
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering as spectral

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# relative and absolute imports for running module and script respectively
# to test functions of this script inside this script itself, without NB
try:
    from ..configs.config_local import Config
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from utils import *
