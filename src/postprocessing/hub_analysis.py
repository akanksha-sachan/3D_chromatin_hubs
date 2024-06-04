import sys
import os
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import seaborn as sns

# add the parent directory of 'src' to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from ..configs.config_local import Config
    from ..preprocessing import *
    from ..module_detection import *
    from ..utils import *
except ImportError:
    from configs.config_local import Config
    from preprocessing import *
    from module_detection import *
    from utils import *

## overlapping node annotations
def overlap_genes_to_nodeset():
    """ overlap TSS bin loci to nodeset locis and store their split per cluster"""
    #search nodeset_attrs starts to overlap gene TSS and gene body bins
    pass

def overlap_sub_compartments_to_nodeset():
    """ overlap subcompartments to nodeset locis ; works ideally with inter-chr nodesets"""
    #search nodeset_attrs starts to overlap subcompartments
    pass

## hub physical properties
def cluster_size_distribution():
        """calculate cluster size distribution split between number of clusters"""
        pass