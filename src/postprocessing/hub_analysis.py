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
def overlap_TSS_cluster_node():
    """ overlap TSS bin loci to nodeset locis and store their split per cluster"""
    #search nodeset_attrs starts to overlap gene TSS and gene body bins
    pass
