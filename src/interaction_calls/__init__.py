# __init__.py for extracting edges from HiC and nodes from annotations
# edge_extraction + validation of edges, node creation has all types of nodes (currently gene)
#sub-module: input is .hic/.cool files to process and extract edges from; bin annotations of genes and REs coming from atac 

from ._3Dinteraction_calls import *
from .node_annotation import *
from .interaction_validation import *
