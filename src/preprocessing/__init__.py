# __init__.py for preprocessing sub-module
## sub-module input is - is .hic/.rna/.atac files (could even be BAM/.fa files) :: output is - processed data in the format for downstream feature selection
# data formats available for bulk HiC: .fa gives BAM which gives .pairs > .hic (Juicer), .mcool (Cooler), .matrix (HiCPro) all have multi-res binned interactions 
# adapt FAN-C code for making input compatible with juicer and cooler visualization and loop calling tools (hiccups)


from .genes_processing import *
from ._2Dcontacts_processing import *
from ._1Dchromatin_processing import *
from .interaction_validation import *
