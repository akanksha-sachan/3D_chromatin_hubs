######### CONTACT MATRIX (.hic/.mcool) PREPROCESSING #########
#use cooler for extracting contact matrix from .mcool file
#use pybind11 straw for extracting contact matrix from .hic file

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import hicstraw
from scipy.sparse import csr_matrix

#### .hic file pre-processing into sparse scipy and dense numpy contact matrices ####

def csr_contact_matrix(norm, hicfile, chr1loc, chr2loc, unit,
                       binsize, is_synapse=False):
    '''
    Extract a specific region of the contact matrix from .hic in scipy's CSR sparse format
    limitations: simplistic storage in rows and columns, not genomic intervals and their indexing
    '''

    tri_list = hicstraw.straw(norm, hicfile, chr1loc,
                           chr2loc, unit, binsize, is_synapse)
    row = [r//binsize for r in tri_list[0]]
    col = [c//binsize for c in tri_list[1]]
    value = tri_list[2]
    N = max(col) + 1

    # re-scale KR matrix to ICE-matrix range
    M = csr_matrix((value, (row, col)), shape=(N, N), dtype=float)
    margs = np.array(M.sum(axis=0)).ravel() + \
        np.array(M.sum(axis=1)).ravel() - M.diagonal(0)
    margs[np.isnan(margs)] = 0
    scale = margs[margs != 0].mean()
    row, col = M.nonzero()
    value = M.data / scale
    M = csr_matrix((value, (row, col)), shape=(N, N), dtype=float)

    return M

def straw_contact_matrix(hicfile, chrom1, chrom2, resolution, gr1, gr2, gc1, gc2, data_type="observed", normalization="KR"):
    """
    Creates a hic object from the HicFile class in straw. Extracts a dense numpy matrix from the hic object.
    """
    hic = hicstraw.HiCFile(hicfile)
    #hic.getChromosomes()
    #hic.getGenomeID()
    #hic.getResolutions()
    #data_type = 'observed' or 'oe'
    # if hic.getResolutions() has 10000, 5000
    mtx = hic.getMatrixZoomData(chrom1, chrom2, data_type, normalization, "BP", resolution) #MatrixZoomData object 

    numpy_matrix = mtx.getRecordsAsMatrix(gr1, gr2, gc1, gc2) #dense matrix loading of a specific genomic region
    hic.close()

    return numpy_matrix

