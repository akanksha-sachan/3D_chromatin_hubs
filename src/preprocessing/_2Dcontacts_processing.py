######### CONTACT MATRIX (.hic/.mcool) QUERYING for Observed and OE counts #########

import io
import json
import os
import struct
import subprocess

import cooler
import hicstraw
import numba as njit
import numpy as np
from memory_profiler import profile
from scipy.sparse import csr_matrix

# import utils from current working code dir
#from ..src.utils.helper import *


###### .mcool/.hic file querying  ######

class Query:
    """
    Querying .mcool/.hic files for observed and OE counts
    """

    def __init__(self, config):
        self.config = config
        self.cooler = cooler.Cooler(config["cool_file"])
        self.hic = hicstraw.HiC(
            config["hic_file"])
        self.temp_dir = config["temp_dir"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)

    def get_chromsizes(self):
        """
        Get chromosomes from .mcool/.hic file
        """
        return self.cooler.chromnames

    def get_bins(self, chrom):
        """
        Get bins from .mcool/.hic file
        """
        bins = self.cooler.bins().fetch(chrom)
        return bins

    def get_observed(self, chrom):
        """
        Get observed counts from .mcool/.hic file
        """
        obs = self.cooler.matrix(balance=False).fetch(chrom)
        return obs

    def get_expected_hic(self, chrom):
        """
        Get expected counts from .hic file
        """
        exp = self.cooler.matrix(balance=True).fetch(chrom)
        return exp

    def get_oe_hic(self, chrom):
        """
        Get observed and expected counts from .mcool/.hic file
        """
        obs = self.get_observed(chrom)
        exp = self.get_expected_hic(chrom)
        return obs, exp

    def get_interactions_sparse(self, chrom):
        """
        Get observed and expected counts from .mcool/.hic file as sparse matrices
        """
        obs = self.get_observed(chrom)
        exp = self.get_expected_hic(chrom)
        return csr_matrix(obs), csr_matrix(exp)


###### Processing utils (norm functions) ######

class Process:
    """
    Post querying processing of observed/oe counts
    """

    def __init__(self, config):
        self.config = config
        self.temp_dir = config["temp_dir"]
        self.resolutions = config["resolutions"]
        if not os.path.exists(self.temp_dir):
            os.mkdir(self.temp_dir)

    def get_oe_cool(self, M, resolution, threshold=1):
        """ 
        The O/E matrix is calculated as the log2 ratio of the raw contact matrix to the expected contact matrix.
        The expected contact matrix is calculated by filling in the average value of the diagonals of the raw contact matrix.
        Remove the NaN bins before calculating O/E so that interpolated edges aren't used
        """
        #process M to mask out the centeromeric bins (NANs currently)
        M = np.nan_to_num(M)
        
        #construct expected matrix
        E = np.zeros_like(M).astype(float)
        sums = []
        for i in range(M.shape[0]):
            contacts = np.diag(M,i)
            #using on non zero diagonal elements to get the denominator
            non_zero_indices = np.nonzero(contacts)[0]
            if len(non_zero_indices) > 0:
                #chr wide expected, not factorized by number of chrs for genome-wide comparison
                expected = contacts.sum() / len(non_zero_indices) 
                
            else:
                expected = 0
            sums.append(expected)
            #uniform distribution of contacts across diagonals assumed
            x_diag,y_diag = np.diag_indices(M.shape[0]-i)
            x,y = x_diag,y_diag+i
            E[x,y] = expected
        E += E.T
        eps=1e-5
        E = np.nan_to_num(E) + eps
        OE = M / E 
        OE[OE == 0] = 1 #to avoid neg inf in log
        OE = np.log(OE) #log transform the OE to get equal-sized bins, as the expected values need to be log binned
        # threshold elements based on M (TODO: decide on an adaptive threshold)
        OE_filtered = np.where(OE > threshold, OE, 0)

        #save OE and expected matrices in temp dir as .npy
        np.save(f"{self.temp_dir}/OE_{resolution}.npy", OE_filtered)
        np.save(f"{self.temp_dir}/E_{resolution}.npy", E)

    def sqrt_norm(self, matrix):
        """
        Square root normalization of a matrix
        """
        coverage = (np.sqrt(np.sum(matrix, axis=-1)))
        with np.errstate(divide='ignore', invalid='ignore'):
            matrix = matrix / coverage.reshape((-1, 1)); matrix = matrix / coverage.reshape((1, -1))
        matrix[np.isnan(matrix)] = 0.0; matrix[np.isinf(matrix)] = 0.0
        return matrix

    def pearson(self, matrix):
        """
        Pearson correlation coefficient between two matrices
        """
        return np.corrcoef(matrix)

    def knight_ruiz_alg(self, A, tol=1e-8, f1=False):
        """
        Knight and Ruiz algorithm for matrix balancing
        """
        n = A.shape[0]
        e = np.ones((n, 1), dtype=np.float64)
        res = []

        Delta = 3
        delta = 0.1
        x0 = np.copy(e)
        g = 0.9

        etamax = eta = 0.1
        stop_tol = tol * 0.5
        x = np.copy(x0)

        rt = tol**2.0
        v = x * (A.dot(x))
        rk = 1.0 - v
        rho_km1 = ((rk.transpose()).dot(rk))[0, 0]
        rho_km2 = rho_km1
        rout = rold = rho_km1

        MVP = 0  # we'll count matrix vector products
        i = 0  # outer iteration count

        if f1:
            print("it in. it res\n")
        k = 0
        while rout > rt:  # outer iteration
            i += 1

            if i > 30:
                break

            k = 0
            y = np.copy(e)
            innertol = max(eta**2.0 * rout, rt)

            while rho_km1 > innertol:  # inner iteration by CG
                k += 1
                if k == 1:
                    Z = rk / v
                    p = np.copy(Z)
                    rho_km1 = (rk.transpose()).dot(Z)
                else:
                    beta = rho_km1 / rho_km2
                    p = Z + beta * p

                if k > 10:
                    break

                w = x * A.dot(x * p) + v * p
                alpha = rho_km1 / (((p.transpose()).dot(w))[0, 0])
                ap = alpha * p
                ynew = y + ap

                if np.amin(ynew) <= delta:

                    if delta == 0:
                        break

                    ind = np.where(ap < 0.0)[0]
                    gamma = np.amin((delta - y[ind]) / ap[ind])
                    y += gamma * ap
                    break

                if np.amax(ynew) >= Delta:
                    ind = np.where(ynew > Delta)[0]
                    gamma = np.amin((Delta - y[ind]) / ap[ind])
                    y += gamma * ap
                    break

                y = np.copy(ynew)
                rk -= alpha * w
                rho_km2 = rho_km1
                Z = rk / v
                rho_km1 = ((rk.transpose()).dot(Z))[0, 0]
            x *= y
            v = x * (A.dot(x))
            rk = 1.0 - v
            rho_km1 = ((rk.transpose()).dot(rk))[0, 0]
            rout = rho_km1
            MVP += k + 1

            rat = rout / rold
            rold = rout
            res_norm = rout**0.5
            eta_o = eta
            eta = g * rat
            if g * eta_o**2.0 > 0.1:
                eta = max(eta, g * eta_o**2.0)
            eta = max(min(eta, etamax), stop_tol / res_norm)
            if f1:
                print(
                    "{:03d} {:06d} {:.3f} {:.e} {:.e} \n".format(
                        i, k, res_norm, rt, rout
                    )
                )
                res.append(res_norm)
        if f1:
            output = "Matrix - vector products = %06i\n" % MVP
            print(output)

        return [x, i, k]
    
print("hello world")