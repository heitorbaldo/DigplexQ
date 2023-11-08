'''
Spectrum-Based Simplicial Measures
'''

import math
import numpy as np
import networkx as nx
from scipy.linalg import expm, logm
from hodgelaplacians import HodgeLaplacians

from digplexq.directed_q_analysis import *
from digplexq.utils import *

__all__ = [
    "hodge_q_laplacian",
    "spectral_n_density",
    "spectral_n_entropy",
    "von_Neuman_n_entropy",
    "KL_divergence",
    "spectral_n_distance",
    "eigenvector_q_centrality",
    #"simplicial_pagerank_centrality",
]


def hodge_q_laplacian(MaxSimp, q):
    '''Returns the Hodge q-Laplacian.
    Parameters
    ---------
    MaxSimp: (array) maximal simplices.
    q: (integer)
    Notes
    -----
    Based on the packege hodgelaplacians ().
    '''
    if q > len(MaxSimp):
        return np.array([[0]])
    
    MaxSimp_conn = connect_array(MaxSimp)
    hl = HodgeLaplacians(MaxSimp_conn)
    Lq = hl.getHodgeLaplacian(q)
    Lq_arr = Lq.toarray()
    return Lq_arr


def spectral_n_density(L, b=1):
    '''Returns spectral density.
    Parameters
    L: Hodge n-Laplacian.
    '''
    Exp = expm(-b*L)
    rho = (Exp)/np.trace(Exp)
    return rho


def spectral_n_entropy(L):
    '''Returns spectral n-entropy Hn.
    Parameters
    ---------
    L: Hodge n-Laplacian.
    '''
    L = L.toarray()
    EigenVal = np.linalg.eig(L)[0]
    
    EVN = []
    for lambda_i in EigenVal:
        EVN.append(round(lambda_i/sum(EigenVal), 6))
    
    H = 0
    for mu_i in EVN:
        if mu_i != 0:
            H += mu_i*math.log2(mu_i)
    return -round(H, 6)


def von_Neuman_n_entropy(L, b=1):
    '''Returns von Neuman entropy
    Parameters
    ---------
    L: Hodge n-Laplacian.
    b: parameter
    '''
    rho = spectral_n_density(L, b)
    S = -np.trace(rho*logm(rho))
    return S


def KL_divergence(L1, L2, b):
    '''Returns Kullback-Liebler divergence.
    Parameters
    ---------
    L1: Hodge n-Laplacian.
    L2: Hodge n-Laplacian.
    ---------
    '''
    rho1 = spectral_n_density(L1, b)
    rho2 = spectral_n_density(L2, b)
    KL = np.trace(rho1*(logm(rho1) - logm(rho2)))
    return KL


def spectral_n_distance(L1, L2):
    '''Returns Kullback-Liebler divergence.
    Parameters
    ---------
    L1: Hodge n-Laplacian.
    L2: Hodge n-Laplacian.
    '''
    Smin = []
    Smax = []
    
    Eigen1 = np.linalg.eig(L1)[0]
    Eigen2 = np.linalg.eig(L2)[0]
    Len1 = len(Eigen1)
    Len2 = len(Eigen2)
    Min = min(Len1, Len2)
    Max = max(Len1, Len2)
    for j in range(Min):
        Smin.append(abs(Eigen1[j] - Eigen2[j]))
    for k in range(Min, Max):
        if Len1 > Len2:
            Smax.append(abs(Eigen1[k]))
        if Len1 < Len2:
            Smax.append(abs(Eigen2[k]))
    S = sum(Smin) + sum(Smax)
    N = sum(Eigen1) + sum(Eigen2)
    return S/N


def eigenvector_q_centrality(Hq, weight=None):
    '''Returns
    Parameters
    ---------
    M: adjacency matrix of the directed q-graph.
    Rerturn: 
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    
    if nx.is_empty(Gq) == True:
        return 0
    
    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        ec = nx.eigenvector_centrality(Gq, weight=weight)
        Max_ec = max(dict_to_array(ec))
        return round(Max_ec, 4)
    
    #If Hq is not weakly connected:
    else:
        EC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            ec = nx.eigenvector_centrality(Sq, weight=weight)
            EC.append(max(dict_to_array(ec)))
        Max_ec = max(EC)
        return round(Max_ec, 4)

