'''
Spectrum-Based Simplicial Measures
'''

import math
import numpy as np
import networkx as nx
from scipy.linalg import expm, logm, sqrtm
from hodgelaplacians import HodgeLaplacians

from digplexq.directed_q_analysis import *
from digplexq.utils import *

__all__ = [
    "hodge_q_laplacian",
    "q_energy",
    "q_spectral_density",
    "q_spectral_entropy",
    "q_vonNeuman_entropy",
    "q_KL_divergence",
    "q_spectral_distance",
    "q_eigenvector_centrality",
    "q_pagerank_centrality",
]


def hodge_q_laplacian(MaxSimp, q):
    '''Returns the Hodge q-Laplacian.
    Parameters
    ---------
    MaxSimp: (array) maximal simplices.
    q: (integer).
    Notes
    -----
    Based on the package hodgelaplacians.
    '''

    n = len(MaxSimp)
    maxdim = len(MaxSimp[n-1][0])

    if q >= maxdim-1:
        return np.array([[0]])
    else:
        MaxSimp_conn = connect_array(MaxSimp)
        hl = HodgeLaplacians(MaxSimp_conn)
        Lq = hl.getHodgeLaplacian(q)
        Lq_arr = Lq.toarray()
        return Lq_arr


def q_energy(Hq):
    '''Returns the q-energy of a digraph.
    Parameters
    ----------
    Hq: (array) q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    Bq = to_binary(Hq)
    BqtBq = np.matmul(Bq.transpose(), Bq)
    sqrt_BqtBq = sqrtm(BqtBq)
    Eq = np.trace(sqrt_BqtBq)
    real_Eq = np.real(Eq)

    return round(real_Eq, 5)


def q_spectral_density(L, b=1):
    '''Returns the q-spectral density of a (di)graph.
    Parameters
    ----------
    L: (array) Hodge q-Laplacian.
    '''
    if np.all(L==0) == True:
        return 0

    Exp = expm(-b*L)
    rho = (Exp)/np.trace(Exp)
    return rho


def q_spectral_entropy(L):
    '''Returns the q-spectral entropy of a (di)graph.
    Parameters
    ---------
    L: (array) Hodge q-Laplacian.
    '''
    if np.all(L==0) == True:
        return 0

    EigenVal = np.linalg.eig(L)[0]

    EVN = []
    for lambda_i in EigenVal:
        EVN.append(round(lambda_i/sum(EigenVal), 6))

    H = 0
    for mu_i in EVN:
        if mu_i != 0:
            H += mu_i*math.log2(mu_i)

    Hq = -np.real(H)
    return round(Hq, 5)


def q_vonNeuman_entropy(L, b=1):
    '''Returns the q-von-Neuman entropy of a (di)graph.
    Parameters
    ---------
    L: (array) Hodge q-Laplacian.
    b: parameter.
    '''
    if np.all(L==0) == True:
        return 0

    rho = q_spectral_density(L, b)
    S = -np.trace(rho*logm(rho))
    return round(S, 5)


def q_KL_divergence(L1, L2, b=1):
    '''Returns the q-Kullback-Liebler divergence.
    Parameters
    ---------
    L1: (array) Hodge q-Laplacian.
    L2: (array) Hodge q-Laplacian.
    ---------
    '''
    rho1 = q_spectral_density(L1, b)
    rho2 = q_spectral_density(L2, b)
    KL = np.trace(rho1*(logm(rho1) - logm(rho2)))
    return round(KL, 5)


def q_spectral_distance(L1, L2):
    '''Returns the q-spectral distance between two (di)graphs.
    Parameters
    ---------
    L1: Hodge q-Laplacian.
    L2: Hodge q-Laplacian.
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

    dist = S/N
    return round(np.real(dist), 5)


def q_eigenvector_centrality(Hq, weight=None):
    '''Returns the q-eigenvector centrality of a digraph.
    Parameters
    ---------
    Hq: (array) q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    ec = nx.eigenvector_centrality_numpy(Gq, weight=weight, max_iter=100, tol=0)
    Max_ec = max(dict_to_array(ec))
    return round(Max_ec, 5)


def q_pagerank_centrality(Hq):
    '''Returns the maximum value of the pagerank centrality.
    Parameters
    ---------
    Hq: (array) q-adjacency matrix.
    Rerturn:
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    prc = nx.pagerank(Gq, alpha=0.85, personalization=None, max_iter=100, tol=1e-06, nstart=None, weight='weight', dangling=None)
    Max_prc = max(dict_to_array(prc))
    return round(Max_prc, 5)
