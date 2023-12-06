'''
Structure-Based simplicial measures.
Most of the measures are the same code present on Networkx.
'''

import math
import numpy as np
import pytest
import statistics as stats
import networkx as nx
from scipy.linalg import expm, logm

from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *
from digplexq.utils import *

__all__ = [
    "average_shortest_q_walk_length",
    "q_eccentricity",
    "q_diameter",
    "q_radius",
    "q_density",
    "in_q_degree",
    "out_q_degree",
    "in_q_degree_centrality",
    "out_q_degree_centrality",
    "upper_q_degree_centrality",
    "q_closeness_centrality",
    "q_harmonic_centrality",
    "q_betweenness_centrality",
    "q_katz_centrality",
    "global_q_reaching_centrality",
    "q_efficiency",
    "global_q_efficiency",
    "average_q_clustering_coefficient",
    "in_q_degree_rich_club_coefficient",
    "out_q_degree_rich_club_coefficient",
    "communicability",
    "q_communicability",
    "q_returnability",
    "q_structural_entropy",
    "in_q_degree_distribution_entropy",
    "out_q_degree_distribution_entropy",
]


#----- Distance-based Measures -----

#Average shortest q-walk length
def average_shortest_q_walk_length(Hq, weight=None):
    '''Returns
    Parameters
    ----------
    Hq: q-dajacency matrix.
    weight: (string) if None; if 'reciprocal'.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        AV = nx.average_shortest_path_length(Gq, weight=weight)
        return AV

    #If Hq is not weakly connected:
    else:
        AV = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            AV.append(nx.average_shortest_path_length(Sq, weight=weight))
        Max_av = max(AV)
        return round(Max_av, 4)
    return Hq


def q_eccentricity(Hq):
    '''Returns
    Parameters
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is strongly connected:
    if nx.is_strongly_connected(Gq) == True:
        ecc = nx.eccentricity(Gq)
        return ecc

    #If Hq is not strongly connected:
    else:
        return math.inf


def q_diameter(Hq):
    '''Returns the diameter of Gq.
    Parameters
    ----------
    Hq: q-adjacency matrix of the directed q-graph.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is strongly connected:
    if nx.is_strongly_connected(Gq) == True:
        Diam = nx.diameter(Gq)
        return Diam

    #If Hq is not strongly connected:
    else:
        return math.inf


def q_radius(Hq):
    '''Returns the radio of a digraph.
    Parameters
    ----------
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is strongly connected:
    if nx.is_strongly_connected(Gq) == True:
        rad = nx.radius(Gq)
        return rad

    #If Hq is not strongly connected:
    else:
        return math.inf


def q_density(Hq):
    '''Returns the
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
        s_den = nx.density(Gq)
        return s_den

    #If Hq is not weakly connected:
    else:
        SDen = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            s_den = nx.density(Sq)
            SDen.append(s_den)
        Max_den = max(SDen)
        return round(Max_den, 4)



#----- Degree Centrality -----

def in_q_degree(Hq):
    '''Returns the in-q-degree of a q-digraph.
    M: adjacency matrix of the directed q-graph.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        inDg = []
        for i in range(len(Hq)):
            inDg.append(Gq.in_degree(i))
        Max_in_dc = max(inDg)
        return Max_in_dc

    #If Hq is not weakly connected:
    else:
        iDC = []
        inDg = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            for i in range(len(Sq)):
                inDg.append(Sq.in_degree(i))
            iDC.append(max(inDg))
        Max_in_dc = max(iDC)
        return round(Max_in_dc, 4)


def out_q_degree(Hq):
    '''Returns the out-q-degree of a q-digraph.
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
        outDg = []
        for i in range(len(Hq)):
            outDg.append(Gq.out_degree(i))
        Max_out_dc = max(outDg)
        return Max_out_dc

    #If Hq is not weakly connected:
    else:
        oDC = []
        outDg = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            for i in range(len(Sq)):
                outDg.append(Gq.out_degree(i))
            oDC.append(max(outDg))
        Max_out_dc = max(oDC)
        return round(Max_out_dc, 4)


def in_q_degree_centrality(Hq):
    '''Returns the in-q-degree centrality of a q-digraph.
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
        in_dc = nx.in_degree_centrality(Gq)
        Max_in_dc = max(dict_to_array(in_dc))
        return Max_in_dc

    #If Hq is not weakly connected:
    else:
        iDC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            in_dc = nx.in_degree_centrality(Sq)
            iDC.append(max(dict_to_array(in_dc)))
        Max_in_dc = max(iDC)
        return round(Max_in_dc, 4)


def out_q_degree_centrality(Hq):
    '''Returns the out-q-degree centrality of a q-digraph.
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
        out_dc = nx.out_degree_centrality(Gq)
        Max_out_dc = max(dict_to_array(out_dc))
        return Max_out_dc

    #If Hq is not weakly connected:
    else:
        oDC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            out_dc = nx.out_degree_centrality(Sq)
            oDC.append(max(dict_to_array(out_dc)))
        Max_out_dc = max(oDC)
        return round(Max_out_dc, 4)


def upper_q_degree_centrality(M, sigma):
    '''Returns
    Parameters
    '''
    CinDeg = (2**q)/(math.comb(f-(q+1), q+h))
    return CinDeg



#----- Closeness Centrality -----
def q_closeness_centrality(Hq, wf_improved=False):
    '''Returns the maximum BC value among the maximum BC value of
    each weakly connected component.
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
        dscc = nx.closeness_centrality(Gq, wf_improved=wf_improved)
        Max_dscc = max(dict_to_array(dscc))
        return Max_dscc

    #If Hq is not weakly connected:
    else:
        CC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            dscc = nx.closeness_centrality(Sq, wf_improved=wf_improved)
            CC.append(max(dict_to_array(dscc)))
        Max_dscc = max(CC)
        return round(Max_dscc, 4)


#----- Harmonic Centrality -----
def q_harmonic_centrality(Hq):
    '''Returns
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
        shc = nx.harmonic_centrality(Gq)
        Max_shc = max(dict_to_array(shc))
        return Max_shc

    #If Hq is not weakly connected:
    else:
        HC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            shc = nx.harmonic_centrality(Sq)
            HC.append(max(dict_to_array(shc)))
        Max_shc = max(HC)
        return round(Max_shc, 4)



def q_betweenness_centrality(Hq, normalized=True, weight=None):
    '''Returns the maximum BC value among the maximum BC value of
    each weakly connected component.
    Parameters
    ----------
    M: adjacency matrix of the directed q-graph.
    weight: (string) if None; if 'reciprocal'.
    Rerturn:
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        dsbc = nx.betweenness_centrality(Gq, normalized=normalized, weight=weight)
        Max_dsbc = max(dict_to_array(dsbc))
        return round(Max_dsbc, 4)

    #If Hq is not weakly connected:
    else:
        BCC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            dsbc = nx.betweenness_centrality(Sq, normalized=normalized, weight=weight)
            BCC.append(max(dict_to_array(dsbc)))
        Max_dsbc = max(BCC)
        return round(Max_dsbc, 4)



def q_katz_centrality(Hq, alpha=0.1, beta=1.0, normalized=True, weight=None):
    '''Returns the maximum Katz centrality.
    Parameters
    ---------
    Hq:
    weight: (string) if None; if 'reciprocal'.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    a1 = 1/max(np.linalg.eigvals(Hq))
    if a1 < 0.1:
        alpha = np.real(a1)/3
    else:
        alpha = 0.1

    katz = nx.katz_centrality_numpy(Gq, alpha=alpha, beta=beta, normalized=normalized, weight=weight)
    
    Max_katz = max(dict_to_array(katz))
    return round(Max_katz, 4)


#global_reaching_centrality(G, weight=None, normalized=True)
def global_q_reaching_centrality(Hq, normalized=False):
    '''Returns
    Parameters
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0

    gsrc = nx.global_reaching_centrality(Gq, normalized=normalized)
    return round(gsrc, 5)



#----- Efficiency and Global Efficiency -----

def q_efficiency(Hq, i):
    '''Returns
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    n = len(Hq)
    nE_i = 0
    
    for j in range(n):
        if j != i:
            try:
                SPL = nx.shortest_path_length(Gq, source=i, target=j, weight=None, method='dijkstra')
                if SPL != 0:
                    nE_i += 1/SPL
                else:
                    pass
            except nx.NetworkXNoPath:
                pass
        else:
            pass
            
    E_i = nE_i/(n-1)
    return round(E_i, 5)


def global_q_efficiency(Hq):
    '''Returns
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    n = len(Hq)
    nGE = 0
    for i in range(n):
        nGE += q_efficiency(Hq, i)
    
    GE = nGE/n
    return round(GE, 5)



#----- Segregation Measures -----

def average_q_clustering_coefficient(Hq, weight=None):
    '''Returns
    Parameters
    ----------
    weight: (string) if None; if 'reciprocal'.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    ACC = nx.average_clustering(Gq, nodes=None, weight=weight, count_zeros=True)
    return round(ACC, 4)



def in_q_degree_rich_club_coefficient(M, k=6):
    '''Returns the in-degree rich-club coefficient phi(k) = E^{in}_k/N^{in}_k(N^{in}_k-1).
    Parameters
    ----------
    M: q-adjacency matrix.
    k: integer.
    '''
    if isinstance(M, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(M, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    N_k = 0
    n = len(M)
    nodes = []
    
    for i in range(n):
        indeg = Gq.in_degree(i)
        if indeg > k:
            N_k += 1
            nodes.append(i)
        else:
            pass
    
    if N_k > 1: 
        H = Gq.subgraph(nodes)
        E_K = Gq.number_of_edges()
        RCC = E_K/(N_k*(N_k - 1))
        return round(RCC, 5)
    else:
        return 0

    
def out_q_degree_rich_club_coefficient(M, k=6):
    '''Returns the out-degree rich-club coefficient phi(k) = E^{out}_k/N^{out}_k(N^{out}_k-1).
    Parameters
    ----------
    M: q-adjacency matrix.
    k: integer.
    '''
    if isinstance(M, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(M, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    N_k = 0
    n = len(M)
    nodes = []
    
    for i in range(n):
        outdeg = Gq.out_degree(i)
        if outdeg > k:
            N_k += 1
            nodes.append(i)
        else:
            pass
    
    if N_k > 1: 
        H = Gq.subgraph(nodes)
        E_K = Gq.number_of_edges()
        RCC = E_K/(N_k*(N_k - 1))
        return round(RCC, 5)
    else:
        return 0



#----- Communicability -----

def communicability(A, i, j):
    '''Returns the communicability of the nodes i and j.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gr = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(Gr) == True:
        return 0
    
    G = expm(A)
    G_ij = G[i,j]
    return round(G_ij, 5)


def q_communicability(Hq):
    '''Returns the communicability of the nodes i and j.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    V = []
    for i in range(len(Hq)):
        for j in range(len(Hq)):
            Gij = communicability(Hq, i, j)
            V.append(Gij)
    
    Max_Gij = max(V)
    return round(Max_Gij, 5)



#----- Returnability -----

def q_returnability(A, normalized=True):
    '''Returns the simplicial q-returnability.
    Parameters
    ---------
    A: (NumPy matrix) q-adjacency matrix.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gr = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(Gr) == True:
        return 0

    Exp = expm(A)
    K = np.trace(Exp) - len(A)

    if normalized == True:
        A_wde = remove_double_edges(A)
        G = nx.Graph(A_wde)
        U = nx.to_numpy_matrix(G)
        
        Exp_und = Exp = expm(U)
        
        K_und = np.trace(Exp_und) - len(U)
        Kn = K/K_und
        return round(Kn, 5)
    else:
        return round(K, 5)

    

#----- Entropies -----

def q_structural_entropy(Hq):
    '''Returns the q-structural entropy of a digraph.
    Parameters
    ----------
    Hq: q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    n = len(Hq)
    sum_Q = 0
    H = 0
    sum_Qi = []
    
    for i in range(n):
        Q_i = 0
        for j in range(n):
            Q_i += communicability(Hq, i, j)
        sum_Qi.append(Q_i)
    
    for i in range(n):
        for j in range(n):
            sum_Q += communicability(Hq, i, j)
    
    if sum_Q == 0:
        return 0
    
    for k in range(n):
        H += (sum_Qi[k]/sum_Q)*math.log2(sum_Qi[k]/sum_Q)
    
    return round(-H, 5)


def in_q_degree_distribution_entropy(Hq):
    '''Returns the in-q-degree distribution entropy of a digraph.
    Parameters
    ----------
    Hq: q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    D = []
    H = 0
    n = len(Hq)
    
    for k in range(n):
        delta_k = 0
        for i in range(n):
            if Gq.in_degree(i) == k:
                delta_k += 1
            else:
                pass
        if delta_k != 0:
            D.append(delta_k)
        else:
            pass
    
    for k in range(len(D)):
        H += (D[k]/n)*math.log2(D[k]/n)
    
    return round(-H, 5)


def out_q_degree_distribution_entropy(Hq):
    '''Returns the out-q-degree distribution entropy of a digraph.
    Parameters
    ----------
    Hq: q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())

    if nx.is_empty(Gq) == True:
        return 0
    
    D = []
    H = 0
    n = len(Hq)
    
    for k in range(n):
        delta_k = 0
        for i in range(n):
            if Gq.out_degree(i) == k:
                delta_k += 1
            else:
                pass
        if delta_k != 0:
            D.append(delta_k)
        else:
            pass
    
    for k in range(len(D)):
        H += (D[k]/n)*math.log2(D[k]/n)
    
    return round(-H, 5)