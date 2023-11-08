'''
Structure-Based simplicial measures.
Most of the measures are the same code present on Networkx.
'''

import math
import numpy as np
import statistics as stats
import networkx as nx
from scipy.linalg import expm, logm, sinm, cosm

from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *


#__all__ = [
#    "average_shortest_q_walk_length",
#    "q_eccentricity",
#    "q_diameter",
#    "q_radius",
#    "q_hubness",
#    "simplicial_density",
#    "in_q_degree",
#    "out_q_degree",
#    "upper_q_degree",
#    "in_q_degree_centrality",
#    "out_q_degree_centrality",
#    "upper_degree_centrality",
#    ""q_harmonic_centrality",",
#    "q_closeness_centrality",
#    "q_betweenness_centrality",
#    "q_katz_centrality",
#    "q_reaching_centrality",
#    "q_efficiency",
#    "q_vulnerability,
#    "q_clustering_coefficient",
#    "q_modularity_coefficient",
#    "q_rich_club_coefficient",
#    "q_communicability",
#    "q_returnability",
#    "q_entropy",
#]


#----- Distance-based Measures -----


#Average shortest q-walk length
def average_shortest_q_walk_length(Hq):
    '''Returns
    Parameters
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    
    if nx.is_empty(Gq) == True:
        return 0
    
    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        AV = nx.average_shortest_path_length(Gq)
        return AV
    
    #If Hq is not weakly connected:
    else:
        AV = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            AV.append(nx.average_shortest_path_length(Sq))
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
    M: adjacency matrix of the directed q-graph.
    Rerturn: 
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
        rad = nx.radius(Gq)
        return rad
    
    #If Hq is not strongly connected:
    else:
        return math.inf


def q_hubness(Hq):
    '''Returns
    Parameters
    '''
    return Hq


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
    
    #If Hq is weakly connected:
    if nx.is_weakly_connected(Gq) == True:
        gsrc = nx.global_reaching_centrality(Gq, normalized=normalized)
        return round(gsrc, 4)
    
    #If Hq is not weakly connected:
    else:
        GRC = []
        wcc = adjacency_matrices_wcc(Hq)
        for c in wcc:
            Sq = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            gsrc = nx.global_reaching_centrality(Sq, normalized=normalized)
            GRC.append(gsrc)
        Max_gsrc = max(GRC)
        return round(Max_gsrc, 4)


#local_reaching_centrality(G, v, paths=None, weight=None, normalized=True)
def local_q_reaching_centrality(Hq, simplex, normalized=False):
    '''Returns
    Parameters
    '''
    return Hq
        

#katz_centrality(G, alpha=0.1, beta=1.0, max_iter=1000, tol=1e-06, nstart=None, normalized=True, weight=None)
def q_katz_centrality(Hq, normalized=True, weight=None):
    '''Returns the maximum Katz centrality.
    Parameters
    ---------
    Hq: 
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    
    if nx.is_empty(Gq) == True:
        return 0
    
    a1 = 1/max(np.linalg.eigvals(Hq))
    if a1 < 0.1:
        alpha = a1/3
    else:
        alpha = 0.1
    
    katz = nx.katz_centrality(Gq, alpha=alpha, normalized=normalized, weight=weight)
    Max_katz = max(dict_to_array(katz))
    return round(Max_katz, 4)

       
    
#flow_hierarchy / flow_hierarchy(G)
def q_flow_hierarchy(Hq):
    '''Returns
    Parameters
    '''
    return Hq



#----- Efficiency and Vulnerability -----

#
def global_q_efficiency(Hq):
    '''Returns
    This function compute the amount of
    directed bicliques
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    
    if nx.is_empty(Gq) == True:
        return 0
    
    return Hq


#
def q_vulnerability(Hq):
    '''Returns
    '''
    return Hq




#----- Segregation Measures -----

#
def q_clustering_coefficient(Hq):
    '''Returns
    '''
    return Hq
        

#
def q_modularity_coefficient(Hq):
    '''Returns the 
    '''
    return Hq


#rich_club_coefficient(G, normalized=True, Q=100, seed=None)
def q_rich_club_coefficient(Hq):
    '''Returns
    '''
    return Hq



#----- Communicability -----


def q_communicability(Hq):
    '''Returns
    '''
    return Hq


#----- Returnability -----

def q_returnability(Hq, normalized=True):
    '''Returns the simplicial q-returnability.
    Parameters
    ---------
    Hq: (NumPy matrix) q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    Gq = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    
    if nx.is_empty(Gq) == True:
        return 0
    
    Exp = expm(Hq)
    K_q = np.trace(Exp) - len(Hq)
    
    if normalized == True:
        return round(K_q, 4)
        
    else:
        return round(K_q, 4)


#----- Entropies -----

def structural_q_entropy(Hq):
    '''Returns
    '''
    return Hq




    