'''
Substructure enumeration.
'''

import numpy as np
import networkx as nx

from digplexq.persistent_homology import *
from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *

__all_ = [
    "enumerate_double_edges",
    "enumerate_directed_n_cycles",
    "enumerate_edges_elementary_directed_quase_clique",
    "enumerate_invariant_elementary_n_paths",
    "count_weakly_q_connected_components",
    "count_n_paths",
    "count_directed_n_cliques",
    "count_elmentary_directed_quasi_cliques",
    "count_double_edges",
    "count_directed_n_cycles",
    "zscore_motif",
    "first_flag_topological_structure_vector",
    "second_flag_topological_structure_vector",
    "third_flag_topological_structure_vector",
    "fourth_flag_topological_structure_vector",
    "fifth_path_topological_structure_vector",
    "first_path_topological_structure_vector",
    "second_path_topological_structure_vector",
    "third_path_topological_structure_vector",
    "fourth_path_topological_structure_vector",
    "fifth_path_topological_structure_vector",
]


#--------------------------------

 
def enumerate_double_edges(M):
    '''Returns all double edges.
    Parameters
    ---------
    '''
    W = []
    for i in range(len(M)):
        for j in range(len(M)):
            if i != j and M[i,j] != 0 and M[j,i] != 0:
                W.append([i,j,i])
            else:
                pass
    return W


def enumerate_directed_n_cycles(M, n):
    '''Returns all directed n-cycles.
    Parameters
    ---------
    M: adjacency matrix.
    n: lenght.
    '''
    DC = []
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    C = list(nx.simple_cycles(G))
    for i in range(len(C)):
        if len(C[i]) == n:
            DC.append(C[i])
        else:
            pass
    return DC
    
    
def enumerate_edges_elementary_directed_quase_clique(path):
    '''Returns all edges in the elementary DQC.
    Parameters
    ---------
    path: (array) elementary n-path.
    '''
    path_edges = []
    for i in range(len(path)-2):
        for j in range(1,3):
            path_edges.append([path[i], path[i+j]])
    
    path_edges.append([path[len(path)-2], path[len(path)-1]])
    return path_edges


def enumerate_invariant_elementary_n_paths(PC, n):
    '''Returns all invariant elementary n-paths.
    M: (NumPy array) a path complex.
    n: (integer) 
    '''
    if n < 1:
        raise ValueError("n must be an integer greater than or equal to 1.")
       
    if n > len(PC):
        return []
    
    #PC = PathComplex(M)
    Elem_Inv_Paths = []
    for path in PC[n-1]:
        DQC_edges = []
        EdgesDQC = enumerate_edges_elementary_directed_quase_clique(path[0])
        for edge in EdgesDQC:
            edge2 = [edge]
            if edge2 in PC[0]:
                DQC_edges.append(edge)
        if len(DQC_edges) == len(EdgesDQC):
            Elem_Inv_Paths.append(path[0])
    return Elem_Inv_Paths



#------------------------------

def count_weakly_q_connected_components(Hq):
    '''Returns all weakly q-connected components.
    Hq: q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a square matrix.")
    
    Struct_Vec = []
    G = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    wcc = list(nx.weakly_connected_components(G))
    
    for i in range(len(wcc)):
        Struct_Vec.append(i+1)
        
    return Struct_Vec
 
    
def count_strongly_q_connected_components(Hq):
    '''Returns all strongly q-connected components.
    Hq: q-adjacency matrix.
    '''
    if isinstance(Hq, np.ndarray) == False:
        raise TypeError("Input must be a square matrix.")
    
    Struct_Vec = []
    G = nx.from_numpy_matrix(Hq, create_using=nx.DiGraph())
    scc = list(nx.strongly_connected_components(G))
    
    for i in range(len(scc)):
        Struct_Vec.append(i+1)
        
    return Struct_Vec


#Amount of n-paths in a path complex:
def count_n_paths(PC):
    '''Returns the number of n-paths.
    Parameters
    ----------
    PC: (NumPy array) path complex.
    '''
    Qtd = []
    for i in range(len(PC)):
        Qtd.append(len(PC[i]))
    return Qtd


#Amount of directed n-cliques in a DFC:
def count_directed_n_cliques(DFC):
    '''Returns the number of directed n-cliques.
    Parameters
    '''
    Qtd = []
    for i in range(len(DFC)):
        Qtd.append(len(DFC[i]))
    return Qtd


def count_elmentary_directed_quasi_cliques(PC, n):
    '''Returns the mount of invariant elementary n-paths in a path complex.
    PC: path complex.
    '''
    qtd_n_dqc = len(enumerate_invariant_elementary_n_paths(PC, n))
    return qtd_n_dqc



def count_edges(M):
    '''Returns all 
    '''
    n = 0
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i,j] > 0:
                n+=1
    return n


def count_double_edges(M):
    '''Returns the number of double-edges.
    '''
    Qtd = len(enumerate_double_edges(M))
    return Qtd
        
    
def count_directed_n_cycles(M):
    '''Returns the number of directed cycles.
    '''
    Qtd = len(enumerate_directed_n_cycles(M))
    return Qtd   


def zscore_motif(M):
    '''Returns all 
    '''
    return M


#------------------------------



def first_flag_topological_structure_vector(DFC):
    '''Returns a vector which entries are the quantities of k-dimensional directed cliques.
    Parameters
    ----------
    DFC: (array) directed flag complex.
    '''
    count_dir_cliques = count_directed_n_cliques(DFC)
    return count_dir_cliques


def second_flag_topological_structure_vector(Hq):
    '''Returns a vector.
    Parameters
    ----------
    Hq: q-adjacency matrix.
    '''
    wcc = count_weakly_q_connected_components(Hq)
    return wcc


def third_flag_topological_structure_vector(Hq):
    '''Returns a vector.
    Parameters
    ----------
    Hq: q-adjacency matrix.
    '''
    scc = count_strongly_q_connected_components(Hq)
    return scc


def fourth_flag_topological_structure_vector(M):
    '''Returns a vector. 
    Parameters
    ----------
    DFC: (array) directed flag complex.
    '''
    v = betti_numbers(M, k=0)    
    return v


def fifth_flag_topological_structure_vector(M):
    '''Returns a vector.
    Parameters
    ----------
    '''
    bcl = barcode_length(M)
    return bcl


#------------------------------

def first_path_topological_structure_vector(PC):
    '''Returns a vector.
    Parameters
    ----------
    PC: (array) path complex.
    '''
    n = len(PC)
    Struct_Vec = []
    for i in range(1,n):
        qtd_i_dqc = len(enumerate_invariant_elementary_n_paths(PC, i))
        Struct_Vec.append(qtd_i_dqc)
    return Struct_Vec


def second_path_topological_structure_vector(PC):
    '''Returns a vector.
    '''
    count_paths = count_n_paths(PC)
    return count_paths


def third_path_topological_structure_vector(M):
    '''Returns a vector. 
    '''
    return M


def fourth_path_topological_structure_vector(M):
    '''Returns a vector. 
    '''
    return M


def fifth_path_topological_structure_vector(M):
    '''Returns a vector. 
    '''
    return M


