'''
Substructure enumeration.
'''

import numpy as np
import networkx as nx

from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *

#__all_ = [
    #"enumerate_double_edges",
    #"enumerate_directed_n_cycles",
    #"enumerate_edges_elementary_directed_quase_clique",
    #"enumerate_invariant_elementary_n_paths",
    #"count_weakly_q_connected_components",
    #"count_n_paths",
    #"count_directed_n_cliques",
    #"count_directed_quasi_cliques",
    #"count_squares",
    #"count_double_edges",
    #"count_n_cycles",
    #"first_topological_structure_vector",
#]


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
    '''
    if n < 1:
        raise ValueError("n must be an integer greater than or equal to 1.")
    
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


#Amount of squares
def count_squares(M):
    '''Returns all 
    '''
    Qtd = []
    for i in range(len(DFC)):
        Qtd.append(len(DFC[i]))
    return Qtd


#Amount of edges
def count_edges(M):
    '''Returns all 
    '''
    n = 0
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i,j] > 0:
                n+=1
    return n

#Amount of double edges
def count_double_edges(M):
    '''Returns all 
    '''
    Qtd = []
    for i in range(len(DFC)):
        Qtd.append(len(DFC[i]))
    return Qtd
        
    
#Amount of cycles
def count_n_cycles(M):
    '''Returns all 
    '''
    Qtd = []
    for i in range(len(DFC)):
        Qtd.append(len(DFC[i]))
    return Qtd   



#------------------------------

def first_structure_vector(M):
    '''Returns all 
    '''
    return M


def second_structure_vector(M):
    '''Returns all 
    '''
    return M


def first_topological_structure_vector(M):
    '''Returns all 
    '''
    return M


def second_topological_structure_vector(PC, n):
    '''Returns a list of all 
    '''
    Struct_Vec = []
    for i in range(n):
        qtd_i_dqc = len(enumerate_invariant_elementary_n_paths(PC, i))
        Struct_Vec.append(qtd_i_dqc)
    return Struct_Vec

