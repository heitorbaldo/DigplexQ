'''
Functions for substructure enumeration.
'''

import numpy as np
import networkx as nx

from digplexq.persistent_homology import *
from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *

__all_ = [
    "union_k",
    "intersection_k",
    "enumerate_double_edges",
    "enumerate_directed_n_cycles",
    "enumerate_edges_elementary_directed_quase_clique",
    "enumerate_invariant_elementary_n_paths",
    "f_count",
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
]



def union_k(X, Y, k):
    '''Returns the union of the k-simplices of X and Y.
    '''
    if X == [] and Y == []:
        return []
    
    Diff = []
    n1 = len(X)
    n2 = len(Y)
    
    if k >= n1 and k >= n2:
        return []
        
    if k > n1 and k < n2:
        return Y[k]
    
    if k < n1 and k > n2:
        return X[k]
    
    if k < n1 and k < n2:
        for simplex in Y[k]:
            if simplex not in X[k]:
                Diff.append(simplex)
        Union_k = Diff + X[k]
        return Union_k


def intersection_k(X, Y, k):
    '''Returns the common k-simplices between X and Y.
    '''
    if X == [] or Y == []:
        return []
    
    Inter = []
    n1 = len(X)
    n2 = len(Y)
    
    if k >= n1 and k >= n2:
        return []
        
    if k > n1 and k < n2:
        return Y[k]
    
    if k < n1 and k > n2:
        return X[k]
    
    if k < n1 and k < n2:
        for simplex in Y[k]:
            if simplex in X[k]:
                Inter.append(simplex)
        return Inter


def enumerate_double_edges(A, q=None):
    '''Returns all double edges.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    W = []
    for i in range(len(M)):
        for j in range(len(M)):
            if i != j and M[i,j] != 0 and M[j,i] != 0:
                W.append([i,j,i])
            else:
                pass
    return W


def enumerate_directed_n_cycles(A, n, q=None):
    '''Returns all directed n-cycles.
    Parameters
    ---------
    A: adjacency matrix.
    n: length.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    DC = []
    C = list(nx.simple_cycles(G))
    for i in range(len(C)):
        if len(C[i]) == n:
            DC.append(C[i])
        else:
            pass
    return DC
    
    
def enumerate_edges_elementary_directed_quase_clique(path):
    '''Returns all edges in the elementary directed quasi-clique.
    
    Parameters
    ---------
    path: (array) Elementary n-path.
    '''
    path_edges = []
    for i in range(len(path)-2):
        for j in range(1,3):
            path_edges.append([path[i], path[i+j]])
    
    path_edges.append([path[len(path)-2], path[len(path)-1]])
    return path_edges


def enumerate_invariant_elementary_n_paths(PC, n):
    '''Returns all invariant elementary n-paths.
    
    Parameters
    ----------
    PC: (array) Path complex.
    n: (integer) Path length.
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

def f_count(X, Y, i=1):
    '''Returns the vector whose entries are the number of directed k-cliques in the digraph.
    Parameters
    ----------
    X: (array) directed flag complex.
    Y: (array) directed flag complex.
    '''
    f = []
    n1 = len(X)
    n2 = len(Y)
       
    if i == 1:
        for k in range(n1):
            counter = 0
            for simplex in union_k(X, Y, k):
                if simplex in X[k]:
                    counter += 1
            f.append(counter)
    
    if i == 2:
        for k in range(n2):
            counter = 0
            for simplex in union_k(X, Y, k):
                if simplex in Y[k]:
                    counter += 1
            f.append(counter)
        
    return f


def count_weakly_q_connected_components(A, q=None):
    '''Returns all weakly q-connected components.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "weakly_connected_components()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    Struct_Vec = []
    wcc = list(nx.weakly_connected_components(G))
    
    for i in range(len(wcc)):
        Struct_Vec.append(i+1)
        
    return Struct_Vec
 
    
def count_strongly_q_connected_components(A, q=None):
    '''Returns all strongly q-connected components.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "strongly_connected_components()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    Struct_Vec = []
    scc = list(nx.strongly_connected_components(G))
    
    for i in range(len(scc)):
        Struct_Vec.append(i+1)
        
    return Struct_Vec


#Amount of n-paths in a path complex:
def count_n_paths(PC):
    '''Returns the number of n-paths.
    
    Parameters
    ----------
    PC: (array) Path complex.
    '''
    Qtd = []
    for i in range(len(PC)):
        Qtd.append(len(PC[i]))
    return Qtd


#Amount of directed n-cliques in a DFC:
def count_directed_n_cliques(DFC):
    '''Returns the number of directed n-cliques.
    
    Parameters
    ----------
    DFC (array) Directed flag complex.
    '''
    Qtd = []
    for i in range(len(DFC)):
        Qtd.append(len(DFC[i]))
    return Qtd


def count_elmentary_directed_quasi_cliques(PC, n):
    '''Returns the mount of invariant elementary n-paths in a path complex.
    
    Parameters
    ----------
    PC (array) Path complex.
    '''
    qtd_n_dqc = len(enumerate_invariant_elementary_n_paths(PC, n))
    return qtd_n_dqc



def count_edges(A):
    '''Returns the number of edges.
    Parameters
    ----------
    A: (array) Adjacency matrix.
    '''
    n = 0
    for i in range(len(A)):
        for j in range(len(A)):
            if M[i,j] > 0:
                n+=1
    return n


def count_double_edges(A, q=None):
    '''Returns the number of double-edges.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    '''
    Qtd = len(enumerate_double_edges(A, q))
    return Qtd
        
    
def count_directed_n_cycles(A, q=None):
    '''Returns the number of directed cycles.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    '''
    Qtd = len(enumerate_directed_n_cycles(A, q))
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


def second_flag_topological_structure_vector(A, q=None):
    '''Returns a vector whose entries are the number of weakly connected components of the digraph.
    Parameters
    ----------
    A: (array) Adjacency matrix.
    '''
    wcc = count_weakly_q_connected_components(A, q)
    return wcc


def third_flag_topological_structure_vector(A, q=None):
    '''Returns a vector whose entries are the number of strongly connected components of the digraph.
    Parameters
    ----------
    A: (array) Adjacency matrix.
    '''
    scc = count_strongly_q_connected_components(A, q)
    return scc


def fourth_flag_topological_structure_vector(A, k=0):
    '''Returns the vector whose entries are the Betti numbers associated to the filtration.
    Parameters
    ----------
    A: adjacency matrix.
    k: (integer) k-th homology group.
    '''
    v_k = betti_numbers(A, k=k)
    return v_k


def fifth_flag_topological_structure_vector(A, k=1):
    '''Returns the vector whose entries are the length of each barcode (li = di - bi).
    Parameters
    ----------
    A: adjacency matrix.
    k: (integer) k-th homology group.
    '''
    bc_k = barcode_length(A, k=k)
    return bc_k


#------------------------------

def first_path_topological_structure_vector(PC):
    '''Returns a vector whose entries are the number of invariant elementary directed 
    quasi-clique present in the digraph.
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
    '''Returns a vector whose entries are the length of each n-path in the digraph.
    Parameters
    ----------
    PC: (array) path complex.
    '''
    count_paths = count_n_paths(PC)
    return count_paths



