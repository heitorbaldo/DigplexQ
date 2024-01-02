'''
Digraph-Based Complexes. 
'''

import numpy as np
import networkx as nx
import warnings
warnings.filterwarnings("ignore")

from digplexq.utils import *

__all__ = [
    "retrieve_k_cliques",
    "intersection",
    "k_simplices",
    "ordered_simplex",
    "retrieve_n_paths",
    "n_paths",
    "split_path",
    "remove_repeated_paths",
    "FlagTournaplex",
    "DirectedFlagComplex",
    "PathComplex",
]


#----------------------------

def retrieve_k_cliques(M, k):
    """Returns all k-cliques of the underlying undirected graph.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    k: (integer) Maximum dimension of the cliques.
       
    Notes
    -----
    Based on the clique enumerator algorithm [1].
    [1] Clique Enumerator algorithm, Zhang, et al. (2005).
    """        
    W = []
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    H = G.to_undirected()
    C = list(nx.enumerate_all_cliques(H))
    for i in range(len(C)):
        if len(C[i]) == k+1:
            W.append(C[i])
        else:
            pass
    return W


def intersection(cl1, cl2):
    """Returns all vertices belonging to the 
    intersection of two cliques.
    
    Parameters
    ----------
    cl1: (array) Clique 1.
    cl2: (array) Clique 2.
    """
    intersec = [i for i in cl1 if i in cl2]
    return intersec


def k_simplices(DFC, k):
    '''Returns all k-dimensional simplices.
    
    Parameters
    ----------
    DFC: (array) Directed Flag Complex.
    k: (integer) Dimension of the simplices.
    '''
    kSimp = []
    for simplex in DFC:
        if len(simplex) == k:
            kSimp.append(simplex)
    return kSimp


def ordered_simplex(M, simplex):
    '''Returns the ordered simplex.
    
    Parameters
    ----------
    M: (array) Adjacency matrix
    simplex: (array) Simplex.
    '''
    P=[]
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    H = G.subgraph(simplex)
    for i in simplex:
        for j in simplex:
            if i != j:
                path = list(nx.all_simple_paths(H, source=i, target=j))
                if path != []:
                    P.append(path)
    for paths in P:
        for k in range(len(paths)):
            if len(paths[k]) == len(simplex):
                ord_simplex = paths[k]
    return ord_simplex


#-------------------------

def retrieve_n_paths(M, s, t, n):
    W = []
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    P = list(nx.all_simple_paths(G, source=s, target=t, cutoff=n))
    for i in range(len(P)):
        if len(P[i]) == n:
            W.append(P[i])
        else:
            pass
    return W


def n_paths(M, n):
    P = []
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    for i in range(len(M)):
        for j in range(len(M)):
            P.append(retrieve_n_paths(M, i, j, n))
            
    P = [P[k] for k in range(len(P)) if len(P[k]) != 0]
    return P


def split_path(p):
    for i in range(len(p)):
        if len(p[i]) > 1:
            for path in p[i]:
                p.append([path])
    return p


def remove_repeated_paths(paths):
    l = len(paths)
    P = []
    for j in range(l):
        if len(paths[j]) == 1:
            P.append(paths[j])
    return P



#----- Flag Tournaplex -----

def FlagTournaplex(M):
    """Returns the flag tournaplex of a given digraph.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    """
    if isinstance(M, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    FT=[]
    for j in range(len(M)):
        if retrieve_k_cliques(M, j) != []:
            FT.append(retrieve_k_cliques(M, j))
    return FT



#----- Directed Flag Complex -----

def DirectedFlagComplex(M, split='by_dimension_without_nodes'):
    """Returns the directed flag complex of a given digraph.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    split: (string) 'by_dimension', 'nodes', 'none'.
    """
    if isinstance(M, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    if np.all(M==0) == True:
        return []
    
    Len = []
    DFC = []
    nodes = []
    Ordered_DFC = []
    With_Nodes = []
    Final_DFC = []
    G = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
    H = G.to_undirected()
    Cliques = list(nx.enumerate_all_cliques(H))
    for clique in Cliques:
        D = nx.from_numpy_matrix(M, create_using=nx.DiGraph())
        S = G.subgraph(clique)
        if nx.is_directed_acyclic_graph(S) == True:
                DFC.append(clique)
    
    for simplex in DFC:
        if len(simplex) > 1:
            o_simplex = ordered_simplex(M, simplex)
            Ordered_DFC.append(o_simplex)
    
    for simplex in Ordered_DFC:
        Len.append(len(simplex))
    
    if Len != []:
        for k in range(2, max(Len)+1):
            Final_DFC.append(k_simplices(Ordered_DFC, k))
        
        for clique in Cliques:   
            if len(clique) == 1:
                nodes.append(clique)

        if split == 'by_dimension_with_nodes':
            With_Nodes.append(nodes)
            for k in range(2, max(Len)+1):
                With_Nodes.append(k_simplices(Ordered_DFC, k))
            return With_Nodes
    
        if split == 'by_dimension_without_nodes':
            return Final_DFC
    
        if split == 'none_without_nodes':
            return Ordered_DFC
    else:
        return []


    
#----- Path Complex -----
    
def PathComplex(M, n):
    '''Returns the path complex of a digraph.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    n: (integer) maximum (n-1)-path length.
    '''
    if isinstance(M, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")
        
    if isinstance(n, int) == False:
        raise TypeError("n must be an integer.")
        
    if n < 2:
        raise ValueError("n must be an integer greater or equal to 2.")
        
    if np.all(M==0) == True:
        return []
    
    P = []
    PC = []
    for i in range(2, n+1):
        if n_paths(M, i) != []:
            P.append(n_paths(M, i))
            
    for k in range(len(P)):
        Pk = remove_repeated_paths(split_path(P[k]))
        PC.append(Pk)
        
    return PC

    
       
#----- Dowker Complex -----
#TODO

#----- Journey Complex -----
#TODO
    


