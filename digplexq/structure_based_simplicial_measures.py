'''
Structure-Based simplicial measures.
Most of the measures implemented here are based on the respective measures from the Networkx library.
'''

import math
import numpy as np
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
    "q_communicability",
    "q_communicability_max",
    "q_returnability",
    "q_structural_entropy",
    "in_q_degree_distribution_entropy",
    "out_q_degree_distribution_entropy",
]


#----- Distance-based Measures -----

def average_shortest_q_walk_length(A, q=None):
    '''Returns the average shortest q-walk length of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "average_shortest_path_length()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    #If G is weakly connected:
    if nx.is_weakly_connected(G) == True:
        AV = nx.average_shortest_path_length(G, weight=None)
        return AV

    #If G is not weakly connected:
    else:
        AV = []
        wcc = adjacency_matrices_wcc(A)
        for c in wcc:
            S = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
            AV.append(nx.average_shortest_path_length(S, weight=None))
        Max_av = max(AV)
        return round(Max_av, 4)


def q_eccentricity(A, q=None):
    '''Returns the q-eccentricity of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "eccentricity()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    #If G is strongly connected:
    if nx.is_strongly_connected(G) == True:
        ecc = nx.eccentricity(G)
        return ecc

    #If G is not strongly connected:
    else:
        return math.inf


def q_diameter(A, q=None):
    '''Returns the q-diameter of a q-digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "diameter()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    #If G is strongly connected:
    if nx.is_strongly_connected(G) == True:
        Diam = nx.diameter(G)
        return Diam

    #If G is not strongly connected:
    else:
        return math.inf


def q_radius(A, q=None):
    '''Returns the q-radius of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "radius()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    #If G is strongly connected:
    if nx.is_strongly_connected(G) == True:
        rad = nx.radius(G)
        return rad

    #If G is not strongly connected:
    else:
        return math.inf


def q_density(A, q=None):
    '''Returns the q-density of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "density()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    s_den = nx.density(G)
    return s_den



def in_q_degree(A, q=None, i=None):
    '''Returns the in-q-degree of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    i: (integer) Node's label. If None, it returns the maximum
    in-degree of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "in_degree()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    if i != None:
        indeg = G.in_degree(i)
        return indeg  
    
    inDg = []
    for i in range(len(A)):
        inDg.append(G.in_degree(i))
    Max_in_dc = max(inDg)
    return Max_in_dc


def out_q_degree(A, i=None, q=None):
    '''Returns the out-q-degree of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    i: (integer) Node's label. If None, it returns the maximum
    out-degree of the digraph.
    q: (integer) Level of clique organization of the graph.

    Notes
    ----------
    Based on the Networkx's function "out_degree()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    if i != None:
        outdeg = G.out_degree(i)
        return outdeg  

    outDg = []
    for i in range(len(A)):
        outDg.append(G.out_degree(i))
    Max_out_dc = max(outDg)
    return Max_out_dc


    
#----- Measures of Centrality -----
    
def in_q_degree_centrality(A, q=None, results="max"):
    '''Returns the in-q-degree centrality of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    in-degree centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "in_degree_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if results == "nodes":
        in_dc = nx.in_degree_centrality(G)
        return in_dc
    
    if results == "max":
        in_dc = nx.in_degree_centrality(G)
        Max_in_dc = max(dict_to_array(in_dc))
        return Max_in_dc


def out_q_degree_centrality(A, q=None, results="max"):
    '''Returns the out-q-degree centrality of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    out-degree centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "out_degree_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if results == "nodes":
        out_dc = nx.out_degree_centrality(G)
        return out_dc
    
    if results == "max":
        out_dc = nx.out_degree_centrality(G)
        Max_out_dc = max(dict_to_array(out_dc))
        return Max_out_dc


def q_closeness_centrality(A, q=None, results="nodes", wf_improved=False):
    '''Returns the q-closeness centrality of the nodes of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    closeness centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "closeness_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    if results == "nodes":
        dscc = nx.closeness_centrality(G, wf_improved=wf_improved)
        return dscc
    
    if results == "max":
        #If G is weakly connected:
        if nx.is_weakly_connected(G) == True:
            dscc = nx.closeness_centrality(G, wf_improved=wf_improved)
            Max_dscc = max(dict_to_array(dscc))
            return Max_dscc

        #If G is not weakly connected:
        else:
            CC = []
            wcc = adjacency_matrices_wcc(A)
            for c in wcc:
                S = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
                dscc = nx.closeness_centrality(S, wf_improved=wf_improved)
                CC.append(max(dict_to_array(dscc)))
            Max_dscc = max(CC)
            return round(Max_dscc, 4)


def q_harmonic_centrality(A, q=None, results="max"):
    '''Returns the q-harmonic centrality of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    harmonic centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "harmonic_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    if results == "nodes":
        shc = nx.harmonic_centrality(G)
        return shc
    
    if results == "max":
        #If G is weakly connected:
        if nx.is_weakly_connected(G) == True:
            shc = nx.harmonic_centrality(G)
            Max_shc = max(dict_to_array(shc))
            return Max_shc

        #If G is not weakly connected:
        else:
            HC = []
            wcc = adjacency_matrices_wcc(A)
            for c in wcc:
                S = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
                shc = nx.harmonic_centrality(S)
                HC.append(max(dict_to_array(shc)))
            Max_shc = max(HC)
            return round(Max_shc, 4)



def q_betweenness_centrality(A, q=None, results="max", normalized=False):
    '''Returns the q-betweenness centrality of the nodes of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    betweenness centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "betweenness_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if results == "nodes":
        dsbc = nx.betweenness_centrality(G, normalized=normalized, weight=None)
        return dsbc
    
    if results == "max":
        #If G is weakly connected:
        if nx.is_weakly_connected(G) == True:
            dsbc = nx.betweenness_centrality(G, normalized=normalized, weight=None)
            Max_dsbc = max(dict_to_array(dsbc))
            return round(Max_dsbc, 4)

        #If G is not weakly connected:
        else:
            BCC = []
            wcc = adjacency_matrices_wcc(A)
            for c in wcc:
                S = nx.from_numpy_matrix(c, create_using=nx.DiGraph())
                dsbc = nx.betweenness_centrality(S, normalized=normalized, weight=None)
                BCC.append(max(dict_to_array(dsbc)))
            Max_dsbc = max(BCC)
            return round(Max_dsbc, 4)


def q_katz_centrality(A, q=None, results="max", alpha=0.1, beta=1.0, normalized=False):
    '''Returns the q-Katz centrality of the nodes of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    results: (string) If "max", it returns the maximum
    Katz centrality of the digraph.
    
    Notes
    ----------
    Based on the Networkx's function "katz_centrality_numpy()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    a1 = 1/max(np.linalg.eigvals(A))
    if a1 < 0.1:
        alpha = np.real(a1)/3
    else:
        alpha = 0.1

    if results == "nodes":
        katz = nx.katz_centrality_numpy(G, alpha=alpha, beta=beta, normalized=normalized, weight=None)
        return katz
    
    if results == "max":
        katz = nx.katz_centrality_numpy(G, alpha=alpha, beta=beta, normalized=normalized, weight=None)
        Max_katz = max(dict_to_array(katz))
        return round(Max_katz, 4)


def global_q_reaching_centrality(A, q=None, normalized=True):
    '''Returns the q-reaching centrality of a digraph.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "global_reaching_centrality()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    gsrc = nx.global_reaching_centrality(G, normalized=normalized)
    return round(gsrc, 5)




#----- Efficiency and Global Efficiency -----

def q_efficiency(A, i, q=None):
    '''Returns the q-efficiency of a digraph's node.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    i: (integer) Node's label.
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
    
    n = len(A)
    nE_i = 0
    
    for j in range(n):
        if j != i:
            try:
                SPL = nx.shortest_path_length(G, source=i, target=j, weight=None, method='dijkstra')
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


def global_q_efficiency(A, q=None):
    '''Returns the global q-efficiency of a digraph.
    
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
    
    n = len(A)
    nGE = 0
    for i in range(n):
        nGE += q_efficiency(A, i)
    
    GE = nGE/n
    return round(GE, 5)



#----- Segregation Measures -----

def average_q_clustering_coefficient(A, q=None):
    '''Returns the average q-clustering coefficient of a q-digraph.

    Parameters
    ----------
    A: (array) Adjacency matrix.
    q: (integer) Level of clique organization of the graph.
    
    Notes
    ----------
    Based on the Networkx's function "average_clustering()".
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    ACC = nx.average_clustering(G, nodes=None, weight=None, count_zeros=True)
    return round(ACC, 4)



def in_q_degree_rich_club_coefficient(A, k=6, q=None):
    '''Returns the in-q-degree rich-club coefficient phi(k) = E^{in}_k/N^{in}_k(N^{in}_k-1).
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    k: (integer) Parameter.
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
    
    N_k = 0
    n = len(A)
    nodes = []
    
    for i in range(n):
        indeg = G.in_degree(i)
        if indeg > k:
            N_k += 1
            nodes.append(i)
        else:
            pass
    
    if N_k > 1: 
        H = G.subgraph(nodes)
        E_K = G.number_of_edges()
        RCC = E_K/(N_k*(N_k - 1))
        return round(RCC, 5)
    else:
        return 0

    
def out_q_degree_rich_club_coefficient(A, k=6, q=None):
    '''Returns the out-q-degree rich-club coefficient phi(k) = E^{out}_k/N^{out}_k(N^{out}_k-1).

    Parameters
    ----------
    A: (array) Adjacency matrix.
    k: (integer) Parameter.
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
    
    N_k = 0
    n = len(A)
    nodes = []
    
    for i in range(n):
        outdeg = G.out_degree(i)
        if outdeg > k:
            N_k += 1
            nodes.append(i)
        else:
            pass
    
    if N_k > 1: 
        H = G.subgraph(nodes)
        E_K = G.number_of_edges()
        RCC = E_K/(N_k*(N_k - 1))
        return round(RCC, 5)
    else:
        return 0



#----- Communicability -----

def q_communicability(A, i, j, q=None):
    '''Returns the q-communicability of the nodes i and j.
    
    Parameters
    ----------
    A: Adjacency matrix.
    i: (integer) Node's label.
    j: (integer) Node's label.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a NumPy square matrix.")

    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())

    if nx.is_empty(G) == True:
        return 0
    
    if q != None:
        A = fast_q_adjacency_matrix(A, q)
        G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    
    G = expm(A)
    G_ij = G[i,j]
    return round(G_ij, 5)


def q_communicability_max(A, q=None):
    '''Returns the maximum q-communicability among all the nodes i and j.
    
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
    
    V = []
    for i in range(len(A)):
        for j in range(len(A)):
            Gij = q_communicability(A, i, j, q)
            V.append(Gij)
    
    Max_Gij = max(V)
    return round(Max_Gij, 5)



#----- Returnability -----

def q_returnability(A, q=None, normalized=True):
    '''Returns the q-returnability of a q-digraph.

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

def q_structural_entropy(A, q=None):
    '''Returns the q-structural entropy of a digraph.
    
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
    
    n = len(A)
    sum_Q = 0
    H = 0
    sum_Qi = []
    
    for i in range(n):
        Q_i = 0
        for j in range(n):
            Q_i += q_communicability(A, i, j, q)
        sum_Qi.append(Q_i)
    
    for i in range(n):
        for j in range(n):
            sum_Q += q_communicability(A, i, j, q)
    
    if sum_Q == 0:
        return 0
    
    for k in range(n):
        H += (sum_Qi[k]/sum_Q)*math.log2(sum_Qi[k]/sum_Q)
    
    return round(-H, 5)


def in_q_degree_distribution_entropy(A, q=None):
    '''Returns the in-q-degree distribution entropy of a digraph.
    
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
    
    D = []
    H = 0
    n = len(A)
    
    for k in range(n):
        delta_k = 0
        for i in range(n):
            if G.in_degree(i) == k:
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


def out_q_degree_distribution_entropy(A, q=None):
    '''Returns the out-q-degree distribution entropy of a digraph.
    
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
    
    D = []
    H = 0
    n = len(A)
    
    for k in range(n):
        delta_k = 0
        for i in range(n):
            if G.out_degree(i) == k:
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