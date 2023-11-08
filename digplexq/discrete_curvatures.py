'''
Discrete Curvatures: Ollivier-Ricci and Forman-Ricci Curvatures.
'''

import numpy as np
import math as m
from digplexq.simplicial_weights import *
from digplexq.utils import *
from digplexq.digraph_based_complexes import *

__all__ = [
    "incident_edges_simplex",
    "incident_edges_node",
    "forman_ricci_curvature",
]



def incident_edges_simplex(DFC, simplex):
    '''Returns the set of incident edges.
    Parameters
    '''
    inc_edges = []
    for simp in DFC[0]:
        for i in range(len(simplex)):
            if simp[0] == simplex[i] or simp[1] == simplex[i]:
                inc_edges.append(simp)
            else:
                pass
    return inc_edges


def incident_edges_node(DFC, simplex, k):
    '''Returns the set of incident edges.
    Parameters
    ---------
    DFC:
    simplex:
    k: node
    '''
    inc_edges_node = []
    IE = incident_edges_simplex(DFC, simplex)
    
    for edge in IE:
        if edge[0] == k or edge[1] == k:
                inc_edges.append(edge)
        else:
            pass
    return inc_edges_node


def forman_ricci_curvature(M, DFC, simplex, weight_func='max'):
    '''Returns the Forman-Ricci curvature of a weighted simplex or weighted path.
    Parameters
    ---------
    M: NumPy matrix.
    simplex: NumPy array.
    Returns
    ---------
    '''
    S1_k = S2_k = 0
    w_simplex = simplex_weight(M, simplex, weight_func)
    w_nodes = node_weights_simplex(M, simplex, weight_func)
    sum_node_weights = sum(w_nodes)
    IE = incident_edges_simplex(DFC, simplex)

    for k in range(len(simplex)):
        S1_k += w_nodes[k]/w_simplex
        
    for k in range(len(simplex)):
        IE_k = incident_edges_node(DFC, simplex, k)
        for edge in IE_k:
            S2_k += w_nodes[k] / m.sqrt(w_simplex*simplex_weight(M, edge, weight_func))
    
    FRC = w_simplex*(S1_k - S2_k)
    return FRC
   


