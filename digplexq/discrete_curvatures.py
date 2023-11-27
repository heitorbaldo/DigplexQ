'''
Discrete Curvatures: Ollivier-Ricci and Forman-Ricci Curvatures.
'''

import numpy as np
import math

from digplexq.digraph_based_complexes import *
from digplexq.simplicial_weights import *
from digplexq.utils import *


__all__ = [
    "incoming_edges",
    "outgoing_edges",
    "incoming_edges_vertex",
    "outgoing_edges_vertex",
    "q_forman_ricci_curvature",
    "in_q_forman_ricci_curvature",
    "out_q_forman_ricci_curvature",
]


def incoming_edges(DFC, edge_curv):
    '''Returns the set of incoming edges.
    Parameters
    ----------
    '''
    in_edges = []
    for edge in DFC[0]:
        if edge[1] == edge_curv[0]:
                in_edges.append(edge)
        else:
            pass
    return in_edges


def outgoing_edges(DFC, edge_curv):
    '''Returns the set of outgoing edges.
    Parameters
    ----------
    '''
    out_edges = []
    for edge in DFC[0]:
        if edge[0] == edge_curv[1]:
                out_edges.append(edge)
        else:
            pass
    return out_edges

def incoming_edges_vertex(DFC, v):
    '''Returns the set of incoming edges.
    Parameters
    ----------
    '''
    in_edges = []
    for edge in DFC[0]:
        if edge[1] == v:
                in_edges.append(edge)
        else:
            pass
    return in_edges


def outgoing_edges_vertex(DFC, v):
    '''Returns the set of outgoing edges.
    Parameters
    ----------
    '''
    out_edges = []
    for edge in DFC[0]:
        if edge[0] == v:
                out_edges.append(edge)
        else:
            pass
    return out_edges
   

def q_forman_ricci_curvature(M, DFC, edge_curv, weight_func='max'):
    '''Returns the Forman-Ricci curvature of a weighted directed edge.
    Parameters
    ----------
    M: NumPy matrix.
    edge: NumPy array.
    '''
    S0_k = 0
    S1_k = 0
    S2_k = 0
    
    w_edge = simplex_weight(M, edge_curv, weight_func)
    w_nodes = node_weights_simplex(M, edge_curv, weight_func)
    
    IE = incoming_edges(DFC, edge_curv)
    OE = outgoing_edges(DFC, edge_curv)

    for k in range(2):
        S0_k += w_nodes[k]/w_edge
        
    for edge in IE:
        S1_k += w_nodes[0] / math.sqrt(w_edge*simplex_weight(M, edge, weight_func))
        
    for edge in OE:
        S2_k += w_nodes[1] / math.sqrt(w_edge*simplex_weight(M, edge, weight_func))
    
    FRC = w_edge*(S0_k - S1_k - S2_k)
    return round(FRC, 5)
   

def in_q_forman_ricci_curvature(M, DFC, v, weight_func='max'):
    '''Returns the in-q-Forman-Ricci curvature of a vertex.
    Parameters
    ----------
    M: NumPy matrix.
    v: integer.
    '''
    inFRC = 0
    IE = incoming_edges_vertex(DFC, v)
    
    for edge in IE:    
        inFRC += q_forman_ricci_curvature(M, DFC, edge, weight_func='max')
    
    return round(inFRC, 5)


def out_q_forman_ricci_curvature(M, DFC, v, weight_func='max'):
    '''Returns the in-q-Forman-Ricci curvature of a vertex.
    Parameters
    ----------
    M: NumPy matrix.
    v: integer.
    '''
    outFRC = 0
    IO = outgoing_edges_vertex(DFC, v)
    
    for edge in IO:    
        outFRC += q_forman_ricci_curvature(M, DFC, edge, weight_func='max')
    
    return round(outFRC, 5)

