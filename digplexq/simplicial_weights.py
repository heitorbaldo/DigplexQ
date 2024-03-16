'''
Node weights, simplex weights, path weights.
'''

import numpy as np


__all__ = [
    "node_weights_max",
    "node_weights_max_in_out",
    "node_weights_simplex",
    "simplex_weight",
]
        
    
def node_weights_max(M):
    '''Returns the node weights W(i) = max(W_max, W_diff).
    
    Parameters
    ----------
    M: (array) adjacency matrix.
    '''
    W_diff = []
    W_max = []
    Weights = []
    for i in range(len(M)):
        sw1 = sw2 = 0
        for j in range(len(M)):
            if i != j:
                sw1 += M[j,i] 
                sw2 += M[i,j]
        W_diff.append(abs(sw1-sw2))
     
    for i in range(len(M)):
        wmax_node = []
        for j in range(len(M)):
            wmax_node.append(M[i,j])
            wmax_node.append(M[j,i])
        W_max.append(max(wmax_node))
        
    for k in range(len(M)):
        Weights.append(max(W_max[k], W_diff[k]))
    return Weights


def node_weights_max_in_out(M):
    '''Returns the node weights W(i) = max(w_out_deg(i), w_in_deg(i)).
    
    Parameters
    ----------
    M: (array) adjacency matrix.
    '''
    if np.all(M==0) == True:
        return 0
    
    W_max = []
    W_max_in = []
    W_max_out = []
    for i in range(len(M)):
        sw1 = sw2 = 0
        for j in range(len(M)):
            if i != j:
                sw1 += M[j,i] 
                sw2 += M[i,j]
        W_max_in.append(sw1)
        W_max_out.append(sw2)
    
    for k in range(len(W_max_in)):
        W_max.append(max(W_max_out[k], W_max_in[k]))
    return  W_max


def node_weights_simplex(M, simplex, weight_func='max'):
    '''Returns the weights of the nodes that belong to the simplex.
    
    Parameters
    ----------
    M: (array) adjacency matrix.
    simplex: (array) simplex.
    weight_func: edge-to-node weight function.
    '''
    if np.all(M==0) == True:
        return 0
    
    w_max = []
    w_in_out = []
    for i in range(len(simplex)):
        w_max.append(node_weights_max(M)[simplex[i]])
        w_in_out.append(node_weights_max_in_out(M)[simplex[i]])
        
    if weight_func == 'max':
        return  w_max
    if weight_func == 'max_in_out':
        return  w_in_out


def simplex_weight(M, simplex, weight_func='max_in_out'):
    '''Returns the multiplication of all the weights of the nodes that belong to the simplex.
    
    Parameters
    ----------
    M: (array) adjacency matrix.
    simplex: (array) simplex.
    weight_func: edge-to-node weight function.
    '''
    if simplex == []:
        return 0
    
    w_max = w_in_out = 1
    for i in range(len(simplex)):
        w_max *= node_weights_max(M)[simplex[i]]
        w_in_out *= node_weights_max_in_out(M)[simplex[i]]
       
    w_in_edge = node_weights_max_in_out(M)[simplex[0]]
        
    if weight_func == 'max':
        return round(w_max, 2)
    if weight_func == 'max_in_out':
        return round(w_in_out, 2)
    if weight_func == 'in_edge':
        return round(w_in_edge, 2)


          

          