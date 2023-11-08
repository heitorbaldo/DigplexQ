'''
Similarity Comparison of Static and Dynamic Path Complexes

'''

import numpy as np
import networkx as nx



#----- Structure Vectors -----

def struct_vector(P):
    
    S = 
    
def struct_feature_vec(P):
    S = struct_vector(P)
    Q = 
    W= 
    
    V = [S, Q, W]
    return V


def total_struct_dist(P1, P2):
    
    Tds = struct_feature_vec(P1) - struct_feature_vec(P2)
    
    return Tds


#----- Topological Distances -----

def n_topological_digraph_distance(x, y):1
    '''
    x: structure vector for X.
    y: structure vector for Y.
    '''
    T_n_ftdd = (LA.norm(abs(np.array(x) - np.array(y)), 2))/(LA.norm(abs(np.array(x)), 2)+LA.norm(abs(np.array(y)), 2))
    return round(T_n_ftdd, 6)






    