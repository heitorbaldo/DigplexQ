'''
Distances. Kernels. Similarity Comparison.

'''

import math
import numpy as np
import networkx as nx
from numpy import linalg as LA

from digplexq.persistent_homology import *
from digplexq.substructure_enumeration import *
from digplexq.utils import *


__all_ = [
    "first_topological_distance",
    "second_topological_distance",
    "fourth_topological_distance",
    "fifth_topological_distance",
    "histogram_cosine_kernel",
    "jaccard_kernel",
]


#----- Topological Distances -----

def first_topological_distance(M1, M2, comp="flag"):
    '''Returns the 1st topological distance.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    if comp == "flag":
        DFC1 = DirectedFlagComplex(M1, split='by_dimension_without_nodes')
        DFC2 = DirectedFlagComplex(M2, split='by_dimension_without_nodes')
        v1 = first_flag_topological_structure_vector(DFC1)
        v2 = first_flag_topological_structure_vector(DFC2)
        
    if comp == "path":
        PC1 = PathComplex(M1, 5)
        PC2 = PathComplex(M2, 5)
        v1 = first_path_topological_structure_vector(PC1)
        v2 = first_path_topological_structure_vector(PC2)
      
    if len(v1) != len(v2):
        diff = abs(len(v1)-len(v2))
        if len(v1) > len(v2):  
            for i in range(diff):
                v2.append(0)
        else:
            for i in range(diff):
                v1.append(0)
    
    x = np.array(v1)
    y = np.array(v2)
    
    if np.all(x==0) == True and np.all(y==0) == True:
        return 0
    
    ftd = (LA.norm(abs(x - y), 2))/( LA.norm(abs(x), 2) + LA.norm(abs(y), 2) )
    
    return round(ftd, 5)


def second_topological_distance(M1, M2, comp="flag"):
    '''Returns the 2nd topological distance.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    if comp == "flag":
        DFC1 = DirectedFlagComplex(M1, split='by_dimension_without_nodes')
        DFC2 = DirectedFlagComplex(M2, split='by_dimension_without_nodes')
        v1 = second_flag_topological_structure_vector(DFC1)
        v2 = second_flag_topological_structure_vector(DFC2)
        
    if comp == "path":
        PC1 = PathComplex(M1, 5)
        PC2 = PathComplex(M2, 5)
        v1 = second_path_topological_structure_vector(PC1)
        v2 = second_path_topological_structure_vector(PC2)
      
    if len(v1) != len(v2):
        diff = abs(len(v1)-len(v2))
        if len(v1) > len(v2):  
            for i in range(diff):
                v2.append(0)
        else:
            for i in range(diff):
                v1.append(0)
    
    x = np.array(v1)
    y = np.array(v2)
    
    if np.all(x==0) == True and np.all(y==0) == True:
        return 0
    
    ftd = (LA.norm(abs(x - y), 2))/( LA.norm(abs(x), 2) + LA.norm(abs(y), 2) )
    
    return round(ftd, 5)



def fourth_topological_distance(M1, M2, comp="flag"):
    '''Returns the 4th topological distance.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    if comp == "flag":
        v1 = fourth_flag_topological_structure_vector(M1)
        v2 = fourth_flag_topological_structure_vector(M2)
        
    if comp == "path":
        v1 = fourth_path_topological_structure_vector(M1)
        v2 = fourth_path_topological_structure_vector(M2)
      
    if len(v1) != len(v2):
        diff = abs(len(v1)-len(v2))
        if len(v1) > len(v2):  
            for i in range(diff):
                v2.append(0)
        else:
            for i in range(diff):
                v1.append(0)
    
    x = np.array(v1)
    y = np.array(v2)
    
    if np.all(x==0) == True and np.all(y==0) == True:
        return 0
    
    ftd = (LA.norm(abs(x - y), 2))/( LA.norm(abs(x), 2) + LA.norm(abs(y), 2) )
    
    return round(ftd, 5)


def fifth_topological_distance(M1, M2, comp="flag"):
    '''Returns the 5th topological distance.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    if comp == "flag":
        v1 = fifth_flag_topological_structure_vector(M1)
        v2 = fifth_flag_topological_structure_vector(M2)
        
    if comp == "path":
        v1 = fifth_path_topological_structure_vector(M1)
        v2 = fifth_path_topological_structure_vector(M2)
      
    if len(v1) != len(v2):
        diff = abs(len(v1)-len(v2))
        if len(v1) > len(v2):  
            for i in range(diff):
                v2.append(0)
        else:
            for i in range(diff):
                v1.append(0)
    
    x = np.array(v1)
    y = np.array(v2)
    
    if np.all(x==0) == True and np.all(y==0) == True:
        return 0
    
    ftd = (LA.norm(abs(x - y), 2))/( LA.norm(abs(x), 2) + LA.norm(abs(y), 2) )
    
    return round(ftd, 5)



#----- Simplicial Kernels -----

def histogram_cosine_kernel(M1, M2):
    '''Returns the histogram cosine kernel (HCK) between two directed flag complexes.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''  
    if np.all(M1==0) == True or np.all(M2==0) == True:
        return 1.0
    
    X = DirectedFlagComplex(M1, split='by_dimension_with_nodes')
    Y = DirectedFlagComplex(M2, split='by_dimension_with_nodes')
    
    #v1 = f_count(X, Y, i=1)
    #v2 = f_count(X, Y, i=2)
    
    v1 = first_flag_topological_structure_vector(X)
    v2 = first_flag_topological_structure_vector(Y)
    
    if len(v1) != len(v2):
        diff = abs(len(v1)-len(v2))
        if len(v1) > len(v2):  
            for i in range(diff):
                v2.append(0)
        else:
            for i in range(diff):
                v1.append(0)
    
    x = np.array(v1)
    y = np.array(v2)
    
    if np.all(x==0) == True and np.all(y==0) == True:
        return 0

    dot_product = np.dot(x, y)
    norm1 = np.linalg.norm(x)
    norm2 = np.linalg.norm(y)

    cosine = dot_product / (norm1 * norm2)
    return round(cosine, 6)


def jaccard_kernel(M1, M2):
    '''Returns the Jaccard kernel between two directed flag complexes.
    
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    if np.all(M1==0) == True or np.all(M2==0) == True:
        return 1.0
    
    X = DirectedFlagComplex(M1, split='none_without_nodes')
    Y = DirectedFlagComplex(M2, split='none_without_nodes')
        
    Diff = []
    Inter = []
    
    n1 = len(X)
    n2 = len(Y)
    
    for simplex in Y:
        if simplex not in X:
            Diff.append(simplex)
    
    for simplex in Y:
        if simplex in X:
            Inter.append(simplex)
    
    Union = X + Diff
    J = 1 - len(Inter)/len(Union)
    return round(J, 5)

    