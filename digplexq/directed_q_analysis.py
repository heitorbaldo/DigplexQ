'''
Directed Q-Analysis on Directed Flag Complexes.
'''

import warnings
warnings.filterwarnings("ignore")
import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize

from digplexq.utils import *
from digplexq.simplicial_weights import *
from digplexq.digraph_based_complexes import *


__all__ = [
    "k_faces",
    "face_maps",
    "q_i_j_connectivity",
    "is_in_q_near",
    "is_out_q_near",
    "ComplementElements",
    "MaximalSimplices",
    "q_index",
    "q_adjacency_matrix",
    "weighted_q_adjacency_matrix",
    "adjacency_matrices_wcc",
    "fast_q_adjacency_matrix",
    "MaximalSimplices_WIN",
    "lower_q_index",
    "lower_q_adjacency_matrix",
    "fast_lower_q_adjacency_matrix",
]


#----- Face maps -----

def k_faces(DFC, simplex, k):
    '''Returns all k-faces of a given simplex.
    
    Parameters
    ----------
    DFC: (array) Directed flag complex (DFC = DirectedFlagComplex(M, "by_dimension_with_nodes")).
    simplex: (array) simplex.
    k: (integer) Dimension of the k-faces.
    '''
    if len(simplex)-1 <= k:
        raise ValueError("k must be less than or equal to the dimension of the simplex.")
        
    kFaces = []
    for face in DFC[k]:
        if(intersection(simplex, face) == face):
            kFaces.append(face)  
    return kFaces


def face_maps(DFC, simplex, k):
    '''Returns the indices and the faces.
    
    Parameters
    -------
    DFC: (array) Directed flag complex.
    simplex: (array) Simplex.
    k: (integer) Dimension of the face.
    '''
    ordered_faces = []
    Faces = k_faces(DFC, simplex, k)
    for i in range(len(simplex)):
        for face in Faces:
            if simplex[i] not in face:
                ordered_faces.append(face)
                
    d = {"i":[],"d_i":[]}
    for i in range(len(ordered_faces)):
        d["i"].append(i)
        d["d_i"].append(ordered_faces[i])
    return d



#----- (q,i,j)-connectivity -----

def q_i_j_connectivity(DFC, simplex1, simplex2):
    '''Returns the (q,i,j)-connectivity.
    
    Parameters
    ----------
    DFC: (array) Directed flag complex.
    simplex1, simplex2: (array) Simplices.
    '''
    k1 = len(simplex1)-2
    k2 = len(simplex2)-2
    Index1 = face_maps(DFC, simplex1, k1)["i"] 
    Index2 = face_maps(DFC, simplex2, k2)["i"] 
    Faces1 = face_maps(DFC, simplex1, k1)["d_i"]
    Faces2 = face_maps(DFC, simplex2, k2)["d_i"]
    Q_i_j = []
    
    if (intersection(simplex1, simplex2) == simplex1 or
        intersection(simplex1, simplex2) == simplex2):
        for i in range(len(Faces1)):
            for j in range(len(Faces2)):
                I = intersection(Faces1[i], Faces2[j])
                if I != []:
                    q_i_j = [len(I)-1, Index1[i], Index2[j]]
                    Q_i_j.append(q_i_j)
    else:
        for i in range(len(Faces1)):
            for j in range(len(Faces2)):
                I = intersection(Faces1[i], Faces2[j])
                if I != []:
                    q_i_j = [len(I)-1, Index1[i], Index2[j]]
                    Q_i_j.append(q_i_j)       
    return Q_i_j


def is_in_q_near(DFC_nodes, q, simplex1, simplex2):
    '''Returns True if simplex1 is in-q-near to simplex2.
    
    Parameters
    -------
    DFC_nodes: (array) Directed flag complex.
    q: (integer) Level of the q-connectivity.
    simplex1, simplex2: (array) Simplices.
    '''
    q_i_j = q_i_j_connectivity(DFC_nodes, simplex1, simplex2)
    
    c = []
    if q_i_j != []:
        for conn in q_i_j:
            if conn[0] == q and conn[1] >= conn[2]:
                c.append(1)
            else:
                c.append(0)
        if sum(c) > 0:
            return True
        else:
            return False
    else:
        return False


def is_out_q_near(DFC_nodes, q, simplex1, simplex2):
    '''Returns True if simplex1 is out-q-near to simplex2.
    
    Parameters
    ----------   
    DFC_nodes: (array) Directed flag complex.
    q: (integer) Level of the q-connectivity.
    simplex1, simplex2: (array) Simplices.
    '''
    q_i_j = q_i_j_connectivity(DFC_nodes, simplex1, simplex2)
    
    c = []
    if q_i_j != []:
        for conn in q_i_j:
            if conn[0] == q and conn[1] <= conn[2]:
                c.append(1)
            else:
                c.append(0)
        if sum(c) > 0:
            return True
        else:
            return False
    else:
        return False



#----- Maximal Simplices -----

def ComplementElements(list1, list2):
    Complement = []
    for i in list1:
        if i not in list2:
            Complement.append(i)
    return Complement


def MaximalSimplices(DFC):
    '''Returns all maximal simplices of a directed flag complex.
    
    Parameters
    -------
    DFC: (array) Directed flag complex (DirectedFlagComplex(M, "by_dimension_without_nodes")).
    '''
    if DFC == []:
        return []
    
    Faces = []
    Split_Faces = []
    Maximal = []
    Maximal_None = []
    Len = []
        
    for i in range(0, len(DFC)-1):
        for j in range(len(DFC[i])):
            for k in range(len(DFC[i+1])):
                if(set(DFC[i][j]).issubset(set(DFC[i+1][k])) == True):
                    Faces.append(DFC[i][j])
                    break
    
    for face in Faces:
        Len.append(len(face))
    
    if Len != []:
        for k in range(1, max(Len)+1):
            Split_Faces.append(k_simplices(Faces, k))
    else:
        Split_Faces = Faces
   
    for l in range(len(Split_Faces)):
        Complements = ComplementElements(DFC[l], Split_Faces[l])
        if Complements != []:
            Maximal.append(Complements)
    
    Maximal.append(DFC[len(Split_Faces)])
    return Maximal

    

#----- Directed q-Graph -----    

def q_index(MaxSimp, q):
    '''Returns the index associated to q.
    
    Parameters
    ---------
    MaxSimp: (array) Maximal simplices.
    q: (integer) Level of the q-connectivity.
    '''
    d = len(MaxSimp)
    dim_simp = len(MaxSimp[0][0])
    max_dim_simp = len(MaxSimp[d-1][0])
    diff = max_dim_simp - q
    
    if q+d < dim_simp:
        index = 0
    
    sum_total = 0
    for i in range(d):
        sum_total += len(MaxSimp[i])
    
    if(diff >= 0 and q >= dim_simp and dim_simp == 1):
        sum_partial = 0
        for j in range(q, d):
            sum_partial += len(MaxSimp[j])
        index = sum_total - sum_partial 
        
    elif(diff >= 0 and q >= dim_simp and dim_simp > 1):
        sum_partial = 0
        for j in range(q-dim_simp+1, d):
            sum_partial += len(MaxSimp[j])
        index = sum_total - sum_partial
    else:
        index = 0
    return index


def q_adjacency_matrix(DFC_dim_nodes, q):
    '''Returns the q-adjacency matrix of a directed flag complex.
    
    Parameters
    -------
    DFC_dim_none: (array) Directed flag complex (without nodes).
    DFC_dim_nodes: (array) Directed flag complex (with nodes).
    q: (integer) Level of the q-connectivity.
    '''
    MS = MaximalSimplices(DFC_dim_nodes)
    MaxSimp = connect_array(MS)
    Len = len(MaxSimp)
    
    if MS == []:
        return np.array([[0]])
    
    H = np.zeros((Len, Len))
    for i in range(len(MaxSimp)):
        for j in range(len(MaxSimp)):
            Q = []
            q_i_j = q_i_j_connectivity(DFC_dim_nodes, MaxSimp[i], MaxSimp[j])
            if i != j and q_i_j != []:
                for conn in q_i_j:
                    Q.append(conn[0])
                Qmax = max(Q)
                if Qmax >= q and is_out_q_near(DFC_dim_nodes, Qmax, MaxSimp[i], MaxSimp[j]) == True:
                    H[i, j] = 1
                else:
                    H[i, j] = 0
            else:
                H[i, j] = 0
    
    
    d = len(MS)
    max_dim_simp = len(MS[d-1][0])
    if max_dim_simp < q+1:
        Hq = np.array([[0]])
    else:
        Hq = H[q_index(MS, q):, q_index(MS, q):]  
    return Hq
    
    
       
def weighted_q_adjacency_matrix(M, DFC_dim_nodes, q):
    '''Returns the weighted q-adjacency matrix of a directed flag complex.
    
    Parameters
    -------
    DFC_dim_none: (array) Directed flag complex (without nodes).
    DFC_dim_nodes: (array) Directed flag complex (with nodes).
    q: (integer) Level of the q-connectivity.
    '''
    MS = MaximalSimplices(DFC_dim_nodes)
    MaxSimp = connect_array(MS)
    Len = len(MaxSimp)
    
    if MS == []:
        return np.array([[0]])
    
    H = np.zeros((Len, Len))
    for i in range(len(MaxSimp)):
        for j in range(len(MaxSimp)):
            Q = []
            q_i_j = q_i_j_connectivity(DFC_dim_nodes, MaxSimp[i], MaxSimp[j])
            if i != j and q_i_j != []:
                for conn in q_i_j:
                    Q.append(conn[0])
                Qmax = max(Q)
                if Qmax >= q and is_out_q_near(DFC_dim_nodes, Qmax, MaxSimp[i], MaxSimp[j]) == True:
                    H[i, j] = simplex_weight(M, MaxSimp[i]) #Node-to-edge transformation: h(i,j) = w(i).
                else:
                    H[i, j] = 0
            else:
                H[i, j] = 0
    
    d = len(MS)
    max_dim_simp = len(MS[d-1][0])
    if max_dim_simp < q+1:
        Hq = np.array([[0]])
    else:
        Hq = H[q_index(MS, q):, q_index(MS, q):]
        Hq = normalize(Hq, axis=1, norm='max')  #Row-normalized matrix w.r.t. the maximum norm.
    
    return Hq

    
    
def adjacency_matrices_wcc(A):
    '''Returns the adjacency matrices of the weakly connected components.
    
    Parameters
    ----------
    A: (array) Adjacency matrix.
    '''
    if isinstance(A, np.ndarray) == False:
        raise TypeError("Input must be a square matrix.")
    
    WCC = []
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    wcc = list(nx.weakly_connected_components(G))
    for c in wcc:
        S = G.subgraph(c)
        WCC.append(nx.adjacency_matrix(S).toarray())
    return WCC
      
    

def fast_q_adjacency_matrix(M, q):
    '''Returns the q-adjacency matrix.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    q: (integer) Level of the q-connectivity.
    '''
    DFC_dim_nodes = DirectedFlagComplex(M, "by_dimension_with_nodes")
    Hq = q_adjacency_matrix(DFC_dim_nodes, q)
    return Hq



#----- Lower Directed q-Graph -----  

def MaximalSimplices_WIN(DFC):
    '''Returns all maximal simplices (without isolated nodes) of a directed flag complex.
    
    Parameters
    -------
    DFC: (array) Directed flag complex (DirectedFlagComplex(M, "by_dimension_without_nodes")).
    '''
    if DFC == []:
        return []
    
    Faces = []
    Split_Faces = []
    Maximal = []
    Maximal_None = []
    Len = []
    for i in range(0, len(DFC)-1):
        for j in range(len(DFC[i])):
            for k in range(len(DFC[i+1])):
                if(set(DFC[i][j]).issubset(set(DFC[i+1][k])) == True):
                    Faces.append(DFC[i][j])
                    break
    
    for face in Faces:
        Len.append(len(face))
    
    if Len != []:
        for k in range(2, max(Len)+1):
            Split_Faces.append(k_simplices(Faces, k))
    else:
        Split_Faces = Faces
   
    for l in range(len(Split_Faces)):
        Complements = ComplementElements(DFC[l], Split_Faces[l])
        if Complements != []:
            Maximal.append(Complements)
    
    Maximal.append(DFC[len(Split_Faces)])
    return Maximal

    
def lower_q_index(MaxSimp, q):
    '''Returns the index associated to q.
    
    Parameters
    ---------
    MaxSimp: (array) Maximal simplices.
    q: (integer) Level of the q-connectivity.
    '''
    index = 0
    
    d = len(MaxSimp)
    dim_simp = len(MaxSimp[0][0])
    max_dim_simp = len(MaxSimp[d-1][0])
    diff = max_dim_simp - q
    
    if q+d < dim_simp:
        index = 0  
    
    if(dim_simp <= q+d and diff > 1):
        m = q+2-dim_simp
        for p in range(1,m+1):
            index += len(MaxSimp[p-1])  
    else:
        index = 0    
    return index


def lower_q_adjacency_matrix(DFC_dim_none, DFC_dim_nodes, q):
    '''Returns the (lower) q-adjacency matrix of a directed flag complex.
    
    Parameters
    -------
    DFC_dim_none: (array) Directed flag complex (without nodes).
    DFC_dim_nodes: (array) Directed flag complex (with nodes).
    q: (integer) Level of the q-connectivity.
    '''
    MS = MaximalSimplices_WIN(DFC_dim_none)
    MaxSimp = connect_array(MS)
    Len = len(MaxSimp)
    
    if MS == []:
        return np.array([[0]])
    
    H = np.zeros((Len, Len))
    for i in range(len(MaxSimp)):
        for j in range(len(MaxSimp)):
            Q = []
            q_i_j = q_i_j_connectivity(DFC_dim_nodes, MaxSimp[i], MaxSimp[j])
            if i != j and q_i_j != []:
                for conn in q_i_j:
                    Q.append(conn[0])
                Qmax = max(Q)
                if Qmax >= q and is_out_q_near(DFC_dim_nodes, Qmax, MaxSimp[i], MaxSimp[j]) == True:
                    H[i, j] = 1
                else:
                    H[i, j] = 0
            else:
                H[i, j] = 0
    
    if len(DFC_dim_none) <= q:
        Hq = np.array([[0]])
    
    else:
        Hq = H[lower_q_index(MS, q):, lower_q_index(MS, q):]  
    return Hq



def fast_lower_q_adjacency_matrix(M, q):
    '''Returns the lower q-adjacency matrix.
    
    Parameters
    ----------
    M: (array) Adjacency matrix.
    q: (integer) Level of the q-connectivity.
    '''
    DFC_dim_nodes = DirectedFlagComplex(M, "by_dimension_with_nodes")
    DFC_dim_none = DirectedFlagComplex(M, "by_dimension_without_nodes")
    Hq = lower_q_adjacency_matrix(DFC_dim_none, DFC_dim_nodes, q)
    return Hq
