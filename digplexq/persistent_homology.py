'''
Persistent Homology for Directed Flag Complexes.
The codes in this file are based on the following packages: Persim, Giotto-TDA, Gudhi.
'''

import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform

import gudhi as gd
import persim as pers
from pyflagser import *
from gtda.homology import FlagserPersistence
from gtda.diagrams import PersistenceEntropy, PersistenceImage, BettiCurve, PairwiseDistance


__all__ = [
    "remove_inf",
    "persistence_diagram",
    "birth_death_betti",
    "barcode_length",
    "betti_numbers",
    "wasserstein_distance",
    "bottleneck_distance",
    "betti_distance",
    "landscape_distance",
    "silhouette_distance",
    "persistence_image_distance",
]

#---------------

def remove_inf(diag):
    R = diag
    R_without_inf = []
    for i in range(len(R[0])):
        if R[0][i][1] != math.inf:
            R_without_inf.append(R[0][i])
    return np.array(R_without_inf)


def persistence_diagram(M):
    '''Returns the persistence diagram.
    M: (array)
    Based on Flagser.
    '''
    SM = csr_matrix(M)
    diagram = FlagserPersistence(homology_dimensions=(0,1,2)).fit_transform_plot([SM]);
    return diagram


def birth_death_betti(M):
    '''Returns the triple [birth, death, dimension].
    M: (array)
    Based on Flagser.
    '''
    SM = csr_matrix(M)
    bdb = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM]);
    return bdb


def barcode_length(M):
    '''Returns the length of each barcode.
    '''
    L = []
    SM = csr_matrix(M)
    diag = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM]);
    n = len(diag[0])
    for i in range(n):
        l_i = round(abs(diag[0][i][1] - diag[0][i][0]), 5)
        L.append(l_i)
    
    return L


def betti_numbers(M, k=0):
    '''Returns the k-th Betti numbers associated with the filtration.
    M: adjacency matrix.
    '''
    SM = csr_matrix(M)
    diag = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM]);
    
    b = BettiCurve().fit_transform(diag);
    
    if k == 0:
        v = b[0][0]
    if k == 1:
        v = b[2][0]
    
    return v


def wasserstein_distance(M1, M2):
    '''Returns the Wasserstein distance between two digraphs.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    wd = pers.wasserstein(diagrams1[0], diagrams2[0])
    return round(wd, 5)


def bottleneck_distance(M1, M2):
    '''Returns the bottleneck distance.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    bnd = pers.bottleneck(diagrams1[0], diagrams2[0])
    return round(bnd, 5)


def betti_distance(M1, M2):
    '''Returns the L2 distance between two Betti curves.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    
    diag1, diag2 = map(remove_inf, (diagrams1, diagrams2))
    
    d = PairwiseDistance(metric='betti',order=None, metric_params={'p':2.0})
    bd1 = d.fit(X=[diag1])
    bd2 = bd1.transform([diag2])
    
    if len(bd2[0][0]) == 2:
        return round(bd2[0][0][1], 5)
    else:
        return round(bd2[0][0][0], 5)


def landscape_distance(M1, M2):
    '''Returns the L2 distance between two persistence landscapes.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    
    diag1, diag2 = map(remove_inf, (diagrams1, diagrams2))
    
    d = PairwiseDistance(metric='landscape',order=None, metric_params={'p':2.0})
    ld1 = d.fit(X=[diag1])
    ld2 = ld1.transform([diag2])
    
    if len(ld2[0][0]) == 2:
        return round(ld2[0][0][1], 5)
    else:
        return round(ld2[0][0][0], 5)
    
    
def silhouette_distance(M1, M2):
    '''Returns the L2 distance between two silhouettes.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    
    diag1, diag2 = map(remove_inf, (diagrams1, diagrams2))
    
    d = PairwiseDistance(metric='silhouette',order=None, metric_params={'p':2.0})
    ld1 = d.fit(X=[diag1])
    ld2 = ld1.transform([diag2])
    
    if len(ld2[0][0]) == 2:
        return round(ld2[0][0][1], 5)
    else:
        return round(ld2[0][0][0], 5)

    
def persistence_image_distance(M1, M2):
    '''Returns the L2 distance between two Gaussian-smoothed diagrams.
    Parameters
    ----------
    M1: adjacency matrix.
    M2: adjacency matrix.
    '''
    SM1 = csr_matrix(M1)
    SM2 = csr_matrix(M2)
    diagrams1 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM1]);
    diagrams2 = FlagserPersistence(homology_dimensions=(0,1)).fit_transform([SM2]);
    
    diag1, diag2 = map(remove_inf, (diagrams1, diagrams2))
    
    d = PairwiseDistance(metric='persistence_image',order=None, metric_params={'p':2.0})
    ld1 = d.fit(X=[diag1])
    ld2 = ld1.transform([diag2])
    
    if len(ld2[0][0]) == 2:
        return round(ld2[0][0][1], 5)
    else:
        return round(ld2[0][0][0], 5)
