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
from gtda.graphs import GraphGeodesicDistance
from gtda.homology import VietorisRipsPersistence, SparseRipsPersistence, FlagserPersistence
from gtda.diagrams import PersistenceEntropy, PersistenceImage, BettiCurve


__all__ = [
    "remove_inf",
    "persistence_diagram",
    "wasserstein_distance",
    "bottleneck_distance",
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
    diagram = FlagserPersistence(homology_dimensions=(0,1,2,3,4)).fit_transform_plot([SM]);
    return diagram


def wasserstein_distance(pers_diag_1, pers_diag_2):
    '''Returns
    Parameters
    '''
    wd = pers.wasserstein(pers_diag_1[0], pers_diag_2[0])
    return round(wd, 4)


def bottleneck_distance(pers_diag_1, pers_diag_2):
    '''Returns
    Parameters
    '''
    diag1, diag2 = map(remove_inf, (pers_diag_1, pers_diag_2))
    diag11 = diag1[:, 0:2]
    diag22 = diag2[:, 0:2]
    bnd = gd.bottleneck_distance(diag11, diag22, 0.1)
    return bnd

