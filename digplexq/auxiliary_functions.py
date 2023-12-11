'''
Auxiliary functions
'''

import numpy as np
import networkx as nx
import warnings
warnings.filterwarnings("ignore")


__all__=[
    "avg_std_q_measure",
]


def avg_std_q_measure(M, start, end, func):
    '''Returns the mean and the standard deviation for each level q.
    Parameters
    ---------
    '''
    X = []
    H0 = []
    H1 = []
    H2 = []
    H3 = []
    for k in range(start, end):
        DFC_dim_nodes = DirectedFlagComplex(M.transpose(), "by_dimension_with_nodes")
        DFC_dim_none = DirectedFlagComplex(M.transpose(), "by_dimension_without_nodes")
        R = []
        for q in range(0, 4):
            Hq = q_adjacency_matrix(DFC_dim_none, DFC_dim_nodes, q)
            result = func(Hq)
            R.append(result)
        X.append(R)

    for k in range(len(X)):
        H0.append(X[k][0])
        H1.append(X[k][1])
        H2.append(X[k][2])
        H3.append(X[k][3])

    Y0 = np.array(H0)
    Y1 = np.array(H1)
    Y2 = np.array(H2)
    Y3 = np.array(H3)

    mean0 = round(Y0.mean(), 4)
    std0 = round(Y0.std(), 4)
    mean1 = round(Y1.mean(), 4)
    std1 = round(Y1.std(), 4)
    mean2 = round(Y2.mean(), 4)
    std2 = round(Y2.std(), 4)
    mean3 = round(Y3.mean(), 4)
    std3 = round(Y3.std(), 4)

    Final_Res = [[mean0, std0], [mean1, std1], [mean2, std2], [mean3, std3]]

    return Final_Res
