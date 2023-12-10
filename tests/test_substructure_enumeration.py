'''
Substructure enumeration
'''

import unittest
import numpy as np

from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *
from digplexq.substructures_enumeration import *


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

A1 = np.matrix([[0, 0, 0], [1, 0, 0], [1, 1, 0]]) #contains a 2-clique
A2 = np.matrix([[0, 0, 1], [1, 0, 1], [0, 0, 0]]) #contains a 2-clique
A3 = np.matrix([[0, 0, 0], [1, 0, 1], [1, 0, 0]]) #contains a 2-clique
A4 = np.matrix([[0, 1, 1], [0, 0, 1], [0, 0, 0]]) #contains a 2-clique
A5 = np.matrix([[0, 1, 0], [0, 0, 0], [1, 1, 0]]) #contains a 2-clique
A6 = np.matrix([[0, 0, 0], [1, 0, 1], [1, 0, 0]]) #contains a 2-clique
A7 = np.matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) #does not contain a 2-clique

B1 = np.matrix([[0, 1, 1, 1], [0, 0, 1, 1], [1, 1, 0, 0], [0, 0, 1, 0]]) #contains a 3-clique
B2 = np.matrix([[0, 1, 0, 1], [0, 0, 1, 1], [1, 1, 0, 0], [0, 0, 1, 0]]) #does not contain a 3-clique

C1 = np.matrix([[0, 1, 1], [1, 0, 0], [1, 1, 0]]) #contains 2 double edges
C2 = np.matrix([[0, 1, 1], [1, 0, 0], [0, 1, 0]]) #contains 1 double edges
C3 = np.matrix([[0, 0, 0], [1, 0, 0], [1, 1, 0]]) #does not contain a double edge

D1 = np.matrix([[0, 1, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 0]]) #contains a directed square
D2 = np.matrix([[0, 0, 0, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0, 0, 0, 0]]) #does not contain a directed square


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class TestCountDirectedCliques(unittest.TestCase):
    
    def test_count_directed_n_cliques(self):
        self.assertEqual(count_directed_n_cliques(A1), 1)
        self.assertEqual(count_directed_n_cliques(A2), 1)
        self.assertEqual(count_directed_n_cliques(A3), 1)
        self.assertEqual(count_directed_n_cliques(A4), 1)
        self.assertEqual(count_directed_n_cliques(A5), 1)
        self.assertEqual(count_directed_n_cliques(A6), 1)
        self.assertEqual(count_directed_n_cliques(A7), 0)


class TestCoundDoubleEdges(unittest.TestCase):
   
    def test_count_double_edges(self):
        self.assertEqual(len(double_edges(C1))/2, 2)
        self.assertEqual(len(double_edges(C2))/2, 1)
        self.assertEqual(len(double_edges(C3))/2, 0)
        
        
class TestCountSquares(unittest.TestCase):
   
    def test_count_squares(self):
        self.assertEqual(len(squares(D1)), 1)
        self.assertEqual(len(squares(D2)), 0)


if __name__ == '__main__':
    unittest.main()