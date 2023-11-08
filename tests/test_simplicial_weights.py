'''
unittest / numpy.testing
'''

import unittest
from unittest import TestCase
import numpy as np
import warnings
warnings.filterwarnings("ignore")

from digplexq.simplicial_weights import *


A = np.matrix([[0, 8, 4], [0, 0, 0], [0, 7, 0]])
B = np.matrix([[0, 0, 9], [4, 0, 0], [1, 13, 0]])
C = np.matrix([[0, 1, 9], [4, 0, 0], [1, 13, 0]])

weights_A = [12, 15, 7]
weights_B = [9, 13, 13]
weights_C = [10, 14, 14]

class TestWeights(TestCase):
    
    @staticmethod
    def test_node_weights_max():
        np.testing.assert_array_equal(node_weights_max(A), weights_A)
        np.testing.assert_array_equal(node_weights_max(B), weights_B)
        
    @staticmethod
    def test_simplex_weight_max_in_out():
        np.testing.assert_array_equal(node_weights_max_in_out(C), weights_C)
        
    @staticmethod
    def test_simplex_weight_max():
        np.testing.assert_array_equal(simplex_weight(A, [0,1,2], weight_func='max'), 1260)
        np.testing.assert_array_equal(simplex_weight(B, [0,1,2], weight_func='max'), 1521)
        

if __name__ == '__main__':
    unittest.main()