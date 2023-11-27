'''
Unit Test
'''

import math
import unittest
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from digplexq.structure_based_simplicial_measures import *
from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *
from digplexq.utils import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

M = np.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 1., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])

indd_M = -round(((6/10)*math.log2(6/10) + (3/10)*math.log2(3/10) + (1/10)*math.log2(1/10)), 5)
outdd_M = -round(((5/10)*math.log2(5/10) + (5/10)*math.log2(5/10)), 5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class TestStructureBasedMeasures(unittest.TestCase):
    
    def test_in_q_degree_distribution_entropy(self):
        self.assertEqual(in_q_degree_distribution_entropy(M), indd_M)
        
    def test_out_q_degree_distribution_entropy(self):
        self.assertEqual(out_q_degree_distribution_entropy(M), outdd_M)

        

if __name__ == '__main__':
    unittest.main()