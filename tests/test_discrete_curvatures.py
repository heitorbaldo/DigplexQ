'''
Unit Test
'''

import math
import unittest
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from digplexq.discrete_curvatures import *
from digplexq.digraph_based_complexes import *
from digplexq.directed_q_analysis import *
from digplexq.utils import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

M = np.array([[0., 1., 1., 0., 1., 0.],
       [0., 0., 0., 0., 0., 0.],
       [0., 1., 0., 0., 0., 1.],
       [1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 1.],
       [1., 0., 0., 0., 0., 0.]])

DFC_dim_none = DirectedFlagComplex(M, "by_dimension_without_nodes")

frc_21 = round(4*((2/4+2/4) - (2/math.sqrt(4*3*2))), 5)
infrc_2 = -7.14162
outfrc_1 = 0

frc_30 = round(3*((1/3+3/3) - (3/math.sqrt(3*3*2)) - (3/math.sqrt(3*3*2)) - (3/math.sqrt(3*3*1))), 5)
infrc_3 = 0 
outfrc_3 = -3.24264


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

class TestFormanRicciCurvature(unittest.TestCase):
    
    def test_q_forman_ricci_curvature(self):
        self.assertEqual(q_forman_ricci_curvature(M, DFC_dim_none, [2,1]), frc_21)
        self.assertEqual(q_forman_ricci_curvature(M, DFC_dim_none, [3,0]), frc_30)
        
    def test_in_q_forman_ricci_curvature(self):
        self.assertEqual(in_q_forman_ricci_curvature(M, DFC_dim_none, 2), infrc_2)
        self.assertEqual(in_q_forman_ricci_curvature(M, DFC_dim_none, 3), infrc_3)
        
    def test_out_q_forman_ricci_curvature(self):
        self.assertEqual(out_q_forman_ricci_curvature(M, DFC_dim_none, 1), outfrc_1)
        self.assertEqual(out_q_forman_ricci_curvature(M, DFC_dim_none, 3), outfrc_3)

        

if __name__ == '__main__':
    unittest.main()