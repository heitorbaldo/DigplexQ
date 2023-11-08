'''
Discrete Morse Theory on Digraphs
'''

import numpy as np
from pyqpath.path_weights import *


__all__ = [
       "DiscreteMorse",
       "CriticalPaths"
]


class DiscreteMorse:
    
    def __init__(self):
        pass
    
    def discMorsefunc_node(self, M, path, v):
        lw = LocalWeights()
        f = lw.local_node_weights(M, path)[v]*local_node_weights(M, path)[v+1]
        return f


    #F(e) = sum_{i} f(i):
    def discMorsefunc_path(self, M, k, path):
        F = 0
        if path[0] != path[len(path)-1]:
            for i in range(k-2):
                F += self.discMorsefunc_node(M, path[0:k], i)
        else:
            for i in range(k-2):
                 F += self.discMorsefunc_node(M, path[0:k+1], i)
        return F

    

class CriticalPaths:
    
    def __init__(self):
        pass
    
    def critical_path(self, M, k, path):
        dm = DiscreteMorse()
        if dm.discMorsefunc_path(M, k, path) == dm.discMorsefunc_path(M, k-1, path) or \
            dm.discMorsefunc_path(M, k, path[0:k]) == dm.discMorsefunc_path(M, k, path[0:k+1]):
            print("%s is not critical." % path)
        else:
            print("%s is critical." % path)

          
      
          