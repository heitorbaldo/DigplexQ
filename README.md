
<img src="digplexq-logo.png" alt="drawing" width="350"/>

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]()
[![PyPI version](https://img.shields.io/pypi/v/digplexq)](https://pypi.org/project/digplexq/)

------
DigplexQ is a Python package to perform computations with digraph-based complexes (directed flag complexes and path complexes). It is an "adjacency matrix-centered" package since it was designed so that the user can perform all computations just by entering an adjacency matrix as input.

This package implements several quantitative methods for the analysis of digraph-based complexes, mainly based on concepts from directed Q-Analysis. Only lower q-adjacencies were implemented so far.

* Free software: MIT license
* Documentation: [Documentation](https://github.com/heitorbaldo/DigplexQ/blob/main/documentation_digplexq.pdf)

Installation
--------

```bash
pip3 install digplexq
```

Examples
--------

```python
from digplexq.directed_q_analysis import *
from digplexq.digraph_based_complexes import *
from digplexq.structure_based_simplicial_measures import *
from digplexq.random_digraphs import *
from digplexq.utils import *

M = directed_erdos_renyi_GnM_model(20, 40, weight=False)
M = remove_double_edges(M) #remove double edges.

#Directed flag complex:
DFC = DirectedFlagComplex(M, "by_dimension_with_nodes")

#Maximal directed simplices:
maxsimp = MaximalSimplices(DFC)

#Lower q-Adjacency matrix:
fast_q_adjacency_matrix(M, q=1)

#in-q-degree centrality
in_q_degree_centrality(M, q=1, results="nodes")
```

More examples are available in the [Jupyter Notebook](https://github.com/heitorbaldo/DigplexQ/blob/main/Tutorial_DigplexQ.ipynb).


Issues
--------
For versions prior to v0.0.8, you may encounter the following errors:

```
AttributeError: module 'networkx' has no attribute 'from_numpy_matrix'
```
```
AttributeError: module 'networkx' has no attribute 'to_numpy_matrix'
```

To solve these issues, downgrade networkx or use ```nx.from_numpy_array()``` instead of ```nx.from_numpy_matrix()```, and ```nx.to_numpy_array()``` instead of ```nx.to_numpy_matrix()```.

Warnings
--------
This package only handles digraphs without double edges, so it is recommended to use the ```remove_double_edges()``` function before performing computations. Also, it is not optimized; therefore, it is only efficient for small digraphs.

Dependencies
--------

* [NumPy](https://github.com/numpy/numpy)
* [SciPy](https://scipy.org/)
* [NetworkX](https://github.com/networkx/networkx)
* [gtda](https://giotto-ai.github.io/gtda-docs/0.5.1/library.html)
* [persim](https://persim.scikit-tda.org/en/latest/)
* [hodgelaplacians](https://github.com/tsitsvero/hodgelaplacians)


