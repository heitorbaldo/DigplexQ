.. image:: digplexq-logo.png
   :width: 350px

.. image:: https://colab.research.google.com/assets/colab-badge.svg
        :target: https://colab.research.google.com/github/

.. image:: https://img.shields.io/pypi/v/digplexq.svg
        :target: https://pypi.python.org/pypi/digplexq

.. image:: https://readthedocs.org/projects/digplexq/badge/?version=latest
        :target: https://digplexq.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/heitorbaldo/digplexq/shield.svg
     :target: https://pyup.io/repos/github/heitorbaldo/digplexq/
     :alt: Updates


DigplexQ is a Python package to perform computations with digraph-based complexes (e.g., directed flag complexes and path complexes)

* Free software: MIT license
* Documentation: https://digplexq.readthedocs.io.

Installation
--------

```python
pip install digplexq
```

Examples
--------

```python
import digplexq

M = directed_erdos_renyi_GnM_model(20, 40, weight=False)
M = remove_double_edges(M) #remove double edges.

DFC_dim_none = DirectedFlagComplex(M, "by_dimension_without_nodes")
```

Dependencies
--------

* NumPy
* SciPy
* Networkx
* gtda
* persim
* hodgelaplacians

