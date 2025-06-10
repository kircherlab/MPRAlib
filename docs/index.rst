MPRAlib documentation
=========================

.. image:: https://github.com/kircherlab/MPRAlib/actions/workflows/tests.yml/badge.svg?branch=master
   :target: https://github.com/kircherlab/MPRAlib/actions/workflows/tests.yml
   :alt: Tests

.. image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/kircherlab/MPRAlib/python-coverage-comment-action-data/endpoint.json
   :target: https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html
   :alt: Coverage badge

.. image:: https://badge.fury.io/py/mpralib.svg
   :target: https://badge.fury.io/py/mpralib
   :alt: PyPI version

.. image:: https://img.shields.io/conda/vn/bioconda/mpralib?label=bioconda
   :target: https://bioconda.github.io/recipes/mpralib/README.html
   :alt: Bioconda Version

MPRAlib is a Python library and CLI for processing MPRA (Massively Parallel Reporter Assay) data.

Installation
------------

PyPI
^^^^

.. code-block:: bash

    pip install mpralib

Conda
^^^^^

From the bioconda channel:

.. code-block:: bash

    conda install -c bioconda mpralib

Usage
-----

Command Line Interface
^^^^^^^^^^^^^^^^^^^^^^

Use the ``mpralib`` command to access various functionalities.

Validate a file
~~~~~~~~~~~~~~~

MPRAlib provides a CLI tool for validating MPRA data files against supported schemas.

.. code-block:: bash

    mpralib validate-file <schema> --input <input_file>

- ``<schema>``: One of ``reporter-sequence-design``, ``reporter-barcode-to-element-mapping``, ``reporter-experiment-barcode``, ``reporter-experiment``, ``reporter-element``, ``reporter-variant``, ``reporter-genomic-element``, ``reporter-genomic-variant``
- ``<input_file>``: Path to your data file (e.g., ``.tsv.gz``, ``.bed.gz``)

**Example:**

.. code-block:: bash

    mpralib validate-file reporter-sequence-design --input data/reporter_sequence_design.example.tsv.gz

Python API
^^^^^^^^^^

MPRAlib is primarily intended to be used as a library. Please see our notebook `mpralib.ipynb <https://github.com/kircherlab/MPRAlib/blob/master/examples/mpralib.ipynb>`_ for a more detailed example.

License
-------

MIT License

Links
-----

- `Documentation <https://github.com/kircherlab/MPRAlib>`_
- `Issues <https://github.com/kircherlab/MPRAlib/issues>`_



.. toctree::
   :maxdepth: 4
   :caption: Documentation
   :hidden:

   doc

.. toctree::
   :maxdepth: 4
   :caption: Tutorial
   :hidden:

   tutorial

.. toctree::
   :maxdepth: 4
   :caption: Project Info
   :hidden:

   project


.. toctree::
   :maxdepth: 4
   :caption: API
   :hidden:

   mpralib
