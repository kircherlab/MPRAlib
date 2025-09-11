MPRAlib documentation
=========================

.. image:: https://app.readthedocs.org/projects/mpralib/badge/?version=latest
   :target: https://mpralib.readthedocs.io/latest/?badge=latest
   :alt: Documentation Status

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


Functional analysis
~~~~~~~~~~~~~~~~~~~

Filter barcodes using multiple filters, like setting min/max counts or detect barcode outliers

.. code-block:: bash

    mpralib functional <schema> <inputs>

- ``<schema>``: One of ``activities``, ``compute-correlation``, ``filter``
- ``<inputs>``: Please use ``--help`` for more details on the schema.

**Example:**

.. code-block:: bash

    mpralib functional filter --input data/reporter_experiment_barcode.example.tsv.gz --method max_count --method-values '{"rna_max_count": 500, "dna_max_count": 300}' --output-barcode data/reporter_experiment_barcode.filtered.tsv.gz


Plotting
~~~~~~~~

Generate plots of your data.

.. code-block:: bash

    mpralib plot <schema> --input <input_file> --bc-threshold <bc_threshold> --output <output_file> <other_inputs>

- ``<schema>``: One of ``barcodes-per-oligo``, ``correlation``, ``dna-vs-rna``, ``outlier``
- ``<input_file>``: MPRA experiment or experiment barcode.
- ``<bc_threshold>``: Minimum number of barcodes per oligo to include in the plot.
- ``<output_file>``: Path to save the plot (e.g., ``.png``, ``.pdf``)
- ``<other_inputs>``: Please use ``--help`` for more details on the schema.

**Example:**

.. code-block:: bash

    mpralib plot correlation --input data/reporter_experiment_barcode.example.tsv.gz --oligos --bc-threshold 10 --modality activity --output data/test.png


Combine counts with other outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combines counts data with MPRA sequence design and quantification from other tools like BCalm to create several other output data, like bed file tracks. Can also create a variant map from the MPRA design file as an input for mpralm and BCalm.

.. code-block:: bash

    mpralib combine <schema> <inputs>

- ``<schema>``: One of ``get-counts``, ``get-reporter-elements``, ``get-reporter-genomic-elements``, ``get-reporter-genomic-variants``, ``get-reporter-variants``, ``get-variant-counts``, ``get-variant-map``
- ``<inputs>``: Please use ``--help`` for more details on the schema.

**Example:**

.. code-block:: bash

    mpralib combine get-variant-map --input data/reporter_experiment_barcode.example.tsv.gz --sequence-design data/mpra_sequence_design.example.tsv.gz --output data/variant_map_of_oligo.tsv.gz


Python API
^^^^^^^^^^

MPRAlib is primarily intended to be used as a library. Please see our notebook `mpralib.ipynb <https://github.com/kircherlab/MPRAlib/blob/master/examples/mpralib.ipynb>`_ for a more detailed example.

License
-------

MIT License

.. toctree::
   :maxdepth: 4
   :hidden:

   doc/overview

.. toctree::
   :maxdepth: 4
   :caption: Documentation
   :hidden:

   doc/getting-started
   doc/cli

.. toctree::
   :maxdepth: 4
   :caption: Tutorial
   :hidden:

   tutorial/tutorial

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
