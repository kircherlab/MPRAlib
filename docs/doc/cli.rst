=======================
Command Line Interface
=======================

Use the ``mpralib`` command in your terminal to access various functionalities.

Validate a file
---------------

MPRAlib provides a CLI tool for validating MPRA data files against supported schemas.

.. code-block:: bash

    mpralib validate-file <schema> --input <input_file>

- ``<schema>``: One of ``reporter-sequence-design``, ``reporter-barcode-to-element-mapping``, ``reporter-experiment-barcode``, ``reporter-experiment``, ``reporter-element``, ``reporter-variant``, ``reporter-genomic-element``, ``reporter-genomic-variant``
- ``<input_file>``: Path to your data file (e.g., ``.tsv.gz``, ``.bed.gz``)

**Example:**

.. code-block:: bash

    mpralib validate-file reporter-sequence-design --input data/reporter_sequence_design.example.tsv.gz
