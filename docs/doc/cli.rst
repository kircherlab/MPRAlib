=======================
Command Line Interface
=======================

Use the ``mpralib`` command in your terminal to access various functionalities.

You can list all commands with:

.. code-block:: bash

    mpralib --help


Validate a file
----------------

MPRAlib provides a CLI tool for validating MPRA data files against supported schemas.

.. code-block:: bash

    mpralib validate-file <schema> --input <input_file>

- ``<schema>``: One of ``reporter-sequence-design``, ``reporter-barcode-to-element-mapping``, ``reporter-experiment-barcode``, ``reporter-experiment``, ``reporter-element``, ``reporter-variant``, ``reporter-genomic-element``, ``reporter-genomic-variant``
- ``<input_file>``: Path to your data file (e.g., ``.tsv.gz``, ``.bed.gz``)

**Example:**

.. code-block:: bash

    mpralib validate-file reporter-sequence-design --input data/reporter_sequence_design.example.tsv.gz


Functional analysis
-------------------

Run core analysis utilities on barcode-level input files.

.. code-block:: bash

    mpralib functional <command> <options>

- ``<command>``: One of ``activities``, ``compute-correlation``, ``filter``

``activities``
^^^^^^^^^^^^^^

Generate activity or barcode output files.

Required options:

- ``--input``
- ``--output``

Optional options:

- ``--bc-threshold`` (default: ``1``)
- ``--element-level/--barcode-level`` (default: ``--element-level``)

**Example:**

.. code-block:: bash

    mpralib functional activities --input data/reporter_experiment_barcode.example.tsv.gz --bc-threshold 10 --element-level --output data/reporter_experiment.activity.tsv.gz

``compute-correlation``
^^^^^^^^^^^^^^^^^^^^^^^

Compute pairwise correlations on oligo-level data.

Required options:

- ``--input``

Optional options:

- ``--bc-threshold`` (default: ``1``)
- ``--correlation-method``: ``pearson``, ``spearman``, ``all`` (default: ``all``)
- ``--correlation-on``: ``dna_normalized``, ``rna_normalized``, ``activity``, ``all`` (default: ``all``)

**Example:**

.. code-block:: bash

    mpralib functional compute-correlation --input data/reporter_experiment_barcode.example.tsv.gz --bc-threshold 10 --correlation-method pearson --correlation-on activity

``filter``
^^^^^^^^^^

Filter barcodes using supported filtering methods and optional method-specific parameters.

Required options:

- ``--input``
- ``--method`` (choices are derived from ``BarcodeFilter``)

Optional options:

- ``--method-values`` (JSON string or Python dict)
- ``--bc-threshold`` (default: ``1``)
- ``--output-activity``
- ``--output-barcode``

**Example:**

.. code-block:: bash

    mpralib functional filter --input data/reporter_experiment_barcode.example.tsv.gz --method max_count --method-values '{"rna_max_count": 500, "dna_max_count": 300}' --output-barcode data/reporter_experiment_barcode.filtered.tsv.gz


Plotting
--------

Generate plots of your data.

.. code-block:: bash

    mpralib plot <command> <options>

- ``<command>``: One of ``barcodes-per-oligo``, ``correlation``, ``dna-vs-rna``, ``outlier``

**Example:**

.. code-block:: bash

    mpralib plot correlation --input data/reporter_experiment_barcode.example.tsv.gz --oligos --bc-threshold 10 --modality activity --output data/test.png

Plot command options:

- ``correlation``: ``--input`` (required), ``--output`` (required), ``--oligos/--barcodes``, ``--bc-threshold``, ``--modality``, ``--replicate`` (can be provided multiple times)
- ``dna-vs-rna``: ``--input`` (required), ``--output`` (required), ``--oligos/--barcodes``, ``--bc-threshold``, ``--replicate`` (can be provided multiple times)
- ``barcodes-per-oligo``: ``--input`` (required), ``--output`` (required), ``--replicate`` (can be provided multiple times)
- ``outlier``: ``--input`` (required), ``--output`` (required), ``--bc-threshold``


Combine counts with other outputs
---------------------------------

Combines counts data with MPRA sequence design and quantification from other tools like BCalm to create several other output data, like bed file tracks. Can also create a variant map from the MPRA design file as an input for mpralm and BCalm.

.. code-block:: bash

    mpralib combine <command> <options>

- ``<command>``: One of ``get-counts``, ``get-reporter-elements``, ``get-reporter-genomic-elements``, ``get-reporter-genomic-variants``, ``get-reporter-variants``, ``get-variant-counts``, ``get-variant-map``

**Example:**

.. code-block:: bash

    mpralib combine get-variant-map --input data/reporter_experiment_barcode.example.tsv.gz --sequence-design data/mpra_sequence_design.example.tsv.gz --output data/variant_map_of_oligo.tsv.gz

Combine command highlights:

- ``get-variant-map``: required ``--sequence-design``, ``--output``; optional ``--input``
- ``get-counts``: required ``--input``, ``--sequence-design``, ``--output``; optional ``--bc-threshold``, ``--scaling-factor``, ``--pseudo-count``, ``--oligos/--barcodes``, ``--normalized-counts/--counts``, ``--elements-only/--all-oligos``
- ``get-variant-counts``: required ``--input``, ``--sequence-design``, ``--output``; optional ``--bc-threshold``, ``--normalized-counts/--counts``, ``--scaling-factor``, ``--pseudo-count``, ``--oligos/--barcodes``
- ``get-reporter-elements``: required ``--input``, ``--sequence-design``, ``--statistics``, ``--output-reporter-elements``; optional ``--bc-threshold``
- ``get-reporter-variants``: required ``--input``, ``--sequence-design``, ``--statistics``, ``--output-reporter-variants``; optional ``--bc-threshold``
- ``get-reporter-genomic-elements``: required ``--input``, ``--sequence-design``, ``--statistics``, ``--reference``, ``--output-reporter-genomic-elements``; optional ``--bc-threshold``
- ``get-reporter-genomic-variants``: required ``--input``, ``--sequence-design``, ``--statistics``, ``--reference``, ``--output-reporter-genomic-variants``; optional ``--bc-threshold``

For detailed usage of a subcommand, use:

.. code-block:: bash

    mpralib <group> <command> --help
