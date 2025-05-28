# MPRAlib

![Tests](https://github.com/kircherlab/MPRAlib/actions/workflows/tests.yml/badge.svg?branch=master)
[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/kircherlab/MPRAlib/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)
[![PyPI version](https://badge.fury.io/py/MPRAlib.svg)](https://badge.fury.io/py/MPRAlib)
[![Bioconda Version](https://img.shields.io/conda/vn/bioconda/mpralib?label=bioconda)](https://bioconda.github.io/recipes/mpralib/README.html)

MPRAlib is a Python library and CLI for processing MPRA (Massively Parallel Reporter Assay) data.

## Installation

### PyPI

```bash
pip install MPRAlib
```

### Conda

```bash
conda install -c bioconda mpralib
```

## Usage

### Command Line Interface

Use the `mpralib` command to access various functionalities.

#### Validate a file

MPRAlib provides a CLI tool for validating MPRA data files against supported schemas.

```bash
mpralib validate-file <schema> --input <input_file>
```

- `<schema>`: One of `reporter-sequence-design`, `reporter-barcode-to-element-mapping`, `reporter-experiment-barcode`, `reporter-experiment`, `reporter-element`, `reporter-variant`, `reporter-genomic-element`, `reporter-genomic-variant`
- `<input_file>`: Path to your data file (e.g., `.tsv.gz`, `.bed.gz`)

**Example:**

```bash
mpralib validate-file reporter-sequence-design --input data/reporter_sequence_design.example.tsv.gz
```

### Python API

In general MPRAlib is ment to be used as a library. Please have a look at our notebook [mpralib.ipynb](https://github.com/kircherlab/MPRAlib/blob/master/examples/mpralib.ipynb) for a more detailed example.

## License

MIT License

## Links

- [Documentation](https://github.com/mschubach/MPRAlib)
- [Issues](https://github.com/mschubach/MPRAlib/issues)

