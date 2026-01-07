# MPRAlib

[![Documentation Status](https://readthedocs.org/projects/mpralib/badge/?version=latest)](https://mpralib.readthedocs.io/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/913714670.svg)](https://doi.org/10.5281/zenodo.18173084)
[![GitHub License](https://img.shields.io/github/license/kircherlab/MPRAlib)](https://github.com/kircherlab/MPRAlib/blob/master/LICENSE)
[![GitHub Release](https://img.shields.io/github/v/release/kircherlab/MPRAlib)](https://github.com/kircherlab/MPRAlib/releases/latest)
[![PyPI version](https://badge.fury.io/py/mpralib.svg)](https://badge.fury.io/py/mpralib)
[![Bioconda Version](https://img.shields.io/conda/vn/bioconda/mpralib?label=bioconda)](https://bioconda.github.io/recipes/mpralib/README.html)
[![Tests](https://github.com/kircherlab/MPRAlib/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/kircherlab/MPRAsnakeflow/actions/workflows/main.yml)
[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/kircherlab/MPRAlib/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/kircherlab/MPRAlib/blob/python-coverage-comment-action-data/htmlcov/index.html)
[![GitHub Issues](https://img.shields.io/github/issues/kircherlab/MPRAlib)](https://github.com/kircherlab/MPRAlib/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/kircherlab/MPRAlib)](https://github.com/kircherlab/MPRAlib/pulls)

MPRAlib is a Python library and CLI for processing MPRA (Massively Parallel Reporter Assay) data.

## Citation

If you use MPRAlib in your work please cite out recent preprint:

**Uniform processing and analysis of IGVF massively parallel reporter assay data with MPRAsnakeflow**
Jonathan D. Rosen, Arjun Devadas Vasanthakumari, Kilian Salomon, Nikola de Lange, Pyaree Mohan Dash, Pia Keukeleire, Ali Hassan, Alejandro Barrera, Martin Kircher, Michael I. Love, Max Schubach
*bioRxiv* (2025). [2025.09.25.678548](https://doi.org/10.1101/2025.09.25.678548)

## Installation

### PyPI

```bash
pip install mpralib
```

### Conda

From the bioconda channel

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

- [Documentation](https://mpralib.readthedocs.io)
- [Issues](https://github.com/kircherlab/MPRAlib/issues)

