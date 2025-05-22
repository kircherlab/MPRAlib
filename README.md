# MPRAlib

[![Tests](https://github.com/mschubach/MPRAlib/actions/workflows/test.yml/badge.svg)](https://github.com/mschubach/MPRAlib/actions/workflows/test.yml)
[![Coverage Status](https://coveralls.io/repos/github/mschubach/MPRAlib/badge.svg?branch=main)](https://coveralls.io/github/mschubach/MPRAlib?branch=main)
[![PyPI version](https://badge.fury.io/py/MPRAlib.svg)](https://badge.fury.io/py/MPRAlib)
[![Conda version](https://anaconda.org/conda-forge/mpralib/badges/version.svg)](https://anaconda.org/conda-forge/mpralib)

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

You can also use MPRAlib as a Python library:

```python
from mpralib.utils.file_validation import validate_file, ValidationSchema

result = validate_file(
    ValidationSchema.REPORTER_SEQUENCE_DESIGN,
    "data/reporter_sequence_design.example.tsv.gz"
)
print(result)
```

## License

MIT License

## Links

- [Documentation](https://github.com/mschubach/MPRAlib)
- [Issues](https://github.com/mschubach/MPRAlib/issues)

