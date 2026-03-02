# MPRAlib AI Coding Guidelines

## Project Overview

MPRAlib is a Python library and CLI tool for processing Massively Parallel Reporter Assay (MPRA) data. It provides data validation, normalization, filtering, and analysis workflows for barcode-resolved MPRA experiments, with support for multiple standardized file formats defined by the IGVF consortium.

## Architecture

### Core Data Model (AnnData-based)

The library uses **AnnData** objects as the fundamental data structure:
- **`MPRAData`** (abstract base class) - Foundation for all MPRA data types with shared normalization, filtering, and correlation methods
  - **`MPRABarcodeData`** (line [597](src/mpralib/mpradata.py#L597)) - Barcode-level data with multiple filtering strategies and ability to aggregate to oligo-level
  - **`MPRAOligoData`** (line [1136](src/mpralib/mpradata.py#L1136)) - Aggregated oligo-level data (created from barcode data)

**Key data structure**:
- `obs` (observations) = experiment replicates
- `var` (variables) = barcodes/oligos (indexed by barcode or oligo ID)
- `X` matrix = primary counts (RNA)
- `layers` dictionary = additional modalities (DNA, normalized variants)
- `varm` = oligo metadata (variant info, design info, statistics)

### Data Flow Pipeline

1. **File Loading** → `MPRABarcodeData.from_file()` reads TSV.GZ files into AnnData
2. **Design Metadata** → `add_sequence_design()` enriches oligos with variant/genomic info from TSV design file
3. **Normalization** → CPM normalization (counts per million) with configurable pseudocount (default=1) and scaling (default=1e6)
4. **Filtering** → `apply_barcode_filter()` removes outliers using 8+ strategies (MIN_BCS_PER_OLIGO, GLOBAL, OLIGO_SPECIFIC, LARGE_EXPRESSION, RANDOM, MIN_COUNT, MAX_COUNT)
5. **Aggregation** → `oligo_data` property converts barcode→oligo level by summing counts
6. **Output** → Export to standardized formats (activity, counts, variant maps, genomic annotations)

## Code Organization

- **[src/mpralib/cli.py](src/mpralib/cli.py)** - Click-based CLI with 4 command groups: `validate-file`, `functional`, `combine`, `plot`
- **[src/mpralib/mpradata.py](src/mpralib/mpradata.py)** - Core data classes and statistical methods (normalization, filtering, correlation)
- **[src/mpralib/utils/io.py](src/mpralib/utils/io.py)** - File I/O (compressed file detection, export functions, chromosome mapping)
- **[src/mpralib/utils/file_validation.py](src/mpralib/utils/file_validation.py)** - JSONSchema validation against 8 IGVF schemas
- **[src/mpralib/utils/plot.py](src/mpralib/utils/plot.py)** - Matplotlib/Seaborn visualization functions (uses BIH Charité color scheme)
- **[src/mpralib/schemas/](src/mpralib/schemas/)** - 8 JSON schema files defining standardized file formats

## Branding & Color Scheme

All plots use the **BIH Charité color palette**:
```python
BIH_COLORS = {
    "BLAU": "#003754",        # RGB 0, 55, 84 - Primary blue
    "WEISS": "#FFFFFF",       # RGB 255, 255, 255 - White
    "KORALL": "#EA5451",      # RGB 234, 84, 81 - Coral (alerts/outliers)
    "DUNKELROT": "#AF1821",   # RGB 175, 24, 33 - Dark red (critical values)
    "HELLROSA": "#FFB0AC",    # RGB 255, 176, 172 - Light pink
    "GOLD": "#9D7220",        # RGB 157, 114, 32 - Gold (highlights)
    "MINERAL": "#009AA9",     # RGB 0, 154, 169 - Mineral (secondary accent)
    "LAVENDEL": "#7876B6",    # RGB 120, 118, 182 - Lavender (tertiary accent)
}
```
Default plot palette: BLAU, MINERAL, KORALL, LAVENDEL, GOLD. Use specific colors for semantic meaning (e.g., DUNKELROT for filtering thresholds, KORALL for outliers).

## Critical Patterns

### Enum Pattern for Type Safety
Use enums extensively for fixed options (not string matching):
```python
from mpralib.mpradata import Modality, BarcodeFilter, CountSampling

# Load and filter
data = MPRABarcodeData.from_file("file.tsv.gz")
data.apply_barcode_filter(BarcodeFilter.GLOBAL, {"times_zscore": 3.0})

# Access modalities
activity = data.data.layers["activity"]  # or use Modality.ACTIVITY
```

### Metadata Management
Use AnnData's built-in metadata storage, NOT separate dictionaries:
```python
# Store/retrieve metadata (not data-altering)
self._add_metadata("SCALING", 1e6)
scaling = self._get_metadata("SCALING")

# Metadata is in `data.uns` dictionary
```

### Lazy Evaluation Pattern
Normalized layers are computed on-demand, not pre-computed:
```python
# Will compute only if not cached
dna_norm = data.normalized_dna_counts  # property
rna_norm = data.normalized_rna_counts  # property

# Check if computed
if "dna_normalized" in data.data.layers:
    # Already exists
```

### Property-based Data Access
Use `@property` decorators for computed/lazy fields. Setters may trigger recalculations:
```python
# These trigger layer validation and caching
@property
def normalized_dna_counts(self) -> np.ndarray:
    # Lazy computation + caching
```

## Testing

Test structure in [tests/](tests/):
- **[test_mpradata.py](tests/test_mpradata.py)** - Unit tests for data classes (1185 lines, heavy use of fixtures)
- **[test_cli*.py](tests/)** - CLI integration tests using Click's test runner
- **[test_utils_*.py](tests/)** - I/O and plotting utilities

Run tests:
```bash
pytest --tb=short  # Run all tests with short traceback
pytest tests/test_mpradata.py::test_apply_count_sampling_rna  # Single test
```

Test fixture pattern:
```python
@pytest.fixture
def mpra_data():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    return MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS, var=VAR, layers=layers))
```

## File Format Standards

8 standardized schemas (JSON-based, enforced by validation):
1. `reporter_sequence_design` - Maps oligos to variants/genomic positions
2. `reporter_experiment_barcode` - Raw barcode counts (DNA/RNA)
3. `reporter_barcode_to_element_mapping` - Barcode-oligo mappings
4. `reporter_experiment` - Aggregate oligo data
5. `reporter_element` / `reporter_variant` - Analyzed elements/variants
6. `reporter_genomic_element` / `reporter_genomic_variant` - BED-like genomic annotations

All formats: **TSV.GZ** (tab-separated, gzip-compressed), validated via CLI `mpralib validate-file <schema> --input <file>`

## Common Operations

### Loading and Basic Analysis
```python
from mpralib.mpradata import MPRABarcodeData, BarcodeFilter, Modality

# Load
data = MPRABarcodeData.from_file("counts.tsv.gz")

# Filter
data.apply_barcode_filter(BarcodeFilter.MIN_COUNT, {"min_count": 5})

# Aggregate
oligo_data = data.oligo_data

# Correlations
corr = oligo_data.correlation("pearson", Modality.ACTIVITY)
```

### Adding Metadata
```python
from mpralib.utils.io import read_sequence_design_file

design = read_sequence_design_file("design.tsv")
data.add_sequence_design(design, "design.tsv")
print(data.variant_map)  # Variant-to-oligo mapping

# Add custom labels
labels = pd.Series({"oligo1": "control", "oligo2": "variant", "oligo3": "control"})
data.add_labels(labels, label_column="category")
```

### Exporting Results
```python
from mpralib.utils.io import export_activity_file, export_barcode_file

export_activity_file(oligo_data, "activity.tsv")
export_barcode_file(data, "barcodes.tsv")
```

## Exception Handling

Use custom exceptions defined in [src/mpralib/exception.py](src/mpralib/exception.py):
- `MPRAlibException` - Base exception
- `IOException` - File/IO errors
- `SequenceDesignException` - Design file column format issues

Example:
```python
from mpralib.exception import MPRAlibException

try:
    data.add_sequence_design(design, path)
except SequenceDesignException as e:
    # Handle column format issue
```

## Development Priorities

1. **Data integrity** - All layer operations must maintain AnnData structure; test with fixtures
2. **Backward compatibility** - Schema versions are immutable; add new schemas if formats change
3. **Computational efficiency** - Use NumPy/masked arrays for large datasets; avoid redundant normalization
4. **Validation first** - Validate input files against JSON schemas before processing
5. **Clear CLI contracts** - Commands are user-facing; preserve option names and semantics
