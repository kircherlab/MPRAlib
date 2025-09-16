import copy
import gzip
import os
import tempfile
import pandas as pd
import numpy as np
from numpy.typing import NDArray
from typing import Optional
import pytest
import anndata as ad
from mpralib.utils.io import export_barcode_file, read_sequence_design_file, export_counts_file, chromosome_map
from mpralib.exception import SequenceDesignException, MPRAlibException
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, MPRAData


class DummyMPRAData(MPRAData):
    def __init__(
        self,
        replicates: list[str],
        oligos: list[str],
        barcodes: list[str],
        dna_counts: NDArray[np.int32],
        rna_counts: NDArray[np.int32],
        barcode_threshold: int,
        barcode_counts: Optional[NDArray[np.int32]] = None,
    ):
        layers = {"rna": rna_counts, "dna": dna_counts}
        obs = pd.DataFrame(index=replicates)
        var = pd.DataFrame({"oligo": oligos}, index=barcodes)
        super().__init__(ad.AnnData(X=rna_counts, obs=obs, var=var, layers=layers), barcode_threshold)
        if barcode_counts is not None:
            self.barcode_counts = barcode_counts


class DummyMPRABarcodeData(DummyMPRAData, MPRABarcodeData):
    pass


class DummyMPRAOligoData(DummyMPRAData, MPRAOligoData):
    pass


@pytest.fixture
def barcode_data():
    replicates = ["rep1", "rep2"]
    barcodes = ["bc1", "bc2", "bc3"]
    oligos = ["ol1", "ol2", "ol3"]
    dna_counts = np.array([[10, 0, 5], [20, 1, 0]])
    rna_counts = np.array([[5, 0, 2], [10, 1, 0]])
    barcode_threshold = 1
    return DummyMPRABarcodeData(replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold)


@pytest.fixture
def oligo_data():
    replicates = ["rep1", "rep2"]
    barcodes = ["bc1", "bc2", "bc3"]
    oligos = ["ol1", "ol2", "ol3"]
    dna_counts = np.array([[10, 0, 5], [20, 1, 0]], dtype=np.int32)
    rna_counts = np.array([[5, 0, 2], [10, 1, 0]], dtype=np.int32)
    barcode_counts = np.array([[1, 0, 1], [1, 0, 1]], dtype=np.int32)
    barcode_threshold = 1
    return DummyMPRAOligoData(replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold, barcode_counts)


def test_export_counts_file_barcode(tmp_path, barcode_data):
    out_path = tmp_path / "counts.tsv"
    export_counts_file(barcode_data, str(out_path), normalized=False)
    df = pd.read_csv(out_path, sep="\t", na_values=["X"])
    # Only rows with observed=True and barcode_counts >= threshold should remain
    assert "oligo_name" in df.columns
    assert set(df["oligo_name"]) == {"ol1", "ol2", "ol3"}
    assert "dna_count_rep1" in df.columns
    assert "rna_count_rep2" in df.columns
    assert "rna_count_rep2" in df.columns
    assert df.loc[df["oligo_name"] == "ol3", "dna_count_rep1"].values[0] == 5
    # Check that zeros are written as empty
    assert np.isnan(df.loc[df["oligo_name"] == "ol3", "dna_count_rep2"].values[0])


def test_export_counts_file_oligo(tmp_path, oligo_data):
    out_path = tmp_path / "counts.tsv"
    export_counts_file(oligo_data, str(out_path), normalized=False)
    df = pd.read_csv(out_path, sep="\t")
    # Only rows with observed=True and barcode_counts >= threshold should remain
    assert "ID" in df.columns
    assert set(df["ID"]) == {"ol1", "ol3"}
    assert "dna_count_rep1" in df.columns
    assert "rna_count_rep2" in df.columns


def test_export_counts_file_normalized(tmp_path, barcode_data):
    out_path = tmp_path / "counts_norm.tsv"
    export_counts_file(barcode_data, str(out_path), normalized=True)
    df = pd.read_csv(out_path, sep="\t")
    # Check float formatting
    assert all(df.filter(like="dna_count").apply(lambda x: isinstance(x, float) or pd.isna(x)))


def test_export_counts_file_with_filter(tmp_path, barcode_data):
    out_path = tmp_path / "counts_filt.tsv"
    # Mask out the first entry
    filter_mask = np.asarray([[True, False, False], [True, False, False]])
    export_counts_file(barcode_data, str(out_path), normalized=False, filter=filter_mask)
    df = pd.read_csv(out_path, sep="\t")
    # Only ol3 should remain (ol2 fails barcode_threshold, ol1 is masked)
    assert set(df["oligo_name"]) == {"ol2", "ol3"}


def test_export_counts_file_invalid_type(tmp_path):
    class Dummy:
        pass

    dummy = Dummy()
    with pytest.raises(MPRAlibException):
        export_counts_file(dummy, str(tmp_path / "fail.tsv"))  # type: ignore


VALID_TSV = """
name\tsequence\tcategory\tclass\tsource\tref\tchr\tstart\tend\tstrand\tvariant_class\tvariant_pos\tSPDI\tallele\tinfo
oligo1\tATGCATGC\tvariant\ttest\tsrc1\tref1\tchr1\t100\t200\t+\t["SNV"]\t[123]\t["NC_00007.14:117548628:TTTTTTT:TTTTTTTTT"]\t["ref"]\tinfo1
oligo2\tGGCATGCA\telement\tvariant positive control\tsrc2\tref2\tchr2\t300\t400\t-\tNA\tNA\tNA\tNA\tinfo2
oligo3\tATGCATGC\tvariant\ttest\tsrc1\tref1\tchr1\t100\t200\t+\t["SNV"]\t[123]\t["NC_00007.14:117548628:TTTTTTT:TTTTTTTTT"]\t["alt"]\tinfo1
"""

INVALID_SEQ_TSV = """
name\tsequence\tcategory\tclass\tsource\tref\tchr\tstart\tend\tstrand\tvariant_class\tvariant_pos\tSPDI\tallele\tinfo
oligo1\tATGCXTGC\tvariant\ttest\tsrc1\tref1\tchr1\t100\t200\t+\t["SNV"]\t[123]\t["NC_00007.14:117548628:TTTTTTT:TTTTTTTTT"]\t["alt"]\tinfo1
"""

INVALID_CATEGORY_TSV = """
name\tsequence\tcategory\tclass\tsource\tref\tchr\tstart\tend\tstrand\tvariant_class\tvariant_pos\tSPDI\tallele\tinfo
oligo1\tATGCATGC\tbadcat\ttest\tsrc1\tref1\tchr1\t100\t200\t+\tNA\tNA\tNA\tNA\tinfo1
"""

INVALID_CLASS_TSV = """
name\tsequence\tcategory\tclass\tsource\tref\tchr\tstart\tend\tstrand\tvariant_class\tvariant_pos\tSPDI\tallele\tinfo
oligo1\tATGCATGC\tvariant\tbadclass\tsrc1\tref1\tchr1\t100\t200\t+\t["SNV"]\t[123]\t["NC_00007.14:117548628:TTTTTTT:TTTTTTTTT"]\t["alt"]\tinfo1
"""


def write_temp_tsv(content):
    fd, path = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w") as f:
        f.write(content)
    return path


def test_read_sequence_design_file_valid():
    path = write_temp_tsv(VALID_TSV)
    df = read_sequence_design_file(path)
    assert isinstance(df, pd.DataFrame)
    assert "sequence" in df.columns
    assert df.shape[0] == 3
    assert set(df.index) == {"oligo1", "oligo2", "oligo3"}
    os.remove(path)


def test_read_sequence_design_file_invalid_sequence():
    path = write_temp_tsv(INVALID_SEQ_TSV)
    with pytest.raises(SequenceDesignException) as excinfo:
        read_sequence_design_file(path)
    assert "sequence" in str(excinfo.value)
    os.remove(path)


def test_read_sequence_design_file_invalid_category():
    path = write_temp_tsv(INVALID_CATEGORY_TSV)
    with pytest.raises(SequenceDesignException) as excinfo:
        read_sequence_design_file(path)
    assert "category" in str(excinfo.value)
    os.remove(path)


def test_read_sequence_design_file_invalid_class():
    path = write_temp_tsv(INVALID_CLASS_TSV)
    with pytest.raises(SequenceDesignException) as excinfo:
        read_sequence_design_file(path)
    assert "class" in str(excinfo.value)
    os.remove(path)


def test_chromosome_map_columns():
    df = chromosome_map()
    expected_columns = ["ucsc", "assembly", "genbank", "refseq", "release"]
    assert list(df.columns) == expected_columns


def test_chromosome_map_release_values():
    df = chromosome_map()
    assert set(df["release"].unique()) == {"GRCh37", "GRCh38"}


def test_chromosome_map_non_empty():
    df = chromosome_map()
    assert not df.empty
    assert df.shape[0] > 0


def test_chromosome_map_no_duplicate_rows():
    df = chromosome_map()
    assert df.duplicated().sum() == 0


OBS = pd.DataFrame(index=["rep1", "rep2", "rep3"])
VAR = pd.DataFrame(
    {"oligo": ["oligo1", "oligo1", "oligo2", "oligo3", "oligo3"]},
    index=["barcode1", "barcode2", "barcode3", "barcode4", "barcode5"],
)
COUNTS_DNA = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
COUNTS_RNA = np.array([[1, 2, 4, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])

FILTER = np.array(
    [
        [False, True, False],
        [False, False, False],
        [True, False, False],
        [False, False, True],
        [False, False, True],
    ]
)


@pytest.fixture
def mpra_data():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    return MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))


@pytest.fixture
def mpra_data_with_bc_filter(mpra_data):
    data = copy.deepcopy(mpra_data)
    data.var_filter = FILTER
    return data


@pytest.fixture
def files():
    input_file = os.path.join(os.path.dirname(__file__), "data", "reporter_experiment_barcode.input.tsv.gz")
    output_file_activity = tempfile.NamedTemporaryFile(delete=False).name
    output_file_barcode = tempfile.NamedTemporaryFile(delete=False).name
    yield {"input": input_file, "output_activity": output_file_activity, "output_barcode": output_file_barcode}
    os.remove(output_file_activity)
    os.remove(output_file_barcode)


def test_export_barcode_file(mpra_data, files):
    export_barcode_file(mpra_data, files["output_barcode"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "io", "reporter_experiment_barcode.simple.out.tsv.gz"
    )

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content


def test_export_barcode_file_var_filter(mpra_data_with_bc_filter, files):
    export_barcode_file(mpra_data_with_bc_filter, files["output_barcode"])
    assert os.path.exists(files["output_barcode"])

    with open(files["output_barcode"], "r") as f:
        output_content = f.read()
    expected_output_file = os.path.join(
        os.path.dirname(__file__), "data", "io", "reporter_experiment_barcode.simple_with_bc_filter.out.tsv.gz"
    )

    with gzip.open(expected_output_file, "rt") as f:
        expected_content = f.read()
    assert output_content == expected_content
