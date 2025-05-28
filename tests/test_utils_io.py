import os
import tempfile
import pandas as pd
import numpy as np
import pytest
import anndata as ad
from mpralib.utils.io import read_sequence_design_file, export_counts_file
from mpralib.exception import SequenceDesignException, MPRAlibException
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, MPRAData


class DummyMPRAData(MPRAData):
    def __init__(self, replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold, barcode_counts=None):
        layers = {"rna": rna_counts, "dna": dna_counts}
        obs = pd.DataFrame(index=replicates)
        var = pd.DataFrame({"oligo": oligos}, index=barcodes)
        super().__init__(ad.AnnData(X=rna_counts, obs=obs, var=var, layers=layers), barcode_threshold)
        if barcode_counts:
            self.barcode_counts = pd.DataFrame(
                barcode_counts,
                index=self.obs_names,
                columns=self.var_names,
            )


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
    dna_counts = np.array([[10, 0, 5], [20, 1, 0]])
    rna_counts = np.array([[5, 0, 2], [10, 1, 0]])
    barcode_counts = [[1, 0, 1], [1, 0, 1]]
    barcode_threshold = 1
    return DummyMPRAOligoData(replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold, barcode_counts)


def test_export_counts_file_barcode(tmp_path, barcode_data):
    out_path = tmp_path / "counts.tsv"
    export_counts_file(barcode_data, str(out_path), normalized=False)
    df = pd.read_csv(out_path, sep="\t", na_values=["X"])
    # Only rows with observed=True and barcode_counts >= threshold should remain
    assert "name" in df.columns
    assert set(df["name"]) == {"ol1", "ol2", "ol3"}
    assert "dna_count_rep1" in df.columns
    assert "rna_count_rep2" in df.columns
    assert "rna_count_rep2" in df.columns
    assert df.loc[df["name"] == "ol3", "dna_count_rep1"].values[0] == 5
    # Check that zeros are written as empty
    assert np.isnan(df.loc[df["name"] == "ol3", "dna_count_rep2"].values[0])


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
    assert np.all(df.filter(like="dna_count").applymap(lambda x: isinstance(x, float) or np.isnan(x)).values)  # type: ignore


def test_export_counts_file_with_filter(tmp_path, barcode_data):
    out_path = tmp_path / "counts_filt.tsv"
    # Mask out the first entry
    filter_mask = np.asarray([[True, False, False], [True, False, False]])
    export_counts_file(barcode_data, str(out_path), normalized=False, filter=filter_mask)
    df = pd.read_csv(out_path, sep="\t")
    # Only ol3 should remain (ol2 fails barcode_threshold, ol1 is masked)
    assert set(df["name"]) == {"ol2", "ol3"}


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
