import numpy as np
import pandas as pd
import pytest
import anndata as ad
import seaborn as sns
import copy
from mpralib.utils.plot import dna_vs_rna, correlation, barcodes_per_oligo
from mpralib.mpradata import MPRAData, MPRABarcodeData, MPRAOligoData, Modality


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
    return copy.deepcopy(DummyMPRABarcodeData(replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold))


@pytest.fixture
def oligo_data():
    replicates = ["rep1", "rep2"]
    barcodes = ["bc1", "bc2", "bc3"]
    oligos = ["ol1", "ol2", "ol3"]
    dna_counts = np.array([[10, 0, 5], [20, 1, 0]])
    rna_counts = np.array([[5, 0, 2], [10, 1, 0]])
    barcode_counts = [[1, 0, 1], [1, 0, 1]]
    barcode_threshold = 1
    return copy.deepcopy(
        DummyMPRAOligoData(replicates, oligos, barcodes, dna_counts, rna_counts, barcode_threshold, barcode_counts)
    )


@pytest.mark.parametrize(
    "layer",
    [
        Modality.DNA,
        Modality.RNA,
        Modality.RNA_NORMALIZED,
        Modality.DNA_NORMALIZED,
        Modality.ACTIVITY,
    ],
)
def test_correlation_returns_pairgrid(oligo_data, layer):
    g = correlation(oligo_data, layer)
    assert isinstance(g, sns.PairGrid)
    # Should have as many axes as replicates
    n = len(oligo_data.obs_names)
    assert g.axes.shape == (n, n)


def test_correlation_with_replicates_subset(oligo_data):
    g = correlation(oligo_data, Modality.DNA, replicates=["rep1", "rep2"])
    assert isinstance(g, sns.PairGrid)
    assert g.axes.shape == (2, 2)
    # Check that the columns are correct
    colnames = [ax.get_xlabel() for ax in g.axes[0]]
    assert any("Replicate" in c for c in colnames)


def test_correlation_masking(oligo_data):
    # Set barcode_counts below threshold for one oligo
    oligo_data.barcode_counts[0, 2] = 0
    g = correlation(oligo_data, Modality.DNA)
    # Should still return a PairGrid
    assert isinstance(g, sns.PairGrid)


@pytest.mark.parametrize(
    "layer",
    [
        Modality.DNA,
        Modality.RNA,
        Modality.RNA_NORMALIZED,
        Modality.DNA_NORMALIZED,
        Modality.ACTIVITY,
    ],
)
def test_correlation_returns_pairgrid_barcode(barcode_data, layer):
    g = correlation(barcode_data, layer)
    assert isinstance(g, sns.PairGrid)
    # Should have as many axes as replicates
    n = len(barcode_data.obs_names)
    assert g.axes.shape == (n, n)


def test_correlation_with_replicates_subset_barcode(barcode_data):
    g = correlation(barcode_data, Modality.DNA, replicates=["rep1", "rep2"])
    assert isinstance(g, sns.PairGrid)
    assert g.axes.shape == (2, 2)
    # Check that the columns are correct
    colnames = [ax.get_xlabel() for ax in g.axes[0]]
    assert any("Replicate" in c for c in colnames)


def test_dna_vs_rna_oligo_returns_jointgrid(oligo_data):
    g = dna_vs_rna(oligo_data)
    assert isinstance(g, sns.axisgrid.JointGrid)
    # Check axis labels
    assert g.ax_joint.get_xlabel() == "DNA [log10]"
    assert g.ax_joint.get_ylabel() == "RNA [log10]"


def test_dna_vs_rna_barcode_returns_jointgrid(barcode_data):
    g = dna_vs_rna(barcode_data)
    assert isinstance(g, sns.axisgrid.JointGrid)
    # Check axis labels
    assert g.ax_joint.get_xlabel() == "DNA [log10]"
    assert g.ax_joint.get_ylabel() == "RNA [log10]"


def test_dna_vs_rna_oligo_with_replicates(oligo_data):
    g = dna_vs_rna(oligo_data, replicates=["rep2"])
    assert isinstance(g, sns.axisgrid.JointGrid)
    # Check axis labels
    assert g.ax_joint.get_xlabel() == "DNA [log10]"
    assert g.ax_joint.get_ylabel() == "RNA [log10]"


def test_barcodes_per_oligo(oligo_data):
    g = barcodes_per_oligo(oligo_data)
    assert isinstance(g, sns.axisgrid.FacetGrid)
    assert g.axes.shape == (1, 2)


def test_barcodes_per_oligo_with_replicates(oligo_data):
    g = barcodes_per_oligo(oligo_data, replicates=["rep2"])
    assert isinstance(g, sns.axisgrid.FacetGrid)
    assert g.axes.shape == (1, 1)
