import numpy as np
import pandas as pd
import anndata as ad
import copy
import pytest
from mpralib.mpradata import MPRABarcodeData, CountSampling, BarcodeFilter, MPRAData, Modality, MPRAOligoData
from mpralib.mpradata import MPRAlibException


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


def test_apply_count_sampling_rna(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA, proportion=0.5)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    assert np.all(rna_sampling <= np.asarray(mpra_data.data.layers["rna"]))
    assert np.all(rna_sampling >= 0)


def test_apply_count_sampling_dna(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.DNA, proportion=0.5)
    dna_sampling = np.asarray(mpra_data.data.layers["dna_sampling"])
    assert np.all(dna_sampling <= mpra_data.data.layers["dna"])
    assert np.all(dna_sampling >= 0)


def test_apply_count_sampling_rna_and_dna(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=0.5)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    dna_sampling = np.asarray(mpra_data.data.layers["dna_sampling"])
    assert np.all(rna_sampling <= mpra_data.data.layers["rna"])
    assert np.all(rna_sampling >= 0)
    assert np.all(dna_sampling <= np.asarray(mpra_data.data.layers["dna"]))
    assert np.all(dna_sampling >= 0)


def test_apply_count_sampling_rna_total(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA, total=10)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    assert np.sum(rna_sampling) <= 30


def test_apply_count_sampling_total(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, total=10)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    dna_sampling = np.asarray(mpra_data.data.layers["dna_sampling"])
    assert np.sum(rna_sampling) <= 30
    assert np.sum(dna_sampling) <= 30


def test_apply_count_sampling_max_value_rna(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA, max_value=2)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    assert np.all(rna_sampling <= 2)


def test_apply_count_sampling_max_value(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, max_value=2)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    dna_sampling = np.asarray(mpra_data.data.layers["dna_sampling"])
    assert np.all(rna_sampling <= 2)
    assert np.all(dna_sampling <= 2)


def test_apply_count_sampling_aggregate_over_replicates(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, total=10, aggregate_over_replicates=True)
    rna_sampling = np.asarray(mpra_data.data.layers["rna_sampling"])
    dna_sampling = np.asarray(mpra_data.data.layers["dna_sampling"])
    assert np.sum(rna_sampling) <= 11
    assert np.sum(dna_sampling) <= 11


def test_barcode_counts(mpra_data):
    supporting_barcodes = mpra_data.barcode_counts
    expected_barcodes = np.array([[2, 2, 1, 2, 2], [2, 2, 1, 2, 2], [2, 2, 1, 2, 2]])
    np.testing.assert_array_equal(supporting_barcodes, expected_barcodes)


def test_barcode_counts_with_filter(mpra_data_with_bc_filter):
    supporting_barcodes = mpra_data_with_bc_filter.barcode_counts
    expected_barcodes = np.array([[2, 2, 1, 2, 2], [2, 2, 1, 2, 2], [2, 2, 1, 2, 2]])
    np.testing.assert_array_equal(supporting_barcodes, expected_barcodes)


def test_raw_dna_counts(mpra_data):
    expected_dna_counts = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
    np.testing.assert_array_equal(mpra_data.raw_dna_counts, expected_dna_counts)


def test_raw_dna_counts_with_modification(mpra_data):
    mpra_data.data.layers["dna"] = np.array([[10, 20, 30, 10, 20], [40, 50, 60, 40, 50], [70, 80, 90, 100, 1000]])
    expected_dna_counts = np.array([[10, 20, 30, 10, 20], [40, 50, 60, 40, 50], [70, 80, 90, 100, 1000]])
    np.testing.assert_array_equal(mpra_data.raw_dna_counts, expected_dna_counts)


def test_filtered_dna_counts(mpra_data, mpra_data_with_bc_filter):
    expected_dna_counts = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
    np.testing.assert_array_equal(mpra_data.dna_counts, expected_dna_counts)
    expected_filtered_dna_counts = np.array([[1, 2, 0, 1, 2], [0, 5, 6, 4, 5], [7, 8, 9, 0, 0]])
    np.testing.assert_array_equal(mpra_data_with_bc_filter.dna_counts, expected_filtered_dna_counts)


def test_dna_counts_with_sampling(mpra_data):
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.DNA, proportion=0.5)
    dna_sampling = mpra_data.data.layers["dna_sampling"]
    np.testing.assert_array_equal(mpra_data.dna_counts, dna_sampling)


def test_dna_counts_with_filter(mpra_data, mpra_data_with_bc_filter):
    mpra_data.apply_count_sampling(CountSampling.DNA, max_value=2)
    expected_filtered_dna_counts = np.array([[1, 2, 2, 1, 2], [2, 2, 2, 2, 2], [2, 2, 2, 2, 2]])
    np.testing.assert_array_equal(mpra_data.dna_counts, expected_filtered_dna_counts)
    mpra_data_with_bc_filter.apply_count_sampling(CountSampling.DNA, max_value=2)
    expected_filtered_dna_counts = np.array([[1, 2, 0, 1, 2], [0, 2, 2, 2, 2], [2, 2, 2, 0, 0]])
    np.testing.assert_array_equal(mpra_data_with_bc_filter.dna_counts, expected_filtered_dna_counts)


@pytest.fixture
def mpra_data_barcode():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    return MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))


def test_apply_barcode_filter_min_count(mpra_data_barcode):
    mpra_data_barcode.apply_barcode_filter(BarcodeFilter.MIN_COUNT, params={"rna_min_count": 4, "dna_min_count": 3})
    expected_filter = np.array(
        [
            [True, False, False],
            [True, False, False],
            [False, False, False],
            [True, False, False],
            [True, False, False],
        ]
    )
    np.testing.assert_array_equal(mpra_data_barcode.var_filter, expected_filter)


def test_barcode_filter_other_to_min_max_counts(mpra_data_barcode):

    result = mpra_data_barcode._barcode_filter_min_max_counts(
        BarcodeFilter.LARGE_EXPRESSION, counts=mpra_data_barcode.raw_dna_counts, count_threshold=1
    )
    expected_filter = np.array(
        [
            [False, False, False],
            [False, False, False],
            [False, False, False],
            [False, False, False],
            [False, False, False],
        ]
    )
    np.testing.assert_array_equal(result, expected_filter)


def test_apply_barcode_filter_max_count(mpra_data_barcode):
    mpra_data_barcode.apply_barcode_filter(BarcodeFilter.MAX_COUNT, params={"rna_max_count": 9, "dna_max_count": 100})
    expected_filter = np.array(
        [
            [False, False, False],
            [False, False, False],
            [False, False, False],
            [False, False, True],
            [False, False, True],
        ]
    )
    np.testing.assert_array_equal(mpra_data_barcode.var_filter, expected_filter)

    mpra_data_barcode.var_filter = None
    mpra_data_barcode.apply_barcode_filter(BarcodeFilter.MAX_COUNT, params={"dna_max_count": 99})
    expected_filter = np.array(
        [
            [False, False, False],
            [False, False, False],
            [False, False, False],
            [False, False, False],
            [False, False, True],
        ]
    )
    np.testing.assert_array_equal(mpra_data_barcode.var_filter, expected_filter)  # type: ignore


@pytest.fixture
def mpra_data_norm():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))
    data.scaling = 10
    return data


@pytest.fixture
def mpra_data_norm_with_bc_filter(mpra_data_norm):
    data = copy.deepcopy(mpra_data_norm)
    data.var_filter = FILTER
    return data


def test_normalize_counts(mpra_data_norm):
    mpra_data_norm._normalize()
    dna_normalized = mpra_data_norm.normalized_dna_counts
    expected_dna_normalized = np.array(
        [
            [1.428, 2.142, 2.857, 1.428, 2.142],
            [1.724, 2.068, 2.413, 1.724, 2.068],
            [0.575, 0.647, 0.719, 0.791, 7.266],
        ]
    )
    np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)
    expected_rna_normalized = np.array(
        [[1.333, 2.0, 3.333, 1.333, 2.0], [1.724, 2.069, 2.414, 1.724, 2.069], [0.576, 0.647, 0.719, 0.791, 7.266]]
    )
    rna_normalized = mpra_data_norm.normalized_rna_counts
    np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)


def test_normalize_without_pseudocount(mpra_data_norm):
    mpra_data = copy.deepcopy(mpra_data_norm)
    mpra_data.pseudo_count = 0
    mpra_data._normalize()
    dna_normalized = np.asarray(mpra_data.data.layers["dna_normalized"])
    expected_dna_normalized = np.array(
        [[1.111, 2.222, 3.333, 1.111, 2.222], [1.667, 2.083, 2.5, 1.667, 2.083], [0.522, 0.597, 0.672, 0.746, 7.463]]
    )
    np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)
    expected_rna_normalized = np.array(
        [[1.0, 2.0, 4.0, 1.0, 2.0], [1.667, 2.083, 2.5, 1.667, 2.083], [0.522, 0.597, 0.672, 0.746, 7.463]]
    )
    rna_normalized = np.asarray(mpra_data.data.layers["rna_normalized"])
    np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)


def test_normalize_counts_with_bc_filter(mpra_data_norm_with_bc_filter):
    mpra_data_norm_with_bc_filter._normalize()

    dna_normalized = np.asarray(mpra_data_norm_with_bc_filter.data.layers["dna_normalized"])
    expected_normalized = np.array(
        [
            [1.428, 2.142, 2.857, 1.428, 2.142],
            [1.724, 2.068, 2.413, 1.724, 2.068],
            [0.575, 0.647, 0.719, 0.791, 7.266],
        ]
    )
    np.testing.assert_almost_equal(dna_normalized, expected_normalized, decimal=3)

    dna_normalized = mpra_data_norm_with_bc_filter.normalized_dna_counts
    expected_normalized = np.array(
        [
            [1.428, 2.142, 0.0, 1.428, 2.142],
            [0.0, 2.068, 2.413, 1.724, 2.068],
            [0.575, 0.647, 0.719, 0.0, 0.0],
        ]
    )
    np.testing.assert_almost_equal(dna_normalized, expected_normalized, decimal=3)

    rna_normalized = np.asarray(mpra_data_norm_with_bc_filter.data.layers["rna_normalized"])
    expected_normalized = np.array(
        [[1.333, 2.0, 3.333, 1.333, 2.0], [1.724, 2.069, 2.414, 1.724, 2.069], [0.576, 0.647, 0.719, 0.791, 7.266]]
    )
    np.testing.assert_almost_equal(rna_normalized, expected_normalized, decimal=3)
    rna_normalized = mpra_data_norm_with_bc_filter.normalized_rna_counts
    expected_normalized = np.array(
        [[1.333, 2.0, 0.0, 1.333, 2.0], [0.0, 2.069, 2.414, 1.724, 2.069], [0.576, 0.647, 0.719, 0.0, 0.0]]
    )
    np.testing.assert_almost_equal(rna_normalized, expected_normalized, decimal=3)


@pytest.fixture
def mpra_oligo_data():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    mpra_barcode_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))
    data = mpra_barcode_data.oligo_data
    data.scaling = 10
    return data


@pytest.fixture
def mpra_oligo_data_with_bc_filter():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    mpra_barcode_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))
    mpra_barcode_data.var_filter = FILTER
    data = mpra_barcode_data.oligo_data
    data.scaling = 10
    return data


def test_oligo_normalize_counts(mpra_oligo_data):
    dna_normalized = mpra_oligo_data.normalized_dna_counts
    expected_dna_normalized = np.array([[1.786, 2.857, 1.786], [1.897, 2.414, 1.897], [0.612, 0.719, 4.029]])
    np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)
    expected_rna_normalized = np.array([[1.667, 3.333, 1.667], [1.897, 2.414, 1.897], [0.612, 0.719, 4.029]])
    rna_normalized = mpra_oligo_data.normalized_rna_counts
    np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)


def test_oligo_normalize_without_pseudocount(mpra_oligo_data):
    mpra_data = copy.deepcopy(mpra_oligo_data)
    mpra_data.pseudo_count = 0
    dna_normalized = mpra_data.normalized_dna_counts
    expected_dna_normalized = np.array([[1.667, 3.333, 1.667], [1.875, 2.5, 1.875], [0.56, 0.672, 4.104]])
    np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)
    expected_rna_normalized = np.array([[1.5, 4.0, 1.5], [1.875, 2.5, 1.875], [0.56, 0.672, 4.104]])
    rna_normalized = mpra_data.normalized_rna_counts
    np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)


def test_oligo_normalize_counts_with_bc_filter(mpra_oligo_data_with_bc_filter):
    dna_normalized = mpra_oligo_data_with_bc_filter.normalized_dna_counts
    expected_normalized = np.array([[2.5, 0.0, 2.5], [2.5, 2.917, 2.292], [3.148, 3.704, 0.0]])
    np.testing.assert_almost_equal(dna_normalized, expected_normalized, decimal=3)
    rna_normalized = mpra_oligo_data_with_bc_filter.normalized_rna_counts
    np.testing.assert_almost_equal(rna_normalized, expected_normalized, decimal=3)


@pytest.fixture
def mpra_corr_data():
    layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
    return MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers)).oligo_data


def test_correlation(mpra_corr_data):
    mpra_corr_data._compute_correlation(mpra_corr_data.activity, "log2FoldChange")
    assert "pearson_correlation_log2FoldChange" in mpra_corr_data.data.obsp
    assert "spearman_correlation_log2FoldChange" in mpra_corr_data.data.obsp


def test_pearson_correlation(mpra_corr_data):
    x = mpra_corr_data.correlation(method="pearson", count_type=Modality.ACTIVITY)
    y = mpra_corr_data.correlation(method="pearson", count_type=Modality.RNA_NORMALIZED)
    z = mpra_corr_data.correlation(method="pearson", count_type=Modality.DNA_NORMALIZED)
    np.testing.assert_equal(x, np.array([[1.0, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]]))
    np.testing.assert_almost_equal(
        y, np.array([[1.0, 1.0, -0.475752], [1.0, 1.0, -0.475752], [-0.475752, -0.475752, 1.0]]), decimal=3
    )
    np.testing.assert_almost_equal(z, np.array([[1.0, 1.0, -0.476], [1.0, 1.0, -0.476], [-0.476, -0.476, 1.0]]), decimal=3)


def test_spearman_correlation(mpra_corr_data):
    x = mpra_corr_data.correlation(method="spearman", count_type=Modality.ACTIVITY)
    y = mpra_corr_data.correlation(method="spearman", count_type=Modality.RNA_NORMALIZED)
    z = mpra_corr_data.correlation(method="spearman", count_type=Modality.DNA_NORMALIZED)
    np.testing.assert_equal(x, np.array([[1.0, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]]))
    np.testing.assert_almost_equal(y, np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), decimal=3)
    np.testing.assert_equal(z, np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))


@pytest.fixture
def mpra_complexity_data():
    counts_dna = np.array([[0, 2, 0, 1, 2], [4, 5, 6, 4, 5], [7, 0, 9, 10, 0]])
    counts_rna = np.array([[1, 2, 0, 1, 2], [4, 5, 6, 4, 5], [7, 0, 9, 10, 0]])
    layers = {"rna": counts_rna, "dna": counts_dna}
    return MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))


def test_lincoln_complexity(mpra_complexity_data):
    complexity = mpra_complexity_data.complexity()
    np.testing.assert_equal(complexity, np.array([[4, 5, 6], [5, 5, 5], [6, 5, 3]]))


def test_chapman_complexity(mpra_complexity_data):
    complexity = mpra_complexity_data.complexity(method="chapman")
    np.testing.assert_equal(complexity, np.array([[4, 5, 5], [5, 5, 5], [5, 5, 3]]))


def test_fail_complexity(mpra_complexity_data):
    with pytest.raises(ValueError):
        mpra_complexity_data.complexity(method="unknown")


def test_read_and_write(tmp_path, mpra_data):
    out_path = tmp_path / "bc_data.h5ad"
    mpra_data.write(out_path)

    data = MPRABarcodeData.read(out_path)
    assert isinstance(data, MPRABarcodeData)
    assert data.data.shape == mpra_data.data.shape
    assert data.data.layers.keys() == mpra_data.data.layers.keys()
    assert np.all(data.rna_counts == mpra_data.rna_counts)
    assert np.all(data.activity == mpra_data.activity)
    assert data.pseudo_count == mpra_data.pseudo_count
    assert data.scaling == mpra_data.scaling


def test_read_and_write_oligo(tmp_path, mpra_oligo_data):
    out_path = tmp_path / "oligo_data.h5ad"
    mpra_oligo_data.write(out_path)

    data = MPRAOligoData.read(out_path)
    assert isinstance(data, MPRAOligoData)
    assert data.data.shape == mpra_oligo_data.data.shape
    assert data.data.layers.keys() == mpra_oligo_data.data.layers.keys()
    assert np.all(data.rna_counts == mpra_oligo_data.rna_counts)
    assert np.all(data.activity == mpra_oligo_data.activity)


def test_read_and_write_with_modifications(tmp_path, mpra_data):
    out_path = tmp_path / "bc_data_mod.h5ad"
    mpra_data.scaling = 10.0
    mpra_data.pseudo_count = 0
    mpra_data.write(out_path)

    data = MPRABarcodeData.read(out_path)
    assert isinstance(data, MPRABarcodeData)
    assert data.scaling == 10.0
    assert data.pseudo_count == 0


def test_modality_from_string():
    assert Modality.from_string("DNA") == Modality.DNA
    assert Modality.from_string("dna") == Modality.DNA
    assert Modality.from_string("RNA") == Modality.RNA
    assert Modality.from_string("rna_normalized") == Modality.RNA_NORMALIZED
    assert Modality.from_string("ACTIVITY") == Modality.ACTIVITY
    with pytest.raises(ValueError):
        Modality.from_string("invalid_modality")


def test_barcodefilter_from_string():
    assert BarcodeFilter.from_string("MIN_BCS_PER_OLIGO") == BarcodeFilter.MIN_BCS_PER_OLIGO
    assert BarcodeFilter.from_string("min_bcs_per_oligo") == BarcodeFilter.MIN_BCS_PER_OLIGO
    assert BarcodeFilter.from_string("GLOBAL") == BarcodeFilter.GLOBAL
    assert BarcodeFilter.from_string("global") == BarcodeFilter.GLOBAL
    assert BarcodeFilter.from_string("OLIGO_SPECIFIC") == BarcodeFilter.OLIGO_SPECIFIC
    assert BarcodeFilter.from_string("LARGE_EXPRESSION") == BarcodeFilter.LARGE_EXPRESSION
    assert BarcodeFilter.from_string("random") == BarcodeFilter.RANDOM
    assert BarcodeFilter.from_string("MIN_COUNT") == BarcodeFilter.MIN_COUNT
    assert BarcodeFilter.from_string("max_count") == BarcodeFilter.MAX_COUNT
    with pytest.raises(ValueError):
        BarcodeFilter.from_string("not_a_filter")


def test_barcode_filter_min_bcs_per_oligo_basic(mpra_data):
    # Should flag all barcodes that have less than the minimum number of barcodes per oligo
    # In the test data, no barcode is a global outlier, so all should be False
    mask = mpra_data._barcode_filter_min_bcs_per_oligo()
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_min_bcs_per_oligo_filter_oligo2(mpra_data):
    # Should flag all barcodes that have less than the minimum number of barcodes per oligo
    # In the test data, oligo2 has only 1 barcode, so it should be flagged when threshold is 2
    mask = mpra_data._barcode_filter_min_bcs_per_oligo(threshold=2)
    expected = np.zeros_like(mask, dtype=bool)
    expected[2, :] = True
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_basic(mpra_data):
    # Should flag barcodes with RNA counts that are global outliers (z-score > 3)
    # In the test data, no barcode is a global outlier, so all should be False
    mask = mpra_data._barcode_filter_global_outliers()
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_basic_with_bc_threshold(mpra_data):
    # Should flag barcodes with RNA counts that are global outliers (z-score > 3)
    # In the test data, no barcode is a global outlier, so all should be False
    mpra_data.barcode_threshold = 2
    mask = mpra_data._barcode_filter_global_outliers(apply_bc_threshold=True)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_basic_with_bc_threshold_aggregated(mpra_data):
    # Should flag barcodes with RNA counts that are global outliers (z-score > 3)
    # In the test data, no barcode is a global outlier, so all should be False
    mpra_data.data.layers["rna"][0, 0] = 0
    mpra_data.data.layers["dna"][0, 0] = 0
    mpra_data.barcode_threshold = 2
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=0.1, apply_bc_threshold=True)
    expected = np.ones_like(mask, dtype=bool)
    expected[[0, 1], 0] = False  # not applied to oligo1 of rep1 (less than 2 barcodes observed)
    expected[2, :] = False  # not applied to oligo2 (1 barcode)
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=0.1, apply_bc_threshold=True, aggregated_bc_threshold=True)
    expected = np.ones_like(mask, dtype=bool)
    expected[0, 0] = False  # not applied to bcarcode1 of rep1 because not observed
    expected[2, :] = False  # not applied to oligo2 (1 barcode)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_basic_with_bc_threshold_with_outlier(mpra_data):
    # Should flag barcodes with RNA counts that are global outliers (z-score > 3)
    # In the test data, no barcode is a global outlier, so all should be False
    mpra_data.barcode_threshold = 2
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=1.7, apply_bc_threshold=True)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=1.4, apply_bc_threshold=True)
    expected[4, 2] = True
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_basic_with_outlier(mpra_data):
    # Should flag barcodes with RNA counts that are global outliers (z-score > 3)
    # In the test data, no barcode is a global outlier, so all should be False
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=1.7)
    expected = np.zeros_like(mask, dtype=bool)
    expected[4, 2] = True
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_with_outliers(mpra_data):
    # Introduce a global outlier in RNA counts for a barcode
    mpra_data.data.layers["rna"][1, 0] = 100  # Make barcode2 in rep1 a strong outlier
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=1.7)
    # Only barcode2 in rep1 should be flagged as True
    expected = np.zeros_like(mask, dtype=bool)
    expected[0, 1] = True
    expected[4, 2] = True
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_all_rna_zero(mpra_data):
    # Set all RNA counts to zero, should not raise error and all should be False
    mpra_data.data.layers["rna"][:] = 0
    mask = mpra_data._barcode_filter_global_outliers(times_zscore=1.0)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_global_outliers_all_zero_true(mpra_data):
    # Set all RNA counts to zero, should not raise error and all should be False
    mpra_data.data.layers["rna"][:] = 0
    mpra_data.data.layers["dna"][:] = 0
    mask = mpra_data._barcode_filter_global_outliers()
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_oligo_specific_outliers_basic(mpra_data):
    # No barcode should be flagged as outlier with high z-score threshold
    mask = mpra_data._barcode_filter_oligo_specific_outliers(times_zscore=3.0)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_oligo_specific_outliers_with_outlier(mpra_data):
    # change to only 2 oligos becasue we want to have different zscores (at least 3 barcodes needed)
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mask = mpra_data._barcode_filter_oligo_specific_outliers(times_zscore=1.15)
    expected = np.zeros_like(mask, dtype=bool)
    expected[4, 2] = True  # barcode5, rep3
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_oligo_specific_outliers_with_outliers(mpra_data):
    # change to only 2 oligos becasue we want to have different zscores (at least 3 barcodes needed)
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mask = mpra_data._barcode_filter_oligo_specific_outliers(times_zscore=1.0)
    expected = np.zeros_like(mask, dtype=bool)
    expected[4, 2] = True  # barcode5, rep3
    expected[2, 0] = True  # barcode3, rep1
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_oligo_specific_outliers_all_zero(mpra_data):
    # All RNA counts zero, should not raise error and all should be False
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.data.layers["rna"][:] = 0
    mask = mpra_data._barcode_filter_oligo_specific_outliers(times_zscore=1.0)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_oligo_specific_outliers_all_zero_true(mpra_data):
    # All RNA counts zero, should not raise error and all should be False
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.data.layers["rna"][:] = 0
    mpra_data.data.layers["dna"][:] = 0
    mask = mpra_data._barcode_filter_oligo_specific_outliers(times_zscore=1.0)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_large_expression_outliers_basic(mpra_data):
    # No barcode should be flagged as outlier with high z-score threshold
    mask = mpra_data._barcode_filter_large_expression_outliers()
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_large_expression_outliers_with_outlier(mpra_data):
    # change to only 2 oligos becasue we want to have different zscores (at least 3 barcodes needed)
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mask = mpra_data._barcode_filter_large_expression_outliers(times_activity=0.12)
    expected = np.zeros_like(mask, dtype=bool)
    expected[2, :] = True  # barcode3, all replicates
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_large_expression_outliers_with_outliers(mpra_data):
    # change to only 2 oligos becasue we want to have different zscores (at least 3 barcodes needed)
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.data.layers["rna"][0, 0] = 100  # Make barcode1 in rep1 a strong outlier
    mask = mpra_data._barcode_filter_large_expression_outliers(times_activity=1.14)
    expected = np.zeros_like(mask, dtype=bool)
    expected[0, :] = True  # barcode1, all replicates
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_large_expression_outliers_all_zero(mpra_data):
    # All RNA counts zero, should not raise error and all should be False
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.data.layers["rna"][:] = 0
    mask = mpra_data._barcode_filter_large_expression_outliers(times_activity=0.6)
    expected = np.zeros_like(mask, dtype=bool)
    expected[3, :] = True  # barcode3, all replicates
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_large_expression_outliers_all_zero_true(mpra_data):
    # All RNA counts zero, should not raise error and all should be False
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.data.layers["rna"][:] = 0
    mpra_data.data.layers["dna"][:] = 0
    mask = mpra_data._barcode_filter_large_expression_outliers(times_activity=0.6)
    expected = np.zeros_like(mask, dtype=bool)
    np.testing.assert_array_equal(mask, expected)


def test_barcode_filter_combine(mpra_data):
    # change to only 2 oligos because we want to have different zscores (at least 3 barcodes needed)
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo3", "oligo3", "oligo3"]
    mpra_data.apply_barcode_filter(BarcodeFilter.MIN_BCS_PER_OLIGO, {"threshold": 3})
    mpra_data.apply_barcode_filter(BarcodeFilter.OLIGO_SPECIFIC, {"times_zscore": 1.0})
    mask = mpra_data.var_filter
    expected = np.zeros_like(mask, dtype=bool)
    expected[0, :] = True  # barcode1, all replicates
    expected[1, :] = True  # barcode2, all replicates
    expected[4, 2] = True  # barcode5, rep3
    expected[2, 0] = True  # barcode3, rep1
    np.testing.assert_array_equal(mask, expected)


def test_mpra_data_from_file():

    assert MPRAData.from_file("test") is None

    assert MPRAData.from_file("ffasfasdfas") is None


def test_mprabarcode_from_file(tmp_path):
    # Create a mock barcode-level MPRA TSV file
    file_path = tmp_path / "test_bc.tsv"
    df = pd.DataFrame(
        {
            "barcode": ["barcode1", "barcode2", "barcode3"],
            "oligo_name": ["oligo1", "oligo1", "oligo2"],
            "dna_count_rep1": [10, 20, 30],
            "rna_count_rep1": [1, 2, 3],
            "dna_count_rep2": [40, 50, 60],
            "rna_count_rep2": [4, 5, 6],
        }
    )
    df_out = df[["barcode", "oligo_name", "dna_count_rep1", "rna_count_rep1", "dna_count_rep2", "rna_count_rep2"]]
    df_out.to_csv(file_path, sep="\t", index=False)

    # Read using from_file
    data = MPRABarcodeData.from_file(str(file_path))
    assert isinstance(data, MPRABarcodeData)
    assert data.data.shape == (2, 3)  # 2 replicates, 3 barcodes
    assert "rna" in data.data.layers
    assert "dna" in data.data.layers
    assert all(data.data.var["oligo"].isin(["oligo1", "oligo2"]))
    assert data.data.uns["file_path"] == str(file_path)
    assert data.data.uns["normalized"] is False


def test_mpraoligo_from_file(tmp_path):
    # Create a mock AnnData file for oligo-level data
    X = np.array([[1, 2], [3, 4]])
    adata = ad.AnnData(X)
    adata.layers["rna"] = X
    adata.layers["dna"] = X + 10
    adata.var_names = ["oligoA", "oligoB"]
    adata.obs_names = ["rep1", "rep2"]
    file_path = tmp_path / "test_oligo.h5ad"
    adata.write(file_path)

    data = MPRAOligoData.from_file(str(file_path))
    assert isinstance(data, MPRAOligoData)
    assert data.data.shape == (2, 2)
    assert "rna" in data.data.layers
    assert "dna" in data.data.layers
    assert list(data.data.var_names) == ["oligoA", "oligoB"]
    assert list(data.data.obs_names) == ["rep1", "rep2"]


def test_drop_correlation_removes_correlation_layers(mpra_data_norm):
    # Compute correlation to create correlation layers
    mpra_data_norm._normalize()
    mpra_data_norm._compute_correlation(mpra_data_norm.activity, "activity")
    mpra_data_norm._compute_correlation(mpra_data_norm.normalized_rna_counts, "rna_normalized")
    mpra_data_norm._compute_correlation(mpra_data_norm.normalized_dna_counts, "dna_normalized")

    # Ensure correlation layers exist
    for layer in ["activity", "rna_normalized", "dna_normalized"]:
        for method in ["pearson", "spearman"]:
            assert f"{method}_correlation_{layer}" in mpra_data_norm.data.obsp
            assert f"{method}_correlation_{layer}_pvalue" in mpra_data_norm.data.obsp
        assert mpra_data_norm._get_metadata(f"correlation_{layer}") is True

    # Drop correlation
    mpra_data_norm._drop_correlation()

    # Ensure correlation layers are removed and metadata is updated
    for layer in ["activity", "rna_normalized", "dna_normalized"]:
        for method in ["pearson", "spearman"]:
            assert f"{method}_correlation_{layer}" not in mpra_data_norm.data.obsp
            assert f"{method}_correlation_{layer}_pvalue" not in mpra_data_norm.data.obsp
        assert mpra_data_norm._get_metadata(f"correlation_{layer}") is False


def test_drop_correlation_no_error_if_layers_missing(mpra_data_norm):
    # Ensure no correlation layers exist
    for layer in ["activity", "rna_normalized", "dna_normalized"]:
        for method in ["pearson", "spearman"]:
            mpra_data_norm.data.obsp.pop(f"{method}_correlation_{layer}", None)
            mpra_data_norm.data.obsp.pop(f"{method}_correlation_{layer}_pvalue", None)
        mpra_data_norm._add_metadata(f"correlation_{layer}", False)

    # Should not raise error
    mpra_data_norm._drop_correlation()
    # Metadata should remain False
    for layer in ["activity", "rna_normalized", "dna_normalized"]:
        assert mpra_data_norm._get_metadata(f"correlation_{layer}") is False


def test_data_property_getter(mpra_data):
    # Should return the underlying AnnData object
    assert isinstance(mpra_data.data, ad.AnnData)
    # Should have expected shape and layers
    assert mpra_data.data.shape == (3, 5)
    assert "rna" in mpra_data.data.layers
    assert "dna" in mpra_data.data.layers


def test_data_property_setter(mpra_data):
    # Create a new AnnData object and set it
    new_counts = np.ones((3, 5), dtype=np.int32)
    new_obs = pd.DataFrame(index=["repA", "repB", "repC"])
    new_var = pd.DataFrame({"oligo": ["oligoA"] * 5}, index=[f"barcode{i}" for i in range(1, 6)])
    new_layers = {"rna": new_counts.copy(), "dna": new_counts.copy()}
    new_adata = ad.AnnData(X=new_counts.copy(), obs=new_obs.copy(), var=new_var.copy(), layers=new_layers)
    mpra_data.data = new_adata
    assert mpra_data.data is new_adata
    assert mpra_data.data.shape == (3, 5)
    assert np.all(mpra_data.data.layers["rna"] == 1)
    assert np.all(mpra_data.data.layers["dna"] == 1)


@pytest.mark.parametrize(
    "has_sampling",
    [False, True],
)
@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_rna_counts_var_filter_and_sampling(mpra_data, has_sampling, modality):
    # Set up a filter that masks some barcodes
    filter_mask = np.zeros_like(mpra_data.var_filter, dtype=bool)
    filter_mask[1, 0] = True  # mask barcode2, rep1
    filter_mask[3, 2] = True  # mask barcode4, rep3
    mpra_data.var_filter = filter_mask

    if has_sampling:
        # Simulate a sampling layer
        if modality == "rna":
            sampled = mpra_data.raw_rna_counts.copy()
        else:
            sampled = mpra_data.raw_dna_counts.copy()
        sampled[0, 0] = 0
        sampled[2, 1] = 0
        mpra_data.data.layers[f"{modality}_sampling"] = sampled
        expected = sampled * ~filter_mask.T
    else:
        if modality == "rna":
            expected = mpra_data.raw_rna_counts * ~filter_mask.T
        else:
            expected = mpra_data.raw_dna_counts * ~filter_mask.T

    if modality == "rna":
        result = mpra_data.rna_counts
    else:
        result = mpra_data.dna_counts
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_rna_counts_no_var_filter(mpra_data, modality):
    mpra_data.var_filter = None
    if modality == "rna":
        expected = mpra_data.raw_rna_counts
        result = mpra_data.rna_counts
    else:
        expected = mpra_data.raw_dna_counts
        result = mpra_data.dna_counts
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_rna_counts_all_filtered(mpra_data, modality):
    mpra_data.var_filter = np.ones_like(mpra_data.var_filter, dtype=bool)
    if modality == "rna":
        expected = np.zeros_like(mpra_data.raw_rna_counts)
        result = mpra_data.rna_counts
    else:
        expected = np.zeros_like(mpra_data.raw_dna_counts)
        result = mpra_data.dna_counts
    np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_rna_counts_with_sampling_and_all_filtered(mpra_data, modality):
    mpra_data.var_filter = np.ones_like(mpra_data.var_filter, dtype=bool)
    if modality == "rna":
        sampled = mpra_data.raw_rna_counts.copy()
    else:
        sampled = mpra_data.raw_dna_counts.copy()
    sampled[0, 0] = 99
    mpra_data.data.layers[f"{modality}_sampling"] = sampled
    expected = np.zeros_like(sampled)
    if modality == "rna":
        result = mpra_data.rna_counts
    else:
        result = mpra_data.dna_counts
    np.testing.assert_array_equal(result, expected)


def test_mpra_data_data_property(mpra_data):
    assert isinstance(mpra_data.data, ad.AnnData)
    new_data = copy.deepcopy(mpra_data.data)
    mpra_data.data = new_data
    assert mpra_data.data is new_data


def test_mpra_data_obs_var_names(mpra_data):
    assert list(mpra_data.obs_names) == ["rep1", "rep2", "rep3"]
    assert list(mpra_data.var_names) == ["barcode1", "barcode2", "barcode3", "barcode4", "barcode5"]
    assert mpra_data.n_obs == 3
    assert mpra_data.n_vars == 5


def test_mpra_data_oligos(mpra_data):
    oligos = mpra_data.oligos
    assert isinstance(oligos, pd.Series)
    assert set(oligos) == {"oligo1", "oligo2", "oligo3"}


def test_mpra_data_raw_counts(mpra_data):
    assert np.array_equal(mpra_data.raw_rna_counts, COUNTS_RNA)
    assert np.array_equal(mpra_data.raw_dna_counts, COUNTS_DNA)


def test_mpra_data_total_counts(mpra_data):
    # Should match sum over axis 1
    assert np.array_equal(mpra_data.total_rna_counts, COUNTS_RNA.sum(axis=1))
    assert np.array_equal(mpra_data.total_dna_counts, COUNTS_DNA.sum(axis=1))


def test_mpra_data_drop_total_counts(mpra_data):
    mpra_data.total_rna_counts  # populate
    mpra_data.total_dna_counts
    mpra_data.drop_total_counts()
    assert "rna_counts" not in mpra_data.data.obs
    assert "dna_counts" not in mpra_data.data.obs


def test_mpra_data_var_filter_setter_and_getter(mpra_data):
    mpra_data.var_filter = FILTER
    assert np.array_equal(mpra_data.var_filter, FILTER)
    mpra_data.var_filter = None
    assert np.all(mpra_data.var_filter == False)  # noqa: E712


def test_mpra_data_add_sequence_design(mpra_data):
    # Create a minimal sequence design DataFrame
    df = pd.DataFrame(
        {
            "category": ["cat1", "cat2", "cat3"],
            "class": ["class1", "class2", "class3"],
            "ref": ["ref1", "ref2", "ref3"],
            "chr": ["chr1", "chr2", "chr3"],
            "start": [1, 2, 3],
            "end": [10, 20, 30],
            "strand": ["+", "-", "+"],
            "variant_class": ["vc1", "vc2", "vc3"],
            "variant_pos": [100, 200, 300],
            "SPDI": [["spdi1"], ["spdi2"], ["spdi3"]],
            "allele": [["ref"], ["alt"], ["ref"]],
        },
        index=["oligo1", "oligo2", "oligo3"],
    )
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo2", "oligo3", "oligo3"]
    mpra_data.add_sequence_design(df, "dummy_path")
    assert "category" in mpra_data.data.var
    assert mpra_data._get_metadata("sequence_design_file") == "dummy_path"


def test_mpra_data_variant_map(mpra_data):
    # Setup sequence design
    df = pd.DataFrame(
        {
            "category": ["cat1", "cat2", "cat3"],
            "class": ["class1", "class2", "class3"],
            "ref": ["ref1", "ref2", "ref3"],
            "chr": ["chr1", "chr2", "chr3"],
            "start": [1, 2, 3],
            "end": [10, 20, 30],
            "strand": ["+", "-", "+"],
            "variant_class": ["vc1", "vc2", "vc3"],
            "variant_pos": [100, 200, 300],
            "SPDI": [["spdi1"], ["spdi2"], ["spdi3"]],
            "allele": [["ref"], ["alt"], ["ref"]],
        },
        index=["oligo1", "oligo2", "oligo3"],
    )
    mpra_data.data.var["oligo"] = ["oligo1", "oligo1", "oligo2", "oligo3", "oligo3"]
    mpra_data.add_sequence_design(df, "dummy_path")
    variant_map = mpra_data.variant_map
    assert isinstance(variant_map, pd.DataFrame)
    assert "REF" in variant_map.columns
    assert "ALT" in variant_map.columns


@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_mpra_data_drop_normalized(mpra_data, modality):
    mpra_data._normalize()
    assert f"{modality}_normalized" in mpra_data.data.layers
    mpra_data.drop_normalized()
    assert f"{modality}_normalized" not in mpra_data.data.layers
    assert mpra_data._get_metadata("normalized") is False


@pytest.mark.parametrize("modality", ["rna", "dna"])
def test_mpra_data_drop_correlation(mpra_data, modality):
    mpra_data._normalize()
    if modality == "rna":
        mpra_data._compute_correlation(mpra_data.normalized_rna_counts, "rna_normalized")
    else:
        mpra_data._compute_correlation(mpra_data.normalized_dna_counts, "dna_normalized")
    assert mpra_data._get_metadata(f"correlation_{modality}_normalized") is True
    mpra_data._drop_correlation()
    assert mpra_data._get_metadata(f"correlation_{modality}_normalized") is False


def test_mpra_data_correlation_invalid_method(mpra_data):
    mpra_data._normalize()
    with pytest.raises(ValueError):
        mpra_data.correlation(method="invalid", count_type=Modality.RNA_NORMALIZED)
    with pytest.raises(ValueError):
        mpra_data.correlation(method="pearson", count_type=Modality.RNA)


def test_drop_count_sampling_removes_sampling_layers_and_metadata(mpra_data):
    # Apply count sampling to create sampling layers and metadata
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=0.5)
    assert "rna_sampling" in mpra_data.data.layers
    assert "dna_sampling" in mpra_data.data.layers
    assert "count_sampling" in mpra_data.data.uns

    # Drop count sampling
    mpra_data.drop_count_sampling()

    # Sampling layers and metadata should be removed
    assert "rna_sampling" not in mpra_data.data.layers
    assert "dna_sampling" not in mpra_data.data.layers
    assert "count_sampling" not in mpra_data.data.uns


def test_drop_count_sampling_also_drops_normalized_and_barcode_counts(mpra_data):
    # Apply count sampling and normalization
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=0.5)
    mpra_data._normalize()
    mpra_data.barcode_counts = np.ones_like(mpra_data.raw_rna_counts)

    # Confirm layers exist
    assert "rna_normalized" in mpra_data.data.layers
    assert "dna_normalized" in mpra_data.data.layers
    assert "barcode_counts" in mpra_data.data.layers

    # Drop count sampling
    mpra_data.drop_count_sampling()

    # Normalized and barcode counts layers should be removed
    assert "rna_normalized" not in mpra_data.data.layers
    assert "dna_normalized" not in mpra_data.data.layers
    assert "barcode_counts" not in mpra_data.data.layers


def test_drop_count_sampling_removes_total_counts(mpra_data):
    # Apply count sampling to create total counts
    np.random.seed(42)
    mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=0.5)
    _ = mpra_data.total_rna_counts
    _ = mpra_data.total_dna_counts
    assert "rna_counts" in mpra_data.data.obs
    assert "dna_counts" in mpra_data.data.obs

    # Drop count sampling
    mpra_data.drop_count_sampling()

    # Total counts should be removed
    assert "rna_counts" not in mpra_data.data.obs
    assert "dna_counts" not in mpra_data.data.obs


def test_modality_enum():
    assert Modality.DNA.value == "dna"
    assert Modality.RNA.value == "rna"
    assert Modality.DNA_NORMALIZED.value == "dna_normalized"
    assert Modality.RNA_NORMALIZED.value == "rna_normalized"
    assert Modality.ACTIVITY.value == "activity"
    assert Modality.from_string("DNA") == Modality.DNA
    assert Modality.from_string("rna_normalized") == Modality.RNA_NORMALIZED
    with pytest.raises(ValueError):
        Modality.from_string("not_a_modality")


def test_countsampling_enum():
    assert CountSampling.RNA.value == "RNA"
    assert CountSampling.DNA.value == "DNA"
    assert CountSampling.RNA_AND_DNA.value == "RNA_AND_DNA"


def test_barcodefilter_enum():
    assert BarcodeFilter.MIN_BCS_PER_OLIGO.value == "min_bcs_per_oligo"
    assert BarcodeFilter.GLOBAL.value == "global"
    assert BarcodeFilter.OLIGO_SPECIFIC.value == "oligo_specific"
    assert BarcodeFilter.LARGE_EXPRESSION.value == "large_expression"
    assert BarcodeFilter.RANDOM.value == "random"
    assert BarcodeFilter.MIN_COUNT.value == "min_count"
    assert BarcodeFilter.MAX_COUNT.value == "max_count"
    assert BarcodeFilter.from_string("MIN_COUNT") == BarcodeFilter.MIN_COUNT
    with pytest.raises(ValueError):
        BarcodeFilter.from_string("not_a_filter")


def test_mpra_data_scaling_and_pseudocount(mpra_data):
    old_scaling = mpra_data.scaling
    mpra_data.scaling = 12345.0
    assert mpra_data.scaling == 12345.0
    mpra_data.scaling = old_scaling
    old_pseudocount = mpra_data.pseudo_count
    mpra_data.pseudo_count = 7
    assert mpra_data.pseudo_count == 7
    mpra_data.pseudo_count = old_pseudocount


def test_mpra_data_metadata(mpra_data):
    mpra_data._add_metadata("test_key", "test_value")
    assert mpra_data._get_metadata("test_key") == "test_value"


def test_mpra_data_var_filter_setter(mpra_data):
    mpra_data.var_filter = None
    assert np.all(mpra_data.var_filter == False)  # noqa: E712
    mask = np.ones_like(mpra_data.var_filter, dtype=bool)
    mpra_data.var_filter = mask
    assert np.all(mpra_data.var_filter == True)  # noqa: E712


def test_mpra_data_write_and_read(tmp_path, mpra_data):
    out_path = tmp_path / "test_write.h5ad"
    mpra_data.write(out_path)
    loaded = MPRABarcodeData.read(out_path)
    assert isinstance(loaded, MPRABarcodeData)
    assert loaded.data.shape == mpra_data.data.shape


def test_mpraoligodata_from_file(tmp_path, mpra_oligo_data):
    out_path = tmp_path / "oligo_data_test.h5ad"
    mpra_oligo_data.write(out_path)
    loaded = MPRAOligoData.from_file(out_path)
    assert isinstance(loaded, MPRAOligoData)
    assert loaded.data.shape == mpra_oligo_data.data.shape


def test_mpraoligodata_barcode_counts_exception(mpra_oligo_data):
    mpra_oligo_data.data.layers.pop("barcode_counts", None)
    with pytest.raises(MPRAlibException):
        _ = mpra_oligo_data.barcode_counts


def test_mpradata_normalize_layer(mpra_data):
    norm = mpra_data._normalize_layer(mpra_data.dna_counts, mpra_data.total_dna_counts)
    assert norm.shape == mpra_data.dna_counts.shape
    assert np.all(norm >= 0)


def test_mpradata_drop_normalized(mpra_data):
    mpra_data._normalize()
    assert "rna_normalized" in mpra_data.data.layers
    mpra_data.drop_normalized()
    assert "rna_normalized" not in mpra_data.data.layers
    assert "dna_normalized" not in mpra_data.data.layers


def test_mpradata_drop_total_counts(mpra_data):
    mpra_data.total_rna_counts
    mpra_data.total_dna_counts
    mpra_data.drop_total_counts()
    assert "rna_counts" not in mpra_data.data.obs
    assert "dna_counts" not in mpra_data.data.obs


def test_mprabarcode_complexity(mpra_data):
    result = mpra_data.complexity()
    assert result.shape[0] == mpra_data.n_obs
    assert result.shape[1] == mpra_data.n_obs


def test_mprabarcode_apply_count_sampling(mpra_data):
    mpra_data.apply_count_sampling(CountSampling.RNA, proportion=0.5)
    assert "rna_sampling" in mpra_data.data.layers
    mpra_data.apply_count_sampling(CountSampling.DNA, proportion=0.5)
    assert "dna_sampling" in mpra_data.data.layers


def test_mprabarcode_apply_barcode_filter(mpra_data):
    mpra_data.apply_barcode_filter(BarcodeFilter.RANDOM, params={"proportion": 0.5})
    assert mpra_data.var_filter.shape == (mpra_data.n_vars, mpra_data.n_obs)
    mpra_data.apply_barcode_filter(BarcodeFilter.MIN_COUNT, params={"rna_min_count": 1})
    assert mpra_data.var_filter.shape == (mpra_data.n_vars, mpra_data.n_obs)


def test_mprabarcode_drop_count_sampling(mpra_data):
    mpra_data.apply_count_sampling(CountSampling.RNA, proportion=0.5)
    mpra_data.drop_count_sampling()
    assert "rna_sampling" not in mpra_data.data.layers
    assert "dna_sampling" not in mpra_data.data.layers


def test_mprabarcode_drop_barcode_counts(mpra_data):
    mpra_data.barcode_counts = np.ones_like(mpra_data.raw_rna_counts)
    mpra_data.drop_barcode_counts()
    assert "barcode_counts" not in mpra_data.data.layers


def test_mpraoligodata_normalize_layer(mpra_oligo_data):
    counts = np.array([[1, 2, 3], [3, 4, 5], [3, 4, 5]], dtype=np.int32)
    total_counts = np.array([3, 7, 7], dtype=np.int32)
    mpra_oligo_data.barcode_counts = np.ones_like(counts)
    norm = mpra_oligo_data._normalize_layer(counts, total_counts)
    assert norm.shape == counts.shape
    assert np.all(norm >= 0)
