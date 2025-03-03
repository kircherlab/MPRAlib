import unittest
import numpy as np
import pandas as pd
import anndata as ad
import copy
from mpralib.mpradata import MPRABarcodeData, CountSampling, BarcodeFilter


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


class TestMPRAdataSampling(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing

        layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
        self.mpra_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))

        self.mpra_data_with_bc_filter = copy.deepcopy(self.mpra_data)

        self.mpra_data_with_bc_filter.var_filter = pd.DataFrame(
            FILTER,
            index=self.mpra_data.data.var_names,
            columns=self.mpra_data.data.obs_names,
        )

    def test_apply_count_sampling_rna(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA, proportion=0.5)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertTrue(np.all(rna_sampling <= self.mpra_data.data.layers["rna"]))
        self.assertTrue(np.all(rna_sampling >= 0))

    def test_apply_count_sampling_dna(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.DNA, proportion=0.5)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.all(dna_sampling <= self.mpra_data.data.layers["dna"]))
        self.assertTrue(np.all(dna_sampling >= 0))

    def test_apply_count_sampling_rna_and_dna(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, proportion=0.5)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.all(rna_sampling <= self.mpra_data.data.layers["rna"]))
        self.assertTrue(np.all(rna_sampling >= 0))
        self.assertTrue(np.all(dna_sampling <= self.mpra_data.data.layers["dna"]))
        self.assertTrue(np.all(dna_sampling >= 0))

    def test_apply_count_sampling_rna_total(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA, total=10)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertLessEqual(np.sum(rna_sampling), 30)

    def test_apply_count_sampling_total(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, total=10)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertLessEqual(np.sum(rna_sampling), 30)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertLessEqual(np.sum(dna_sampling), 30)

    def test_apply_count_sampling_max_value_rna(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA, max_value=2)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertTrue(np.all(rna_sampling <= 2))

    def test_apply_count_sampling_max_value(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, max_value=2)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertTrue(np.all(rna_sampling <= 2))
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.all(dna_sampling <= 2))

    def test_apply_count_sampling_aggregate_over_replicates(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, total=10, aggregate_over_replicates=True)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertLessEqual(np.sum(rna_sampling), 11)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertLessEqual(np.sum(dna_sampling), 11)

    def test_compute_supporting_barcodes(self):
        # Manually set grouped_data for testing
        supporting_barcodes = self.mpra_data._compute_supporting_barcodes()

        expected_barcodes = np.array([[2, 1, 2], [2, 1, 2], [2, 1, 2]])

        np.testing.assert_array_equal(supporting_barcodes.to_numpy(), expected_barcodes)

    def test_compute_supporting_barcodes_with_filter(self):

        supporting_barcodes = self.mpra_data_with_bc_filter._compute_supporting_barcodes()

        expected_barcodes = np.array([[2, 0, 2], [1, 1, 2], [2, 1, 0]])

        np.testing.assert_array_equal(supporting_barcodes.to_numpy(), expected_barcodes)

    def test_raw_dna_counts(self):
        expected_dna_counts = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
        np.testing.assert_array_equal(self.mpra_data.raw_dna_counts, expected_dna_counts)

    def test_raw_dna_counts_with_modification(self):
        self.mpra_data.data.layers["dna"] = np.array([[10, 20, 30, 10, 20], [40, 50, 60, 40, 50], [70, 80, 90, 100, 1000]])
        expected_dna_counts = np.array([[10, 20, 30, 10, 20], [40, 50, 60, 40, 50], [70, 80, 90, 100, 1000]])
        np.testing.assert_array_equal(self.mpra_data.raw_dna_counts, expected_dna_counts)

    def test_filtered_dna_counts(self):
        # Test without barcode filter
        expected_dna_counts = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
        np.testing.assert_array_equal(self.mpra_data.dna_counts, expected_dna_counts)

        # Test with barcode filter
        expected_filtered_dna_counts = np.array([[1, 2, 0, 1, 2], [0, 5, 6, 4, 5], [7, 8, 9, 0, 0]])
        np.testing.assert_array_equal(self.mpra_data_with_bc_filter.dna_counts, expected_filtered_dna_counts)

    def test_dna_counts_with_sampling(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.DNA, proportion=0.5)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        np.testing.assert_array_equal(self.mpra_data.dna_counts, dna_sampling)

    def test_dna_counts_with_filter(self):
        self.mpra_data.apply_count_sampling(CountSampling.DNA, max_value=2)
        expected_filtered_dna_counts = np.array([[1, 2, 2, 1, 2], [2, 2, 2, 2, 2], [2, 2, 2, 2, 2]])
        np.testing.assert_array_equal(self.mpra_data.dna_counts, expected_filtered_dna_counts)
        self.mpra_data_with_bc_filter.apply_count_sampling(CountSampling.DNA, max_value=2)
        expected_filtered_dna_counts = np.array([[1, 2, 0, 1, 2], [0, 2, 2, 2, 2], [2, 2, 2, 0, 0]])
        np.testing.assert_array_equal(self.mpra_data_with_bc_filter.dna_counts, expected_filtered_dna_counts)


class TestMPRAdataBarcodeFilter(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing
        layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
        self.mpra_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))

    def test_apply_barcode_filter_min_count(self):
        self.mpra_data.apply_barcode_filter(BarcodeFilter.MIN_COUNT, params={"rna_min_count": 4, "dna_min_count": 3})
        expected_filter = np.array(
            [
                [True, False, False],
                [True, False, False],
                [False, False, False],
                [True, False, False],
                [True, False, False],
            ]
        )
        np.testing.assert_array_equal(self.mpra_data.var_filter.values, expected_filter)

    def test_apply_barcode_filter_max_count(self):
        self.mpra_data.apply_barcode_filter(BarcodeFilter.MAX_COUNT, params={"rna_max_count": 9, "dna_max_count": 100})
        expected_filter = np.array(
            [
                [False, False, False],
                [False, False, False],
                [False, False, False],
                [False, False, True],
                [False, False, True],
            ]
        )
        np.testing.assert_array_equal(self.mpra_data.var_filter.values, expected_filter)

        self.mpra_data.var_filter = None

        self.mpra_data.apply_barcode_filter(BarcodeFilter.MAX_COUNT, params={"dna_max_count": 99})
        expected_filter = np.array(
            [
                [False, False, False],
                [False, False, False],
                [False, False, False],
                [False, False, False],
                [False, False, True],
            ]
        )
        np.testing.assert_array_equal(self.mpra_data.var_filter.values, expected_filter)


class TestMPRABarcodeDataNormalization(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing
        layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
        self.mpra_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))

        self.mpra_data.SCALING = 10

        self.mpra_data_with_bc_filter = copy.deepcopy(self.mpra_data)

        self.mpra_data_with_bc_filter.var_filter = pd.DataFrame(
            FILTER,
            index=self.mpra_data.data.var_names,
            columns=self.mpra_data.data.obs_names,
        )

    def test_normalize_counts(self):
        self.mpra_data._normalize()
        dna_normalized = self.mpra_data.normalized_dna_counts
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
        rna_normalized = self.mpra_data.normalized_rna_counts
        np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)

    def test_normalize_without_pseudocount(self):
        mpra_data = copy.deepcopy(self.mpra_data)
        mpra_data.PSEUDOCOUNT = 0
        mpra_data._normalize()
        dna_normalized = mpra_data.data.layers["dna_normalized"]
        expected_dna_normalized = np.array(
            [[1.111, 2.222, 3.333, 1.111, 2.222], [1.667, 2.083, 2.5, 1.667, 2.083], [0.522, 0.597, 0.672, 0.746, 7.463]]
        )
        np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)
        expected_rna_normalized = np.array(
            [[1.0, 2.0, 4.0, 1.0, 2.0], [1.667, 2.083, 2.5, 1.667, 2.083], [0.522, 0.597, 0.672, 0.746, 7.463]]
        )
        rna_normalized = mpra_data.data.layers["rna_normalized"]
        np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)

    def test_normalize_counts_with_bc_filter(self):
        self.mpra_data_with_bc_filter._normalize()
        dna_normalized = self.mpra_data_with_bc_filter.data.layers["dna_normalized"]
        expected_normalized = np.array(
            [[2.0, 3.0, 0.0, 2.0, 3.0], [0.0, 2.5, 2.917, 2.083, 2.5], [2.963, 3.333, 3.704, 0.0, 0.0]]
        )
        np.testing.assert_almost_equal(dna_normalized, expected_normalized, decimal=3)

        rna_normalized = self.mpra_data_with_bc_filter.data.layers["rna_normalized"]
        np.testing.assert_almost_equal(rna_normalized, expected_normalized, decimal=3)


class TestMPRAOligoDataNormalization(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing
        layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
        mpra_barcode_data = MPRABarcodeData(ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers))

        self.mpra_data = mpra_barcode_data.oligo_data
        self.mpra_data.SCALING = 10

        mpra_barcode_data = copy.deepcopy(mpra_barcode_data)

        mpra_barcode_data.var_filter = pd.DataFrame(
            FILTER,
            index=mpra_barcode_data.data.var_names,
            columns=mpra_barcode_data.data.obs_names,
        )

        self.mpra_data_with_bc_filter = mpra_barcode_data.oligo_data
        self.mpra_data_with_bc_filter.SCALING = 10

    def test_normalize_counts(self):
        dna_normalized = self.mpra_data.normalized_dna_counts
        expected_dna_normalized = np.array([[1.786, 2.857, 1.786], [1.897, 2.414, 1.897], [0.612, 0.719, 4.029]])
        np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)

        expected_rna_normalized = np.array([[1.667, 3.333, 1.667], [1.897, 2.414, 1.897], [0.612, 0.719, 4.029]])
        rna_normalized = self.mpra_data.normalized_rna_counts
        np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)

    def test_normalize_without_pseudocount(self):
        mpra_data = copy.deepcopy(self.mpra_data)
        mpra_data.PSEUDOCOUNT = 0

        dna_normalized = mpra_data.normalized_dna_counts
        expected_dna_normalized = np.array([[1.667, 3.333, 1.667], [1.875, 2.5, 1.875], [0.56, 0.672, 4.104]])
        np.testing.assert_almost_equal(dna_normalized, expected_dna_normalized, decimal=3)

        expected_rna_normalized = np.array([[1.5, 4.0, 1.5], [1.875, 2.5, 1.875], [0.56, 0.672, 4.104]])
        rna_normalized = mpra_data.normalized_rna_counts
        np.testing.assert_almost_equal(rna_normalized, expected_rna_normalized, decimal=3)

    def test_normalize_counts_with_bc_filter(self):

        dna_normalized = self.mpra_data_with_bc_filter.normalized_dna_counts
        expected_normalized = np.array([[2.5, np.nan, 2.5], [2.5, 2.917, 2.292], [3.148, 3.704, np.nan]])
        np.testing.assert_almost_equal(dna_normalized, expected_normalized, decimal=3)

        rna_normalized = self.mpra_data_with_bc_filter.normalized_rna_counts
        np.testing.assert_almost_equal(rna_normalized, expected_normalized, decimal=3)


class TestMPRAdataCorrelation(unittest.TestCase):

    def setUp(self):
        layers = {"rna": COUNTS_RNA.copy(), "dna": COUNTS_DNA.copy()}
        self.mpra_data = MPRABarcodeData(
            ad.AnnData(X=COUNTS_RNA.copy(), obs=OBS.copy(), var=VAR.copy(), layers=layers)
        ).oligo_data

    def test_correlation(self):
        self.mpra_data._compute_correlation(self.mpra_data.activity, "log2FoldChange")
        self.assertIn("pearson_correlation_log2FoldChange", self.mpra_data.data.obsp)
        self.assertIn("spearman_correlation_log2FoldChange", self.mpra_data.data.obsp)

    def test_pearson_correlation(self):
        x = self.mpra_data.pearson_correlation_activity
        y = self.mpra_data.pearson_correlation_rna
        z = self.mpra_data.pearson_correlation_dna
        np.testing.assert_equal(x, np.array([[1.0, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]]))
        np.testing.assert_almost_equal(
            y, np.array([[1.0, 1.0, -0.475752], [1.0, 1.0, -0.475752], [-0.475752, -0.475752, 1.0]]), decimal=3
        )
        np.testing.assert_almost_equal(z, np.array([[1.0, 1.0, -0.476], [1.0, 1.0, -0.476], [-0.476, -0.476, 1.0]]), decimal=3)

    def test_spearman_correlation(self):
        x = self.mpra_data.spearman_correlation_activity
        y = self.mpra_data.spearman_correlation_rna
        z = self.mpra_data.spearman_correlation_dna
        np.testing.assert_equal(x, np.array([[1.0, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]]))
        np.testing.assert_almost_equal(y, np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), decimal=3)
        np.testing.assert_equal(z, np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))


if __name__ == "__main__":
    unittest.main()
