import unittest
import numpy as np
import pandas as pd
import anndata as ad
import copy
from mpralib.mpradata import MPRAdata, CountSampling


class TestMPRAdata(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing
        obs = pd.DataFrame(index=["rep1", "rep2", "rep3"])
        var = pd.DataFrame(
            {"oligo": ["oligo1", "oligo1", "oligo2", "oligo3", "oligo3"]},
            index=["barcode1", "barcode2", "barcode3", "barcode4", "barcode5"],
        )
        X = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
        layers = {"rna": X.copy(), "dna": X.copy()}
        self.mpra_data = MPRAdata(ad.AnnData(X=X, obs=obs, var=var, layers=layers))

        self.mpra_data_with_bc_filter = copy.deepcopy(self.mpra_data)

        self.mpra_data_with_bc_filter.barcode_filter = pd.DataFrame(
            np.array(
                [
                    [False, True, False],
                    [False, False, False],
                    [True, False, False],
                    [False, False, True],
                    [False, False, True],
                ]
            ),
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
        self.assertTrue(np.sum(rna_sampling) <= 29)

    def test_apply_count_sampling_total(self):
        np.random.seed(42)
        self.mpra_data.apply_count_sampling(CountSampling.RNA_AND_DNA, total=10)
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertTrue(np.sum(rna_sampling) <= 29)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.sum(dna_sampling) <= 29)

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
        self.assertTrue(np.sum(rna_sampling) <= 11)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.sum(dna_sampling) <= 11)

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


if __name__ == "__main__":
    unittest.main()
