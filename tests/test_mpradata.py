import unittest
import numpy as np
import pandas as pd
import anndata as ad
from mpralib.mpradata import MPRAdata, CountSampling


class TestMPRAdata(unittest.TestCase):

    def setUp(self):
        # Create a sample AnnData object for testing
        obs = pd.DataFrame(index=["rep1", "rep2", "rep3"])
        var = pd.DataFrame(
            index=["barcode1", "barcode2", "barcode3", "barcode4", "barcode5"]
        )
        X = np.array([[1, 2, 3, 1, 2], [4, 5, 6, 4, 5], [7, 8, 9, 10, 100]])
        layers = {"rna": X.copy(), "dna": X.copy()}
        self.adata = ad.AnnData(X=X, obs=obs, var=var, layers=layers)
        self.mpra_data = MPRAdata(self.adata)

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
        np.random.seed(4242)
        self.mpra_data.apply_count_sampling(
            CountSampling.RNA_AND_DNA, total=10, aggregate_over_replicates=True
        )
        rna_sampling = self.mpra_data.data.layers["rna_sampling"]
        self.assertTrue(np.sum(rna_sampling) <= 10)
        dna_sampling = self.mpra_data.data.layers["dna_sampling"]
        self.assertTrue(np.sum(dna_sampling) <= 11)


if __name__ == "__main__":
    unittest.main()
