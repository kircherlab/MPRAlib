import numpy as np
import ast
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, pearsonr
from enum import Enum
import logging


class CountSampling(Enum):
    RNA = "RNA"
    DNA = "DNA"
    RNA_AND_DNA = "RNA_AND_DNA"


class BarcodeFilter(Enum):
    RNA_ZSCORE = "RNA_ZSCORE"
    MAD = "MAD"
    RANDOM = "RANDOM"
    MIN_COUNT = "MIN_COUNT"
    MAX_COUNT = "MAX_COUNT"


class MPRAdata:

    LOGGER = logging.getLogger(__name__)

    SCALING = 1e6
    PSEUDOCOUNT = 1

    def __init__(self, data: ad.AnnData):
        self._data = data
        self._grouped_data = None
        self.barcode_filter = None

    @property
    def data(self) -> ad.AnnData:
        return self._data

    @data.setter
    def data(self, new_data: ad.AnnData) -> None:
        self._data = new_data

    @property
    def grouped_data(self):
        if not self._grouped_data:
            self._group_data()
        return self._grouped_data

    @grouped_data.setter
    def grouped_data(self, new_data: ad.AnnData) -> None:
        if new_data is None and self._grouped_data is not None:
            self.LOGGER.info("Dropping grouped data")
        self._grouped_data = new_data

    @property
    def normalized_dna_counts(self):
        if "dna_normalized" not in self.data.layers:
            self._normalize()
        return self.data.layers["dna_normalized"]

    @property
    def raw_rna_counts(self):
        return self.data.layers["rna"]

    @property
    def rna_counts(self):
        if "rna_sampling" in self.data.layers:
            return self.data.layers["rna_sampling"] * ~self.barcode_filter.T.values
        else:
            return self.raw_rna_counts * ~self.barcode_filter.T.values

    @property
    def raw_dna_counts(self):
        return self.data.layers["dna"]

    @property
    def dna_counts(self):
        if "dna_sampling" in self.data.layers:
            return self.data.layers["dna_sampling"] * ~self.barcode_filter.T.values
        else:
            return self.raw_dna_counts * ~self.barcode_filter.T.values

    @property
    def normalized_rna_counts(self):
        if "rna_normalized" not in self.data.layers:
            self._normalize()
        return self.data.layers["rna_normalized"]

    @property
    def barcodes(self):
        return self.data.var_names

    @property
    def oligos(self):
        return self.data.var["oligo"]

    @property
    def replicates(self):
        return self.data.obs_names

    @property
    def n_replicates(self):
        return self.data.n_obs

    @property
    def n_raw_barcodes(self):
        return self.data.n_vars

    @property
    def barcode_threshold(self) -> int:
        return self.data.uns["barcode_threshold"]

    @barcode_threshold.setter
    def barcode_threshold(self, barcode_threshold: int) -> None:
        self._add_metadata("barcode_threshold", barcode_threshold)
        self.grouped_data = None

    @property
    def barcode_filter(self) -> pd.DataFrame:
        return self.data.varm["barcode_filter"]

    @barcode_filter.setter
    def barcode_filter(self, new_data: pd.DataFrame) -> None:
        if new_data is None:
            self.data.varm["barcode_filter"] = pd.DataFrame(
                np.full((self.data.n_vars, self.data.n_obs), False),
                index=self.barcodes,
                columns=self.replicates,
            )
            if "barcode_filter" in self.data.uns:
                del self.data.uns["barcode_filter"]
        else:
            self.data.varm["barcode_filter"] = new_data

        self.grouped_data = None
        self.drop_normalized()

    def drop_normalized(self):

        self.grouped_data = None

        self.LOGGER.info("Dropping normalized data")

        self.data.layers.pop("rna_normalized", None)
        self.data.layers.pop("dna_normalized", None)
        self.data.uns["normalized"] = False

    @property
    def spearman_correlation(self) -> np.ndarray:
        """
        Calculate and return the Spearman correlation for the grouped data.

        This property checks if the Spearman correlation has already been computed
        and stored in the `grouped_data` attribute. If not, it calls the `_correlation`
        method to compute it. The computed Spearman correlation is then retrieved
        from the `grouped_data.obsp` attribute.

        Returns:
            np.ndarray: The Spearman correlation matrix for the grouped data.
        """
        if not self.grouped_data.uns["correlation"]:
            self._correlation()
        return self.grouped_data.obsp["spearman_correlation"]

    @property
    def pearson_correlation(self) -> np.ndarray:
        """
        Computes and returns the Pearson correlation matrix.

        This property checks if the Pearson correlation matrix is already computed
        and stored in the `grouped_data` attribute. If not, it calls the `_correlation`
        method to compute it. The Pearson correlation matrix is then retrieved from
        the `grouped_data` attribute.

        Returns:
            np.ndarray: The Pearson correlation matrix for the grouped data..
        """
        if not self.grouped_data.uns["correlation"]:
            self._correlation()
        return self.grouped_data.obsp["pearson_correlation"]

    @property
    def variant_map(self) -> pd.DataFrame:
        """
        Generates a DataFrame mapping variant IDs to their reference and alternate alleles.
        Returns:
            pd.DataFrame: A DataFrame with columns 'REF' and 'ALT', indexed by variant IDs ('ID').
                          'REF' contains lists of reference alleles, and 'ALT' contains lists of alternate alleles.
        Raises:
            ValueError: If the metadata file is not loaded in `self.data.uns`.
        """
        if "metadata_file" not in self.data.uns:
            raise ValueError("Metadata file not loaded.")

        spdis = set([item for sublist in self.grouped_data.var["SPDI"].values for item in sublist])

        df = {"ID": [], "REF": [], "ALT": []}
        for spdi in spdis:
            df["ID"].append(spdi)
            spdi_data = self.grouped_data[:, self.grouped_data.var["SPDI"].apply(lambda x: spdi in x)]
            spdi_idx = spdi_data.var["SPDI"].apply(lambda x: x.index(spdi))
            refs = []
            alts = []
            for idx, value in spdi_data.var["allele"].items():
                if "ref" == value[spdi_idx[idx]]:
                    refs.append(self.oligos[idx])
                else:
                    alts.append(self.oligos[idx])
            df["REF"].append(refs)
            df["ALT"].append(alts)

        return pd.DataFrame(df, index=df["ID"]).set_index("ID")

    @property
    def element_dna_counts(self) -> pd.DataFrame:
        return self._element_counts("dna")

    @property
    def element_rna_counts(self) -> pd.DataFrame:
        return self._element_counts("rna")

    def _element_counts(self, layer: str) -> pd.DataFrame:
        df = {"ID": self.grouped_data.var["oligo"]}
        for replicate in self.grouped_data.obs_names:
            df["counts_" + replicate] = self.grouped_data[replicate, :].layers[layer][0]

        df = pd.DataFrame(df).set_index("ID")
        return df[(df.T != 0).all()]

    @property
    def variant_dna_counts(self) -> pd.DataFrame:
        return self._variant_counts("dna")

    @property
    def variant_rna_counts(self) -> pd.DataFrame:
        return self._variant_counts("rna")

    def _variant_counts(self, layer: str) -> pd.DataFrame:
        df = {"ID": []}
        for replicate in self.grouped_data.obs_names:
            df["counts_" + replicate + "_REF"] = []
            df["counts_" + replicate + "_ALT"] = []

        for spdi, row in self.variant_map.iterrows():
            df["ID"].append(spdi)
            counts_ref = self.grouped_data.layers[layer][:, self.grouped_data.var["oligo"].isin(row["REF"])].sum(axis=1)
            counts_alt = self.grouped_data.layers[layer][:, self.grouped_data.var["oligo"].isin(row["ALT"])].sum(axis=1)
            idx = 0
            for replicate in self.grouped_data.obs_names:
                df["counts_" + replicate + "_REF"].append(counts_ref[idx])
                df["counts_" + replicate + "_ALT"].append(counts_alt[idx])
                idx += 1
        df = pd.DataFrame(df).set_index("ID")
        return df[(df.T != 0).all()]

    def add_metadata_file(self, metadata_file):
        metadata = pd.read_csv(metadata_file, sep="\t", header=0, na_values=["NA"]).drop_duplicates()
        metadata["name"] = pd.Categorical(metadata["name"].str.replace(" ", "_"))
        metadata.set_index("name", inplace=True)
        metadata["variant_class"] = metadata["variant_class"].fillna("[]").apply(ast.literal_eval)
        metadata["variant_pos"] = metadata["variant_pos"].fillna("[]").apply(ast.literal_eval)
        metadata["SPDI"] = metadata["SPDI"].fillna("[]").apply(ast.literal_eval)
        metadata["allele"] = metadata["allele"].fillna("[]").apply(ast.literal_eval)

        matched_metadata = metadata.loc[self.oligos]

        self.data.var["category"] = pd.Categorical(matched_metadata["category"])
        self.data.var["class"] = pd.Categorical(matched_metadata["class"])
        self.data.var["ref"] = pd.Categorical(matched_metadata["ref"])
        self.data.var["chr"] = pd.Categorical(matched_metadata["chr"])
        self.data.var["start"] = matched_metadata["start"].values
        self.data.var["end"] = matched_metadata["end"].values
        self.data.var["strand"] = pd.Categorical(matched_metadata["strand"])
        self.data.var["variant_class"] = matched_metadata["variant_class"].values
        self.data.var["variant_pos"] = matched_metadata["variant_pos"].values
        self.data.var["SPDI"] = matched_metadata["SPDI"].values
        self.data.var["allele"] = matched_metadata["allele"].values

        self._add_metadata("metadata_file", metadata_file)

        # need to reset grouped data after adding metadata
        self.grouped_data = None

    def write(self, file_data_path: str, file_grouped_path: str = None) -> None:
        self.data.write(file_data_path)

        if file_grouped_path and self.grouped_data:
            self.grouped_data.write(file_grouped_path)

    @classmethod
    def from_file(cls, file_path: str) -> "MPRAdata":
        """
        Create an instance of the class from a file.
        This method reads data from a specified file, processes it, and returns an instance of the class containing the
        data in an AnnData object.
        Parameters
        ----------
        file_path : str
            The path to the file containing the data.
        Returns
        -------
        cls
            An instance of the class containing the processed data.
        Notes
        -----
        The input file is expected to be a tab-separated values (TSV) file with a header and an index column.
        The method assumes that the RNA and DNA replicate columns are interleaved, starting with DNA.
        The method processes the data to create an AnnData object with RNA and DNA layers, and additional metadata.
        """
        data = pd.read_csv(file_path, sep="\t", header=0, index_col=0)
        data = data.fillna(0)

        replicate_columns_rna = data.columns[2::2]
        replicate_columns_dna = data.columns[1::2]

        anndata_replicate_rna = data.loc[:, replicate_columns_rna].transpose().astype(int)
        anndata_replicate_dna = data.loc[:, replicate_columns_dna].transpose().astype(int)

        anndata_replicate_rna.index = [replicate.split("_")[2] for replicate in replicate_columns_rna]
        anndata_replicate_dna.index = [replicate.split("_")[2] for replicate in replicate_columns_dna]

        adata = ad.AnnData(anndata_replicate_rna)
        adata.layers["rna"] = adata.X.copy()
        adata.layers["dna"] = anndata_replicate_dna

        adata.var["oligo"] = pd.Categorical(data["oligo_name"].values)

        adata.uns["file_path"] = file_path
        adata.uns["date"] = pd.to_datetime("today").strftime("%Y-%m-%d")

        adata.uns["normalized"] = False
        adata.uns["barcode_threshold"] = None

        adata.varm["barcode_filter"] = pd.DataFrame(
            np.full((adata.n_vars, adata.n_obs), False),
            index=adata.var_names,
            columns=adata.obs_names,
        )

        return cls(adata)

    def _normalize(self):

        self.LOGGER.info("Normalizing data")

        self.data.layers["dna_normalized"] = self._normalize_layer(
            self.dna_counts, ~self.barcode_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self.data.layers["rna_normalized"] = self._normalize_layer(
            self.rna_counts, ~self.barcode_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self._add_metadata("normalized", True)

    def _normalize_layer(
        self, counts: np.ndarray, not_barcode_filter: np.ndarray, scaling: float, pseudocount: int
    ) -> np.ndarray:

        # I do a pseudo count when normalizing to avoid division by zero when computing logfold ratios.
        # barcode filter has to be used again because we want to have a zero on filtered values.
        total_counts = np.sum((counts + pseudocount) * not_barcode_filter, axis=1)

        # Avoid division by zero when pseudocount is set to 0
        total_counts[total_counts == 0] = 1
        return ((counts + pseudocount) / total_counts[:, np.newaxis] * scaling) * not_barcode_filter

    def _group_data(self):

        if not self.data.uns["normalized"]:
            self._normalize()

        # Convert the result back to an AnnData object
        self.grouped_data = ad.AnnData(self._sum_counts_by_oligo(self.rna_counts))

        self.grouped_data.layers["rna"] = self.grouped_data.X.copy()
        self.grouped_data.layers["dna"] = self._sum_counts_by_oligo(self.dna_counts)

        self.grouped_data.layers["barcodes"] = self._compute_supporting_barcodes()

        self.grouped_data.layers["rna_normalized"] = (
            self._sum_counts_by_oligo(self.normalized_rna_counts) / self.grouped_data.layers["barcodes"]
        )
        self.grouped_data.layers["dna_normalized"] = (
            self._sum_counts_by_oligo(self.normalized_dna_counts) / self.grouped_data.layers["barcodes"]
        )

        self.grouped_data.obs_names = self.replicates
        self.grouped_data.var_names = self.data.var["oligo"].unique().tolist()

        # Subset of vars using the firs occurence of oligo name
        indices = self.data.var["oligo"].drop_duplicates(keep="first").index
        self.grouped_data.var = self.data.var.loc[indices]

        if self.barcode_threshold is not None:
            mask = self.grouped_data.layers["barcodes"] < self.barcode_threshold
            self.grouped_data.X[mask] = 0
            for layer_name in self.grouped_data.layers.keys():
                if self.grouped_data.layers[layer_name].dtype.kind in {"i", "u"}:
                    self.grouped_data.layers[layer_name][mask] = 0
                else:
                    self.grouped_data.layers[layer_name][mask] = np.nan
        self.grouped_data.obs = self.data.obs

        self.grouped_data.uns = self.data.uns

        self.grouped_data.uns["correlation"] = False

        self._compute_activities()

    def _sum_counts_by_oligo(self, counts):

        grouped = pd.DataFrame(
            counts,
            index=self.replicates,
            columns=self.barcodes,
        ).T.groupby(self.data.var["oligo"], observed=True)

        # Perform an operation on each group, e.g., mean
        return grouped.sum().T

    def _compute_supporting_barcodes(self):

        grouped = pd.DataFrame(
            (self.raw_rna_counts + self.raw_dna_counts) * ~self.barcode_filter.T.values,
            index=self.replicates,
            columns=self.barcodes,
        ).T.groupby(self.oligos, observed=True)

        return grouped.apply(lambda x: (x != 0).sum()).T

    def _compute_activities(self):
        self.grouped_data.layers["log2FoldChange"] = np.log2(
            self.grouped_data.layers["rna_normalized"] / self.grouped_data.layers["dna_normalized"]
        )

    def _barcode_filter_rna_zscore(self, times_zscore=3):

        barcode_mask = pd.DataFrame(
            self.raw_dna_counts + self.raw_rna_counts,
            index=self.replicates,
            columns=self.barcodes,
        ).T.apply(lambda x: (x != 0))

        df_rna = pd.DataFrame(self.raw_rna_counts, index=self.replicates, columns=self.barcodes).T
        grouped = df_rna.where(barcode_mask).groupby(self.oligos, observed=True)

        mask = ((df_rna - grouped.transform("mean")) / grouped.transform("std")).abs() > times_zscore

        return mask

    def _barcode_filter_mad(self, times_mad=3, n_bins=20):

        # sum up DNA and RNA counts across replicates
        DNA_sum = pd.DataFrame(self.raw_dna_counts, index=self.replicates, columns=self.barcodes).T.sum(axis=1)
        RNA_sum = pd.DataFrame(self.raw_rna_counts, index=self.replicates, columns=self.barcodes).T.sum(axis=1)
        df_sums = pd.DataFrame({"DNA_sum": DNA_sum, "RNA_sum": RNA_sum}).fillna(0)
        # removing all barcodes with 0 counts in RNA an more DNA count than number of replicates/observations
        df_sums = df_sums[(df_sums["DNA_sum"] > self.data.n_obs) & (df_sums["RNA_sum"] > 0)]

        # remove all barcodes where oligo has less barcodes as the number of replicates/observations
        df_sums = df_sums.groupby(self.data.var["oligo"], observed=True).filter(lambda x: len(x) >= self.data.n_obs)

        # Calculate ratio, ratio_med, ratio_diff, and mad
        df_sums["ratio"] = np.log2(df_sums["DNA_sum"] / df_sums["RNA_sum"])
        df_sums["ratio_med"] = df_sums.groupby(self.data.var["oligo"], observed=True)["ratio"].transform("median")
        df_sums["ratio_diff"] = df_sums["ratio"] - df_sums["ratio_med"]
        # df_sums['mad'] = (df_sums['ratio'] - df_sums['ratio_med']).abs().mean()

        # Calculate quantiles within  n_bins
        qs = np.unique(np.quantile(np.log10(df_sums["RNA_sum"]), np.arange(0, n_bins) / n_bins))

        # Create bins based on rna_count
        df_sums["bin"] = pd.cut(
            np.log10(df_sums["RNA_sum"]),
            bins=qs,
            include_lowest=True,
            labels=[str(i) for i in range(0, len(qs) - 1)],
        )
        # Filter based on ratio_diff and mad
        df_sums["ratio_diff_med"] = df_sums.groupby("bin", observed=True)["ratio_diff"].transform("median")
        df_sums["ratio_diff_med_dist"] = np.abs(df_sums["ratio_diff"] - df_sums["ratio_diff_med"])
        df_sums["mad"] = df_sums.groupby("bin", observed=True)["ratio_diff_med_dist"].transform("median")

        m = df_sums.ratio_diff > times_mad * df_sums.mad
        df_sums = df_sums[~m]

        return self.barcode_filter.apply(lambda col: col | ~self.barcode_filter.index.isin(df_sums.index))

    def _barcode_filter_random(self, proportion=1.0, total=None, aggegate_over_replicates=True):

        if aggegate_over_replicates and total is None:
            total = self.barcode_filter.shape[0]
        elif total is None:
            total = self.barcode_filter.size

        num_true_cells = int(total * (1.0 - proportion))
        true_indices = np.random.choice(total, num_true_cells, replace=False)

        mask = pd.DataFrame(
            np.full((self.data.n_vars, self.data.n_obs), False),
            index=self.barcodes,
            columns=self.replicates,
        )

        if aggegate_over_replicates:
            mask.iloc[true_indices, :] = True
        else:
            flat_df = mask.values.flatten()

            flat_df[true_indices] = True

            mask = pd.DataFrame(
                flat_df.reshape(self.barcode_filter.shape),
                index=self.barcodes,
                columns=self.replicates,
            )
        return mask

    def _barcode_filter_min_count(self, rna_min_count=None, dna_min_count=None):

        return self._barcode_filter_min_max_count(BarcodeFilter.MIN_COUNT, rna_min_count, dna_min_count)

    def _barcode_filter_max_count(self, rna_max_count=None, dna_max_count=None):

        return self._barcode_filter_min_max_count(BarcodeFilter.MAX_COUNT, rna_max_count, dna_max_count)

    def _barcode_filter_min_max_counts(self, barcode_filter, counts, count_threshold):
        if barcode_filter == BarcodeFilter.MIN_COUNT:
            return (counts < count_threshold).T
        elif barcode_filter == BarcodeFilter.MAX_COUNT:
            return (counts > count_threshold).T

    def _barcode_filter_min_max_count(self, barcode_filter, rna_count=None, dna_count=None):
        mask = pd.DataFrame(
            np.full((self.n_raw_barcodes, self.n_replicates), False),
            index=self.barcodes,
            columns=self.replicates,
        )
        if rna_count is not None:
            mask = mask | self._barcode_filter_min_max_counts(barcode_filter, self.raw_rna_counts, rna_count)
        if dna_count is not None:
            mask = mask | self._barcode_filter_min_max_counts(barcode_filter, self.raw_dna_counts, dna_count)

        return mask

    def apply_barcode_filter(self, barcode_filter: BarcodeFilter, params: dict = {}):
        filter_switch = {
            BarcodeFilter.RNA_ZSCORE: self._barcode_filter_rna_zscore,
            BarcodeFilter.MAD: self._barcode_filter_mad,
            BarcodeFilter.RANDOM: self._barcode_filter_random,
            BarcodeFilter.MIN_COUNT: self._barcode_filter_min_count,
            BarcodeFilter.MAX_COUNT: self._barcode_filter_max_count,
        }

        filter_func = filter_switch.get(barcode_filter)
        if filter_func:
            self.barcode_filter = self.barcode_filter | filter_func(**params)
        else:
            raise ValueError(f"Unsupported barcode filter: {barcode_filter}")

        self._add_metadata("barcode_filter", [barcode_filter.value])

    def drop_count_sampling(self):

        self.drop_normalized()

        self.LOGGER.info("Dropping count sampling")

        del self.data.uns["count_sampling"]
        self.data.layers.pop("rna_sampling", None)
        self.data.layers.pop("dna_sampling", None)

    def _calculate_proportions(self, proportion, total, aggregate_over_replicates, counts, replicates):
        pp = [1.0] * replicates

        if proportion is not None:
            pp = [proportion] * replicates

        if total is not None:
            if aggregate_over_replicates:
                for i, p in enumerate(pp):
                    pp[i] = min(total / np.sum(counts), p)
            else:
                for i, p in enumerate(pp):
                    pp[i] = min(total / np.sum(counts[i, :]), p)
        return pp

    def _sample_individual_counts(self, x, proportion):
        return int(
            np.floor(x * proportion)
            + (0.0 if x != 0 and np.random.rand() > (x * proportion - np.floor(x * proportion)) else 1.0)
        )

    def _apply_sampling(self, layer_name, counts, proportion, total, max_value, aggregate_over_replicates):
        self.data.layers[layer_name] = counts.copy()

        if total is not None or proportion is not None:

            pp = self._calculate_proportions(
                proportion, total, aggregate_over_replicates, self.data.layers[layer_name], self.n_replicates
            )

            vectorized_sample_individual_counts = np.vectorize(self._sample_individual_counts)

            for i, p in enumerate(pp):
                self.data.layers[layer_name][i, :] = vectorized_sample_individual_counts(
                    self.data.layers[layer_name][i, :], proportion=p
                )

        if max_value is not None:
            self.data.layers[layer_name] = np.clip(self.data.layers[layer_name], None, max_value)

    def apply_count_sampling(
        self,
        count_type: CountSampling,
        proportion: float = None,
        total: int = None,
        max_value: int = None,
        aggregate_over_replicates: bool = False,
    ) -> None:

        if count_type == CountSampling.RNA or count_type == CountSampling.RNA_AND_DNA:
            self._apply_sampling("rna_sampling", self.rna_counts, proportion, total, max_value, aggregate_over_replicates)

        if count_type == CountSampling.DNA or count_type == CountSampling.RNA_AND_DNA:
            self._apply_sampling("dna_sampling", self.dna_counts, proportion, total, max_value, aggregate_over_replicates)

        self._add_metadata(
            "count_sampling",
            [
                {
                    count_type.value: {
                        "proportion": proportion,
                        "total": total,
                        "max_value": max_value,
                        "aggregate_over_replicates": aggregate_over_replicates,
                    }
                }
            ],
        )

        self.grouped_data = None

        self.drop_normalized()

    def _add_metadata(self, key, value):
        if isinstance(value, list):
            if key not in self.data.uns:
                self.data.uns[key] = value
            else:
                self.data.uns[key] = self.data.uns[key] + value
        else:
            self.data.uns[key] = value

    def _correlation(self):
        """
        Compute Pearson and Spearman correlation matrices for the log2FoldChange layer of grouped_data.
        This method calculates both Pearson and Spearman correlation coefficients and their corresponding p-values
        for each pair of rows in the log2FoldChange layer of the grouped_data attribute. The results are stored in
        the obsp attribute of grouped_data with the following keys:
        - "pearson_correlation": Pearson correlation coefficients matrix.
        - "pearson_correlation_pvalue": p-values for the Pearson correlation coefficients.
        - "spearman_correlation": Spearman correlation coefficients matrix.
        - "spearman_correlation_pvalue": p-values for the Spearman correlation coefficients.
        Additionally, a flag indicating that the correlation has been computed is stored in the uns attribute of
        grouped_data with the key "correlation".
        Returns:
            None
        """
        num_columns = self.grouped_data.shape[0]
        self.grouped_data.obsp["pearson_correlation"] = np.zeros((num_columns, num_columns))
        self.grouped_data.obsp["pearson_correlation_pvalue"] = np.zeros((num_columns, num_columns))
        self.grouped_data.obsp["spearman_correlation"] = np.zeros((num_columns, num_columns))
        self.grouped_data.obsp["spearman_correlation_pvalue"] = np.zeros((num_columns, num_columns))

        for i in range(num_columns):
            for j in range(i, num_columns):
                mask = ~np.isnan(self.grouped_data.layers["log2FoldChange"][i, :]) & ~np.isnan(
                    self.grouped_data.layers["log2FoldChange"][j, :]
                )
                x = self.grouped_data.layers["log2FoldChange"][i, mask]
                y = self.grouped_data.layers["log2FoldChange"][j, mask]

                # Spearman correlation
                (
                    self.grouped_data.obsp["spearman_correlation"][i, j],
                    self.grouped_data.obsp["spearman_correlation_pvalue"][i, j],
                ) = spearmanr(x, y)
                self.grouped_data.obsp["spearman_correlation"][j, i] = self.grouped_data.obsp["spearman_correlation"][i, j]
                self.grouped_data.obsp["spearman_correlation_pvalue"][j, i] = self.grouped_data.obsp[
                    "spearman_correlation_pvalue"
                ][i, j]

                # Pearson correlation
                (
                    self.grouped_data.obsp["pearson_correlation"][i, j],
                    self.grouped_data.obsp["pearson_correlation_pvalue"][i, j],
                ) = pearsonr(x, y)
                self.grouped_data.obsp["pearson_correlation"][j, i] = self.grouped_data.obsp["pearson_correlation"][i, j]
                self.grouped_data.obsp["pearson_correlation_pvalue"][j, i] = self.grouped_data.obsp[
                    "pearson_correlation_pvalue"
                ][i, j]

        self.grouped_data.uns["correlation"] = True
