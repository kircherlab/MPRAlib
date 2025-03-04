from abc import ABC, abstractmethod
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


class MPRAData(ABC):

    LOGGER = logging.getLogger(__name__)
    LOGGER.setLevel(logging.INFO)

    SCALING = 1e6
    PSEUDOCOUNT = 1

    @classmethod
    @abstractmethod
    def from_file(cls, file_path: str) -> "MPRAData":
        pass

    def __init__(self, data: ad.AnnData, barcode_threshold: int = 0):
        self._data = data
        self.var_filter = None
        self.barcode_threshold = barcode_threshold

    @property
    def data(self) -> ad.AnnData:
        return self._data

    @data.setter
    def data(self, new_data: ad.AnnData) -> None:
        self._data = new_data

    @property
    def var_names(self):
        return self.data.var_names

    @property
    def n_vars(self):
        return self.data.n_vars

    @property
    def obs_names(self):
        return self.data.obs_names

    @property
    def n_obs(self):
        return self.data.n_obs

    @property
    def oligos(self):
        return self.data.var["oligo"]

    @property
    def raw_rna_counts(self):
        return self.data.layers["rna"]

    @property
    def normalized_dna_counts(self):
        if "dna_normalized" not in self.data.layers:
            self._normalize()
        return self.data.layers["dna_normalized"]

    @property
    def rna_counts(self):
        if "rna_sampling" in self.data.layers:
            return self.data.layers["rna_sampling"] * ~self.var_filter.T.values
        else:
            return self.raw_rna_counts * ~self.var_filter.T.values

    @property
    def raw_dna_counts(self):
        return self.data.layers["dna"]

    @property
    def dna_counts(self):
        if "dna_sampling" in self.data.layers:
            return self.data.layers["dna_sampling"] * ~self.var_filter.T.values
        else:
            return self.raw_dna_counts * ~self.var_filter.T.values

    @property
    def normalized_rna_counts(self):
        if "rna_normalized" not in self.data.layers:
            self._normalize()
        return self.data.layers["rna_normalized"]

    @property
    def activity(self):
        if "log2FoldChange" not in self.data.layers:
            self._compute_activities()
        return self.data.layers["log2FoldChange"]

    def _compute_activities(self) -> None:
        ratio = np.divide(
            self.normalized_rna_counts,
            self.normalized_dna_counts,
            where=self.normalized_dna_counts != 0,
        )
        with np.errstate(divide='ignore'):
            log2ratio = np.log2(ratio)
            log2ratio[np.isneginf(log2ratio)] = np.nan
        self.data.layers["log2FoldChange"] = log2ratio * ~self.var_filter.T.values

    @property
    def observed(self):
        return (self.dna_counts + self.rna_counts) > 0

    @property
    def var_filter(self) -> pd.DataFrame:
        return self.data.varm["var_filter"]

    @var_filter.setter
    def var_filter(self, new_data: pd.DataFrame) -> None:
        if new_data is None:
            self.data.varm["var_filter"] = pd.DataFrame(
                np.full((self.data.n_vars, self.data.n_obs), False),
                index=self.var_names,
                columns=self.obs_names,
            )
            if "var_filter" in self.data.uns:
                del self.data.uns["var_filter"]
        else:
            self.data.varm["var_filter"] = new_data

        self.drop_normalized()

    @property
    def barcode_threshold(self) -> int:
        return self._get_metadata("barcode_threshold")

    @barcode_threshold.setter
    def barcode_threshold(self, barcode_threshold: int) -> None:
        self._add_metadata("barcode_threshold", barcode_threshold)
        self._drop_correlation()

    @abstractmethod
    def _normalize(self):
        pass

    def drop_normalized(self):

        self.LOGGER.info("Dropping normalized data")

        self.data.layers.pop("rna_normalized", None)
        self.data.layers.pop("dna_normalized", None)
        self._add_metadata("normalized", False)

    def correlation(self, method="pearson", count_type="activity") -> np.ndarray:
        """
        Calculates and return the correlation for activity or  normalized counts.

        Returns:
            np.ndarray: The Pearson or Spearman correlation matrix.
        """
        if count_type == "dna":
            return self._correlation(method, self.normalized_dna_counts, "dna_normalized")
        elif count_type == "rna":
            return self._correlation(method, self.normalized_rna_counts, "rna_normalized")
        elif count_type == "activity":
            return self._correlation(method, self.activity, "log2FoldChange")
        else:
            raise ValueError(f"Unsupported count type: {count_type}")

    def _correlation(self, method, data, layer):
        if not self._get_metadata(f"correlation_{layer}"):
            self._compute_correlation(data, layer)
        return self.data.obsp[f"{method}_correlation_{layer}"]

    def _compute_correlation(self, data, layer):
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
        num_columns = self.n_obs
        for correlation in ["pearson", "spearman"]:
            self.data.obsp[f"{correlation}_correlation_{layer}"] = np.zeros((num_columns, num_columns))
            self.data.obsp[f"{correlation}_correlation_{layer}_pvalue"] = np.zeros((num_columns, num_columns))

        def compute_correlation(x, y, method) -> tuple:
            if method == "spearman":
                return spearmanr(x, y)
            elif method == "pearson":
                return pearsonr(x, y)
            else:
                raise ValueError(f"Unsupported correlation method: {method}")

        for i in range(num_columns):
            for j in range(i, num_columns):
                mask = ~np.isnan(data[i, :]) & ~np.isnan(data[j, :])
                x = data[i, mask]
                y = data[j, mask]

                for method in ["spearman", "pearson"]:
                    corr, pvalue = compute_correlation(x, y, method)
                    self.data.obsp[f"{method}_correlation_{layer}"][i, j] = corr
                    self.data.obsp[f"{method}_correlation_{layer}_pvalue"][i, j] = pvalue
                    self.data.obsp[f"{method}_correlation_{layer}"][j, i] = corr
                    self.data.obsp[f"{method}_correlation_{layer}_pvalue"][j, i] = pvalue

        self._add_metadata(f"correlation_{layer}", True)

    def _drop_correlation(self) -> None:

        for layer in ["rna_normalized", "dna_normalized", "log2FoldChange"]:
            if self._get_metadata(f"correlation_{layer}"):
                self.LOGGER.info(f"Dropping correlation for layer {layer}")

                self._add_metadata(f"correlation_{layer}", False)
                for method in ["pearson", "spearman"]:
                    del self.data.obsp[f"{method}_correlation_{layer}"]
                    del self.data.obsp[f"{method}_correlation_{layer}_pvalue"]

    def write(self, file_data_path: str) -> None:
        self.data.write(file_data_path)

    def _add_metadata(self, key, value):
        if isinstance(value, list):
            if key not in self.data.uns:
                self.data.uns[key] = value
            else:
                self.data.uns[key] = self.data.uns[key] + value
        else:
            self.data.uns[key] = value

    def _get_metadata(self, key):
        if key in self.data.uns:
            return self.data.uns[key]
        else:
            return None


class MPRABarcodeData(MPRAData):

    @property
    def oligo_data(self) -> ad.AnnData:
        self.LOGGER.info("Computing oligo data")

        return self._oligo_data()

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

        spdis = set([item for sublist in self.oligo_data.var["SPDI"].values for item in sublist])

        df = {"ID": [], "REF": [], "ALT": []}
        for spdi in spdis:
            df["ID"].append(spdi)
            spdi_data = self.oligo_data[:, self.oligo_data.var["SPDI"].apply(lambda x: spdi in x)]
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
        df = {"ID": self.oligo_data.var["oligo"]}
        for replicate in self.oligo_data.obs_names:
            df["counts_" + replicate] = self.oligo_data[replicate, :].layers[layer][0]

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
        for replicate in self.oligo_data.obs_names:
            df["counts_" + replicate + "_REF"] = []
            df["counts_" + replicate + "_ALT"] = []

        for spdi, row in self.variant_map.iterrows():
            df["ID"].append(spdi)
            counts_ref = self.oligo_data.layers[layer][:, self.oligo_data.var["oligo"].isin(row["REF"])].sum(axis=1)
            counts_alt = self.oligo_data.layers[layer][:, self.oligo_data.var["oligo"].isin(row["ALT"])].sum(axis=1)
            idx = 0
            for replicate in self.oligo_data.obs_names:
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
        self.oligo_data = None

    @classmethod
    def from_file(cls, file_path: str) -> "MPRABarcodeData":
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

        adata.varm["var_filter"] = pd.DataFrame(
            np.full((adata.n_vars, adata.n_obs), False),
            index=adata.var_names,
            columns=adata.obs_names,
        )

        return cls(adata)

    def _normalize(self):

        self.drop_normalized()

        self.LOGGER.info("Normalizing data")

        self.data.layers["dna_normalized"] = self._normalize_layer(
            self.dna_counts, self.observed, ~self.var_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self.data.layers["rna_normalized"] = self._normalize_layer(
            self.rna_counts, self.observed, ~self.var_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self._add_metadata("normalized", True)

    def _normalize_layer(
        self, counts: np.ndarray, observed: np.ndarray, not_var_filter: np.ndarray, scaling: float, pseudocount: int
    ) -> np.ndarray:

        # I do a pseudo count when normalizing to avoid division by zero when computing logfold ratios.
        # barcode filter has to be used again because we want to have a zero on filtered values.
        total_counts = np.sum((counts + (pseudocount * observed)) * not_var_filter, axis=1)

        # Avoid division by zero when pseudocount is set to 0
        total_counts[total_counts == 0] = 1
        return ((counts + (pseudocount * observed)) / total_counts[:, np.newaxis] * scaling) * not_var_filter

    def _oligo_data(self) -> "MPRAOligoData":

        # Convert the result back to an AnnData object
        oligo_data = ad.AnnData(self._sum_counts_by_oligo(self.rna_counts))

        oligo_data.layers["rna"] = oligo_data.X.copy()
        oligo_data.layers["dna"] = self._sum_counts_by_oligo(self.dna_counts)

        oligo_data.layers["barcode_counts"] = self._compute_supporting_barcodes()

        oligo_data.obs_names = self.obs_names
        oligo_data.var_names = self.data.var["oligo"].unique().tolist()

        # Subset of vars using the firs occurence of oligo name
        indices = self.data.var["oligo"].drop_duplicates(keep="first").index
        oligo_data.var = self.data.var.loc[indices]

        oligo_data.obs = self.data.obs

        for key, value in self.data.uns.items():
            oligo_data.uns[f"MPRABarcodeData_{key}"] = value

        oligo_data.uns["correlation_log2FoldChange"] = False
        oligo_data.uns["correlation_rna_normalized"] = False
        oligo_data.uns["correlation_dna_normalized"] = False

        return MPRAOligoData(oligo_data, self.barcode_threshold)

    def _sum_counts_by_oligo(self, counts):

        grouped = pd.DataFrame(
            counts,
            index=self.obs_names,
            columns=self.var_names,
        ).T.groupby(self.data.var["oligo"], observed=True)

        # Perform an operation on each group, e.g., mean
        return grouped.sum().T

    def _compute_supporting_barcodes(self):

        grouped = pd.DataFrame(
            self.observed * ~self.var_filter.T.values,
            index=self.obs_names,
            columns=self.var_names,
        ).T.groupby(self.oligos, observed=True)

        return grouped.apply(lambda x: x.sum()).T

    def _barcode_filter_rna_zscore(self, times_zscore=3):

        barcode_mask = pd.DataFrame(
            self.raw_dna_counts + self.raw_rna_counts,
            index=self.obs_names,
            columns=self.var_names,
        ).T.apply(lambda x: (x != 0))

        df_rna = pd.DataFrame(self.raw_rna_counts, index=self.obs_names, columns=self.var_names).T
        grouped = df_rna.where(barcode_mask).groupby(self.oligos, observed=True)

        mask = ((df_rna - grouped.transform("mean")) / grouped.transform("std")).abs() > times_zscore

        return mask

    def _barcode_filter_mad(self, times_mad=3, n_bins=20):

        # sum up DNA and RNA counts across replicates
        DNA_sum = pd.DataFrame(self.raw_dna_counts, index=self.obs_names, columns=self.var_names).T.sum(axis=1)
        RNA_sum = pd.DataFrame(self.raw_rna_counts, index=self.obs_names, columns=self.var_names).T.sum(axis=1)
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

        return self.var_filter.apply(lambda col: col | ~self.var_filter.index.isin(df_sums.index))

    def _barcode_filter_random(self, proportion=1.0, total=None, aggegate_over_replicates=True):

        if aggegate_over_replicates and total is None:
            total = self.var_filter.shape[0]
        elif total is None:
            total = self.var_filter.size

        num_true_cells = int(total * (1.0 - proportion))
        true_indices = np.random.choice(total, num_true_cells, replace=False)

        mask = pd.DataFrame(
            np.full((self.data.n_vars, self.data.n_obs), False),
            index=self.var_names,
            columns=self.obs_names,
        )

        if aggegate_over_replicates:
            mask.iloc[true_indices, :] = True
        else:
            flat_df = mask.values.flatten()

            flat_df[true_indices] = True

            mask = pd.DataFrame(
                flat_df.reshape(self.var_filter.shape),
                index=self.var_names,
                columns=self.obs_names,
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
            np.full((self.n_vars, self.n_obs), False),
            index=self.var_names,
            columns=self.obs_names,
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
            self.var_filter = self.var_filter | filter_func(**params)
        else:
            raise ValueError(f"Unsupported barcode filter: {barcode_filter}")

        self._add_metadata("var_filter", [barcode_filter.value])

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
                proportion, total, aggregate_over_replicates, self.data.layers[layer_name], self.n_obs
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

        self.drop_normalized()


class MPRAOligoData(MPRAData):

    @property
    def barcode_counts(self):
        return self.data.layers["barcode_counts"]

    def correlation(self, method="pearson", count_type="activity") -> np.ndarray:
        """
        Calculates and return the correlation for activity or  normalized counts.

        Returns:
            np.ndarray: The Pearson or Spearman correlation matrix.
        """
        if count_type == "dna":
            filtered = self.normalized_dna_counts.copy()
            layer_name = "dna_normalized"
        elif count_type == "rna":
            layer_name = "rna_normalized"
            filtered = self.normalized_rna_counts.copy()
        elif count_type == "activity":
            filtered = self.activity.copy()
            layer_name = "log2FoldChange"
        else:
            raise ValueError(f"Unsupported count type: {count_type}")

        filtered[self.barcode_counts < self.barcode_threshold] = np.nan
        return self._correlation(method, filtered, layer_name)

    @classmethod
    def from_file(cls, file_path: str) -> "MPRAOligoData":

        return MPRAOligoData(ad.read(file_path))

    def _normalize(self):

        self.drop_normalized()

        self.LOGGER.info("Normalizing data")

        self.data.layers["dna_normalized"] = self._normalize_layer(
            self.dna_counts, ~self.var_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self.data.layers["rna_normalized"] = self._normalize_layer(
            self.rna_counts, ~self.var_filter.T.values, self.SCALING, self.PSEUDOCOUNT
        )
        self._add_metadata("normalized", True)

    def _normalize_layer(self, counts: np.ndarray, not_var_filter: np.ndarray, scaling: float, pseudocount: int) -> np.ndarray:

        # I do a pseudo count when normalizing to avoid division by zero when computing logfold ratios.
        # Pseudocount has also be done per barcode.
        # var filter has to be used again because we want to have a zero on filtered values.
        total_counts = np.sum((counts + (pseudocount * self.barcode_counts)) * not_var_filter, axis=1)

        # Avoid division by zero when pseudocount is set to 0
        total_counts[total_counts == 0] = 1
        scaled_counts = (counts + (pseudocount * self.barcode_counts)) / total_counts[:, np.newaxis] * scaling
        return (
            np.divide(scaled_counts, self.barcode_counts, where=self.barcode_counts != 0, out=np.zeros_like(scaled_counts))
            * not_var_filter
        )
