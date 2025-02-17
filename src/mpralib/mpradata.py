import numpy as np
import ast
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, pearsonr
from enum import Enum


class OutlierFilter(Enum):
    RNA_ZSCORE = "RNA_ZSCORE"
    MAD = "MAD"
    RANDOM = "RANDOM"


class MPRAdata:

    SCALING = 1e6

    def __init__(self, data: ad.AnnData):
        self._data = data
        self._grouped_data = None

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
        self._grouped_data = new_data

    @property
    def rna_counts(self):
        return self.data.layers["rna"]

    @property
    def dna_counts(self):
        return self.data.layers["dna"]

    @property
    def oligos(self):
        return self.data.var["oligo"]

    @property
    def barcode_threshold(self) -> int:
        return self.data.uns["barcode_threshold"]

    @barcode_threshold.setter
    def barcode_threshold(self, barcode_threshold: int) -> None:
        self.data.uns["barcode_threshold"] = barcode_threshold
        self.grouped_data = None

    @property
    def filter(self) -> pd.DataFrame:
        return self.data.varm["filter"]

    @filter.setter
    def filter(self, new_data: pd.DataFrame) -> None:
        if new_data is None:
            self.data.varm["filter"] = pd.DataFrame(
                np.full((self.data.n_vars, self.data.n_obs), False),
                index=self.data.var_names,
                columns=self.data.obs_names,
            )
            self.data.uns["filter"] = []
        else:
            self.data.varm["filter"] = new_data
        self.grouped_data = None
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

        spdis = set(
            [
                item
                for sublist in self.grouped_data.var["SPDI"].values
                for item in sublist
            ]
        )

        df = {"ID": [], "REF": [], "ALT": []}
        for spdi in spdis:
            df["ID"].append(spdi)
            spdi_data = self.grouped_data[
                :, self.grouped_data.var["SPDI"].apply(lambda x: spdi in x)
            ]
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
            counts_ref = self.grouped_data.layers[layer][
                :, self.grouped_data.var["oligo"].isin(row["REF"])
            ].sum(axis=1)
            counts_alt = self.grouped_data.layers[layer][
                :, self.grouped_data.var["oligo"].isin(row["ALT"])
            ].sum(axis=1)
            idx = 0
            for replicate in self.grouped_data.obs_names:
                df["counts_" + replicate + "_REF"].append(counts_ref[idx])
                df["counts_" + replicate + "_ALT"].append(counts_alt[idx])
                idx += 1
        df = pd.DataFrame(df).set_index("ID")
        return df[(df.T != 0).all()]

    def add_metadata_file(self, metadata_file):
        metadata = pd.read_csv(
            metadata_file, sep="\t", header=0, na_values=["NA"]
        ).drop_duplicates()
        metadata["name"] = pd.Categorical(metadata["name"].str.replace(" ", "_"))
        metadata.set_index("name", inplace=True)
        metadata["variant_class"] = (
            metadata["variant_class"].fillna("[]").apply(ast.literal_eval)
        )
        metadata["variant_pos"] = (
            metadata["variant_pos"].fillna("[]").apply(ast.literal_eval)
        )
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

        self.data.uns["metadata_file"] = metadata_file

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

        anndata_replicate_rna = (
            data.loc[:, replicate_columns_rna].transpose().astype(int)
        )
        anndata_replicate_dna = (
            data.loc[:, replicate_columns_dna].transpose().astype(int)
        )

        anndata_replicate_rna.index = [
            replicate.split("_")[2] for replicate in replicate_columns_rna
        ]
        anndata_replicate_dna.index = [
            replicate.split("_")[2] for replicate in replicate_columns_dna
        ]

        adata = ad.AnnData(anndata_replicate_rna)
        adata.layers["rna"] = adata.X.copy()
        adata.layers["dna"] = anndata_replicate_dna

        adata.obs["replicate"] = [
            replicate.split("_")[2] for replicate in replicate_columns_rna
        ]

        adata.var["oligo"] = pd.Categorical(data["oligo_name"].values)

        adata.uns["file_path"] = file_path
        adata.uns["date"] = pd.to_datetime("today").strftime("%Y-%m-%d")

        adata.uns["normalized"] = False
        adata.uns["barcode_threshold"] = None

        adata.varm["filter"] = pd.DataFrame(
            np.full((adata.n_vars, adata.n_obs), False),
            index=adata.var_names,
            columns=adata.obs_names,
        )
        adata.uns["filter"] = []

        return cls(adata)

    def _number_of_barcodes(self):
        bdata = self.data[
            :,
            self.data.var.oligo
            == "A:HNF4A-ChMod_chr10:11917871-11917984__chr10:11917842-11918013_:015",
        ]
        print(bdata.var.rna)
        filter = (bdata.layers["rna"] > 0) & (bdata.layers["dna"] > 0)
        print(filter)
        filtered_adata = bdata[filter, :]
        print(filtered_adata)
        barcode_counts = np.sum(bdata.X != 0, axis=1)
        print(barcode_counts)

        return bdata.X.shape[0]

    def _normalize(self):

        self._normalize_layer("dna")
        self._normalize_layer("rna")

    def _normalize_layer(self, layer_name):
        total_counts = np.sum(
            self.data.layers[layer_name] * ~self.data.varm["filter"].T.values, axis=1
        )
        total_counts[total_counts == 0] = 1
        self.data.layers[layer_name + "_normalized"] = (
            self.data.layers[layer_name] / total_counts[:, np.newaxis] * self.SCALING
        )
        self.data.uns["normalized"] = True

    def _group_data(self):

        if not self.data.uns["normalized"]:
            self._normalize()

        # Convert the result back to an AnnData object
        self.grouped_data = ad.AnnData(self._group_sum_data_layer("rna"))

        self.grouped_data.layers["rna"] = self.grouped_data.X.copy()
        self.grouped_data.layers["dna"] = self._group_sum_data_layer("dna")

        self.grouped_data.layers["barcodes"] = self._compute_supporting_barcodes()

        self.grouped_data.layers["rna_normalized"] = (
            self._group_sum_data_layer("rna_normalized")
            / self.grouped_data.layers["barcodes"]
        )
        self.grouped_data.layers["dna_normalized"] = (
            self._group_sum_data_layer("dna_normalized")
            / self.grouped_data.layers["barcodes"]
        )

        self.grouped_data.obs_names = self.data.obs_names
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

    def _group_sum_data_layer(self, layer_name):

        grouped = pd.DataFrame(
            self.data.layers[layer_name] * ~self.filter.T.values,
            index=self.data.obs_names,
            columns=self.data.var_names,
        ).T.groupby(self.data.var["oligo"], observed=True)

        # Perform an operation on each group, e.g., mean
        return grouped.sum().T

    def _compute_supporting_barcodes(self):

        grouped = pd.DataFrame(
            (self.dna_counts + self.dna_counts) * ~self.filter.T.values,
            index=self.data.obs_names,
            columns=self.data.var_names,
        ).T.groupby(self.oligos, observed=True)

        return grouped.apply(lambda x: (x != 0).sum()).T

    def _compute_activities(self):
        self.grouped_data.layers["log2FoldChange"] = np.log2(
            self.grouped_data.layers["rna_normalized"]
            / self.grouped_data.layers["dna_normalized"]
        )

    def _filter_rna_zscore(self, times_zscore=3):

        barcode_mask = pd.DataFrame(
            self.dna_counts + self.rna_counts,
            index=self.data.obs_names,
            columns=self.data.var_names,
        ).T.apply(lambda x: (x != 0))

        df_rna = pd.DataFrame(
            self.rna_counts, index=self.data.obs_names, columns=self.data.var_names
        ).T
        grouped = df_rna.where(barcode_mask).groupby(self.oligos, observed=True)

        mask = (
            (df_rna - grouped.transform("mean")) / grouped.transform("std")
        ).abs() > times_zscore

        return mask

    def _filter_mad(self, times_mad=3, n_bins=20):

        # sum up DNA and RNA counts across replicates
        DNA_sum = pd.DataFrame(
            self.rna_counts, index=self.data.obs_names, columns=self.data.var_names
        ).T.sum(axis=1)
        RNA_sum = pd.DataFrame(
            self.rna_counts, index=self.data.obs_names, columns=self.data.var_names
        ).T.sum(axis=1)
        df_sums = pd.DataFrame({"DNA_sum": DNA_sum, "RNA_sum": RNA_sum}).fillna(0)
        # removing all barcodes with 0 counts in RNA an more DNA count than number of replicates/observations
        df_sums = df_sums[
            (df_sums["DNA_sum"] > self.data.n_obs) & (df_sums["RNA_sum"] > 0)
        ]

        # remove all barcodes where oligo has less barcodes as the number of replicates/observations
        df_sums = df_sums.groupby(self.data.var["oligo"], observed=True).filter(
            lambda x: len(x) >= self.data.n_obs
        )

        # Calculate ratio, ratio_med, ratio_diff, and mad
        df_sums["ratio"] = np.log2(df_sums["DNA_sum"] / df_sums["RNA_sum"])
        df_sums["ratio_med"] = df_sums.groupby(self.data.var["oligo"], observed=True)[
            "ratio"
        ].transform("median")
        df_sums["ratio_diff"] = df_sums["ratio"] - df_sums["ratio_med"]
        # df_sums['mad'] = (df_sums['ratio'] - df_sums['ratio_med']).abs().mean()

        # Calculate quantiles within  n_bins
        qs = np.unique(
            np.quantile(np.log10(df_sums["RNA_sum"]), np.arange(0, n_bins) / n_bins)
        )

        # Create bins based on rna_count
        df_sums["bin"] = pd.cut(
            np.log10(df_sums["RNA_sum"]),
            bins=qs,
            include_lowest=True,
            labels=[str(i) for i in range(0, len(qs) - 1)],
        )
        # Filter based on ratio_diff and mad
        df_sums["ratio_diff_med"] = df_sums.groupby("bin", observed=True)[
            "ratio_diff"
        ].transform("median")
        df_sums["ratio_diff_med_dist"] = np.abs(
            df_sums["ratio_diff"] - df_sums["ratio_diff_med"]
        )
        df_sums["mad"] = df_sums.groupby("bin", observed=True)[
            "ratio_diff_med_dist"
        ].transform("median")

        m = df_sums.ratio_diff > times_mad * df_sums.mad
        df_sums = df_sums[~m]

        # print(self.data.varm["filter"][~self.data.varm["filter"].index.isin(df_sums.index)])
        return self.filter.apply(
            lambda col: col | ~self.filter.index.isin(df_sums.index)
        )

    def _filter_random(self, proportion=1.0, total=None, aggegate_over_obs=True):

        if aggegate_over_obs and total is None:
            total = self.filter.shape[0]
        elif total is None:
            total = self.filter.size

        num_true_cells = int(total * (1.0 - proportion))
        true_indices = np.random.choice(total, num_true_cells, replace=False)

        mask = pd.DataFrame(
            np.full((self.data.n_vars, self.data.n_obs), False),
            index=self.data.var_names,
            columns=self.data.obs_names,
        )

        if aggegate_over_obs:
            mask.iloc[true_indices, :] = True
        else:
            flat_df = mask.values.flatten()

            flat_df[true_indices] = True

            mask = pd.DataFrame(
                flat_df.reshape(self.filter.shape),
                index=self.data.var_names,
                columns=self.data.obs_names,
            )
        return mask

    def filter_outlier(self, outlier_filter, params):
        filter_switch = {
            OutlierFilter.RNA_ZSCORE: self._filter_rna_zscore,
            OutlierFilter.MAD: self._filter_mad,
            OutlierFilter.RANDOM: self._filter_random,
        }

        filter_func = filter_switch.get(outlier_filter)
        if filter_func:
            self.filter = self.filter | filter_func(**params)
        else:
            raise ValueError(f"Unsupported outlier filter: {outlier_filter}")

        self.data.uns["filter"].append(outlier_filter.value)

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
        self.grouped_data.obsp["pearson_correlation"] = np.zeros(
            (num_columns, num_columns)
        )
        self.grouped_data.obsp["pearson_correlation_pvalue"] = np.zeros(
            (num_columns, num_columns)
        )
        self.grouped_data.obsp["spearman_correlation"] = np.zeros(
            (num_columns, num_columns)
        )
        self.grouped_data.obsp["spearman_correlation_pvalue"] = np.zeros(
            (num_columns, num_columns)
        )

        for i in range(num_columns):
            for j in range(i, num_columns):
                mask = ~np.isnan(
                    self.grouped_data.layers["log2FoldChange"][i, :]
                ) & ~np.isnan(self.grouped_data.layers["log2FoldChange"][j, :])
                x = self.grouped_data.layers["log2FoldChange"][i, mask]
                y = self.grouped_data.layers["log2FoldChange"][j, mask]

                # Spearman correlation
                (
                    self.grouped_data.obsp["spearman_correlation"][i, j],
                    self.grouped_data.obsp["spearman_correlation_pvalue"][i, j],
                ) = spearmanr(x, y)
                self.grouped_data.obsp["spearman_correlation"][j, i] = (
                    self.grouped_data.obsp["spearman_correlation"][i, j]
                )
                self.grouped_data.obsp["spearman_correlation_pvalue"][j, i] = (
                    self.grouped_data.obsp["spearman_correlation_pvalue"][i, j]
                )

                # Pearson correlation
                (
                    self.grouped_data.obsp["pearson_correlation"][i, j],
                    self.grouped_data.obsp["pearson_correlation_pvalue"][i, j],
                ) = pearsonr(x, y)
                self.grouped_data.obsp["pearson_correlation"][j, i] = (
                    self.grouped_data.obsp["pearson_correlation"][i, j]
                )
                self.grouped_data.obsp["pearson_correlation_pvalue"][j, i] = (
                    self.grouped_data.obsp["pearson_correlation_pvalue"][i, j]
                )

        self.grouped_data.uns["correlation"] = True
