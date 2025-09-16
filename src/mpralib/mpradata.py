from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import spearmanr, pearsonr
from enum import Enum
import logging
import os
from mpralib.exception import MPRAlibException
from typing import Any, Callable, Optional
from numpy.typing import NDArray


class Modality(Enum):
    """An enumeration representing different data modalities in MPRA (Massively Parallel Reporter Assay) experiments."""

    DNA = "DNA"
    """str: Represents DNA data modality. """
    RNA = "RNA"
    """str: Represents RNA data modality."""
    DNA_NORMALIZED = "DNA_NORMALIZED"
    """str: Represents normalized DNA data modality."""
    RNA_NORMALIZED = "RNA_NORMALIZED"
    """str: Represents normalized RNA data modality."""
    ACTIVITY = "ACTIVITY"
    """str: Represents activity data modality, typically calculated as the log2 ratio of normalized RNA to DNA counts."""

    def __new__(cls, value):
        obj = object.__new__(cls)
        obj._value_ = str(value).lower()
        return obj

    @classmethod
    def from_string(cls, value: str) -> "Modality":
        """Creates a Modality enum member from a string value.

        Args:
            value (str): The string representation of the enum member.

        Returns:
            The corresponding Modality enum member.

        Raises:
            ValueError: If the provided string does not match any Modality member.
        """
        for member in cls:
            if member.value == value.lower():
                return member
        raise ValueError(f"{value} is not a valid {cls.__name__}")


class CountSampling(Enum):
    """Enumeration representing the types of count sampling available for MPRA data."""

    RNA = "RNA"
    """str: Represents RNA count sampling."""
    DNA = "DNA"
    """str: Represents DNA count sampling."""
    RNA_AND_DNA = "RNA_AND_DNA"
    """str: Represents both RNA and DNA count sampling."""


class BarcodeFilter(Enum):
    """Enumeration of available barcode filtering methods."""

    MIN_BCS_PER_OLIGO = "MIN_BCS_PER_OLIGO"
    """str: Filter barcodes based on a minimum number of barcodes per oligo."""
    GLOBAL = "GLOBAL"
    """str: Filter barcodes based on RNA z-score."""
    OLIGO_SPECIFIC = "OLIGO_SPECIFIC"
    """str: Filter barcodes based on standard deviation per oligo"""
    LARGE_EXPRESSION = "LARGE_EXPRESSION"
    """str: Filter barcodes based on Median Absolute Deviation (MAD)."""
    RANDOM = "RANDOM"
    """str: Randomly filter barcodes."""
    MIN_COUNT = "MIN_COUNT"
    """str: Filter barcodes with counts below a specified minimum."""
    MAX_COUNT = "MAX_COUNT"
    """str: Filter barcodes with counts above a specified maximum."""

    def __new__(cls, value):
        obj = object.__new__(cls)
        obj._value_ = str(value).lower()
        return obj

    @classmethod
    def from_string(cls, value: str) -> "BarcodeFilter":
        """Creates a BarcodeFilter enum member from a string value.

        Args:
            value (str): The string representation of the enum member.

        Returns:
            The corresponding BarcodeFilter enum member.

        Raises:
            ValueError: If the provided string does not match any BarcodeFilter member.
        """
        for member in cls:
            if member.value == value.lower():
                return member
        raise ValueError(f"{value} is not a valid {cls.__name__}")


class MPRAData(ABC):
    """Abstract base class for handling MPRA (Massively Parallel Reporter Assay) data using AnnData objects.

    This class provides a standardized interface and core functionality for managing, normalizing filtering, and analyzing MPRA data, including DNA/RNA counts, barcode handling, activity computation, and correlation analysis. It is designed to be subclassed for specific MPRA data formats.

    Args:
        data (anndata.AnnData): The AnnData object containing MPRA data.
        barcode_threshold (int, optional): Minimum barcode count threshold for filtering. Defaults to 0.

    Attributes:
        _SCALING (float): Default scaling factor for normalization.
        _PSEUDOCOUNT (int): Default pseudocount for normalization.
        _data (anndata.AnnData): The AnnData object containing MPRA data.

    Raises:
        ValueError: If required metadata (e.g., sequence design file) is not loaded.
    """  # noqa: E501

    LOGGER = logging.getLogger(__name__)
    """logging.Logger: Logger for the class."""
    LOGGER.setLevel(logging.INFO)

    _SCALING = 1e6
    _PSEUDOCOUNT = 1

    @classmethod
    @abstractmethod
    def from_file(cls, file_path: str) -> "MPRAData":
        """Create an instance of the class from a file.

        This method reads data from a specified file (reporter experiment barcode or reporter experiment file format), processes it, and returns an instance of the class containing the data in an AnnData object.

        Args:
            file_path (str): Path to the input file containing reporter experiment barcode or reporter experiment data.

        Returns:
            An instance of MPRAData containing the processed data in an AnnData object.

        Raises:
            IOError: If the file cannot be read or parsed.
            ValueError: If the file format is invalid.
        """  # noqa: E501
        pass

    def __init__(self, data: ad.AnnData, barcode_threshold: int = 0):
        self._data = data
        self.var_filter = None
        self.barcode_threshold = barcode_threshold

        # Initialize scaling and pseudo count metadata
        scaling = self._get_metadata("SCALING")
        if scaling is not None:
            self._SCALING = scaling
        else:
            self._add_metadata("SCALING", self._SCALING)

        pseudo_count = self._get_metadata("PSEUDOCOUNT")
        if pseudo_count is not None:
            self._PSEUDOCOUNT = pseudo_count
        else:
            self._add_metadata("PSEUDOCOUNT", self._PSEUDOCOUNT)

    @property
    def scaling(self) -> float:
        """float: Scaling factor for normalization."""
        return self._SCALING

    @scaling.setter
    def scaling(self, new_scaling: float) -> None:
        if new_scaling != self._SCALING:
            self.drop_normalized()
            self._SCALING = new_scaling
            self._add_metadata("SCALING", self._SCALING)

    @property
    def pseudo_count(self) -> int:
        """int: Pseudocount added during normalization to avoid division by zero."""
        return self._PSEUDOCOUNT

    @pseudo_count.setter
    def pseudo_count(self, new_pseudo_count: int) -> None:
        if new_pseudo_count != self._PSEUDOCOUNT:
            self.drop_normalized()
            self._PSEUDOCOUNT = new_pseudo_count
            self._add_metadata("PSEUDOCOUNT", self._PSEUDOCOUNT)

    @property
    def data(self) -> ad.AnnData:
        """ad.AnnData: The underlying AnnData object containing MPRA data."""
        return self._data

    @data.setter
    def data(self, new_data: ad.AnnData) -> None:
        self._data = new_data

    @property
    def var_names(self) -> pd.Index:
        """pd.Index: Returns the variable names (samples) of the dataset."""
        return self.data.var_names

    @property
    def n_vars(self) -> int:
        """int: Returns the number of variables (samples) in the dataset."""
        return self.data.n_vars

    @property
    def obs_names(self) -> pd.Index:
        """pd.Index: Returns the observation names (barcodes) of the dataset."""
        return self.data.obs_names

    @property
    def n_obs(self) -> int:
        """int: Returns the number of observations (barcodes) in the dataset."""
        return self.data.n_obs

    @property
    def oligos(self) -> pd.Series:
        """pd.Series: Returns the oligo names for each variable in the dataset."""
        return self.data.var["oligo"]

    @property
    def raw_rna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the raw RNA counts from the dataset."""
        return np.asarray(self.data.layers["rna"])

    @property
    def normalized_dna_counts(self) -> NDArray[np.float32]:
        """NDArray[np.float32]: Returns the normalized DNA counts from the dataset, applying the variable filter if present."""
        if "dna_normalized" not in self.data.layers:
            self._normalize()
        return np.asarray(self.data.layers["dna_normalized"]) * ~self.var_filter.T

    @property
    def rna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the raw RNA or, if present, sampled RNA counts, applying the variable filter if present."""  # noqa: E501
        if "rna_sampling" in self.data.layers:
            return np.asarray(self.data.layers["rna_sampling"]) * ~self.var_filter.T
        else:
            return self.raw_rna_counts * ~self.var_filter.T

    @property
    def raw_dna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the raw DNA counts from the dataset."""
        return np.asarray(self.data.layers["dna"])

    @property
    def dna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the raw DNA or, if present, sampled DNA counts, applying the variable filter if present."""  # noqa: E501
        if "dna_sampling" in self.data.layers:
            return np.asarray(self.data.layers["dna_sampling"]) * ~self.var_filter.T
        else:
            return self.raw_dna_counts * ~self.var_filter.T

    @property
    def normalized_rna_counts(self) -> NDArray[np.float32]:
        """NDArray[np.float32]: Returns the normalized RNA counts from the dataset, applying the variable filter if present."""
        if "rna_normalized" not in self.data.layers:
            self._normalize()
        return np.asarray(self.data.layers["rna_normalized"]) * ~self.var_filter.T

    @property
    def activity(self) -> NDArray[np.float32]:
        """NDArray[np.float32]: Returns the activity values calculated from normalized RNA and DNA counts, applying the variable filter if present."""  # noqa: E501
        if "activity" not in self.data.layers:
            self._compute_activities()
        return np.asarray(self.data.layers["activity"]) * ~self.var_filter.T

    def _compute_activities(self) -> None:
        ratio = np.divide(
            self.normalized_rna_counts,
            self.normalized_dna_counts,
            where=self.normalized_dna_counts != 0,
        )
        with np.errstate(divide="ignore"):
            log2ratio = np.log2(ratio)
            log2ratio[np.isneginf(log2ratio)] = np.nan
        self.data.layers["activity"] = log2ratio

    @property
    def total_dna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the total DNA counts for each replicate. Usually it are the total raw counts per replicate. Only when sampled data is availabe it returns the sampled counts."""  # noqa: E501
        if "dna_counts" not in self.data.obs:
            if "dna_sampling" in self.data.layers:
                self.data.obs["dna_counts"] = np.sum(np.asarray(self.data.layers["dna_sampling"]), axis=1)
            else:
                self.data.obs["dna_counts"] = np.sum(self.raw_dna_counts, axis=1)
        return np.asarray(self.data.obs["dna_counts"])

    @property
    def total_rna_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the total RNA counts for each replicate. Usually it are the total raw counts per replicate. Only when sampled data is availabe it returns the sampled counts."""  # noqa: E501
        if "rna_counts" not in self.data.obs:
            if "rna_sampling" in self.data.layers:
                self.data.obs["rna_counts"] = np.sum(np.asarray(self.data.layers["rna_sampling"]), axis=1)
            else:
                self.data.obs["rna_counts"] = np.sum(self.raw_rna_counts, axis=1)
        return np.asarray(self.data.obs["rna_counts"])

    def drop_total_counts(self) -> None:
        """Removes total RNA and DNA counts from the dataset."""
        if "rna_counts" in self.data.obs:
            del self.data.obs["rna_counts"]
        if "dna_counts" in self.data.obs:
            del self.data.obs["dna_counts"]

    @property
    def observed(self) -> NDArray[np.bool_]:
        """:class:`NDArray[np.bool_]`: Returns a boolean NumPy array indicating which barcodes (observations) have non-zero counts in either DNA or RNA. Uses sampled counts when available. otherwise raw counts"""  # noqa: E501

        if "dna_sampling" in self.data.layers:
            dna_counts = np.asarray(self.data.layers["dna_sampling"])
        else:
            dna_counts = self.raw_dna_counts
        if "rna_sampling" in self.data.layers:
            rna_counts = np.asarray(self.data.layers["rna_sampling"])
        else:
            rna_counts = self.raw_rna_counts

        return (dna_counts + rna_counts) > 0

    @property
    def var_filter(self) -> NDArray[np.bool_]:
        """:class:`NDArray[np.bool_]`: Returns a boolean NumPy array indicating which variables (samples) are filtered out."""

        return np.asarray(self.data.varm["var_filter"], dtype=np.bool_)

    @var_filter.setter
    def var_filter(self, new_data: Optional[NDArray[np.bool_]]) -> None:
        if new_data is None:
            self.data.varm["var_filter"] = np.full((self.data.n_vars, self.data.n_obs), False)
            if "var_filter" in self.data.uns:
                del self.data.uns["var_filter"]
        else:
            self.data.varm["var_filter"] = new_data

    @property
    @abstractmethod
    def barcode_counts(self) -> NDArray[np.int32]:
        """NDArray[np.int32]: Returns the barcode counts matrix, which is the number of observed and not filtered barcodes for each oligo."""  # noqa: E501
        pass

    @barcode_counts.setter
    def barcode_counts(self, new_data: NDArray[np.int32]) -> None:
        self.data.layers["barcode_counts"] = new_data

    @abstractmethod
    def drop_barcode_counts(self) -> None:
        """Removes or clears the barcode counts data from the current object.

        Raises:
            NotImplementedError: If the method is not yet implemented.
        """
        pass

    @property
    def barcode_threshold(self) -> int:
        """int: Returns the threshold for barcode filtering."""
        threshold = self._get_metadata("barcode_threshold")
        if threshold is None:
            # If no threshold is set, default to 0
            return 0
        else:
            return threshold

    @barcode_threshold.setter
    def barcode_threshold(self, barcode_threshold: int) -> None:
        self._add_metadata("barcode_threshold", barcode_threshold)
        self._drop_correlation()

    @property
    def variant_map(self) -> pd.DataFrame:
        """pd.DataFrame: Returns a DataFrame mapping SPDI IDs to alleles and oligos.

        Raises:
            ValueError: If the sequence design file is not loaded in the metadata.
        """
        if not self._get_metadata("sequence_design_file") and not (
            isinstance(self, MPRAOligoData) and self._get_metadata("MPRABarcodeData_sequence_design_file")
        ):
            raise ValueError("Sequence design file not loaded.")

        oligos = self.data.var["oligo"].repeat(self.data.var["SPDI"].apply(lambda x: len(x)).tolist())

        spdis = np.concatenate(self.data.var["SPDI"].values)

        alleles = np.concatenate(self.data.var["allele"].values)

        df = pd.DataFrame({"ID": spdis, "allele": alleles, "oligo": oligos})
        df["REF"] = df["oligo"].where(df["allele"] == "ref", None)
        df["ALT"] = df["oligo"].where(df["allele"] == "alt", None)
        df = df.groupby("ID").agg({"REF": lambda x: list(filter(None, x)), "ALT": lambda x: list(filter(None, x))})
        df = df[(df["REF"].apply(len) >= 1) & (df["ALT"].apply(len) >= 1)]

        return df

    def _normalize(self) -> None:
        self.drop_normalized()

        self.LOGGER.info("Normalizing data")

        if "dna_sampling" in self.data.layers:
            dna_counts = np.asarray(self.data.layers["dna_sampling"])
        else:
            dna_counts = self.raw_dna_counts
        if "rna_sampling" in self.data.layers:
            rna_counts = np.asarray(self.data.layers["rna_sampling"])
        else:
            rna_counts = self.raw_rna_counts

        self.data.layers["dna_normalized"] = self._normalize_layer(dna_counts, self.total_dna_counts)
        self.data.layers["rna_normalized"] = self._normalize_layer(rna_counts, self.total_rna_counts)
        self._add_metadata("normalized", True)

    @abstractmethod
    def _normalize_layer(
        self,
        counts: NDArray[np.int32],
        total_counts: NDArray[np.int32],
    ) -> NDArray[np.float32]:
        pass

    def drop_normalized(self) -> None:
        """Removes normalized RNA and DNA data layers from the dataset.

        This method deletes the "rna_normalized" and "dna_normalized" layers from the `self.data.layers` attribute, logs the operation, updates the metadata to indicate that normalization is no longer present, and drops any associated correlation data.
        """  # noqa: E501

        self.LOGGER.info("Dropping normalized data")

        self.data.layers.pop("rna_normalized", None)
        self.data.layers.pop("dna_normalized", None)
        self._drop_correlation()
        self._add_metadata("normalized", False)

    def correlation(self, method: str = "pearson", count_type: Modality = Modality.ACTIVITY) -> NDArray[np.float32]:
        """Calculates and return the correlation for activity or normalized counts.

        Returns:
            The Pearson or Spearman correlation matrix.
        """

        if method not in {"pearson", "spearman"}:
            raise ValueError(f"Unsupported correlation method: {method}")

        if count_type == Modality.DNA_NORMALIZED:
            filtered = self.normalized_dna_counts.copy()
            layer_name = str(count_type.value)
        elif count_type == Modality.RNA_NORMALIZED:
            layer_name = str(count_type.value)
            filtered = self.normalized_rna_counts.copy()
        elif count_type == Modality.ACTIVITY:
            filtered = self.activity.copy()
            layer_name = str(count_type.value)
        else:
            raise ValueError(f"Unsupported count type: {count_type}")

        filtered[self.barcode_counts < self.barcode_threshold] = np.nan

        return self._correlation(method, filtered, layer_name)

    def _correlation(self, method: str, data: NDArray[np.float32], layer: str) -> NDArray[np.float32]:
        if not self._get_metadata(f"correlation_{layer}"):
            self._compute_correlation(data, layer)
        correlation = self.data.obsp[f"{method}_correlation_{layer}"]
        correlation = np.asarray(correlation, dtype=np.float32)
        return correlation

    def _compute_correlation(self, data: NDArray[np.float32], layer) -> None:

        correlation_methods = ["pearson", "spearman"]

        def compute_correlation(x, y, method: str) -> tuple:
            if method == "spearman":
                return spearmanr(x, y)
            elif method == "pearson":
                return pearsonr(x, y)
            else:
                raise ValueError(f"Unsupported correlation method: {method}")

        # apply var filter to data
        data[self.var_filter.T] = np.nan

        num_columns = self.n_obs

        correlations = {}
        pvalues = {}
        for method in correlation_methods:
            correlations[method] = np.zeros((num_columns, num_columns), dtype=np.float32)
            pvalues[method] = np.zeros((num_columns, num_columns), dtype=np.float32)

        for i in range(num_columns):
            for j in range(i, num_columns):
                mask = ~np.isnan(data[i, :]) & ~np.isnan(data[j, :])
                x = data[i, mask]
                y = data[j, mask]

                for method in correlation_methods:
                    corr, pvalue = compute_correlation(x, y, method)
                    correlations[method][i, j] = corr
                    pvalues[method][i, j] = pvalue
                    correlations[method][j, i] = corr
                    pvalues[method][j, i] = pvalue

        for method in correlation_methods:
            self.data.obsp[f"{method}_correlation_{layer}"] = correlations[method]
            self.data.obsp[f"{method}_correlation_{layer}_pvalue"] = pvalues[method]

        self._add_metadata(f"correlation_{layer}", True)

    def _drop_correlation(self) -> None:

        for layer in ["rna_normalized", "dna_normalized", "activity"]:
            if self._get_metadata(f"correlation_{layer}"):
                self.LOGGER.info(f"Dropping correlation for layer {layer}")

                self._add_metadata(f"correlation_{layer}", False)
                for method in ["pearson", "spearman"]:
                    del self.data.obsp[f"{method}_correlation_{layer}"]
                    del self.data.obsp[f"{method}_correlation_{layer}_pvalue"]

    def _add_metadata(self, key: str, value: Any) -> None:
        if isinstance(value, list):
            if key not in self.data.uns:
                self.data.uns[key] = value
            else:
                self.data.uns[key] = self.data.uns[key] + value
        else:
            self.data.uns[key] = value

    def _get_metadata(self, key) -> Optional[Any]:
        if key in self.data.uns:
            return self.data.uns[key]
        else:
            return None

    def add_sequence_design(self, df_sequence_design: pd.DataFrame, sequence_design_file_path) -> None:
        """Add sequence design metadata to the object's data.

        Args:
            df_sequence_design (pd.DataFrame): DataFrame containing sequence design information, indexed by oligo identifiers.
            sequence_design_file_path (str): Path to the file from which the sequence design data was loaded to store it into the metadata.
        """  # noqa: E501

        df_matched_metadata = df_sequence_design.loc[self.oligos]

        self.data.var["category"] = pd.Categorical(df_matched_metadata["category"])
        self.data.var["class"] = pd.Categorical(df_matched_metadata["class"])
        self.data.var["ref"] = pd.Categorical(df_matched_metadata["ref"])
        self.data.var["chr"] = pd.Categorical(df_matched_metadata["chr"])
        self.data.var["start"] = df_matched_metadata["start"].values
        self.data.var["end"] = df_matched_metadata["end"].values
        self.data.var["strand"] = pd.Categorical(df_matched_metadata["strand"])
        self.data.var["variant_class"] = df_matched_metadata["variant_class"].values
        self.data.var["variant_pos"] = df_matched_metadata["variant_pos"].values
        self.data.var["SPDI"] = df_matched_metadata["SPDI"].values
        self.data.var["allele"] = df_matched_metadata["allele"].values

        self._add_metadata("sequence_design_file", sequence_design_file_path)

    def write(self, file_data_path: os.PathLike) -> None:
        """Writes the AnnData object to a file.

        Args:
            file_data_path (os.PathLike): The path where the AnnData object will be saved.
        """
        self.data.write(file_data_path)
        self.LOGGER.info(f"Data written to {file_data_path}")

    @classmethod
    def read(cls, file_data_path: str) -> "MPRAData":
        """Reads an AnnData object from a file.

        Args:
            file_data_path (str): The path from which the AnnData object will be read.

        Returns:
            MPRAData: An instance of the class containing the data read from the file.
        """
        data = ad.read_h5ad(file_data_path)
        return cls(data)


class MPRABarcodeData(MPRAData):
    """A class for handling barcode-level MPRA (Massively Parallel Reporter Assay) data, providing methods for data import, normalization, filtering, and aggregation to oligo-level data.

    This class extends `MPRAData` and is designed to work with barcode-resolved MPRA datasets, supporting a variety of barcode filtering strategies, normalization routines, and data transformations. It leverages AnnData for data storage and manipulation.

    Note:
        - Filtering and normalization methods are barcode-aware and can be customized via method parameters.
        - Aggregation to oligo-level data is supported for downstream analysis.
    """  # noqa: E501

    @property
    def barcode_counts(self) -> NDArray[np.int32]:
        if "barcode_counts" not in self.data.layers or self.data.layers["barcode_counts"] is None:
            self.barcode_counts = (
                pd.DataFrame(
                    self.observed,  # FIXME make sure var_filter is applied correctly
                    index=self.obs_names,
                    columns=self.var_names,
                )
                .T.groupby(self.oligos, observed=True)
                .transform("sum")
            ).T.values
        return np.asarray(self.data.layers["barcode_counts"], dtype=np.int32)

    @barcode_counts.setter
    def barcode_counts(self, new_data: NDArray[np.int32]) -> None:
        self.data.layers["barcode_counts"] = new_data

    def drop_barcode_counts(self) -> None:

        self.LOGGER.info("Dropping barcode counts")

        self.data.layers.pop("barcode_counts", None)

    @property
    def oligo_data(self) -> "MPRAOligoData":
        """MPRAOligoData: Returns an instance of `MPRAOligoData` containing aggregated oligo-level data."""
        self.LOGGER.info("Computing oligo data")

        return self._oligo_data()

    @classmethod
    def from_file(cls, file_path: str) -> "MPRABarcodeData":

        data = pd.read_csv(file_path, sep="\t", header=0, index_col=0)
        data = data.fillna(0)

        replicate_columns_rna = data.columns[2::2]
        replicate_columns_dna = data.columns[1::2]

        anndata_replicate_rna = data[replicate_columns_rna].transpose().astype(np.int32)
        anndata_replicate_dna = data[replicate_columns_dna].transpose().astype(np.int32)

        anndata_replicate_rna.index = pd.Index([replicate.split("_")[2] for replicate in replicate_columns_rna])
        anndata_replicate_dna.index = pd.Index([replicate.split("_")[2] for replicate in replicate_columns_dna])

        adata = ad.AnnData(anndata_replicate_rna)
        adata.layers["rna"] = np.array(adata.X, dtype=np.int32)
        adata.layers["dna"] = np.asarray(anndata_replicate_dna.values, dtype=np.int32)

        adata.var["oligo"] = data["oligo_name"].astype("category")

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

    def complexity(self, method="lincoln") -> NDArray[np.int64]:
        """Calculates and returns the complexity of barcodes using the Lincoln-Peterson or Chapman estimation.

        Args:
            method (str): Either "lincoln" or "chapman".

        Returns:
           The Lincoln-Peterson or Chapman estimate.
        """

        if method not in {"lincoln", "chapman"}:
            raise ValueError("Method must be either 'lincoln' or 'chapman'.")

        n_observed = np.sum(self.observed, axis=1)
        num_rows = self.observed.shape[0]
        results = np.zeros((num_rows, num_rows), dtype=np.int64)
        for i in range(num_rows):
            for j in range(i, num_rows):
                n_recap = np.sum(np.sum(np.logical_and(self.observed[i, :], self.observed[j, :])))
                if method == "lincoln":
                    count = (n_observed[i] * n_observed[j]) / n_recap if n_recap > 0 else 0
                elif method == "chapman":
                    count = ((n_observed[i] + 1) * (n_observed[j] + 1) / (n_recap + 1)) - 1

                count = int(np.floor(count))  # type: ignore
                results[i, j] = count
                results[j, i] = count  # symmetric

        return results

    def _normalize_layer(
        self,
        counts: NDArray[np.int32],
        total_counts: NDArray[np.int32],
    ) -> NDArray[np.float32]:

        # I do a pseudo count when normalizing to avoid division by zero when computing logfold ratios.
        # barcode filter is already in the observed matrix
        total_counts = total_counts + np.sum(self.pseudo_count * self.observed, axis=1)

        # Avoid division by zero when pseudocount is set to 0
        total_counts[total_counts == 0] = 1
        return (
            ((counts + (self.pseudo_count * self.observed)) / total_counts[:, np.newaxis] * self.scaling) * self.observed
        ).astype(np.float32)

    def _oligo_data(self) -> "MPRAOligoData":

        # Convert the result back to an AnnData object
        oligo_data = ad.AnnData(self._sum_counts_by_oligo(self.rna_counts))

        oligo_data.layers["rna"] = np.array(oligo_data.X)
        oligo_data.layers["dna"] = self._sum_counts_by_oligo(self.dna_counts)

        oligo_data.layers["barcode_counts"] = self._sum_counts_by_oligo(self.observed * ~self.var_filter.T)

        oligo_data.obs_names = self.obs_names.tolist()
        oligo_data.var_names = self.data.var["oligo"].unique().tolist()

        # Subset of vars using the firs occurence of oligo name
        indices = self.data.var["oligo"].drop_duplicates(keep="first").index
        oligo_data.var = self.data.var.loc[indices]

        oligo_data.obs = self.data.obs

        for key, value in self.data.uns.items():
            oligo_data.uns[f"MPRABarcodeData_{key}"] = value

        oligo_data.uns["SCALING"] = self.scaling
        oligo_data.uns["PSEUDOCOUNT"] = self.pseudo_count
        oligo_data.uns["correlation_activity"] = False
        oligo_data.uns["correlation_rna_normalized"] = False
        oligo_data.uns["correlation_dna_normalized"] = False

        return MPRAOligoData(oligo_data, self.barcode_threshold)

    def _sum_counts_by_oligo(self, counts: NDArray) -> pd.DataFrame:

        grouped = pd.DataFrame(
            counts,
            index=self.obs_names,
            columns=self.var_names,
        ).T.groupby(self.data.var["oligo"], observed=True)

        # Perform an operation on each group, e.g., mean
        return grouped.sum().T

    def _barcode_filter_min_bcs_per_oligo(self, threshold: int = 0) -> NDArray[np.bool_]:

        if threshold <= 0:
            return np.full((self.n_vars, self.n_obs), False, dtype=bool)

        mask = (self.barcode_counts < threshold).T

        return mask

    def _get_barcode_mask_for_outlier_filtering(
        self, apply_bc_threshold: bool, aggregated_bc_threshold: bool
    ) -> NDArray[np.bool_]:
        """Generates a boolean mask for barcode filtering based on outlier criteria.

        This method reads data from a specified file (reporter experiment barcode or reporter experiment file format), processes it, and returns an instance of the class containing the data in an AnnData object.

        Args:
            apply_bc_threshold (bool): Whether not only observed barcodes (dna + rna >0, and not filtered), but also a sufficient number of barcodes per oligo should be present.
            aggregated_bc_threshold (bool): If True, barcode threshold is used for total observed barcodes per oligo across replicates; if False, applies threshold to observed barcodes per oligo per replicate.

        Returns:
            Boolean mask indicating which barcodes pass the filtering criteria.
        """  # noqa: E501

        barcode_mask = self.observed.T

        if apply_bc_threshold:
            if aggregated_bc_threshold:
                aggregated_barcode_counts = (
                    pd.DataFrame(self.observed.any(axis=0).reshape(-1, 1), index=self.var_names)
                    .groupby(self.oligos, observed=True)
                    .transform("sum")
                    .values
                )
                barcode_mask = barcode_mask * np.tile(aggregated_barcode_counts, (1, self.n_obs)) >= self.barcode_threshold
            else:
                barcode_mask = barcode_mask * (self.barcode_counts.T >= self.barcode_threshold)

        return barcode_mask

    def _barcode_filter_global_outliers(
        self, times_zscore: float = 3.0, apply_bc_threshold: bool = False, aggregated_bc_threshold: bool = False
    ) -> NDArray[np.bool_]:

        barcode_mask = self._get_barcode_mask_for_outlier_filtering(apply_bc_threshold, aggregated_bc_threshold)

        df_rna = pd.DataFrame(self.rna_counts, index=self.obs_names, columns=self.var_names).T.where(barcode_mask)
        mask = ((df_rna - df_rna.mean(axis=0)) / df_rna.std(axis=0)).abs() > times_zscore

        return mask.values.astype(bool)

    def _barcode_filter_oligo_specific_outliers(
        self, times_zscore: float = 3.0, apply_bc_threshold: bool = False, aggregated_bc_threshold: bool = False
    ) -> NDArray[np.bool_]:

        barcode_mask = self._get_barcode_mask_for_outlier_filtering(apply_bc_threshold, aggregated_bc_threshold)

        df_rna = pd.DataFrame(self.rna_counts, index=self.obs_names, columns=self.var_names).T.where(barcode_mask)
        grouped = df_rna.groupby(self.oligos, observed=True)
        mask = ((df_rna - grouped.transform("mean")) / grouped.transform("std").fillna(0).replace(0, 1)).abs() > times_zscore

        return mask.values.astype(bool)

    def _barcode_filter_large_expression_outliers(
        self, times_activity: float = 5.0, apply_bc_threshold: bool = False, aggregated_bc_threshold: bool = False
    ) -> NDArray[np.bool_]:

        ratio = np.divide(
            self.normalized_rna_counts.sum(axis=0),
            self.normalized_dna_counts.sum(axis=0),
            where=self.normalized_dna_counts.sum(axis=0) != 0,
        )
        with np.errstate(divide="ignore"):
            log2ratio = np.log2(ratio)
            log2ratio[np.isneginf(log2ratio)] = np.nan

        barcode_mask = (
            self._get_barcode_mask_for_outlier_filtering(apply_bc_threshold, aggregated_bc_threshold)
            .any(axis=1)
            .reshape(-1, 1)
        )

        log2ratio = pd.DataFrame({"ratio": log2ratio}, index=self.var_names).where(barcode_mask)
        log2ratio["oligo_median"] = log2ratio.groupby(self.oligos, observed=True)["ratio"].transform("median")

        expr_diff = log2ratio["ratio"] - log2ratio["oligo_median"]  # calculate difference from oligo median
        mask = (expr_diff > times_activity).values.astype(bool)  # numpy boolean array

        mask_array = np.tile(mask[:, np.newaxis], (1, self.n_obs))
        return mask_array

    def _barcode_filter_mad(self, times_mad=3, n_bins=20) -> NDArray[np.bool_]:
        # sum up DNA and RNA counts across replicates
        DNA_sum = pd.DataFrame(self.raw_dna_counts, index=self.obs_names, columns=self.var_names).T.sum(axis=1)
        RNA_sum = pd.DataFrame(self.raw_rna_counts, index=self.obs_names, columns=self.var_names).T.sum(axis=1)
        df_sums = pd.DataFrame({"DNA_sum": DNA_sum, "RNA_sum": RNA_sum}).fillna(0)
        # removing all barcodes with 0 counts in RNA and more DNA count than number of replicates/observations
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

    def _barcode_filter_random(
        self, proportion: float = 1.0, total: Optional[int] = None, aggegate_over_replicates: bool = True
    ) -> NDArray[np.bool_]:

        if aggegate_over_replicates and total is None:
            total = self.var_filter.shape[0]
        elif total is None:
            total = self.var_filter.size

        num_true_cells = int(total * (1.0 - proportion))  # type: ignore
        true_indices = np.random.choice(total, num_true_cells, replace=False)

        mask = np.full((self.data.n_vars, self.data.n_obs), False)

        if aggegate_over_replicates:
            mask[true_indices, :] = True
        else:
            flat_df = mask.flatten()

            flat_df[true_indices] = True

            mask = flat_df.reshape(self.var_filter.shape)

        return mask

    def _barcode_filter_min_count(
        self, rna_min_count: Optional[int] = None, dna_min_count: Optional[int] = None
    ) -> NDArray[np.bool_]:

        return self._barcode_filter_min_max_count(BarcodeFilter.MIN_COUNT, rna_min_count, dna_min_count)

    def _barcode_filter_max_count(
        self, rna_max_count: Optional[int] = None, dna_max_count: Optional[int] = None
    ) -> NDArray[np.bool_]:

        return self._barcode_filter_min_max_count(BarcodeFilter.MAX_COUNT, rna_max_count, dna_max_count)

    def _barcode_filter_min_max_counts(
        self, barcode_filter: BarcodeFilter, counts: NDArray[np.int32], count_threshold: int
    ) -> NDArray[np.bool_]:
        if barcode_filter == BarcodeFilter.MIN_COUNT:
            return (counts < count_threshold).T
        elif barcode_filter == BarcodeFilter.MAX_COUNT:
            return (counts > count_threshold).T
        else:
            return np.full((self.n_vars, self.n_obs), False)

    def _barcode_filter_min_max_count(
        self, barcode_filter: BarcodeFilter, rna_count: Optional[int] = None, dna_count: Optional[int] = None
    ) -> NDArray[np.bool_]:
        mask = np.full((self.n_vars, self.n_obs), False)
        if rna_count is not None:
            mask = mask | self._barcode_filter_min_max_counts(barcode_filter, self.raw_rna_counts, rna_count)
        if dna_count is not None:
            mask = mask | self._barcode_filter_min_max_counts(barcode_filter, self.raw_dna_counts, dna_count)

        return mask

    def apply_barcode_filter(self, barcode_filter: BarcodeFilter, params: dict = {}) -> None:
        """Applies a specified barcode filter to the dataset using the provided parameters.

        This method selects the appropriate barcode filtering function based on the `barcode_filter` argument and applies it to update the `var_filter` attribute. Supported filters include RNA z-score, MAD, random, minimum count, and maximum count. After applying the filter, metadata is updated to record the applied filter.

        Args:
            barcode_filter (BarcodeFilter): The type of barcode filter to apply.
            params (dict, optional): Additional parameters to pass to the filter function. Defaults to an empty dictionary.

        Raises:
            ValueError: If an unsupported barcode filter is provided.
        """  # noqa: E501

        filter_switch: dict[BarcodeFilter, Callable[..., NDArray[np.bool_]]] = {
            BarcodeFilter.MIN_BCS_PER_OLIGO: self._barcode_filter_min_bcs_per_oligo,
            BarcodeFilter.GLOBAL: self._barcode_filter_global_outliers,
            BarcodeFilter.LARGE_EXPRESSION: self._barcode_filter_large_expression_outliers,
            BarcodeFilter.OLIGO_SPECIFIC: self._barcode_filter_oligo_specific_outliers,
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

    def drop_count_sampling(self) -> None:
        """Removes count sampling data from the dataset.

        This method performs the following actions:
        - Calls `drop_normalized()` to remove any normalized data.
        - Logs the action of dropping count sampling.
        - Deletes the "count_sampling" entry from the `.uns` attribute of the data.
        - Removes "rna_sampling" and "dna_sampling" layers from the data, if they exist.
        """

        self.drop_normalized()
        self.drop_barcode_counts()
        self.drop_total_counts()

        self.LOGGER.info("Dropping count sampling")

        del self.data.uns["count_sampling"]
        self.data.layers.pop("rna_sampling", None)
        self.data.layers.pop("dna_sampling", None)

    def _calculate_proportions(
        self,
        proportion: Optional[float],
        total: Optional[int],
        aggregate_over_replicates: bool,
        counts: NDArray[np.int32],
        replicates: int,
    ):
        pp = [1.0] * replicates

        if proportion is not None:
            pp = [proportion] * replicates

        if total is not None:
            if aggregate_over_replicates:
                for i, p in enumerate(pp):
                    pp[i] = min(float(total) / float(np.sum(counts)), p)
            else:
                for i, p in enumerate(pp):
                    pp[i] = min(float(total) / float(np.sum(counts[i, :])), p)
        return pp

    def _sample_individual_counts(self, x, proportion: float) -> int:
        return int(
            np.floor(x * proportion)
            + (0.0 if x != 0 or np.random.rand() > (x * proportion - np.floor(x * proportion)) else 1.0)
        )

    def _apply_sampling(
        self,
        layer_name: str,
        counts: NDArray[np.int32],
        proportion: Optional[float],
        total: Optional[int],
        max_value: Optional[int],
        aggregate_over_replicates: bool,
    ) -> None:

        sampled_counts = counts.copy()

        if total is not None or proportion is not None:

            pp = self._calculate_proportions(proportion, total, aggregate_over_replicates, sampled_counts, self.n_obs)

            vectorized_sample_individual_counts = np.vectorize(self._sample_individual_counts)

            for i, p in enumerate(pp):
                sampled_counts[i, :] = vectorized_sample_individual_counts(sampled_counts[i, :], proportion=p)

        if max_value is not None:
            sampled_counts = np.clip(sampled_counts, None, max_value)

        self.data.layers[layer_name] = sampled_counts

    def apply_count_sampling(
        self,
        count_type: CountSampling,
        proportion: Optional[float] = None,
        total: Optional[int] = None,
        max_value: Optional[int] = None,
        aggregate_over_replicates: bool = False,
    ) -> None:
        """Applies count sampling to RNA and/or DNA count data according to the specified parameters.

        Args:
            count_type (CountSampling): Specifies which counts to sample. Options are RNA, DNA, or RNA_AND_DNA.
            proportion (Optional[float]): Proportion of counts to sample (between 0 and 1). If None, this parameter is ignored.
            total (Optional[int]): Total number of counts to sample. If None, this parameter is ignored.
            max_value (Optional[int]): Maximum value for sampled counts. If None, this parameter is ignored.
            aggregate_over_replicates (bool): Whether to aggregate counts over replicates before sampling.

        Side Effects:
            - Adds sampling metadata to the object.
            - Drops any normalized data associated with the object.
        """

        if count_type == CountSampling.RNA or count_type == CountSampling.RNA_AND_DNA:
            self._apply_sampling("rna_sampling", self.raw_rna_counts, proportion, total, max_value, aggregate_over_replicates)

        if count_type == CountSampling.DNA or count_type == CountSampling.RNA_AND_DNA:
            self._apply_sampling("dna_sampling", self.raw_dna_counts, proportion, total, max_value, aggregate_over_replicates)

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
        self.drop_total_counts()
        self.drop_normalized()
        self.drop_barcode_counts()


class MPRAOligoData(MPRAData):
    """MPRAOligoData is a subclass of MPRAData designed to handle MPRA (Massively Parallel Reporter Assay) oligo-level data.

    This class provides methods for loading, normalizing, and managing barcode counts and associated data layers for MPRA experiments. Barcode counts must be pre-set before accessing, as they cannot be computed within this class. The normalization process includes pseudocount handling to avoid division by zero and supports per-barcode normalization.

    Raises:
        MPRAlibException: If barcode counts are not set when accessed.
    """  # noqa: E501

    @property
    def barcode_counts(self) -> NDArray[np.int32]:
        if "barcode_counts" not in self.data.layers or self.data.layers["barcode_counts"] is None:
            raise MPRAlibException(
                "Barcode counts are not set in MPRAOligoData and cannot be computed. They have to be pre-set before accessing."
            )
        return np.asarray(self.data.layers["barcode_counts"], dtype=np.int32)

    @barcode_counts.setter
    def barcode_counts(self, new_data: NDArray[np.int32]) -> None:
        self.data.layers["barcode_counts"] = new_data

    def drop_barcode_counts(self):
        pass

    @classmethod
    def from_file(cls, file_path: str) -> "MPRAOligoData":

        return MPRAOligoData(ad.read_h5ad(file_path))

    def _normalize_layer(
        self,
        counts: NDArray[np.int32],
        total_counts: NDArray[np.int32],
    ) -> NDArray[np.float32]:

        # I do a pseudo count when normalizing to avoid division by zero when computing logfold ratios.
        # Pseudocount has also be done per barcode.
        # var filter has to be used again because we want to have a zero on filtered values.
        total_counts = total_counts + np.sum((self.pseudo_count * self.barcode_counts) * self.observed, axis=1)

        # Avoid division by zero when pseudocount is set to 0
        total_counts[total_counts == 0] = 1
        scaled_counts = (counts + (self.pseudo_count * self.barcode_counts)) / total_counts[:, np.newaxis] * self.scaling
        return np.divide(scaled_counts, self.barcode_counts, where=self.barcode_counts != 0, out=np.zeros_like(scaled_counts))
