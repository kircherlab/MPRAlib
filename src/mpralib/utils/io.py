import pandas as pd
import numpy as np
import ast
import os
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, MPRAData
from mpralib.exception import SequenceDesignException, MPRAlibException


def chromosome_map() -> pd.DataFrame:
    base_path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_path, "../data", "hg19.chromAlias.txt")
    df = pd.read_csv(file_path, sep="\t", header=None, comment="#", dtype="category")
    file_path = os.path.join(base_path, "../data", "hg38.chromAlias.txt")
    df = pd.concat(
        [
            df,
            pd.read_csv(file_path, sep="\t", header=None, comment="#", dtype="category"),
        ],
        ignore_index=True,
    )
    df.columns = ["ucsc", "assembly", "genbank", "refseq"]
    return df


def is_compressed_file(filepath):
    """
    Check if a file is compressed (gzip or bgz).

    Parameters:
        filepath (str): Path to the file to check.

    Returns:
        bool: True if the file is compressed, False otherwise.
    """
    return is_gzip_file(filepath) or is_bgzf(filepath)


def is_gzip_file(filepath):
    """
    Check if a file is a gzip-compressed file based on its magic number.

    Parameters:
        filepath (str or Path): Path to the file to check.

    Returns:
        bool: True if the file is gzip-compressed, False otherwise.
    """
    with open(filepath, "rb") as f:
        magic = f.read(2)
    return magic == b"\x1f\x8b"


def is_bgzf(filepath):
    """
    Check if a file is in BGZF (Blocked GNU Zip Format) format.

    BGZF is a variant of the standard gzip format with extra fields that allow for random access.
    This function reads the file header and checks for the BGZF-specific magic numbers and flags.

    Parameters:
        filepath (str): Path to the file to be checked.

    Returns:
        bool: True if the file is in BGZF format, False otherwise.
    """
    with open(filepath, "rb") as f:
        header = f.read(18)
    return (
        len(header) >= 18
        and header[0:2] == b"\x1f\x8b"  # gzip magic number
        and header[3] == 4  # FLG.FEXTRA set
        and header[12:14] == b"BC"  # BGZF extra subfield
    )


def export_activity_file(mpradata: MPRAOligoData, output_file_path: str) -> None:
    """
    Export activity data from an MPRAdata object to a tab-separated values (TSV) file.

    The function processes the grouped data from the MPRAdata object, extracts relevant information
    for each replicate, and writes the data to a TSV file. The output file contains columns for
    replicate, oligo name, DNA counts, RNA counts, normalized DNA counts, normalized RNA counts,
    log2 fold change, and the number of barcodes. Barcode filters, count sampling
    and barcode thresholds are applied.

    Parameters:
        mpradata (MPRAdata): An object containing MPRA (Massively Parallel Reporter Assay) data.
        output_file_path (str): The file path where the output TSV file will be saved.

    Returns:
        None
    """

    output = pd.DataFrame()

    mpradata.activity

    for replicate in mpradata.obs_names:
        replicate_data = mpradata.data[replicate, :]
        replicate_data = replicate_data[
            :, np.asarray(replicate_data.layers["barcode_counts"]) >= np.asarray(mpradata.barcode_threshold)
        ]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var["oligo"],
            "dna_counts": np.asarray(replicate_data.layers["dna"])[0, :],
            "rna_counts": np.asarray(replicate_data.layers["rna"])[0, :],
            "dna_normalized": np.round(np.asarray(replicate_data.layers["dna_normalized"])[0, :], 4),
            "rna_normalized": np.round(np.asarray(replicate_data.layers["rna_normalized"])[0, :], 4),
            "log2FoldChange": np.round(np.asarray(replicate_data.layers["activity"])[0, :], 4),
            "n_bc": np.asarray(replicate_data.layers["barcode_counts"])[0, :],
        }
        output = pd.concat([output, pd.DataFrame(df)], axis=0)

    output.to_csv(output_file_path, sep="\t", index=False)


def export_barcode_file(mpradata: MPRABarcodeData, output_file_path: str) -> None:
    """
    Export barcode count data to a file.

    This function takes an MPRAdata object and exports its barcode count data
    to a specified file path in tab-separated values (TSV) format. The output file
    will contain columns for barcodes, oligo names, and DNA/RNA counts for each replicate.
    Modifides counts (barcode filter/sampling) if applicable will be written.

    Parameters:
        mpradata (MPRAdata): An object containing MPRA data, including barcodes, oligos,
                         DNA counts, RNA counts, and replicates.
    output_file_path (str): The file path where the output TSV file will be saved.

    Returns:
    None
    """

    output = pd.DataFrame({"barcode": mpradata.var_names, "oligo_name": mpradata.oligos})

    mpradata.normalized_rna_counts

    dna_counts = mpradata.dna_counts
    rna_counts = mpradata.rna_counts
    for i, replicate in enumerate(mpradata.obs_names):
        output[f"dna_count_{replicate}"] = dna_counts[i]
        output[f"rna_count_{replicate}"] = rna_counts[i]
    output.replace(0, "", inplace=True)
    output.to_csv(output_file_path, sep="\t", index=False)


def export_counts_file(
    mpradata: MPRAData, output_file_path: str, normalized: bool = False, filter: np.ndarray | None = None
) -> None:
    if isinstance(mpradata, MPRABarcodeData):
        df = {"ID": mpradata.var_names}
    elif isinstance(mpradata, MPRAOligoData):
        df = {"ID": mpradata.oligos}
    else:
        raise MPRAlibException(f"Invalid MPRA data type {type(mpradata)}. Expected MPRAOligoData or MPRABarcodeData.")

    if normalized:
        dna_counts = mpradata.normalized_dna_counts.copy()
        rna_counts = mpradata.normalized_rna_counts.copy()
    else:
        dna_counts = mpradata.dna_counts.copy()
        rna_counts = mpradata.rna_counts.copy()

    not_observed_mask = ~mpradata.observed
    lower_bc_mask = mpradata.barcode_counts < mpradata.barcode_threshold
    mask = [not_observed_mask, lower_bc_mask]
    if filter is not None:
        mask.append(filter)

    dna_counts = np.ma.masked_array(dna_counts, mask=np.any(mask, axis=0))
    rna_counts = np.ma.masked_array(rna_counts, mask=np.any(mask, axis=0))

    for idx, replicate in enumerate(mpradata.obs_names):
        df["dna_count_" + replicate] = dna_counts[idx, :]
        df["rna_count_" + replicate] = rna_counts[idx, :]

    df = pd.DataFrame(df).set_index("ID")
    df = df.map(lambda x: np.nan if isinstance(x, np.ma.core.MaskedConstant) else x)
    nan_columns = df.isna().all(axis=1)

    if isinstance(mpradata, MPRABarcodeData):
        df.insert(0, "name", mpradata.oligos)

    # remove IDs which are all zero
    df = df[~nan_columns]

    if normalized:
        df.to_csv(output_file_path, sep="\t", index=True, na_rep="", float_format="%.6f")
    else:
        df.to_csv(output_file_path, sep="\t", index=True, na_rep="", float_format="%.0f")


def read_sequence_design_file(file_path: str) -> pd.DataFrame:
    """
    Read sequence design from a tab-separated values (TSV) file.

    This function reads metadata from a TSV file and returns it as a pandas DataFrame.
    The metadata file should contain columns for sample ID, replicate, and any additional
    metadata. The sample ID should correspond to the oligo name in the MPRA data object.

    Parameters:
        file_path (str): The file path of the metadata TSV file.

    Returns:
        pd.DataFrame: A DataFrame containing the metadata.
    """

    df = pd.read_csv(
        file_path,
        sep="\t",
        header=0,
        na_values=["NA"],
        usecols=[
            "name",
            "sequence",
            "category",
            "class",
            "source",
            "ref",
            "chr",
            "start",
            "end",
            "strand",
            "variant_class",
            "variant_pos",
            "SPDI",
            "allele",
            "info",
        ],
    ).drop_duplicates()

    # Set specific columns as arrays
    df["variant_class"] = df["variant_class"].fillna("[]").apply(ast.literal_eval)
    df["variant_pos"] = df["variant_pos"].fillna("[]").apply(ast.literal_eval)
    df["SPDI"] = df["SPDI"].fillna("[]").apply(ast.literal_eval)
    df["allele"] = df["allele"].fillna("[]").apply(ast.literal_eval)

    # Set specific columns as categorical or integer types
    df["category"] = pd.Categorical(df["category"])
    df["class"] = pd.Categorical(df["class"])
    df["chr"] = pd.Categorical(df["chr"])
    df["strand"] = pd.Categorical(df["strand"])
    df["start"] = df["start"].astype("Int64")  # Nullable integer type
    df["end"] = df["end"].astype("Int64")  # Nullable integer type
    df["name"] = pd.Categorical(df["name"].str.replace(r"[\s\[\]]", "_", regex=True))

    # oligo name as index
    df.set_index("name", inplace=True)

    # Validate that the 'sequence' column contains only valid DNA characters
    if np.any(~df["sequence"].str.match(r"^[ATGCatgc]+$", na=False)):
        raise SequenceDesignException("sequence", file_path)

    # Validate that the 'class' column contains only 'variant', 'element', 'synthetic' or 'scrambled'
    valid_categories = {"variant", "element", "synthetic", "scrambled"}
    if not set(df["category"].cat.categories).issubset(valid_categories):
        raise SequenceDesignException("category", file_path)

    # Validate that the 'category' column contains only:
    valid_classes = {
        "test",
        "variant positive control",
        "variant negative control",
        "element active control",
        "element inactive control",
    }
    if not set(df["class"].cat.categories).issubset(valid_classes):
        raise SequenceDesignException("class", file_path)

    return df
