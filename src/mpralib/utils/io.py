import pandas as pd
import numpy as np
import ast
import os
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData, MPRAData
from mpralib.exception import SequenceDesignException


def chromosome_map() -> pd.DataFrame:
    base_path = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_path, "data", "hg19.chromAlias.txt")
    df = pd.read_csv(file_path, sep="\t", header=None, comment="#", dtype="category")
    file_path = os.path.join(base_path, "data", "hg38.chromAlias.txt")
    df = pd.concat(
        [
            df,
            pd.read_csv(file_path, sep="\t", header=None, comment="#", dtype="category"),
        ],
        ignore_index=True,
    )
    df.columns = ["ucsc", "assembly", "genbank", "refseq"]
    return df


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
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpradata.barcode_threshold]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var["oligo"],
            "dna_counts": replicate_data.layers["dna"][0, :],
            "rna_counts": replicate_data.layers["rna"][0, :],
            "dna_normalized": np.round(replicate_data.layers["dna_normalized"][0, :], 4),
            "rna_normalized": np.round(replicate_data.layers["rna_normalized"][0, :], 4),
            "log2FoldChange": np.round(replicate_data.layers["activity"][0, :], 4),
            "n_bc": replicate_data.layers["barcode_counts"][0, :],
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


def export_counts_file(mpradata: MPRAData, output_file_path: str) -> None:
    if isinstance(mpradata, MPRAOligoData):
        df = {"ID": mpradata.oligos}
    else:
        df = {"ID": mpradata.var_names, "name": mpradata.oligos}
    dna_counts = mpradata.dna_counts
    dna_counts[mpradata.barcode_counts < mpradata.barcode_threshold] = 0
    rna_counts = mpradata.rna_counts
    rna_counts[mpradata.barcode_counts < mpradata.barcode_threshold] = 0
    for idx, replicate in enumerate(mpradata.obs_names):
        df["dna_count_" + replicate] = dna_counts[idx, :]
        df["rna_count_" + replicate] = rna_counts[idx, :]

    df = pd.DataFrame(df).set_index("ID")
    # remove IDs which are all zero
    df = df[(df.T != 0).all()]

    df.to_csv(output_file_path, sep="\t", index=True)


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
