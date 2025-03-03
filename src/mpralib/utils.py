import pandas as pd
import numpy as np
import os
from mpralib.mpradata import MPRABarcodeData, MPRAOligoData


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
            "log2FoldChange": np.round(replicate_data.layers["log2FoldChange"][0, :], 4),
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
