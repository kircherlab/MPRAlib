import pandas as pd
import numpy as np
import os
from mpralib.mpradata import MPRAdata


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


def export_activity_file(mpradata: MPRAdata, output_file_path: str) -> None:

    output = pd.DataFrame()

    data = mpradata.grouped_data

    for replicate in data.obs_names:
        replicate_data = data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcodes"] != 0]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var["oligo"],
            "dna_counts": replicate_data.layers["dna"][0, :],
            "rna_counts": replicate_data.layers["rna"][0, :],
            "dna_normalized": np.round(replicate_data.layers["dna_normalized"][0, :], 4),
            "rna_normalized": np.round(replicate_data.layers["rna_normalized"][0, :], 4),
            "log2FoldChange": np.round(replicate_data.layers["log2FoldChange"][0, :], 4),
            "n_bc": replicate_data.layers["barcodes"][0, :],
        }
        output = pd.concat([output, pd.DataFrame(df)], axis=0)

    output.to_csv(output_file_path, sep="\t", index=False)


def export_barcode_activity_file(mpradata: MPRAdata, output_file_path: str) -> None:

    output = pd.DataFrame({
            "barcode": mpradata.barcodes,
            "oligo_name": mpradata.oligos
        })

    for replicate in mpradata.replicates:
        replicate_data = mpradata.data[replicate, :]

        output[f"dna_count_{replicate}"] = replicate_data.layers["dna"][0]
        output[f"rna_count_{replicate}"] = replicate_data.layers["rna"][0]

    output.to_csv(output_file_path, sep="\t", index=False)
