import click
import pandas as pd
import numpy as np
import math
import pysam
from sklearn.preprocessing import MinMaxScaler
from mpralib.mpradata import MPRABarcodeData, BarcodeFilter
from mpralib.utils import chromosome_map, export_activity_file, export_barcode_file


@click.group(help="Command line interface of MPRAlib, a library for MPRA data analysis.")
def cli():
    pass


@cli.command(help="Generating element activity or barcode count files.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA result in a barcode format.",
)
@click.option(
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output (element level only).",
)
@click.option(
    "--element-level/--barcode-level",
    "element_level",
    default=True,
    help="Export activity at the element (default) or barcode level.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def activities(input_file, bc_threshold, element_level, output_file):
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    if element_level:
        export_activity_file(mpradata.oligo_data, output_file)
    else:
        export_barcode_file(mpradata, output_file)


@cli.command(help="Generating correlation of lo2 folds on barcode threshold.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA result in a barcode format.",
)
@click.option(
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output (element level only).",
)
def correlation(input_file, bc_threshold):
    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.barcode_threshold = bc_threshold

    print(mpradata.correlation("pearson", "activity"))

    print(mpradata.correlation("pearson", "dna"))

    print(mpradata.correlation("pearson", "rna"))


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--outlier-rna-zscore-times",
    "rna_zscore_times",
    default=3,
    type=float,
    help="Absolute rna z_score is not allowed to be larger than this value.",
)
@click.option(
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def filter_outliers(input_file, rna_zscore_times, bc_threshold, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    # mpradata.apply_barcode_filter(OutlierFilter.MAD, {})

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": rna_zscore_times})

    print(mpradata.spearman_correlation_activity)
    print(mpradata.pearson_correlation_activity)

    data = mpradata.oligo_data

    print(data.layers["barcodes"].sum())
    print((data.layers["barcodes"] == 0).sum())

    output = pd.DataFrame()

    for replicate in data.obs["replicate"]:
        replicate_data = data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcodes"] >= bc_threshold]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var_names.values,
            "dna_counts": replicate_data.layers["dna"][0, :],
            "rna_counts": replicate_data.layers["rna"][0, :],
            "dna_normalized": np.round(replicate_data.layers["dna_normalized"][0, :], 4),
            "rna_normalized": np.round(replicate_data.layers["rna_normalized"][0, :], 4),
            "log2FoldChange": np.round(replicate_data.layers["log2FoldChange"][0, :], 4),
            "n_bc": replicate_data.layers["barcodes"][0, :],
        }
        output = pd.concat([output, pd.DataFrame(df)], axis=0)

    output.to_csv(output_file, sep="\t", index=False)


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def get_variant_map(input_file, metadata_file, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    variant_map = mpradata.variant_map
    for key in ["REF", "ALT"]:
        variant_map[key] = [",".join(i) for i in variant_map[key]]

    variant_map.to_csv(output_file, sep="\t", index=True)


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--output-dna",
    "output_dna_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of dna counts.",
)
@click.option(
    "--output-rna",
    "output_rna_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of rna counts.",
)
@click.option(
    "--output-mpra-data",
    "output_mpra_data_file",
    required=False,
    type=(click.Path(writable=True), click.Path(writable=True)),
    help="Output file of MPRA data object.",
)
def get_element_counts(input_file, metadata_file, output_dna_file, output_rna_file, output_mpra_data_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")

    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    click.echo("After correlate filtering:")
    click.echo(mpradata.pearson_correlation_activity)

    mpradata.element_dna_counts.to_csv(output_dna_file, sep="\t", index=True)
    mpradata.element_rna_counts.to_csv(output_rna_file, sep="\t", index=True)

    if output_mpra_data_file:
        mpradata.write(output_mpra_data_file[0], output_mpra_data_file[1])


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--output-dna",
    "output_dna_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of dna counts.",
)
@click.option(
    "--output-rna",
    "output_rna_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of rna counts.",
)
def get_variant_counts(input_file, metadata_file, output_dna_file, output_rna_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    click.echo("Initial Pearson correlation:")
    click.echo(mpradata.pearson_correlation_activity)

    mpradata.barcode_threshold = 10

    click.echo("After BC-threshold filtering:")
    click.echo(mpradata.pearson_correlation_activity)

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    click.echo("After ZSCORE filtering:")
    click.echo(mpradata.pearson_correlation_activity)

    mpradata.variant_dna_counts.to_csv(output_dna_file, sep="\t", index=True)
    mpradata.variant_rna_counts.to_csv(output_rna_file, sep="\t", index=True)


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--mpralm",
    "mpralm_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA lm file.",
)
@click.option(
    "--output-reporter-elements",
    "output_reporter_elements_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_elements(input_file, metadata_file, mpralm_file, output_reporter_elements_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")

    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    df = pd.read_csv(mpralm_file, sep="\t", header=0)

    indexes_in_order = [
        mpradata.oligo_data.var["oligo"][mpradata.oligo_data.var["oligo"] == ID].index.tolist() for ID in df["ID"]
    ]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    df.index = indexes_in_order

    df = df.join(mpradata.oligo_data.var["oligo"], how="right")

    mpradata.oligo_data.varm["mpralm_element"] = df

    out_df = mpradata.oligo_data.varm["mpralm_element"][["oligo", "logFC", "P.Value", "adj.P.Val"]]
    out_df["inputCount"] = mpradata.oligo_data.layers["dna_normalized"].mean(axis=0)
    out_df["outputCount"] = mpradata.oligo_data.layers["rna_normalized"].mean(axis=0)
    out_df.dropna(inplace=True)
    out_df["minusLog10PValue"] = -np.log10(out_df["P.Value"])
    out_df["minusLog10QValue"] = -np.log10(out_df["adj.P.Val"])
    out_df.rename(columns={"oligo": "oligo_name", "logFC": "log2FoldChange"}, inplace=True)
    out_df[
        [
            "oligo_name",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    ].to_csv(output_reporter_elements_file, sep="\t", index=False)


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--mpralm",
    "mpralm_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA lm file.",
)
@click.option(
    "--output-reporter-variants",
    "output_reporter_variants_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_variants(input_file, metadata_file, mpralm_file, output_reporter_variants_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    spdi_map = mpradata.variant_map

    df = pd.read_csv(mpralm_file, sep="\t", header=0, index_col=0)

    variant_dna_counts = mpradata.variant_dna_counts
    variant_rna_counts = mpradata.variant_rna_counts

    columns_ref = []
    columns_alt = []
    for replicate in mpradata.oligo_data.obs_names:
        columns_ref.append("counts_" + replicate + "_REF")
        columns_alt.append("counts_" + replicate + "_ALT")

    dna_counts = np.array(variant_dna_counts[columns_ref].sum()) + np.array(variant_dna_counts[columns_alt].sum())
    rna_counts = np.array(variant_rna_counts[columns_ref].sum()) + np.array(variant_rna_counts[columns_alt].sum())

    for spdi, row in spdi_map.iterrows():
        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = (
                (variant_dna_counts.loc[spdi][columns_ref] / dna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "inputCountAlt"] = (
                (variant_dna_counts.loc[spdi][columns_alt] / dna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "outputCountRef"] = (
                (variant_rna_counts.loc[spdi][columns_ref] / rna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "outputCountAlt"] = (
                (variant_rna_counts.loc[spdi][columns_alt] / rna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "variantPos"] = int(
                mpradata.oligo_data.var["variant_pos"][mpradata.oligos.isin(row["REF"])].values[0][0]
            )

    df["minusLog10PValue"] = -np.log10(df["P.Value"])
    df["minusLog10QValue"] = -np.log10(df["adj.P.Val"])
    df.rename(
        columns={
            "CI.L": "CI_lower_95",
            "CI.R": "CI_upper_95",
            "logFC": "log2FoldChange",
        },
        inplace=True,
    )
    df["postProbEffect"] = df["B"].apply(lambda x: math.exp(x) / (1 + math.exp(x)))
    df["variant_id"] = df.index
    df["refAllele"] = df["variant_id"].apply(lambda x: x.split(":")[2])
    df["altAllele"] = df["variant_id"].apply(lambda x: x.split(":")[3])
    df["variantPos"] = df["variantPos"].astype(int)

    df[
        [
            "variant_id",
            "log2FoldChange",
            "inputCountRef",
            "outputCountRef",
            "inputCountAlt",
            "outputCountAlt",
            "minusLog10PValue",
            "minusLog10QValue",
            "postProbEffect",
            "CI_lower_95",
            "CI_upper_95",
            "variantPos",
            "refAllele",
            "altAllele",
        ]
    ].to_csv(output_reporter_variants_file, sep="\t", index=False)


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--mpralm",
    "mpralm_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA lm file.",
)
@click.option(
    "--output-reporter-genomic-elements",
    "output_reporter_genomic_elements_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_genomic_elements(input_file, metadata_file, mpralm_file, output_reporter_genomic_elements_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")

    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    df = pd.read_csv(mpralm_file, sep="\t", header=0)

    indexes_in_order = [
        mpradata.oligo_data.var["oligo"][mpradata.oligo_data.var["oligo"] == ID].index.tolist() for ID in df["ID"]
    ]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    df.index = indexes_in_order

    df = df.join(mpradata.oligo_data.var["oligo"], how="right")

    mpradata.oligo_data.varm["mpralm_element"] = df

    out_df = mpradata.oligo_data.varm["mpralm_element"][["oligo", "logFC", "P.Value", "adj.P.Val"]]
    out_df.loc[:, ["inputCount"]] = mpradata.oligo_data.layers["dna_normalized"].mean(axis=0)
    out_df.loc[:, ["outputCount"]] = mpradata.oligo_data.layers["rna_normalized"].mean(axis=0)
    print(out_df)
    out_df["chr"] = mpradata.oligo_data.var["chr"]
    out_df["start"] = mpradata.oligo_data.var["start"]
    out_df["end"] = mpradata.oligo_data.var["end"]
    out_df["strand"] = mpradata.oligo_data.var["strand"]
    out_df.dropna(inplace=True)
    scaler = MinMaxScaler(feature_range=(0, 1000))
    out_df["score"] = scaler.fit_transform(out_df[["logFC"]]).astype(int)
    out_df["minusLog10PValue"] = -np.log10(out_df["P.Value"])
    out_df["minusLog10QValue"] = -np.log10(out_df["adj.P.Val"])
    out_df.rename(columns={"oligo": "name", "logFC": "log2FoldChange"}, inplace=True)
    out_df = out_df[
        [
            "chr",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    ].sort_values(by=["chr", "start", "end"])

    with pysam.BGZFile(output_reporter_genomic_elements_file, "ab") as f:
        f.write(out_df.to_csv(sep="\t", index=False, header=False).encode())


@cli.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--metadata",
    "metadata_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--mpralm",
    "mpralm_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA lm file.",
)
@click.option(
    "--output-reporter-genomic-variants",
    "output_reporter_genomic_variants_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_genomic_variants(input_file, metadata_file, mpralm_file, output_reporter_genomic_variants_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": 3})

    spdi_map = mpradata.variant_map

    df = pd.read_csv(mpralm_file, sep="\t", header=0, index_col=0)

    variant_dna_counts = mpradata.variant_dna_counts
    variant_rna_counts = mpradata.variant_rna_counts

    columns_ref = []
    columns_alt = []
    for replicate in mpradata.oligo_data.obs_names:
        columns_ref.append("counts_" + replicate + "_REF")
        columns_alt.append("counts_" + replicate + "_ALT")

    dna_counts = np.array(variant_dna_counts[columns_ref].sum()) + np.array(variant_dna_counts[columns_alt].sum())
    rna_counts = np.array(variant_rna_counts[columns_ref].sum()) + np.array(variant_rna_counts[columns_alt].sum())

    for spdi, row in spdi_map.iterrows():
        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = (
                (variant_dna_counts.loc[spdi][columns_ref] / dna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "inputCountAlt"] = (
                (variant_dna_counts.loc[spdi][columns_alt] / dna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "outputCountRef"] = (
                (variant_rna_counts.loc[spdi][columns_ref] / rna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "outputCountAlt"] = (
                (variant_rna_counts.loc[spdi][columns_alt] / rna_counts) * MPRABarcodeData.SCALING
            ).mean()
            df.loc[spdi, "variantPos"] = int(
                mpradata.oligo_data.var["variant_pos"][mpradata.oligos.isin(row["REF"])].values[0][0]
            )
            df.loc[spdi, "strand"] = mpradata.oligo_data.var["strand"][mpradata.oligos.isin(row["REF"])].values[0]

    df["variantPos"] = df["variantPos"].astype(int)
    df["minusLog10PValue"] = -np.log10(df["P.Value"])
    df["minusLog10QValue"] = -np.log10(df["adj.P.Val"])
    df.rename(
        columns={
            "CI.L": "CI_lower_95",
            "CI.R": "CI_upper_95",
            "logFC": "log2FoldChange",
        },
        inplace=True,
    )
    df["postProbEffect"] = df["B"].apply(lambda x: math.exp(x) / (1 + math.exp(x)))
    df["variant_id"] = df.index
    df["refAllele"] = df["variant_id"].apply(lambda x: x.split(":")[2])
    df["altAllele"] = df["variant_id"].apply(lambda x: x.split(":")[3])
    df["start"] = df["variant_id"].apply(lambda x: x.split(":")[1]).astype(int)
    df["end"] = df["start"] + df["refAllele"].apply(lambda x: len(x)).astype(int)

    map = chromosome_map()
    df["chr"] = df["variant_id"].apply(lambda x: map[map["refseq"] == x.split(":")[0]].loc[:, "ucsc"].values[0])

    df.dropna(inplace=True)

    scaler = MinMaxScaler(feature_range=(0, 1000))
    df["score"] = scaler.fit_transform(df[["log2FoldChange"]]).astype(int)

    df = df[
        [
            "chr",
            "start",
            "end",
            "variant_id",
            "score",
            "strand",
            "log2FoldChange",
            "inputCountRef",
            "outputCountRef",
            "inputCountAlt",
            "outputCountAlt",
            "minusLog10PValue",
            "minusLog10QValue",
            "postProbEffect",
            "CI_lower_95",
            "CI_upper_95",
            "variantPos",
            "refAllele",
            "altAllele",
        ]
    ].sort_values(by=["chr", "start", "end"])

    with pysam.BGZFile(output_reporter_genomic_variants_file, "ab") as f:
        f.write(df.to_csv(sep="\t", index=False, header=False).encode())


if __name__ == "__main__":
    cli()
