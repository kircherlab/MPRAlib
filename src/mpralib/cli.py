import click
import pandas as pd
import numpy as np
import math
import pysam
from sklearn.preprocessing import MinMaxScaler
from mpralib.mpradata import MPRABarcodeData, BarcodeFilter
from mpralib.utils import chromosome_map, export_activity_file, export_barcode_file, export_counts_file


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


@cli.command(help="Compute pairwise correlations under a certain barcode threshold.")
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
    "--correlation-method",
    "correlation_method",
    required=False,
    default="all",
    type=click.Choice(["pearson", "spearman", "all"]),
    help="Computing pearson, spearman or both.",
)
@click.option(
    "--correlation-on",
    "correlation_on",
    required=False,
    default="all",
    type=click.Choice(["dna", "rna", "activity", "all"]),
    help="Using a barcode threshold for output (element level only).",
)
def correlation(input_file, bc_threshold, correlation_on, correlation_method):
    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.barcode_threshold = bc_threshold

    if correlation_on == "all":
        correlation_on = ["dna", "rna", "activity"]
    else:
        correlation_on = [correlation_on]

    if correlation_method == "all":
        correlation_method = ["pearson", "spearman"]
    else:
        correlation_method = [correlation_method]

    for method in correlation_method:
        for on in correlation_on:
            click.echo(f"{method} correlation on {on}: {mpradata.correlation(method, on).flatten()[[1, 2, 5]]}")


@cli.command(help="Filter out outliers based on RNA z-score and copute correlations before and afterwards.")
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
    required=False,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def filter_outliers(input_file, rna_zscore_times, bc_threshold, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    oligo_data = mpradata.oligo_data

    click.echo(
        f"Pearson correlation log2FoldChange BEFORE outlier removal: {
            oligo_data.correlation("pearson", "activity").flatten()[[1, 2, 5]]
            }"
    )

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": rna_zscore_times})

    oligo_data = mpradata.oligo_data
    click.echo(
        f"Pearson correlation log2FoldChange AFTER outlier removal: {
            oligo_data.correlation("pearson", "activity").flatten()[[1, 2, 5]]
            }"
    )
    if output_file:
        export_activity_file(oligo_data, output_file)


@cli.command(help="Using a metadata file to generate a oligo to variant mapping.")
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

    mpradata = mpradata.oligo_data

    mpradata.add_metadata_file(metadata_file)

    variant_map = mpradata.variant_map
    for key in ["REF", "ALT"]:
        variant_map[key] = [",".join(i) for i in variant_map[key]]

    variant_map.to_csv(output_file, sep="\t", index=True)


@cli.command(help="Return DNA and RNA counts for all oligos or only thosed tagges as elements/ref snvs.")
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
    "--oligos/--barcodes",
    "use_oligos",
    required=False,
    type=click.BOOL,
    default=True,
    help="Return counst per oligo or per barcode.",
)
@click.option(
    "--elements-only/--all-oligos",
    "elements_only",
    required=False,
    type=click.BOOL,
    default=False,
    help="Only return count data for elements and ref sequence of variants.",
)
@click.option(
    "--output",
    "output_file",
    required=False,
    type=click.Path(writable=True),
    help="Output file of all non zero counts.",
)
def get_counts(input_file, metadata_file, use_oligos, elements_only, output_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    if use_oligos:
        mpradata = mpradata.oligo_data

    mpradata.add_metadata_file(metadata_file)

    click.echo(
        f"Pearson correlation log2FoldChange, all oligos: {
            mpradata.correlation("pearson", "activity").flatten()[[1, 2, 5]]
            }"
    )

    if elements_only:

        mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")
        mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    click.echo(
        f"Pearson correlation log2FoldChange, element oligos: {
            mpradata.correlation("pearson", "activity").flatten()[[1, 2, 5]]
            }"
    )

    # TODO: Implement a barcode_threhold_filter silimar to var_filter marks barcodes/oligos with lessthan a certain number of counts
    if output_file:
        export_counts_file(mpradata, output_file)


@cli.command(help="Write out DNA and RNA counts for REF and ALT oligos.")
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
    help="Output file of all non zero counts, divided by ref and alt.",
)
def get_variant_counts(input_file, metadata_file, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.add_metadata_file(metadata_file)

    mask = mpradata.data.var["category"] == "variant"
    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    df = {"ID": []}
    for replicate in mpradata.obs_names:
        df["dna_counts_" + replicate + "_REF"] = []
        df["dna_counts_" + replicate + "_ALT"] = []
        df["rna_counts_" + replicate + "_REF"] = []
        df["rna_counts_" + replicate + "_ALT"] = []

    for spdi, row in mpradata.variant_map.iterrows():
        df["ID"].append(spdi)
        mask_ref = mpradata.oligos.isin(row["REF"])
        mask_alt = mpradata.oligos.isin(row["ALT"])
        dna_counts_ref = mpradata.dna_counts[:, mask_ref].sum(axis=1)
        dna_counts_alt = mpradata.dna_counts[:, mask_alt].sum(axis=1)
        rna_counts_ref = mpradata.rna_counts[:, mask_ref].sum(axis=1)
        rna_counts_alt = mpradata.rna_counts[:, mask_alt].sum(axis=1)
        idx = 0
        for idx, replicate in enumerate(mpradata.obs_names):
            df["dna_counts_" + replicate + "_REF"].append(dna_counts_ref[idx])
            df["dna_counts_" + replicate + "_ALT"].append(dna_counts_alt[idx])
            df["rna_counts_" + replicate + "_REF"].append(rna_counts_ref[idx])
            df["rna_counts_" + replicate + "_ALT"].append(rna_counts_alt[idx])
    df = pd.DataFrame(df).set_index("ID")

    # remove IDs which are all zero
    df = df[(df.T != 0).all()]

    df.to_csv(output_file, sep="\t", index=True)


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
