import click
import pandas as pd
import numpy as np
import math
import pysam
from sklearn.preprocessing import MinMaxScaler
from mpralib.mpradata import MPRABarcodeData, BarcodeFilter
from mpralib.utils import chromosome_map, export_activity_file, export_barcode_file, export_counts_file

pd.options.mode.copy_on_write = True


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

    print("Generating oligo data...")
    mpradata = mpradata.oligo_data

    print("Adding metadata file...")
    mpradata.add_metadata_file(metadata_file)

    print("Generating variant map...")
    variant_map = mpradata.variant_map

    print("Prepair output file...")
    for key in ["REF", "ALT"]:
        # TODO: what happens a variant has multiple alt alleles?
        # Right now joined by comma. But maybe I hould use one row per ref/alt pai?
        variant_map[key] = [",".join(i) for i in variant_map[key]]

    variant_map = variant_map[variant_map.apply(lambda x: all(x != ""), axis=1)]

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
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output.",
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
def get_counts(input_file, metadata_file, bc_threshold, use_oligos, elements_only, output_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

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
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output.",
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
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of all non zero counts, divided by ref and alt.",
)
def get_variant_counts(input_file, metadata_file, bc_threshold, use_oligos, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = bc_threshold

    mask = mpradata.data.var["category"] == "variant"
    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    if use_oligos:
        # TODO adapt to Barcode thresholds
        mpradata = mpradata.oligo_data

        variant_map = mpradata.variant_map

        df = {"variant_id": []}
        for replicate in mpradata.obs_names:
            for rna_or_dna in ["dna", "rna"]:
                df[f"{rna_or_dna}_count_{replicate}_REF"] = []
                df[f"{rna_or_dna}_count_{replicate}_ALT"] = []

        dna_counts = mpradata.dna_counts.copy()
        rna_counts = mpradata.rna_counts.copy()

        for spdi, row in variant_map.iterrows():

            mask_ref = mpradata.oligos.isin(row["REF"])
            mask_alt = mpradata.oligos.isin(row["ALT"])

            df["variant_id"].append(spdi)

            dna_counts_ref = dna_counts[:, mask_ref].sum(axis=1)
            rna_counts_ref = rna_counts[:, mask_ref].sum(axis=1)

            dna_counts_alt = dna_counts[:, mask_alt].sum(axis=1)
            rna_counts_alt = rna_counts[:, mask_alt].sum(axis=1)

            for idx, replicate in enumerate(mpradata.obs_names):
                df[f"dna_count_{replicate}_REF"].append(dna_counts_ref[idx])
                df[f"rna_count_{replicate}_REF"].append(rna_counts_ref[idx])
                df[f"dna_count_{replicate}_ALT"].append(dna_counts_alt[idx])
                df[f"rna_count_{replicate}_ALT"].append(rna_counts_alt[idx])

        df = pd.DataFrame(df).set_index("variant_id")
        # remove IDs which are all zero
        df = df[(df.T != 0).all()].dropna()
        df.to_csv(output_file, sep="\t", index=True)
    else:
        export_counts_file(mpradata, output_file)


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
    "--elements-only/--all-oligos",
    "elements_only",
    required=False,
    type=click.BOOL,
    default=False,
    help="Only return count data for elements and ref sequence of variants.",
)
@click.option(
    "--statistics",
    "statistics_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file of mpralm or BCalm results.",
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
    "--output-reporter-elements",
    "output_reporter_elements_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_elements(
    input_file, metadata_file, statistics_file, elements_only, bc_threshold, output_reporter_elements_file
):

    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = bc_threshold

    if elements_only:

        mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")
        mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    df = pd.read_csv(statistics_file, sep="\t", header=0)

    indexes_in_order = [mpradata.oligos[mpradata.oligos == ID].index.tolist() for ID in df["ID"]]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    df.index = indexes_in_order

    df = df.join(mpradata.oligos, how="right")

    mpradata.data.varm["mpralm_element"] = df

    out_df = mpradata.data.varm["mpralm_element"][["oligo", "logFC", "P.Value", "adj.P.Val"]]
    out_df.loc[:, "inputCount"] = mpradata.normalized_dna_counts.mean(axis=0)
    out_df.loc[:, "outputCount"] = mpradata.normalized_rna_counts.mean(axis=0)
    out_df.dropna(inplace=True)
    out_df.loc[out_df["P.Value"] == 0, "P.Value"] = np.finfo(float).eps
    out_df.loc[out_df["adj.P.Val"] == 0, "adj.P.Val"] = np.finfo(float).eps
    out_df.loc[:, "minusLog10PValue"] = np.abs(np.log10(out_df.loc[:, "P.Value"]))
    out_df.loc[:, "minusLog10QValue"] = np.abs(np.log10(out_df.loc[:, "adj.P.Val"]))
    out_df = out_df.rename(columns={"oligo": "oligo_name", "logFC": "log2FoldChange"})
    out_df[
        [
            "oligo_name",
            "log2FoldChange",
            "inputCount",
            "outputCount",
            "minusLog10PValue",
            "minusLog10QValue",
        ]
    ].to_csv(output_reporter_elements_file, sep="\t", index=False, float_format="%.4f")


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
    "--statistics",
    "statistics_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file of the result from mpralm or BCalm.",
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
    "--output-reporter-variants",
    "output_reporter_variants_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_variants(input_file, metadata_file, statistics_file, bc_threshold, output_reporter_variants_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = bc_threshold

    mask = mpradata.data.var["category"] == "variant"
    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    mpradata = mpradata.oligo_data

    variant_map = mpradata.variant_map

    df = pd.read_csv(statistics_file, sep="\t", header=0, index_col=0)

    dna_counts = mpradata.dna_counts.copy()
    rna_counts = mpradata.rna_counts.copy()

    rna_counts_sum = rna_counts.sum(axis=1)
    dna_counts_sum = dna_counts.sum(axis=1)

    for spdi, row in variant_map.iterrows():

        mask_ref = mpradata.oligos.isin(row["REF"])
        mask_alt = mpradata.oligos.isin(row["ALT"])

        dna_counts_ref = dna_counts[:, mask_ref].sum(axis=1)
        rna_counts_ref = rna_counts[:, mask_ref].sum(axis=1)

        dna_counts_alt = dna_counts[:, mask_alt].sum(axis=1)
        rna_counts_alt = rna_counts[:, mask_alt].sum(axis=1)

        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = ((dna_counts_ref / dna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "inputCountAlt"] = ((dna_counts_alt / dna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "outputCountRef"] = ((rna_counts_ref / rna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "outputCountAlt"] = ((rna_counts_alt / rna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "variantPos"] = int(mpradata.data.var["variant_pos"][mpradata.oligos.isin(row["REF"])].values[0][0])

    df.loc[df["P.Value"] == 0, "P.Value"] = np.finfo(float).eps
    df.loc[df["adj.P.Val"] == 0, "adj.P.Val"] = np.finfo(float).eps
    df["minusLog10PValue"] = np.abs(np.log10(df["P.Value"]))
    df["minusLog10QValue"] = np.abs(np.log10(df["adj.P.Val"]))
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
    ].to_csv(output_reporter_variants_file, sep="\t", index=False, float_format="%.4f")


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
    "--statistics",
    "statistics_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA lm file.",
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
    "--output-reporter-genomic-elements",
    "output_reporter_genomic_elements_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_genomic_elements(
    input_file, metadata_file, statistics_file, bc_threshold, output_reporter_genomic_elements_file
):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = bc_threshold

    mpradata = mpradata.oligo_data

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")
    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    df = pd.read_csv(statistics_file, sep="\t", header=0)

    indexes_in_order = [mpradata.oligos[mpradata.oligos == ID].index.tolist() for ID in df["ID"]]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    df.index = indexes_in_order

    df = df.join(mpradata.oligos, how="right")

    mpradata.data.varm["mpralm_element"] = df

    out_df = mpradata.data.varm["mpralm_element"][["oligo", "logFC", "P.Value", "adj.P.Val"]]

    out_df.loc[:, ["inputCount"]] = mpradata.normalized_dna_counts.mean(axis=0)
    out_df.loc[:, ["outputCount"]] = mpradata.normalized_rna_counts.mean(axis=0)

    out_df["chr"] = mpradata.data.var["chr"]
    out_df["start"] = mpradata.data.var["start"]
    out_df["end"] = mpradata.data.var["end"]
    out_df["strand"] = mpradata.data.var["strand"]

    # apply ref and element mask
    out_df = out_df[mask]

    out_df.dropna(inplace=True)
    scaler = MinMaxScaler(feature_range=(0, 1000))
    out_df["score"] = scaler.fit_transform(np.abs(out_df[["logFC"]])).astype(int)
    out_df.loc[out_df["P.Value"] == 0, "P.Value"] = np.finfo(float).eps
    out_df.loc[out_df["adj.P.Val"] == 0, "adj.P.Val"] = np.finfo(float).eps
    out_df["minusLog10PValue"] = np.abs(np.log10(out_df["P.Value"]))
    out_df["minusLog10QValue"] = np.abs(np.log10(out_df["adj.P.Val"]))
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
        f.write(out_df.to_csv(sep="\t", index=False, header=False, float_format="%.4f").encode())


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
    "--statistics",
    "statistics_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file of the result from mpralm or BCalm.",
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
    "--output-reporter-genomic-variants",
    "output_reporter_genomic_variants_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_genomic_variants(
    input_file, metadata_file, statistics_file, bc_threshold, output_reporter_genomic_variants_file
):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = bc_threshold

    mask = mpradata.data.var["category"] == "variant"
    mpradata.var_filter = mpradata.var_filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1)

    mpradata = mpradata.oligo_data

    variant_map = mpradata.variant_map

    df = pd.read_csv(statistics_file, sep="\t", header=0, index_col=0)

    dna_counts = mpradata.dna_counts.copy()
    rna_counts = mpradata.rna_counts.copy()

    rna_counts_sum = rna_counts.sum(axis=1)
    dna_counts_sum = dna_counts.sum(axis=1)

    for spdi, row in variant_map.iterrows():

        mask_ref = mpradata.oligos.isin(row["REF"])
        mask_alt = mpradata.oligos.isin(row["ALT"])

        dna_counts_ref = dna_counts[:, mask_ref].sum(axis=1)
        rna_counts_ref = rna_counts[:, mask_ref].sum(axis=1)

        dna_counts_alt = dna_counts[:, mask_alt].sum(axis=1)
        rna_counts_alt = rna_counts[:, mask_alt].sum(axis=1)

        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = ((dna_counts_ref / dna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "inputCountAlt"] = ((dna_counts_alt / dna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "outputCountRef"] = ((rna_counts_ref / rna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "outputCountAlt"] = ((rna_counts_alt / rna_counts_sum) * MPRABarcodeData.SCALING).mean()
            df.loc[spdi, "variantPos"] = int(mpradata.data.var["variant_pos"][mpradata.oligos.isin(row["REF"])].values[0][0])
            df.loc[spdi, "strand"] = mpradata.data.var["strand"][mpradata.oligos.isin(row["REF"])].values[0]

    df["variantPos"] = df["variantPos"].astype(int)
    df.loc[df["P.Value"] == 0, "P.Value"] = np.finfo(float).eps
    df.loc[df["adj.P.Val"] == 0, "adj.P.Val"] = np.finfo(float).eps
    df["minusLog10PValue"] = np.abs(np.log10(df["P.Value"]))
    df["minusLog10QValue"] = np.abs(np.log10(df["adj.P.Val"]))
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
    df["score"] = scaler.fit_transform(np.abs(df[["log2FoldChange"]])).astype(int)

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
        f.write(df.to_csv(sep="\t", index=False, header=False, float_format="%.4f").encode())


if __name__ == "__main__":
    cli()
