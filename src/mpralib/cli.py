import click
import pandas as pd
import numpy as np
import math
import pysam
from sklearn.preprocessing import MinMaxScaler
from mpralib.mpradata import MPRABarcodeData, BarcodeFilter, Modality
from mpralib.utils.io import (
    chromosome_map,
    export_activity_file,
    export_barcode_file,
    export_counts_file,
    read_sequence_design_file,
)
import mpralib.utils.plot as plt
from mpralib.utils.file_validation import validate_tsv_with_schema, ValidationSchema

pd.options.mode.copy_on_write = True


@click.group(help="Command line interface of MPRAlib, a library for MPRA data analysis.")
def cli():
    pass


@cli.group(help="Validate standardized MPRA reporter formats.")
def validate_file():
    pass


@validate_file.command(help="Validate MPRA reporter sequence design file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Sequence Design file to validate.",
)
def reporter_sequence_design(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_SEQUENCE_DESIGN):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate MPRA Reporter Barcode to Element Mapping file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Barcode to Element Mapping file to validate.",
)
def reporter_barcode_to_element_mapping(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_BARCODE_TO_ELEMENT_MAPPING):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Experiment Barcode file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Experiment Barcode file to validate.",
)
def reporter_experiment_barcode(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_EXPERIMENT_BARCODE):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Experiment file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Experiment file to validate.",
)
def reporter_experiment(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_EXPERIMENT):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Element file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Element file to validate.",
)
def reporter_element(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_ELEMENT):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Variant file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Variant file to validate.",
)
def reporter_variant(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_VARIANT):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Genomic Element file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Genomic Element file to validate.",
)
def reporter_genomic_element(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_GENOMIC_ELEMENT):
        raise click.ClickException("Validation failed. Please check the input file.")


@validate_file.command(help="Validate Reporter Genomic Variant file.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="MPRA Reporter Genomic Variant file to validate.",
)
def reporter_genomic_variant(input_file):
    if not validate_tsv_with_schema(input_file, ValidationSchema.REPORTER_GENOMIC_VARIANT):
        raise click.ClickException("Validation failed. Please check the input file.")


@cli.group(help="General functionality.")
def functional():
    pass


@functional.command(help="Generating element activity or barcode count files.")
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


@functional.command(help="Compute pairwise correlations under a certain barcode threshold.")
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
    type=click.Choice(["dna_normalized", "rna_normalized", "activity", "all"]),
    help="Using a barcode threshold for output (element level only).",
)
def compute_correlation(input_file, bc_threshold, correlation_on, correlation_method):
    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.barcode_threshold = bc_threshold

    if correlation_on == "all":
        correlation_on = [Modality.DNA_NORMALIZED, Modality.RNA_NORMALIZED, Modality.ACTIVITY]
    else:
        correlation_on = [Modality.from_string(correlation_on)]

    if correlation_method == "all":
        correlation_method = ["pearson", "spearman"]
    else:
        correlation_method = [correlation_method]

    for method in correlation_method:
        for on in correlation_on:
            click.echo(f"{method} correlation on {on}: {mpradata.correlation(method, on).flatten()[[1, 2, 5]]}")


functional.command(help="Filter out outliers based on RNA z-score and copute correlations before and afterwards.")


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
        f"""Pearson correlation log2FoldChange BEFORE outlier removal: {
            oligo_data.correlation('pearson', Modality.ACTIVITY).flatten()[[1, 2, 5]]
        }"""
    )

    mpradata.apply_barcode_filter(BarcodeFilter.RNA_ZSCORE, {"times_zscore": rna_zscore_times})

    oligo_data = mpradata.oligo_data
    click.echo(
        f"""Pearson correlation log2FoldChange AFTER outlier removal: {
            oligo_data.correlation('pearson', Modality.ACTIVITY).flatten()[[1, 2, 5]]
        }"""
    )
    if output_file:
        export_activity_file(oligo_data, output_file)


@cli.group(help="MPRA sequence design file functionality.")
def sequence_design():
    pass


@sequence_design.command(help="Using a metadata file to generate a oligo to variant mapping.")
@click.option(
    "--input",
    "input_file",
    required=False,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results (barcode output file). "
    "If set only map of present oligos in the barcode count file will be generated.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def get_variant_map(input_file, sequence_design_file, output_file):

    df_sequence_design = read_sequence_design_file(sequence_design_file)

    if input_file:
        print("Read barcode count file...")
        mpradata = MPRABarcodeData.from_file(input_file)

        print("Generating oligo data...")
        mpradata = mpradata.oligo_data

        print("Adding metadata file...")
        mpradata.add_sequence_design(df_sequence_design, sequence_design_file)

        print(mpradata.data.var["SPDI"])
        print(type(mpradata.data.var["SPDI"]))
        print(mpradata.data.var["SPDI"].values)

        print("Generating variant map...")
        variant_map = mpradata.variant_map
    else:
        variant_map = pd.DataFrame(
            {
                "ID": np.concatenate(df_sequence_design["SPDI"].values),
                "allele": np.concatenate(df_sequence_design["allele"].values),
                "oligo": df_sequence_design.index.values.repeat(df_sequence_design["SPDI"].apply(lambda x: len(x))),
            }
        )
        variant_map["REF"] = variant_map.apply(lambda row: row["oligo"] if row["allele"] == "ref" else None, axis=1)
        variant_map["ALT"] = variant_map.apply(lambda row: row["oligo"] if row["allele"] == "alt" else None, axis=1)
        variant_map = variant_map.groupby("ID").agg(
            {"REF": lambda x: list(filter(None, x)), "ALT": lambda x: list(filter(None, x))}
        )

    print("Prepair output file...")
    for key in ["REF", "ALT"]:
        # TODO: what happens a variant has multiple alt alleles?
        # Right now joined by comma. But maybe I hould use one row per ref/alt pai?
        variant_map.loc[:, key] = variant_map[key].apply(lambda x: ",".join(x))

    variant_map = variant_map[variant_map.apply(lambda x: all(x != ""), axis=1)]

    variant_map.to_csv(output_file, sep="\t", index=True)


@sequence_design.command(help="Return DNA and RNA counts for all oligos or only thosed tagges as elements/ref snvs.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
    "--normalized-counts/--counts",
    "normalized_counts",
    required=False,
    type=click.BOOL,
    default=False,
    help="Getting counts or normalized counts.",
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
    required=True,
    type=click.Path(writable=True),
    help="Output file of all non zero counts.",
)
def get_counts(input_file, sequence_design_file, bc_threshold, normalized_counts, use_oligos, elements_only, output_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    if use_oligos:
        mpradata = mpradata.oligo_data

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    element_mask = None
    if elements_only:

        element_mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (
            mpradata.data.var["category"] == "element"
        )
        element_mask = ~np.repeat(np.array(element_mask)[np.newaxis, :], mpradata.n_obs, axis=0)

    export_counts_file(mpradata, output_file, normalized=normalized_counts, filter=element_mask)


@sequence_design.command(help="Write out DNA and RNA counts for REF and ALT oligos.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
    "--normalized-counts/--counts",
    "normalized_counts",
    required=False,
    type=click.BOOL,
    default=False,
    help="Getting counts or normalized counts.",
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
def get_variant_counts(input_file, sequence_design_file, bc_threshold, normalized_counts, use_oligos, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    mpradata.barcode_threshold = bc_threshold

    if use_oligos:
        mpradata = mpradata.oligo_data

        variant_map = mpradata.variant_map

        df = {"variant_id": []}
        for replicate in mpradata.obs_names:
            for rna_or_dna in ["dna", "rna"]:
                df[f"{rna_or_dna}_count_{replicate}_REF"] = []
                df[f"{rna_or_dna}_count_{replicate}_ALT"] = []

        if normalized_counts:
            dna_counts = mpradata.normalized_dna_counts.copy()
            rna_counts = mpradata.normalized_rna_counts.copy()
        else:
            dna_counts = mpradata.dna_counts.copy()
            rna_counts = mpradata.rna_counts.copy()

        not_observed_mask = ~mpradata.observed
        lower_bc_mask = mpradata.barcode_counts < mpradata.barcode_threshold
        not_variant_mask = ~np.repeat(
            np.array(mpradata.data.var["category"] == "variant")[np.newaxis, :], mpradata.n_obs, axis=0
        )

        dna_counts = np.ma.masked_array(dna_counts, mask=np.any([not_observed_mask, lower_bc_mask, not_variant_mask], axis=0))
        rna_counts = np.ma.masked_array(rna_counts, mask=np.any([not_observed_mask, lower_bc_mask, not_variant_mask], axis=0))

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
        # remove IDs which are all zero or just on ref or alt
        df = df.map(lambda x: np.nan if isinstance(x, np.ma.core.MaskedConstant) else x)
        ref_columns = [f"dna_count_{replicate}_REF" for replicate in mpradata.obs_names] + [
            f"rna_count_{replicate}_REF" for replicate in mpradata.obs_names
        ]
        alt_columns = [f"dna_count_{replicate}_ALT" for replicate in mpradata.obs_names] + [
            f"rna_count_{replicate}_ALT" for replicate in mpradata.obs_names
        ]
        nan_columns = df.loc[:, ref_columns].isna().all(axis=1) | df.loc[:, alt_columns].isna().all(axis=1)
        df = df[~nan_columns]

        if normalized_counts:
            df.to_csv(output_file, sep="\t", index=True, na_rep="", float_format="%.6f")
        else:
            df.to_csv(output_file, sep="\t", index=True, na_rep="", float_format="%.0f")
    else:
        export_counts_file(
            mpradata,
            output_file,
            normalized=normalized_counts,
            filter=~np.repeat(np.array(mpradata.data.var["category"] == "variant")[np.newaxis, :], mpradata.n_obs, axis=0),
        )


@sequence_design.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
def get_reporter_elements(input_file, sequence_design_file, statistics_file, bc_threshold, output_reporter_elements_file):

    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    mpradata.barcode_threshold = bc_threshold

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


@sequence_design.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
def get_reporter_variants(input_file, sequence_design_file, statistics_file, bc_threshold, output_reporter_variants_file):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    mpradata.barcode_threshold = bc_threshold

    mpradata = mpradata.oligo_data

    variant_map = mpradata.variant_map

    df = pd.read_csv(statistics_file, sep="\t", header=0, index_col=0, na_values="NA").dropna()

    dna_counts = mpradata.normalized_dna_counts.copy()
    rna_counts = mpradata.normalized_rna_counts.copy()

    for spdi, row in variant_map.iterrows():

        mask_ref = mpradata.oligos.isin(row["REF"])
        mask_alt = mpradata.oligos.isin(row["ALT"])

        dna_counts_ref = dna_counts[:, mask_ref].sum(axis=1)
        rna_counts_ref = rna_counts[:, mask_ref].sum(axis=1)

        dna_counts_alt = dna_counts[:, mask_alt].sum(axis=1)
        rna_counts_alt = rna_counts[:, mask_alt].sum(axis=1)

        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = dna_counts_ref.mean()
            df.loc[spdi, "inputCountAlt"] = dna_counts_alt.mean()
            df.loc[spdi, "outputCountRef"] = rna_counts_ref.mean()
            df.loc[spdi, "outputCountAlt"] = rna_counts_alt.mean()
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


@sequence_design.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
    "--reference",
    "reference",
    required=True,
    type=str,
    help="Using only this reference as denoted inas ref in the sequence design file.",
)
@click.option(
    "--output-reporter-genomic-elements",
    "output_reporter_genomic_elements_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of MPRA data object.",
)
def get_reporter_genomic_elements(
    input_file, sequence_design_file, statistics_file, bc_threshold, reference, output_reporter_genomic_elements_file
):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    mpradata.barcode_threshold = bc_threshold

    mpradata = mpradata.oligo_data

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")
    mask = mask & (mpradata.data.var["ref"] == reference)

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


@sequence_design.command()
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--sequence-design",
    "sequence_design_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Sequence design file.",
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
    input_file, sequence_design_file, statistics_file, bc_threshold, output_reporter_genomic_variants_file
):

    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.add_sequence_design(read_sequence_design_file(sequence_design_file), sequence_design_file)

    mpradata.barcode_threshold = bc_threshold

    mpradata = mpradata.oligo_data

    variant_map = mpradata.variant_map

    df = pd.read_csv(statistics_file, sep="\t", header=0, index_col=0, na_values="NA").dropna()

    dna_counts = mpradata.normalized_dna_counts.copy()
    rna_counts = mpradata.normalized_rna_counts.copy()

    for spdi, row in variant_map.iterrows():

        mask_ref = mpradata.oligos.isin(row["REF"])
        mask_alt = mpradata.oligos.isin(row["ALT"])

        dna_counts_ref = dna_counts[:, mask_ref].sum(axis=1)
        rna_counts_ref = rna_counts[:, mask_ref].sum(axis=1)

        dna_counts_alt = dna_counts[:, mask_alt].sum(axis=1)
        rna_counts_alt = rna_counts[:, mask_alt].sum(axis=1)

        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = dna_counts_ref.mean()
            df.loc[spdi, "inputCountAlt"] = dna_counts_alt.mean()
            df.loc[spdi, "outputCountRef"] = rna_counts_ref.mean()
            df.loc[spdi, "outputCountAlt"] = rna_counts_alt.mean()
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


@cli.group(help="Plotting functions.")
def plot():
    pass


@plot.command()
@click.option(
    "--input",
    "input_file",
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
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output.",
)
@click.option(
    "--modality",
    "modality",
    required=False,
    default="activity",
    type=click.Choice(["dna_normalized", "rna_normalized", "activity"]),
    help="What modality should be plotted.",
)
@click.option(
    "--replicate",
    "replicates",
    required=False,
    multiple=True,
    type=str,
    help="Copare only these two replicates. Ottheriwse all.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output plot file.",
)
def correlation(input_file, use_oligos, bc_threshold, modality, replicates, output_file):
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    if use_oligos:
        mpradata = mpradata.oligo_data

    modality = Modality.from_string(modality)

    if replicates:
        fig = plt.correlation(mpradata, modality, replicates)
    else:
        fig = plt.correlation(mpradata, modality)

    fig.savefig(output_file)


@plot.command(help="Plotting the DNA vs RNA counts (log10, median on  multiple replicates).")
@click.option(
    "--input",
    "input_file",
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
    "--bc-threshold",
    "bc_threshold",
    required=False,
    default=1,
    type=int,
    help="Using a barcode threshold for output.",
)
@click.option(
    "--replicate",
    "replicates",
    required=False,
    multiple=True,
    type=str,
    help="Copare only these two replicates. Ottheriwse all.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output plot file.",
)
def dna_vs_rna(input_file, use_oligos, bc_threshold, replicates, output_file):
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    if use_oligos:
        mpradata = mpradata.oligo_data

    if replicates:
        fig = plt.dna_vs_rna(mpradata, replicates)
    else:
        fig = plt.dna_vs_rna(mpradata)

    fig.savefig(output_file)


@plot.command(help="Histogramm of barcodes per oligo.")
@click.option(
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, readable=True),
    help="Input file path of MPRA results.",
)
@click.option(
    "--replicate",
    "replicates",
    required=False,
    multiple=True,
    type=str,
    help="Show only these replicates. Otheriwse all.",
)
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output plot file.",
)
def barcodes_per_oligo(input_file, replicates, output_file):
    mpradata = MPRABarcodeData.from_file(input_file).oligo_data

    if replicates:
        fig = plt.barcodes_per_oligo(mpradata, replicates)
    else:
        fig = plt.barcodes_per_oligo(mpradata)

    fig.savefig(output_file)


@plot.command(help="Outlier plot.")
@click.option(
    "--input",
    "input_file",
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
def outlier(input_file, bc_threshold):
    mpradata = MPRABarcodeData.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    plt.barcodes_outlier(mpradata)

    # fig.savefig(output_file)


if __name__ == "__main__":
    cli()
