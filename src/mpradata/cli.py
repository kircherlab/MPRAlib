import click
import pandas as pd
import numpy as np
import math
from mpradata import MPRAdata
from mpradata import OutlierFilter


@click.group()
def cli():
    pass


@cli.command()
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
@click.option(
    "--output",
    "output_file",
    required=True,
    type=click.Path(writable=True),
    help="Output file of results.",
)
def activities(input_file, bc_threshold, output_file):
    """Reads a file and generates an MPRAdata object."""
    mpradata = MPRAdata.from_file(input_file)

    data = mpradata.get_grouped_data()

    output = pd.DataFrame()
    
    for replicate in data.obs['replicate']:
        replicate_data = data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers['barcodes'] >= bc_threshold]
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
    
    output.to_csv(output_file, sep='\t', index=False)


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
    mpradata = MPRAdata.from_file(input_file)

    mpradata.barcode_threshold = bc_threshold

    # mpradata.filter_outlier(OutlierFilter.MAD, {})

    mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": rna_zscore_times})
    
    print(mpradata.spearman_correlation)
    print(mpradata.pearson_correlation)

    data = mpradata.grouped_data

    print(data.layers['barcodes'].sum())
    print((data.layers['barcodes'] == 0).sum())

    output = pd.DataFrame()
    
    for replicate in data.obs['replicate']:
        replicate_data = data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers['barcodes'] >= bc_threshold]
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
    
    output.to_csv(output_file, sep='\t', index=False)


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
    mpradata = MPRAdata.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    variant_map = mpradata.variant_map
    for key in ["REF", "ALT"]:
        variant_map[key] = [",".join(i) for i in variant_map[key]]

    variant_map.to_csv(output_file, sep='\t', index=True)


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

    mpradata = MPRAdata.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": 3})

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")

    mpradata.filter = (mpradata.filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1))

    click.echo("After correlate filtering:")
    click.echo(mpradata.pearson_correlation)

    mpradata.element_dna_counts.to_csv(output_dna_file, sep='\t', index=True)
    mpradata.element_rna_counts.to_csv(output_rna_file, sep='\t', index=True)

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
    mpradata = MPRAdata.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    click.echo("Initial Pearson correlation:")
    click.echo(mpradata.pearson_correlation)

    mpradata.barcode_threshold = 10

    click.echo("After BC-threshold filtering:")
    click.echo(mpradata.pearson_correlation)

    mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": 3})

    click.echo("After ZSCORE filtering:")
    click.echo(mpradata.pearson_correlation)

    mpradata.variant_dna_counts.to_csv(output_dna_file, sep='\t', index=True)
    mpradata.variant_rna_counts.to_csv(output_rna_file, sep='\t', index=True)


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

    mpradata = MPRAdata.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": 3})

    mask = mpradata.data.var["allele"].apply(lambda x: "ref" in x).values | (mpradata.data.var["category"] == "element")

    mpradata.filter = (mpradata.filter | ~np.repeat(np.array(mask)[:, np.newaxis], 3, axis=1))

    df = pd.read_csv(mpralm_file, sep='\t', header=0)
    
    indexes_in_order = [mpradata.grouped_data.var["oligo"][
        mpradata.grouped_data.var["oligo"] == ID].index.tolist() for ID in df["ID"]
    ]
    indexes_in_order = [index for sublist in indexes_in_order for index in sublist]
    df.index = indexes_in_order

    df = df.join(mpradata.grouped_data.var["oligo"], how="right")

    mpradata.grouped_data.varm["mpralm_element"] = df

    out_df = mpradata.grouped_data.varm["mpralm"][["oligo", "logFC", "P.Value", "adj.P.Val"]]
    out_df["inputCount"] = mpradata.grouped_data.layers["dna_normalized"].mean(axis=0)
    out_df["outputCount"] = mpradata.grouped_data.layers["rna_normalized"].mean(axis=0)
    out_df.dropna(inplace=True)
    out_df["minusLog10PValue"] = -np.log10(out_df["P.Value"])
    out_df["minusLog10QValue"] = -np.log10(out_df["adj.P.Val"])
    out_df.rename(columns={"oligo": "oligo_name", "logFC": "log2FoldChange"}, inplace=True)
    out_df[['oligo_name', 'log2FoldChange', 'inputCount', 'outputCount', 'minusLog10PValue', 'minusLog10QValue']].to_csv(
        output_reporter_elements_file, sep='\t', index=False
    )


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

    mpradata = MPRAdata.from_file(input_file)

    mpradata.add_metadata_file(metadata_file)

    mpradata.barcode_threshold = 10

    mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": 3})

    spdi_map = mpradata.variant_map

    df = pd.read_csv(mpralm_file, sep='\t', header=0, index_col=0)

    variant_dna_counts = mpradata.variant_dna_counts
    variant_rna_counts = mpradata.variant_rna_counts

    columns_ref = []
    columns_alt = []
    for replicate in mpradata.grouped_data.obs_names:
        columns_ref.append("counts_" + replicate + "_REF")
        columns_alt.append("counts_" + replicate + "_ALT")

    dna_counts = np.array(variant_dna_counts[columns_ref].sum()) + np.array(variant_dna_counts[columns_alt].sum())
    rna_counts = np.array(variant_rna_counts[columns_ref].sum()) + np.array(variant_rna_counts[columns_alt].sum())
    
    for spdi, row in spdi_map.iterrows():
        if spdi in df.index:
            df.loc[spdi, "inputCountRef"] = ((variant_dna_counts.loc[spdi][columns_ref] / dna_counts) * MPRAdata.SCALING).mean()
            df.loc[spdi, "inputCountAlt"] = ((variant_dna_counts.loc[spdi][columns_alt] / dna_counts) * MPRAdata.SCALING).mean()
            df.loc[spdi, "outputCountRef"] = ((variant_rna_counts.loc[spdi][columns_ref] / rna_counts) * MPRAdata.SCALING).mean()
            df.loc[spdi, "outputCountAlt"] = ((variant_rna_counts.loc[spdi][columns_alt] / rna_counts) * MPRAdata.SCALING).mean()
            df.loc[spdi, "variantPos"] = int(mpradata.grouped_data.var["variant_pos"][mpradata.oligos.isin(row["REF"])].values[0][0])
    
    df["minusLog10PValue"] = -np.log10(df["P.Value"])
    df["minusLog10QValue"] = -np.log10(df["adj.P.Val"])
    df.rename(columns={"CI.L": "CI_lower_95", "CI.R": "CI_upper_95", "logFC": "log2FoldChange"}, inplace=True)
    df["postProbEffect"] = df["B"].apply(lambda x: math.exp(x) / (1 + math.exp(x)))
    df["variant_id"] = df.index
    df['refAllele'] = df["variant_id"].apply(lambda x: x.split(":")[2])
    df['altAllele'] = df["variant_id"].apply(lambda x: x.split(":")[3])
    df['variantPos'] = df['variantPos'].astype(int)

    df[['variant_id', 'log2FoldChange', 'inputCountRef', 'outputCountRef', 'inputCountAlt', 'outputCountAlt', 'minusLog10PValue', 'minusLog10QValue', 'postProbEffect', 'CI_lower_95', 'CI_upper_95', 'variantPos', 'refAllele', 'altAllele']].to_csv(
        output_reporter_variants_file, sep='\t', index=False
    )


if __name__ == '__main__':
    cli()
