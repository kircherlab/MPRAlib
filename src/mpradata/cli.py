import click
import pandas as pd
import numpy as np
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


if __name__ == '__main__':
    cli()
