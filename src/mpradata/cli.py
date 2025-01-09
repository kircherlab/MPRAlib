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
        replicate_data = data[replicate,:]
        replicate_data = replicate_data[:,replicate_data.layers['barcodes'] >= bc_threshold]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var_names.values,
            "dna_counts": replicate_data.layers["dna"][0,:],
            "rna_counts": replicate_data.layers["rna"][0,:],
            "dna_normalized": np.round(replicate_data.layers["dna_normalized"][0,:], 4),
            "rna_normalized": np.round(replicate_data.layers["rna_normalized"][0,:], 4),
            "log2FoldChange": np.round(replicate_data.layers["log2FoldChange"][0,:], 4),
            "n_bc": replicate_data.layers["barcodes"][0,:],
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
    default = 3,
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


    # mpradata.filter_outlier(OutlierFilter.MAD, {})

    # mpradata.filter_outlier(OutlierFilter.RNA_ZSCORE, {"times_zscore": rna_zscore_times})

    
    
    data = mpradata.grouped_data

    print(data.layers['barcodes'].sum())

    output = pd.DataFrame()
    
    for replicate in data.obs['replicate']:
        replicate_data = data[replicate,:]
        replicate_data = replicate_data[:,replicate_data.layers['barcodes'] >= bc_threshold]
        df = {
            "replicate": np.repeat(replicate, replicate_data.var_names.size),
            "oligo_name": replicate_data.var_names.values,
            "dna_counts": replicate_data.layers["dna"][0,:],
            "rna_counts": replicate_data.layers["rna"][0,:],
            "dna_normalized": np.round(replicate_data.layers["dna_normalized"][0,:], 4),
            "rna_normalized": np.round(replicate_data.layers["rna_normalized"][0,:], 4),
            "log2FoldChange": np.round(replicate_data.layers["log2FoldChange"][0,:], 4),
            "n_bc": replicate_data.layers["barcodes"][0,:],
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

    data = mpradata.get_grouped_data()

    spdis = set([item for sublist in data.var['SPDI'].values for item in sublist])

    df = {
        "SPDI": [],
        "ref": [],
        "alt": []
    }
    for spdi in spdis:
        df["SPDI"].append(spdi)
        spdi_data = data[:, data.var['SPDI'].apply(lambda x: spdi in x)]
        spdi_idx = spdi_data.var['SPDI'].apply(lambda x: x.index(spdi))
        refs = []
        alts = []
        for idx, value in spdi_data.var["allele"].items():
            if "ref" == value[spdi_idx[idx]]:
                refs.append(idx)
            else:
                alts.append(idx)
        df["ref"].append(refs)
        df["alt"].append(alts)
    for key in ["ref", "alt"]:
        df[key] = [",".join(i) for i in df[key] ]
    
    output = pd.DataFrame(df)
    output.to_csv(output_file, sep='\t', index=False)
    
        

if __name__ == '__main__':
    cli()