import click
import pandas as pd
import numpy as np
from mpralib.mpradata import MPRABarcodeData
from Bio import SeqIO


@click.command()
@click.option(
    "--assignment", "assignment_file", type=click.Path(exists=True), required=True, help="Path to the assignment file."
)
@click.option("--barcode", "barcode_file", type=click.Path(exists=True), required=True, help="Path to the barcode count file.")
@click.option("--design", "design_file", type=click.Path(exists=True), required=True, help="Path to the design file.")
@click.option("--output", "output_file", type=click.Path(), help="Path to the output file for the metrics.")
def cli(assignment_file, barcode_file, design_file, output_file):

    click.echo("===============Assignment===============")
    # median_assigned_barocdes
    assignment_df = pd.read_csv(assignment_file, sep="\t", header=None)

    assignment_grouped = assignment_df.groupby(1).size()
    median_assigned_barcodes = int(assignment_grouped.median())

    metrics_data = {
        "Workflow": ["Assignment"],
        "Metric": ["median_assigned_barocdes"],
        "Value": [median_assigned_barcodes],
        "Description": ["Median number of barcodes assigned to each oligo."],
    }

    median_assigned_barcodes_df = pd.DataFrame(metrics_data)

    # fraction_assigned_oligos
    assigned_oligos = len(assignment_grouped)
    total_oligos = count_fasta_entries(design_file)

    fraction_assigned_oligos = assigned_oligos / total_oligos
    fraction_assigned_oligos_df = pd.DataFrame(
        {
            "Workflow": ["Assignment"],
            "Metric": ["fraction_assigned_oligos"],
            "Value": [round(fraction_assigned_oligos, 4)],
            "Description": ["Fraction of oligos that have at least one barcode assigned."],
        }
    )

    click.echo("===============Counts===============")

    mpra_barcode_data = MPRABarcodeData.from_file(barcode_file)
    mpra_oligo_data = mpra_barcode_data.oligo_data

    # pearson_correlation
    pearson_correlation_df = pearson_correlation(mpra_oligo_data)

    # median_barcodes_passing_filtering
    median_barcodes_passing_filtering_df = median_barcodes_passing_filtering(mpra_oligo_data)

    # median_rna_read_count
    median_rna_read_count_df = median_rna_read_count(mpra_oligo_data)

    # pct_oligos_passing
    pct_oligos_passing_filtering_df = pct_oligos_passing(mpra_oligo_data, assigned_oligos)

    assignment_metrics_df = pd.concat(
        [
            median_assigned_barcodes_df,
            fraction_assigned_oligos_df,
            pearson_correlation_df,
            median_barcodes_passing_filtering_df,
            median_rna_read_count_df,
            pct_oligos_passing_filtering_df,
        ],
        ignore_index=True,
    )
    print(assignment_metrics_df)

    if output_file:
        assignment_metrics_df.to_csv(output_file, sep='\t', index=False)


def pct_oligos_passing(mpra_oligo_data, assigned_oligos):
    mpra_oligo_data.barcode_threshold = 10
    n_oligos_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_oligos_replicate += [len(replicate_data.var["oligo"])]

    pct_oligos_passing_filtering_data = {
        "Workflow": ["Experiment", "Experiment"],
        "Metric": [
            "pct_oligos_passing",
            "pct_oligos_passing",
        ],
        "Value": [
            round(np.median(n_oligos_replicate) / assigned_oligos * 100, 2),
            round(np.min(n_oligos_replicate) / assigned_oligos * 100, 2),
        ],
        "Description": [
            f"Median of replicates: % of oligos that have at least {mpra_oligo_data.barcode_threshold} barcodes wih counts.",
            f"Min of replicates: % of oligos that have at least {mpra_oligo_data.barcode_threshold} barcodes wih counts.",
        ],
    }

    return pd.DataFrame(pct_oligos_passing_filtering_data)


def pearson_correlation(mpra_oligo_data):
    mpra_oligo_data.barcode_threshold = 10

    correlation_data = {
        "Workflow": ["Experiment", "Experiment"],
        "Metric": ["pearson_correlation", "pearson_correlation"],
        "Value": [
            np.median(mpra_oligo_data.correlation().flatten()[[1, 2, 5]]).round(3),
            mpra_oligo_data.correlation().flatten()[[1, 2, 5]].min().round(3),
        ],
        "Description": [
            f"Median or replicates: Pearson correlation log2FoldChange (BC threshold >= {mpra_oligo_data.barcode_threshold})",
            f"Min or replicates: Pearson correlation log2FoldChange (BC threshold >= {mpra_oligo_data.barcode_threshold})",
        ],
    }

    return pd.DataFrame(correlation_data)


def count_fasta_entries(fasta_file):
    count = 0
    for record in SeqIO.parse(fasta_file, "fasta"):
        count += 1
    return count


def median_barcodes_passing_filtering(mpra_oligo_data):

    mpra_oligo_data.barcode_threshold = 1
    n_barcodes_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_barcodes_replicate += [np.median(replicate_data.layers["barcode_counts"])]

    all = mpra_oligo_data.barcode_counts[mpra_oligo_data.barcode_counts >= mpra_oligo_data.barcode_threshold].flatten()

    metrics_data = {
        "Workflow": ["Experiment", "Experiment", "Experiment"],
        "Metric": [
            "median_barcodes_passing_filtering",
            "median_barcodes_passing_filtering",
            "median_barcodes_passing_filtering",
        ],
        "Value": [int(np.median(n_barcodes_replicate)), int(np.min(n_barcodes_replicate)), int(np.median(all[all != 0]))],
        "Description": [
            f"Median or replicates: Median barcodes per oligo (BC threshold >= {mpra_oligo_data.barcode_threshold})",
            f"Min or replicates: Median barcodes per oligo (BC threshold >= {mpra_oligo_data.barcode_threshold})",
            f"Median on all (flatten) Barcode counts where (BC threshold >= {mpra_oligo_data.barcode_threshold})",
        ],
    }

    return pd.DataFrame(metrics_data)


def median_rna_read_count(mpra_oligo_data):
    mpra_oligo_data.barcode_threshold = 1
    n_rna_replicate = []
    for replicate in mpra_oligo_data.obs_names:
        replicate_data = mpra_oligo_data.data[replicate, :]
        replicate_data = replicate_data[:, replicate_data.layers["barcode_counts"] >= mpra_oligo_data.barcode_threshold]
        n_rna_replicate += [np.median(replicate_data.layers["rna"])]

    all = mpra_oligo_data.rna_counts[mpra_oligo_data.barcode_counts >= mpra_oligo_data.barcode_threshold].flatten()

    metrics_data = {
        "Workflow": ["Experiment", "Experiment", "Experiment"],
        "Metric": [
            "median_rna_read_count",
            "median_rna_read_count",
            "median_rna_read_count",
        ],
        "Value": [int(np.median(n_rna_replicate)), int(np.min(n_rna_replicate)), int(np.median(all[all != 0]))],
        "Description": [
            f"Median or replicates: Median RNA counts per oligo (BC threshold >= {mpra_oligo_data.barcode_threshold})",
            f"Min or replicates: Median RNA counts per oligo (BC threshold >= {mpra_oligo_data.barcode_threshold})",
            f"Median on all (flatten) RNA counts where BC threshold >= {mpra_oligo_data.barcode_threshold}",
        ],
    }

    return pd.DataFrame(metrics_data)


if __name__ == "__main__":
    cli()
