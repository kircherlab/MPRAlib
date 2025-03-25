import numpy as np
import pandas as pd
import seaborn as sns
from mpralib.mpradata import MPRAData, Modality, MPRAOligoData

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
custom_palette = sns.color_palette(["#72ACBF", "#BF2675", "#2ecc71", "#f1c40f", "#9b59b6"])
sns.set_theme(style="whitegrid", rc=custom_params, palette=custom_palette)


def correlation(data: MPRAData, layer: Modality, replicates=None) -> sns.JointGrid:

    counts = None
    if layer == Modality.DNA:
        counts = data.dna_counts.copy()
    elif layer == Modality.RNA:
        counts = data.rna_counts.copy()
    elif layer == Modality.RNA_NORMALIZED:
        counts = data.normalized_rna_counts.copy()
    elif layer == Modality.DNA_NORMALIZED:
        counts = data.normalized_dna_counts.copy()
    elif layer == Modality.ACTIVITY:
        counts = data.activity.copy()

    counts = np.ma.masked_array(counts, mask=[data.barcode_counts < data.barcode_threshold])

    if replicates:
        idx = [data.obs_names.get_loc(rep) for rep in replicates]
        counts = pd.DataFrame(counts[idx].T, columns=[f"Replicate {i}" for i in data.obs_names[idx]], index=data.var_names)
    else:
        counts = pd.DataFrame(counts.T, columns=[f"Replicate {i}" for i in data.obs_names], index=data.var_names)

    g = sns.PairGrid(counts)
    g.map_upper(sns.scatterplot, s=1)
    g.map_diag(sns.kdeplot)
    g.map_lower(sns.kdeplot, fill=True)
    g.figure.suptitle("Correlation Plot")
    return g


def dna_vs_rna(data: MPRAData, replicates=None) -> sns.JointGrid:

    counts_dna = data.normalized_dna_counts.copy()
    counts_rna = data.normalized_rna_counts.copy()

    mask = [data.barcode_counts < data.barcode_threshold]

    counts_dna = np.ma.masked_array(counts_dna, mask=mask)
    counts_rna = np.ma.masked_array(counts_rna, mask=mask)

    if replicates:
        idx = [data.obs_names.get_loc(rep) for rep in replicates]
        counts_dna = counts_dna[idx]
        counts_rna = counts_rna[idx]

    median_dna = np.median(counts_dna, axis=0)
    median_rna = np.median(counts_rna, axis=0)

    median_dna = np.ma.masked_equal(median_dna, 0)
    median_rna = np.ma.masked_equal(median_rna, 0)

    df = pd.DataFrame({"DNA [log10]": np.log10(median_dna), "RNA [log10]": np.log10(median_rna)})

    g = sns.jointplot(data=df, x="DNA [log10]", y="RNA [log10]", kind="scatter", s=3)
    g.ax_joint.plot(
        [df["DNA [log10]"].min(), df["DNA [log10]"].max()],
        [df["RNA [log10]"].min(), df["RNA [log10]"].max()],
        linestyle="--",
        color="red",
    )
    g.figure.suptitle("Median normalized counts across replicates")
    g.figure.subplots_adjust(top=0.95)
    return g


def barcodes_per_oligo(data: MPRAOligoData, replicates=None) -> sns.JointGrid:

    bc_counts = data.barcode_counts.copy()

    intercept_median = np.median(bc_counts, axis=1)
    intercept_mean = np.mean(bc_counts, axis=1)

    if replicates:
        idx = [data.obs_names.get_loc(rep) for rep in replicates]
        bc_counts = pd.DataFrame(
            bc_counts[idx].T, columns=[f"Replicate {i}" for i in data.obs_names[idx]], index=data.var_names
        )
    else:
        replicates = data.obs_names
        bc_counts = pd.DataFrame(bc_counts.T, columns=[f"Replicate {i}" for i in data.obs_names], index=data.var_names)

    bc_counts = bc_counts.melt(var_name="replicate", value_name="n_bc")

    g = sns.FacetGrid(bc_counts, col="replicate")

    # Histogram plot
    g.map(sns.histplot, "n_bc")

    i = 0
    for ax in g.axes.flatten():
        ax.axvline(x=intercept_median[i], color="red", linestyle="--")
        ax.axvline(x=intercept_mean[i], color="blue", linestyle="--")
        ax.set_title(f"Replicate {replicates[i]}")
        ax.text(ax.get_xlim()[1] * 0.5, ax.get_ylim()[1] * 0.9, f"{intercept_median[i]:.0f}", color="red", ha="left")
        ax.text(ax.get_xlim()[1] * 0.5, ax.get_ylim()[1] * 0.8, f"{intercept_mean[i]:.2f}", color="blue", ha="left")
        i = i + 1
    g.set_axis_labels("Frequency", "arcodes per oligo")

    return g
