import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpralib.mpradata import MPRAData, Modality, MPRAOligoData, MPRABarcodeData

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
custom_palette = sns.color_palette(["#72ACBF", "#BF2675", "#2ecc71", "#f1c40f", "#9b59b6"])
sns.set_theme(style="white", rc=custom_params, palette=custom_palette)


def correlation(data: MPRAData, layer: Modality, replicates=None) -> sns.PairGrid:

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
        idx = np.array([data.obs_names.get_loc(rep) for rep in replicates])
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
        color="#BF2675",
    )
    g.figure.suptitle("Median normalized counts across replicates")
    g.figure.subplots_adjust(top=0.95)
    return g


def barcodes_per_oligo(data: MPRAOligoData, replicates=None) -> sns.FacetGrid:

    bc_counts = data.barcode_counts.copy()

    intercept_median = np.median(bc_counts, axis=1)
    intercept_mean = np.mean(bc_counts, axis=1)

    if replicates:
        idx = np.array([data.obs_names.get_loc(rep) for rep in replicates])
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
    g.set_axis_labels("Frequency", "Barcodes per oligo")

    return g


def barcodes_outlier(data: MPRABarcodeData):

    # counts_dna = data.normalized_dna_counts.copy()
    # counts_rna = data.normalized_rna_counts.copy()
    counts_dna = data.dna_counts.copy()
    counts_rna = data.rna_counts.copy()

    counts_dna_sum = counts_dna.sum(axis=0)
    counts_rna_sum = counts_rna.sum(axis=0)

    mask = np.any(
        [
            np.any(data.barcode_counts < data.barcode_threshold, axis=0),
            np.any([counts_dna_sum <= 10, counts_rna_sum <= 0], axis=0),
            ~np.any(data.observed, axis=0),
        ],
        axis=0,
    )

    print(np.shape(mask))
    print(np.sum(mask))

    counts_dna_sum = np.ma.masked_array(counts_dna_sum, mask=mask)
    counts_rna_sum = np.ma.masked_array(counts_rna_sum, mask=mask)

    counts_ratio = np.log2(counts_rna_sum / counts_dna_sum)
    counts = pd.DataFrame(
        {
            "dna": counts_dna_sum,
            "ratio": counts_ratio,
            "oligo": data.oligos,
        },
    )

    counts = counts[~mask]

    print(len(counts))

    counts["ratio_med"] = counts.groupby("oligo")["ratio"].transform("median")
    counts["ratio_diff"] = counts["ratio"] - counts["ratio_med"]

    nbin = 20
    qs = np.quantile(np.log10(counts["dna"]), np.arange(0, nbin + 1) / nbin)
    counts["bin"] = pd.cut(
        np.log10(counts["dna"]), bins=qs, include_lowest=True, labels=[str(i) for i in range(0, len(qs) - 1)]
    )

    stats = counts.groupby("bin").agg(mean_diff=("ratio_diff", "mean"), sd_diff=("ratio_diff", "std")).reset_index()

    # Get the last n categories of ordered categories
    # Plotting
    plt.figure(figsize=(10, 6))
    sampled_counts = counts.sample(n=min(10000, len(counts)))
    if isinstance(sampled_counts, pd.Series):
        sampled_counts = sampled_counts.to_frame().T
    sns.scatterplot(data=sampled_counts, x="bin", y="ratio_diff", alpha=0.3)
    filtered_counts = counts[(counts["ratio_diff"] > 5) & (counts["bin"].isin(counts["bin"].cat.categories[9:]))]
    if isinstance(filtered_counts, pd.Series):
        filtered_counts = filtered_counts.to_frame().T
    sns.scatterplot(
        data=filtered_counts,
        x="bin",
        y="ratio_diff",
        alpha=1,
        color="red",
    )
    plt.errorbar(x=stats["bin"], y=stats["mean_diff"], yerr=2 * stats["sd_diff"], fmt="o", color="dodgerblue", linewidth=1)
    plt.axhline(y=5, color="red")
    plt.xlabel("Bin")
    plt.ylabel("Ratio Difference")
    plt.title("Ratio Difference by Bin")
    plt.show()
