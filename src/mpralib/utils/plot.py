
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpralib.mpradata import MPRAData, CountType

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
custom_palette = sns.color_palette(["#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6"])
sns.set_theme(style="whitegrid", rc=custom_params, palette=custom_palette)


def correlation(data: MPRAData, layer: CountType, x, y) -> Figure:
    fig, ax = plt.subplots()

    idx_x = data.obs_names.get_loc(x)
    idx_y = data.obs_names.get_loc(y)

    counts = None
    if layer == CountType.DNA:
        counts = data.dna_counts.copy()
    elif layer == CountType.RNA:
        counts = data.rna_counts.copy()
    elif layer == CountType.RNA_NORMALIZED:
        counts = data.normalized_rna_counts.copy()
    elif layer == CountType.DNA_NORMALIZED:
        counts = data.normalized_dna_counts.copy()
    elif layer == CountType.ACTIVITY:
        counts = data.activity.copy()

    counts = np.ma.masked_array(counts, mask=[data.barcode_counts < data.barcode_threshold])

    sns.scatterplot(x=counts[idx_x], y=counts[idx_y], s=5, color=".15")
    sns.histplot(x=counts[idx_x], y=counts[idx_y], bins=50, pthresh=.1, cmap="mako")
    sns.kdeplot(x=counts[idx_x], y=counts[idx_y], levels=5, color="w", linewidths=1)

    ax.set_title("Correlation Plot")
    ax.set_xlabel(f'Replicate {x}')
    ax.set_ylabel(f'Replicate {y}')

    # Return the figure object
    return fig
