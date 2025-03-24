
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpralib.mpradata import MPRAData, CountType

custom_params = {"axes.spines.right": False, "axes.spines.top": False}
custom_palette = sns.color_palette(["#72ACBF", "#BF2675", "#2ecc71", "#f1c40f", "#9b59b6"])
sns.set_theme(style="whitegrid", rc=custom_params, palette=custom_palette)


def correlation(data: MPRAData, layer: CountType, x=None, y=None) -> Figure:

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

    cmap = sns.light_palette('#BF2675', as_cmap=True)

    if x is None or y is None:
        g = sns.PairGrid(
            pd.DataFrame(counts.T, columns=[f"Replicate {i}" for i in data.obs_names], index=data.var_names)
        )

        g.map_upper(sns.scatterplot, s=1)
        g.map_diag(sns.kdeplot)
        g.map_lower(sns.kdeplot, fill=True)
        g.figure.suptitle("Correlation Plot")
        return g

    else:
        fig, ax = plt.subplots()

        idx_x = data.obs_names.get_loc(x)
        idx_y = data.obs_names.get_loc(y)

        sns.scatterplot(x=counts[idx_x], y=counts[idx_y], s=5)
        sns.histplot(x=counts[idx_x], y=counts[idx_y], bins=50, pthresh=.1, cmap=cmap)

        ax.set_title("Correlation Plot")
        ax.set_xlabel(f'Replicate {x}')
        ax.set_ylabel(f'Replicate {y}')

        # Return the figure object
        return fig
