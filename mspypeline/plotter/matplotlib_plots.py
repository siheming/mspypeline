import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def save_volcano_results(plot_data):
    plot_data.to_csv(os.path.join(self.file_dir_volcano,
                                  f"volcano_plot_data_{g1}_vs_{g2}_full_{col_mapping[col].replace(' ', '_')}.csv"))
    # save significant values
    plot_data[plot_data[col] < 0.05].to_csv(
        os.path.join(self.file_dir_volcano,
                     f"volcano_plot_data_{g1}_vs_{g2}_significant_{col_mapping[col].replace(' ', '_')}.csv"))
    # calculate mean intensity for unique genes for plotting
    unique_g1 = v1[exclusive_1].mean(axis=1).rename(f"{df_to_use} mean intensity")
    unique_g2 = v2[exclusive_2].mean(axis=1).rename(f"{df_to_use} mean intensity")
    # save unique values
    unique_g1.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g1}.csv"),
                     header=True)
    unique_g2.to_csv(os.path.join(self.file_dir_volcano, f"volcano_plot_data_{g1}_vs_{g2}_unique_{g2}.csv"),
                     header=True)

    def get_volcano_significances(fchange, pval, fchange_threshold):
        if pval > 0.05 or abs(fchange) < np.log2(fchange_threshold):
            return "ns"
        elif fchange >= 0:
            return "up"
        elif fchange < 0:
            return "down"
        else:
            raise ValueError("heisenbug")

    significance_to_color = {"ns": "gray", "up": "red", "down": "blue"}
    significance_to_label = {"ns": "non-significant", "up": f"upregulated in {g2}", "down": f"upregulated in {g1}"}

    # plot
    fig = plt.figure(figsize=(7, 7))

    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 8, 1])
    ax_unique_down = plt.subplot(gs[0])
    ax = plt.subplot(gs[1])
    ax_unique_up = plt.subplot(gs[2])

    # hide the spines between ax and ax2
    ax_unique_down.spines['right'].set_visible(False)
    ax_unique_up.spines['left'].set_visible(False)
    ax_unique_down.yaxis.tick_left()
    ax_unique_up.yaxis.tick_right()
    ax_unique_up.yaxis.set_label_position("right")
    # hide the xticks
    ax_unique_down.tick_params(which='both', bottom=False, labelbottom=False)
    ax_unique_up.tick_params(which='both', bottom=False, labelbottom=False)

    # non sign gray, left side significant blue right side red
    significances_plot = [get_volcano_significances(log_fold_change, p_val, fchange_threshold)
                          for log_fold_change, p_val in zip(plot_data["logFC"], plot_data[col])]
    for regulation in significance_to_color:
        mask = [x == regulation for x in significances_plot]
        ax.scatter(plot_data["logFC"][mask], -np.log10(plot_data[col])[mask], s=scatter_size,
                   color=significance_to_color[regulation], label=f"{sum(mask)} {significance_to_label[regulation]}")
    ymin, ymax = ax.get_ybound()
    xmin, xmax = ax.get_xbound()
    m = max(abs(xmin), xmax)
    ax.set_xlim(left=-1 * m, right=m)
    xmin, xmax = ax.get_xbound()
    x_offset = (np.log2(fchange_threshold) / xmax) / 2
    # add line at significance threshold
    axline_kwargs = dict(linestyle="--", color="black", alpha=0.5, linewidth=1)
    if any(plot_data[col] < 0.05):
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0, xmax=0.5 - x_offset)
        ax.axhline(-np.log10(0.05), **axline_kwargs, xmin=0.5 + x_offset, xmax=1)

    # add lines for minimum fold change threshold
    y_percentage = (-np.log10(0.05) + abs(ymin)) / (ymax + abs(ymin))
    if fchange_threshold > 0:
        ax.axvline(-np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
        ax.axvline(np.log2(fchange_threshold), **axline_kwargs, ymin=y_percentage, ymax=1)
    # plot unique values with mean intensity at over maximum
    ax_unique_down.scatter([0] * len(unique_g1), unique_g1, s=scatter_size, color="dodgerblue",
                           label=f"{len(unique_g1)} unique in {g1}")
    ax_unique_up.scatter([0] * len(unique_g2), unique_g2, s=scatter_size, color="coral",
                         label=f"{len(unique_g2)} unique in {g2}")
    # adjust bounds
    ymin_down, ymax_down = ax_unique_down.get_ybound()
    ymin_up, ymax_up = ax_unique_up.get_ybound()
    ax_unique_down.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))
    ax_unique_up.set_ylim(bottom=min(ymin_down, ymin_up), top=max(ymax_down, ymax_up))

    # figure stuff
    if show_suptitle:
        fig.suptitle(f"{g1} vs {g2}")
    ax.set_xlabel(r"$Log_2$ Fold Change")
    ax.set_ylabel(r"-$Log_{10}$" + f" {col_mapping[col]}")
    ax_unique_down.set_ylabel(self.intensity_label_names[df_to_use])
    ax_unique_up.set_ylabel(self.intensity_label_names[df_to_use])
    fig.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False)
    res_path = os.path.join(self.file_dir_volcano,
                            f"volcano_{g1}_{g2}_no_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(res_path, dpi=200, bbox_inches="tight")
    significant_upregulated = plot_data[
        (plot_data["logFC"] > np.log2(fchange_threshold)) & (plot_data[col] < 0.05)].sort_values(by=[col],
                                                                                                 ascending=True).head(
        10)
    significant_downregulated = plot_data[
        (plot_data["logFC"] < -np.log2(fchange_threshold)) & (plot_data[col] < 0.05)].sort_values(by=[col],
                                                                                                  ascending=True).head(
        10)
    significant = pd.concat([significant_upregulated, significant_downregulated])
    texts = []
    for log_fold_change, p_val, gene_name in zip(significant["logFC"], significant[col], significant.index):
        texts.append(ax.text(log_fold_change, -np.log10(p_val), gene_name, ha="center", va="center", fontsize=8))
    adjust_text(texts, arrowprops=dict(width=0.15, headwidth=0, color='gray', alpha=0.6), ax=ax)
    res_path = os.path.join(self.file_dir_volcano,
                            f"volcano_{g1}_{g2}_annotation_{col_mapping[col].replace(' ', '_')}" + FIG_FORMAT)
    fig.savefig(res_path, dpi=200, bbox_inches="tight")
    # TODO scatter plot of significant genes
