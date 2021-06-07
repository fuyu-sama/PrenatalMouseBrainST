#! venv/bin/python
# -*- encoding: utf-8 -*-

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#  Codes are far away from bugs with Buddha's bless
#

# %% environment config
from pathlib import Path

import numpy as np
import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import ranksums
from statsmodels.stat import multitest

session_info.show()
plt.rcParams.update({"font.size": 24})

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    "E175B": "V10M17-085-E175B",
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
colors = [
    "#FAEBD7",
    "#00FFFF",
    "#FFD700",
    "#0000FF",
    "#FF8C00",
    "#EE82EE",
    "#9ACD32",
    "#5F9EA0",
    "#7FFF00",
    "#7FFFD4",
    "#6495ED",
    "#008B8B",
    "#B8860B",
    "#C0C0C0",
    "#000080",
    "#D8BFD8",
    "#00CED1",
    "#9400D3",
    "#8E804B",
    "#0089A7",
    "#CB1B45",
    "#FFB6C1",
    "#00FF00",
    "#800000",
    "#376B6D",
    "#D8BFD8",
    "#F5F5F5",
    "#D2691E",
]
timepoints = {
    "E135": [["E135A"], ["E135B"]],
    "E155": [["E155A"], ["E155B"]],
    "E165": [["E165A"], ["E165B"]],
    "E175": [["E175A1", "E175A2"], ["E175B"]],
    "P0": [["P0A1"], ["P0A2"]],
}


# %% helper func
def quantile_normalization(df: pd.DataFrame) -> pd.DataFrame:
    rank_mean = df.stack().groupby(
        df.rank(method="first").stack().astype(int)).mean()
    return df.rank(method="min").stack().astype(int).map(rank_mean).unstack()


# %% read data
count_full_df = pd.DataFrame()
for idx in idx_full:
    count_path = Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/logcpm/scale_df/logcpm/{idx}-logcpm.csv")
    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    print(f"{idx}:\t{count_df.shape}")
    count_full_df = pd.concat(
        [count_full_df, count_df],
        axis=0,
        verify_integrity=True,
    )
count_full_df = count_full_df.dropna(axis=1)
print(f"Full:\t{count_full_df.shape}")
for idx in idx_full:
    count_path = Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/logcpm/scale_df/logcpm/{idx}-logcpm-inter.csv",
    )
    count_df = count_full_df.reindex(
        index=[i for i in count_full_df.index if idx in i])
    count_df.T.to_csv(count_path)
count_full_df.T.to_csv(
    Path.joinpath(
        Path.home(),
        "workspace/mouse-brain-full/logcpm/scale_df/logcpm/full-logcpm-inter.csv"
    ))
# count_full_df = quantile_normalization(count_full_df.T).T

# %% draw
q = 0.9
count_full_sub = pd.DataFrame()
count_e175 = pd.DataFrame()
for timepoint in timepoints:
    idx = timepoints[timepoint]
    count_sub_0 = count_full_df.reindex(
        index=[i for i in count_full_df.index if i.split("_")[0] in idx[0]])
    count_sub_1 = count_full_df.reindex(
        index=[i for i in count_full_df.index if i.split("_")[0] in idx[1]])
    count_mean_0 = count_sub_0.mean()
    count_mean_1 = count_sub_1.mean()

    # subset by distance
    # d = abs((count_mean_1 - count_mean_0) / (math.sqrt(2)))
    # c = d.where(d <= d.quantile(q), np.nan)
    # c = c.where(np.isnan(c), 1)
    # c = c.replace(np.nan, 0)

    # subset by n_sigma
    # d = abs((count_mean_1 - count_mean_0) / (math.sqrt(2)))
    # c = d.where(d <= (d.mean() + q * d.std()), np.nan)
    # c = c.where(d >= (d.mean() - q * d.std()), np.nan)
    # c = c.where(np.isnan(c), 1)
    # c = c.replace(np.nan, 0)

    # subset by statistic test
    d = pd.Series(dtype=np.float64)
    for gene in count_sub_0.columns:
        d[gene] = ranksums(count_sub_0[gene], count_sub_1[gene])[1]
    d = multitest.fdrcorrection(d, alpha=1e-8)[0]
    c = np.where(d, 0, 1)
    c = pd.Series(c, index=count_sub_0.columns)

    fig, ax_plot = plt.subplots(1, figsize=(10, 10))
    ax_plot.scatter(
        count_mean_0,
        count_mean_1,
        c=c,
        cmap=ListedColormap(["#a3a3c2", "#cc6600"]),
        s=4,
    )
    ax_plot.set_xlabel(" ".join(idx[0]))
    ax_plot.set_ylabel(" ".join(idx[1]))
    r_all = np.corrcoef(
        count_mean_0,
        count_mean_1,
    )[0, 1]
    r_part = np.corrcoef(
        count_mean_0[c > 0],
        count_mean_1[c > 0],
    )[0, 1]
    ax_plot.text(0, 5, f"Pearson for all genes: {r_all:.3f}")
    ax_plot.text(0, 4.7, f"Pearson for selected genes: {r_part:.3f}")
    ax_plot.text(0, 4.4, f"Number of selected genes: {int(c.sum())}")
    ax_plot.plot(range(6), range(6), "--r")
    ax_plot.text(4.7, 5, "y = x", c="r")
    ax_plot.set_title("logcpm ranksums")
    fig.savefig(
        Path.joinpath(
            Path.home(),
            f"workspace/mouse-brain-full/results/correlation/{timepoint}-correlation.jpg",
        ))
