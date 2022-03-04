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
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
from matplotlib import pyplot as plt
from sklearn.mixture import GaussianMixture

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
plt.rcParams.update({"font.size": 16})

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
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
colors = [
    "#FAEBD7", "#00FFFF", "#FFD700", "#0000FF", "#FF8C00", "#EE82EE",
    "#9ACD32", "#5F9EA0", "#7FFF00", "#7FFFD4", "#6495ED", "#008B8B",
    "#B8860B", "#C0C0C0", "#000080", "#D8BFD8", "#00CED1", "#9400D3",
    "#8E804B", "#0089A7", "#CB1B45", "#FFB6C1", "#00FF00", "#800000",
    "#376B6D", "#D8BFD8", "#F5F5F5", "#D2691E"
]
try:
    scale_method = sys.argv[1]
    idx = sys.argv[2]
except IndexError:
    scale_method = "cpm"
    idx = "E165A"


def exp_trans(x, gamma):
    if gamma == 0:
        return x
    return (np.exp(gamma * x) - 1) / gamma


# %% read
for idx in idx_full:
    global_moran = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"results/global_moran/{idx}-{scale_method}-8.csv",
        ),
        index_col=0,
        header=0,
    )

    for n_components in [2, 3]:
        i_values = global_moran["I_value"].dropna()
        i_values = exp_trans(i_values.to_frame(), -6)
        gmm = GaussianMixture(n_components=n_components)
        results = pd.Series(gmm.fit_predict(i_values), index=i_values.index)

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_yticks([])
        ax.hist(i_values.to_numpy(), bins=100, density=True)
        x = np.linspace(i_values.min(), i_values.max(), 5000)
        ax.plot(x, np.exp(gmm.score_samples(x)), lw=2, label="GMM")
        for i in range(n_components):
            ax.plot(
                x,
                stats.norm.pdf(
                    x,
                    gmm.means_[i, 0],
                    gmm.covariances_[i, 0]**(1 / 2),
                ) * gmm.weights_[i],
                lw=2,
                ls="--",
                label=f"Gaussian {i}, weight {gmm.weights_[i]:.3f}",
            )
        ax.set_title(f"{idx} GMM PDF")
        plt.legend()
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/3/{idx}-{n_components}-pdf.jpg",
            ))
        results.to_csv(
            Path.joinpath(WORKDIR, f"results/3/{idx}-{n_components}.csv"))
