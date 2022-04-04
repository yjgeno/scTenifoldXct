import pytest

import itertools

import pandas as pd
import numpy as np

from scTenifoldXct.core import null_test


def generate_fake_df_nn(n_ligand=3000, n_receptors=3000, n_cands=200):
    gene_names = [f"GENE{i}" for i in range(max(n_ligand, n_receptors))]
    iteration = itertools.product(gene_names, gene_names)
    inds, ligands, receptors = [], [], []
    for i, j in iteration:
        inds.append(f"{i}_{j}")
        ligands.append(i)
        receptors.append(j)
    df = pd.DataFrame({"ligand": ligands,
                       "receptor": receptors,
                       "dist": np.random.chisquare(1, (n_ligand * n_receptors,)),
                       "correspondence": np.random.lognormal(0, 4, size=(n_ligand * n_receptors,))},
                      index=inds)
    return df, np.random.choice(df.index, size=(n_cands,), replace=False)


@pytest.mark.parametrize("df_nn,candidates", [
    generate_fake_df_nn(3000, 3000, 200),
    generate_fake_df_nn(1000, 1000, 200),
])
@pytest.mark.parametrize("filter_zeros", [True])
def test_null_test(df_nn, candidates, filter_zeros):
    null_test(df_nn=df_nn, candidates=candidates, filter_zeros=filter_zeros)


def test_chi2_test(xct_skin):
    xct_skin.train_nn(n_steps= 1000, lr = 0.001)
    xct_skin.chi2_test(dof=3, pval=0.05, cal_FDR=True, plot_result=True)