import pytest
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy
import scipy.sparse

from scTenifoldXct.core import scTenifoldXct


@pytest.fixture(scope="session")
def ada_skin():
    data_path = Path(__file__).parent.parent / "./data/LS.h5ad"
    ada = sc.read_h5ad(data_path)
    data = scipy.sparse.csr_matrix.toarray(ada.X)
    counts = np.asarray(np.expm1(data), dtype=int)
    ada.layers['raw'] = counts
    ada.layers['log1p'] = data
    HVG_i = np.argsort(np.asarray(ada.var['vst.variance.standardized']))[-3000:]
    return ada[:, HVG_i]


@pytest.fixture(scope="session")
def xct_skin(ada_skin):
    return scTenifoldXct(data=ada_skin,
                         cell_names=['Inflam. FIB', 'Inflam. DC'],
                         obs_label="ident",
                         species="human",
                         rebuild_GRN=True,
                         GRN_file_dir='./skin_net',
                         verbose = True)


# small dataset
@pytest.fixture(scope="session")
def xct_paul15():
    ada = sc.datasets.paul15()[:, :100]  # raw counts
    ada.layers['raw'] = np.asarray(ada.X, dtype=int)
    sc.pp.log1p(ada)
    ada.layers['log1p'] = ada.X.copy()
    return scTenifoldXct(data=ada, cell_names=['14Mo', '15Mo'],
                         obs_label="paul15_clusters",
                         rebuild_GRN=True, GRN_file_dir='./Net_for_Test',
                         query_DB=None, verbose=True)
