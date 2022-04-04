import pytest

from scTenifoldXct.nn import ManifoldAlignmentNet


@pytest.fixture(scope="session")
def df_nn_skin(xct_skin):
    xct_skin.train_nn(n_steps= 1000, lr = 0.001)
    return xct_skin.null_test(pct=0.025, plot_result=True)


@pytest.fixture(scope="session")
def df_nn_paul15(xct_paul15):
    xct_paul15.train_nn(n_steps= 1000, lr = 0.001)
    return xct_paul15.null_test(pct=0.025, plot_result=True)