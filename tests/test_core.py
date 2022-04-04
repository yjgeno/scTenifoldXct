import pytest
from scTenifoldXct.core import GRN


def test_xct_obj_attrs(xct_skin):
    assert isinstance(xct_skin.net_A, GRN), "net_A should be a GRN object"
    assert isinstance(xct_skin.net_B, GRN), "net_B should be a GRN object"
    assert xct_skin.net_B.shape == xct_skin.net_A.shape
