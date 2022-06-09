from pathlib import Path
from anndata import (
    AnnData,
    read_h5ad,
    read_csv,
    read_excel,
    read_hdf,
    read_loom,
    read_mtx,
    read_text,
)
import scanpy as sc
import pandas as pd
from typing import Union
from scipy import sparse

# https://scanpy.readthedocs.io/en/stable/api.html#reading
file_attrs = ["h5ad", "csv", "xlsx", "h5", "loom", "mtx", "txt"]
read_fun_keys = ["h5ad", "csv", "excel", "hdf", "loom", "mtx", "text"]


def read_counts(counts_path: str, as_sparse: bool = True, **kwargs):
    file_attr = counts_path.split(".")[-1]
    if Path(counts_path).is_file() and file_attr in file_attrs:
        print(f"loading counts from {counts_path}...")
        read_fun_key = read_fun_keys[file_attrs.index(file_attr)]
        read_fun = getattr(sc, f"read_{read_fun_key}")  # define sc.read_{} function
        adata = read_fun(counts_path, **kwargs).transpose() # transpose counts file
    else:
        raise ValueError("incorrect file path given to counts")
    if as_sparse:
        print("making counts sparse...")
        adata.X = (
            sparse.csr_matrix(adata.X) if not sparse.issparse(adata.X) else adata.X
        )
    return adata


def build_adata(
    counts_path: str,
    meta_gene_path: Union[None, str],
    meta_cell_path: Union[None, str],
    sep="\t",
    header = None,
    log_normalize: bool = True,
    **kwargs,
):
    """
    Args:
        counts_path (str): Path to counts.
        meta_gene (Union[None, str]): Path to metadata of variables.
        meta_cell (Union[None, str]): Path to metadata of Observations
        sep (str, optional): The delimiter for metadata. Defaults to "\t".
        log_normalize (bool, optional): Whether log-normalize the counts. Defaults to True.

    Returns:
        AnnData
    """
    adata = read_counts(counts_path, **kwargs)
   
    if meta_gene_path is not None and Path(meta_gene_path).is_file():
        try:
            print("adding metadata for variables...")
            adata.var = pd.read_csv(meta_gene_path, header=header, sep=sep)
        except Exception:
            raise ValueError("incorrect file path given to meta_gene")
    if meta_cell_path is not None and Path(meta_cell_path).is_file():
        try:
            print("adding metadata for observations...")
            adata.obs = pd.read_csv(meta_cell_path, header=header, sep=sep)
        except Exception:
            raise ValueError("incorrect file path given to meta_cell")

    if log_normalize:
        print("normalizing counts and save it to adata.layers[\"log1p\"]")
        adata.layers["raw"] = adata.X
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["log1p"] = adata.X  # save log normalized counts to adata.X

    return adata


if __name__ == "__main__":
    from os import path

    data_dir = path.join(
        path.dirname(path.dirname(path.abspath(__file__))), "tutorials", "data", "filtered_gene_bc_matrices"
    )
    counts_path = str(path.join(data_dir, "matrix.mtx"))
    gene_path = str(path.join(data_dir, "genes.tsv"))
    cell_path = str(path.join(data_dir, "barcodes.tsv"))
    adata = build_adata(counts_path, gene_path, cell_path)
    print(adata)