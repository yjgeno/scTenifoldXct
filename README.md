# scTenifoldXct
Manifold learning to detect cell-cell interactions.
<br/>
<p align="center">
    <img src="LS_git.jpeg" alt="drawing" width="300"/>
</p>
<br/>

### Install Dependencies
We suggest first intall dependencies of scTenifoldXct with `conda`:
```shell
conda env create -f environment.yml
conda activate scTenifold
```
<br/>

### Install scTenifoldXct
Install scTenifoldXct with `pip`:
```shell
pip install git+https://github.com/cailab-tamu/scTenifoldXct.git 
```

or install it manually from source:
```shell
git clone https://github.com/cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
pip install .
```
<br/>

### Usages

#### Quick Start
The following code runs scTenifoldXct on an example dataset located in the tutorials.
```python
import scanpy as sc
import scTenifoldXct as st

adata = sc.read_h5ad('data/adata_short_example.h5ad') # load data
xct = st.scTenifoldXct(data = adata, # an AnnData 
                    source_celltype = 'Inflam. FIB', # sender cell type
                    target_celltype = 'Inflam. DC', # receiver cell type
                    obs_label = 'ident', # colname in adata.obs indicating cell types
                    rebuild_GRN = True, # whether to build GRNs
                    GRN_file_dir = 'Net_example_dev',  # folder path to GRNs
                    verbose = True, # whether to verbose the processing
                    n_cpus = -1) # CPU multiprocessing, -1 to use all
emb = xct.get_embeds(train = True) # Manifold alignment to project data to low-dimensional embeddings
xct_pairs = xct.null_test() # non-parametric test to get significant interactions
print(xct_pairs)
```

### Tutorial
We have included two tutorial notebooks on scTenifoldXct usage and results visualization.

Single-sample interaction analysis:<br> https://github.com/cailab-tamu/scTenifoldXct/blob/master/tutorials/tutorial-short_example.ipynb <br>
Differential interaction analysis across samples:<br> https://github.com/cailab-tamu/scTenifoldXct/blob/master/tutorials/tutorial-merge_short_example.ipynb
