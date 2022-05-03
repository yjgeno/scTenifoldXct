# scTenifoldXct
Manifold learning to detect cell-cell interactions.

### Installation with pip (coming later)
```shell
pip install scTenifoldXct 
```

### Installation from source
```shell
git clone https://github.com/cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
pip install .
```

### Dependencies Installation with `conda`:
```shell
conda env create -f environment.yml
conda activate scTenifold
```

### Usages

#### Quick Start
The following code runs scTenifoldXct on an example dataset located in the tutorials.
```python
import scanpy as sc
import scTenifoldXct as st

adata = sc.read_h5ad('data/adata_short_example.h5ad')
xct = st.scTenifoldXct(data = adata, 
                    cell_names = ['Inflam. FIB', 'Inflam. DC'],
                    obs_label = "ident",
                    rebuild_GRN = True, 
                    GRN_file_dir = './Net_example_dev',  
                    verbose = True,
                    n_cpus = -1)
emb = xct.get_embeds(train = True)
xct_pairs = xct.null_test()
print(xct_pairs)
```

### Tutorial
We have included two tutorial notebooks on scTenifoldXct usage and results visualization.

Single-sample interaction analysis:<br> https://github.com/cailab-tamu/scTenifoldXct/blob/master/tutorials/tutorial-short_example.ipynb <br>
Differential interaction analysis across samples:<br> https://github.com/cailab-tamu/scTenifoldXct/blob/master/tutorials/tutorial-merge_short_example.ipynb
