# scTenifoldXct
Manifold learning to detect cell-cell interactions.

### Installation
```shell
git clone git@github.com:cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
```
Also, scTenifoldXct relies on scTenifoldpy (https://github.com/qwerty239qwe/scTenifoldpy) to construct networks, please install it by
```shell
pip install scTenifoldpy
```

### Dependencies Installation
```shell
conda create -n scTenifoldXct
conda activate scTenifoldXct
pip install -r requirements.txt
```

### Usages
```python
import scTenifoldXct
import scanpy
import numpy as np

adata = sc.datasets.paul15()[:, :100]  # raw counts
adata.layers['raw'] = np.asarray(adata.X, dtype=int)
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()

xct_obj = scTenifoldXct(data=adata, 
                        cell_names=['14Mo', '15Mo'],
                        obs_label="paul15_clusters",
                        rebuild_GRN=True, 
                        GRN_file_dir='./Net_for_Test', 
                        query_DB=None, verbose=True)
```
