# scTenifoldXct
Manifold learning to detect cell-cell interactions.

### Dependencies Installation
```shell
conda create -n scTenifold 
conda activate scTenifold
pip install -r requirements.txt
```

### Installation
```shell
git clone git@github.com:cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
```

### Usages
```python
import scTenifoldXct

xct = scTenifoldXct(data = adata, 
                    cell_names = ['Inflam. FIB', 'Inflam. DC'],
                    obs_label = "ident",
                    rebuild_GRN = True, 
                    GRN_file_dir = './Net_example',  
                    verbose = True)
```
