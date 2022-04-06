# scTenifoldXct
Manifold learning to detect cell-cell interactions.

### Installation
```shell
git clone https://github.com/cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
```

### Dependencies Installation
```shell
conda create -n scTenifold 
conda activate scTenifold
pip install -r requirements.txt
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
