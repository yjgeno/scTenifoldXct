from typing import List
import itertools
import os
from os import PathLike
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from anndata._core.views import ArrayView
import anndata
import scipy
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from scTenifold import cal_pcNet


warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

try:
    from scTenifoldXct.pcNet import pcNet
except ImportError:
    print('module \'pcNet\' not found')

class Xct_metrics():
    '''require adata with layer 'raw' (counts) and 'log1p' (normalized), cell labels in obs 'ident' '''

    # __slots__ = ('_specis', 'genes', 'LRs', '_genes_index_DB', 'TFs')
    def __init__(self, adata, specis = 'Human'):
        if not ('raw' and 'log1p' in adata.layers.keys()):
            raise IndexError('require adata with count and normalized layers named \'raw\' and \'log1p\'')
        else:
            self._specis = specis
            self.genes = adata.var_names
            self.LRs = self.LR_DB()
            self._genes_index_DB = self.get_index(DB = self.LRs)
            self.TFs = self.TF_DB()

    def LR_DB(self):
        '''load omnipath DB for L-R pairs'''
        LR = pd.read_csv('https://raw.githubusercontent.com/yjgeno/Xct/dev/DB/omnipath_intercell_toUse_v2.csv')
        LR_toUse = LR[['genesymbol_intercell_source', 'genesymbol_intercell_target']]
        LR_toUse.columns = ['ligand', 'receptor']
        if self._specis == 'Mouse':
            for col in LR_toUse.columns:
                LR_toUse[col] = LR_toUse[col].str.capitalize()
        elif self._specis == 'Human':
            pass
        else:
            raise KeyError('current ligand/receptor DB only supports \'Mouse\' and \'Human\'')
        del LR

        return LR_toUse

    def TF_DB(self):
        '''load TFome DB for TFs'''
        TFs = pd.read_csv('https://raw.githubusercontent.com/yjgeno/Xct/dev/DB/41587_2020_742_MOESM3_ESM.csv', header=None, names=['TF_symbol'])
        if self._specis == 'Mouse':
            TFs['TF_symbol'] = TFs['TF_symbol'].str.capitalize()
        elif self._specis == 'Human':
            pass
        else:
            raise KeyError('current transcript factor DB only supports \'Mouse\' and \'Human\'')
        return TFs

    def subset(self):
        '''subset adata var with only DB L and R'''
        genes = np.ravel(self.LRs.to_numpy())
        genes = np.unique(genes[genes != None])
        genes_use = self.genes.intersection(genes)
        return [list(self.genes).index(g) for g in genes_use]  # index in orig adata

    def get_index(self, DB):
        '''original index of DB L-R pairs in adata var'''
        g_LRs = DB.iloc[:, :2].to_numpy()  # L-R
        gene_list = [None] + list(self.genes)

        gene_index = np.zeros(len(np.ravel(g_LRs)), dtype = int)
        for g in gene_list:
            g_index = np.asarray(np.where(np.isin(np.ravel(g_LRs), g)))
            if g_index.size == 0:
                continue
            else:
                for i in g_index:
                    gene_index[i] = gene_list.index(g)
        genes_index_DB = np.array(gene_index).reshape(g_LRs.shape)  # gene index refer to subset adata var + 1

        return genes_index_DB

    def get_metric(self, adata: ArrayView, verbose = False):  # require normalized data
        '''compute metrics for each gene'''
        data_norm = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X.copy()  # adata.layers['log1p']
        if verbose:
            print('(cell, feature):', data_norm.shape)

        if (data_norm % 1 != 0).any():  # check space: True for log (float), False for counts (int)
            mean = np.mean(data_norm, axis = 0  )  # .toarray()
            var = np.var(data_norm, axis = 0  )  # .toarray()
            # mean[mean == 0] = 1e-12
            # dispersion = var / mean
            # cv = np.sqrt(var) / mean

            return mean, var  # , dispersion, cv
        else:
            raise ValueError("require log data")


class Xct(Xct_metrics):

    def __init__(self, adata, CellA, CellB, specis = 'Human', build_GRN = False, save_GRN = False, pcNet_name = 'pcNet', queryDB = None, verbose = False):
        '''build_GRN: if True to build GRN thru pcNet, if False to load built GRN files;
            save_GRN: save constructed 2 pcNet;
            pcNet_name: name of GRN (.csv) files, read/write;
            queryDB: 3 modes to construct correspondence w: None, 'comb', 'pairs'
            '''
        super().__init__(adata, specis = specis)
        self._cell_names = CellA, CellB
        self._metric_names = ['mean', 'var']

        if not ('ident' in adata.obs.keys()):
            raise IndexError('require adata with cell labels saved in \'ident\'')
        else:
            ada_A = adata[adata.obs['ident'] == CellA, :] # view
            ada_B = adata[adata.obs['ident'] == CellB, :] # view
        self._cell_numbers = ada_A.shape[0], ada_B.shape[0]
        self.genes_names = list(ada_A.var_names.astype(str)), list(ada_B.var_names.astype(str))
        self._X = ada_A.X, ada_B.X # input array for nn projection
        if verbose:
            print \
                (f'init an Xct object for interactions from {self._cell_names[0]} ({self._cell_numbers[0]}) to {self._cell_names[1]} ({self._cell_numbers[1]})...')

        self._metric_A = self.get_metric(ada_A)
        self._metric_B = self.get_metric(ada_B)

        self.ref = self.fill_metric()
        print(self.ref)

        self.genes_index = self.get_index(DB = self.ref)
        print(self.genes_index)

        pcNet_path_A = f'./data/{pcNet_name}_A.csv'
        pcNet_path_B = f'./data/{pcNet_name}_B.csv'
        if build_GRN:
            if verbose:
                print('building GRN...')
            if os.path.isfile(pcNet_path_A):
                self._net_A = np.genfromtxt(pcNet_path_A, delimiter="\t")
            else:
                self._net_A = pcNet(ada_A.X, nComp=5, symmetric=True)
            if verbose:
                print('GRN of Cell A has been built, start building GRN of Cell B...')

            if os.path.isfile(pcNet_path_B):
                self._net_B = np.genfromtxt(pcNet_path_B, delimiter="\t")
            else:
                self._net_B = pcNet(ada_B.X, nComp=5, symmetric=True)
            if verbose:
                print('GRN of Cell B has been built, building correspondence...')

            if save_GRN:
                os.makedirs('./data', exist_ok = True) # create dir 'data'
                np.savetxt(pcNet_path_A, self._net_A, delimiter="\t")
                np.savetxt(pcNet_path_B, self._net_B, delimiter="\t")
        else:
            try:
                if verbose:
                    print('loading GRNs...')
                self._net_A = np.genfromtxt(pcNet_path_A, delimiter="\t")
                self._net_B = np.genfromtxt(pcNet_path_B, delimiter="\t")
            except ImportError:
                print('require pcNet_name where csv files saved in with tab as delimiter')
        if verbose:
            print('building correspondence...')
        self._w = self._build_w(alpha = 0.55, queryDB = queryDB, scale_w = True)
        if verbose:
            print('init completed.\n')
        del ada_A, ada_B

    def __str__(self):
        info = f'Xct object for interactions from {self._cell_names[0]} ({self._cell_numbers[0]}) to {self._cell_names[1]} ({self._cell_numbers[1]})'
        if '_w' in dir(self):
            return info + f'\n# of genes = {len(self.genes_names[0])} X {len(self.genes_names[1])} \nCorrespondence = {self._w.shape[0]} X {self._w.shape[1]}'
        else:
            return info

    def fill_metric(self, ref_obj = None, verbose = False):
        '''fill the corresponding metrics for genes of selected pairs (L-R candidates)'''
        if ref_obj is None:
            genes_index = self._genes_index_DB
        else:
            genes_index = ref_obj.genes_index

        index_L = genes_index[:, 0]
        index_R = genes_index[:, 1]

        df = pd.DataFrame()

        for metric_A, metric_B, metric in zip(self._metric_A, self._metric_B, self._metric_names):
            filled_L = np.array([0 if i == 0 else metric_A[ i -1] for i in index_L], dtype = float)
            filled_R = np.array([0 if i == 0 else metric_B[ i -1] for i in index_R], dtype = float)

            filled = np.concatenate((filled_L[:, None], filled_R[:, None]), axis=1)
            result = pd.DataFrame(data = filled, columns = [f'{metric}_L', f'{metric}_R'])
            df = pd.concat([df, result], axis=1)

        if ref_obj is None:
            df = pd.concat([self.LRs, df], axis=1) # concat 1:1 since sharing same index
            mask1 = (df['mean_L'] > 0) & (df['mean_R'] > 0) # filter 0 (none or zero expression) of LR
            df = df[mask1]

        else:
            ref_DB = self.LRs.iloc[ref_obj.ref.index, :].reset_index(drop = True, inplace = False)  # match index
            df = pd.concat([ref_DB, df], axis=1)
            df.set_index(pd.Index(ref_obj.ref.index), inplace = True)

        if verbose:
            print('Selected {} LR pairs'.format(df.shape[0]))

        return df

    def _build_w(self, alpha, queryDB = None, scale_w = True):
        '''build w: 3 modes, default None will not query the DB and use all pair-wise corresponding scores'''
        # (1-a)*u^2 + a*var
        metric_A_temp = (( 1 -alpha )* np.square(self._metric_A[0]) + alpha* (self._metric_A[1]))[:, None]
        metric_B_temp = (( 1 -alpha )* np.square(self._metric_B[0]) + alpha* (self._metric_B[1]))[None, :]
        # print(metric_A_temp.shape, metric_B_temp.shape)
        w12 = metric_A_temp @metric_B_temp
        w12_orig = w12.copy()

        def zero_out_w(w, LR_idx):
            lig_idx = np.ravel(np.asarray(LR_idx[:, 0]))
            lig_idx = list(np.unique(lig_idx[lig_idx != 0]) - 1)
            rec_idx = np.ravel(np.asarray(LR_idx[:, 1]))
            rec_idx = list(np.unique(rec_idx[rec_idx != 0]) - 1)
            # reverse select and zeros LR that not in idx list
            mask_lig = np.ones(w.shape[0], dtype=np.bool)
            mask_lig[lig_idx] = 0
            mask_rec = np.ones(w.shape[1], dtype=np.bool)
            mask_rec[rec_idx] = 0

            w[mask_lig, :] = 0
            w[:, mask_rec] = 0
            assert np.count_nonzero(w) == len(lig_idx) * len(rec_idx)

            return w

        if queryDB not in [None, 'comb', 'pairs']:
            raise KeyError('queryDB using the keyword None, \'comb\' or \'pairs\'')
        elif queryDB is None:
            pass
        elif queryDB == 'comb':
            # ada.var index of LR genes (the intersect of DB and object genes, no pair relationship maintained)
            LR_idx_toUse = self._genes_index_DB
            w12 = zero_out_w(w12, LR_idx_toUse)
        elif queryDB == 'pairs':
            # maintain L-R pair relationship, both > 0
            LR_idx_toUse = self._genes_index_DB[(self._genes_index_DB[:, 0] > 0) & (self._genes_index_DB[:, 1] > 0)]
            w12 = zero_out_w(w12, LR_idx_toUse)

        if scale_w:
            mu = 1
            w12 = mu * ((self._net_A +1).sum() + (self._net_B +1).sum()) / \
                        (2 * w12_orig.sum()) * w12  # scale factor using w12_orig

        w = np.block([[self._net_A + 1, w12],
                      [w12.T, self._net_B + 1]])

        return w

    def add_names_to_nets(self):
        '''for graph visualization'''
        self._net_A = pd.DataFrame(self._net_A, columns=self.genes_names[0], index=self.genes_names[0])
        self._net_B = pd.DataFrame(self._net_B, columns=self.genes_names[1], index=self.genes_names[1])
        print('completed.')


def get_candidates(df_filtered):
    '''selected L-R candidates'''
    candidates = [a+'_'+b for a, b in zip(np.asarray(df_filtered['ligand'],dtype=str),
                                          np.asarray(df_filtered['receptor'],dtype=str))]
    return candidates


def get_counts_np(*Xct_objects):
    '''return a list of counts in numpy array, gene by cell'''
    # if not all(isinstance(obj, Xct) for obj in Xct_objects):
    #     raise TypeError('input Xct object(s)')
    # else:
    counts_np = list(itertools.chain(*([obj._X[0].T, obj._X[1].T] for obj in Xct_objects))) #gene by cell
    counts_np = [counts.toarray() if scipy.sparse.issparse(counts) else counts for counts in counts_np]
    return counts_np # a list


def build_W(*Xct_objects):
    '''build a cross-object corresponding matrix for further differential analysis'''
    W12 = np.zeros((Xct_objects[0]._w.shape[0], Xct_objects[1]._w.shape[1]), float)

    mu = 0.9
    scaled_diag = mu * ((Xct_objects[0]._w).sum() + (Xct_objects[1]._w).sum()) / (4 * len(W12))
    np.fill_diagonal(W12, scaled_diag)

    W = np.block([[Xct_objects[0]._w, W12],
                  [W12.T, Xct_objects[1]._w]])

    return W


def _pair_distance(Xct_object, projections, dist_metric='euclidean'):
    '''distances of each directional pair in latent space'''
    X = projections[:len(projections) // 2, :]
    Y = projections[len(projections) // 2:, :]
    dist = scipy.spatial.distance.cdist(X, Y, metric=dist_metric)
    dist_df = pd.DataFrame(dist, index=Xct_object.genes_names[0], columns=Xct_object.genes_names[1])

    return dist_df


def nn_aligned_dist(Xct_object, projections, dist_metric='euclidean', rank=False):
    '''output info of each pair'''
    print(f'computing pair-wise {dist_metric} distances...')
    dist_df = _pair_distance(Xct_object, projections, dist_metric=dist_metric)

    dist_df = pd.DataFrame(dist_df.stack())  # multi_index
    dist_df.reset_index(level=[0, 1], inplace=True)
    dist_df.columns = ['ligand', 'receptor', 'dist']
    dist_df.index = dist_df['ligand'] + '_' + dist_df['receptor']

    print('adding column \'correspondence\'...')
    w12 = Xct_object._w[:Xct_object._net_A.shape[0], Xct_object._net_A.shape[1]:]
    dist_df['correspondence'] = w12.reshape((w12.size,))

    if rank:
        print('adding column \'rank\'...')
        dist_df = dist_df.sort_values(by=['dist'])
        dist_df['rank'] = np.arange(len(dist_df)) + 1

    return dist_df