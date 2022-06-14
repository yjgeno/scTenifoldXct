import numpy as np
import pandas as pd
import scipy
from scipy import sparse

from .core import scTenifoldXct
from .nn import ManifoldAlignmentNet
from .stat import chi2_diff_test
from .visualization import plot_pcNet_method


class merge_scTenifoldXct:
    def __init__(self, 
                *Xcts: scTenifoldXct, 
                verbose: bool = True):
        self.Xcts = Xcts
        self._merge_candidates = list(set(self.Xcts[0]._candidates).union(set(self.Xcts[1]._candidates)))
        self.verbose = verbose
        self.n_dim = 3
        self.mu = 0.9
        # cal big W
        if self.verbose:
            print("merging samples and building correspondence...")
        self._W, self.W12_shape = self._build_W()

        self._nn_trainer = ManifoldAlignmentNet(self._get_data_arrs(),
                                                w=self._W,
                                                n_dim=self.n_dim,
                                                layers=None)
        if self.verbose:
            print("merge_scTenifoldXct init completed")

    def _get_data_arrs(self):  
        '''return a list of counts in numpy array, gene by cell'''
        data_arr_A = [cell_data.X.T.toarray() if scipy.sparse.issparse(cell_data.X) else cell_data.X.T   # gene by cell
                    for _, cell_data in self.Xcts[0]._cell_data_dic.items()]
        data_arr_B = [cell_data.X.T.toarray() if scipy.sparse.issparse(cell_data.X) else cell_data.X.T   # gene by cell
                    for _, cell_data in self.Xcts[1]._cell_data_dic.items()]
        return data_arr_A + data_arr_B  # a list

    def _build_W(self):
        '''build a cross-object corresponding matrix for further differential analysis'''
        W12 = np.zeros((self.Xcts[0]._w.shape[0], self.Xcts[1]._w.shape[1]), float)
        scaled_diag = self.mu * ((self.Xcts[0]._w).sum() + (self.Xcts[1]._w).sum()) / (4 * len(W12)) 
        np.fill_diagonal(W12, scaled_diag)
        # W12 = W12.todok()
        W = sparse.vstack([sparse.hstack([self.Xcts[0]._w, W12]),
            sparse.hstack([W12.T, self.Xcts[1]._w])])      
        return W, W12.shape

    @property
    def trainer(self):
        return self._nn_trainer

    @property
    def W(self):
        return self._W

    @property
    def merge_candidates(self):
        return self._merge_candidates

    def get_embeds(self,
                train = True,
                n_steps=1000,
                lr=0.001,
                verbose=False,
                plot_losses: bool = False,
                losses_file_name: str = None,
                **optim_kwargs
                ):
        if train:
            merge_projections = self._nn_trainer.train(n_steps=n_steps, lr=lr, verbose=verbose, **optim_kwargs)
        else:
            merge_projections = self._nn_trainer.reload_embeds()
        if plot_losses:
            self._nn_trainer.plot_losses(losses_file_name)
       
        return merge_projections

    def plot_losses(self, **kwargs):
        self._nn_trainer.plot_losses(**kwargs)


    def nn_aligned_diff(self, 
                    merge_projections,
                    dist_metric: str = "euclidean",
                    rank: bool = False
                    ):
        '''pair-wise difference of aligned distance'''
        projections_split = np.array_split(merge_projections, 2)
        self._aligned_result_A = self.Xcts[0]._nn_trainer.nn_aligned_dist(projections_split[0],
                                                                gene_names_x=self.Xcts[0]._genes[self.Xcts[0]._cell_names[0]],
                                                                gene_names_y=self.Xcts[0]._genes[self.Xcts[0]._cell_names[1]],
                                                                w12_shape=self.Xcts[0].w12_shape,
                                                                dist_metric=dist_metric,
                                                                rank=rank,
                                                                verbose=self.verbose)
        self._aligned_result_B = self.Xcts[1]._nn_trainer.nn_aligned_dist(projections_split[1],
                                                                gene_names_x=self.Xcts[1]._genes[self.Xcts[1]._cell_names[0]],
                                                                gene_names_y=self.Xcts[1]._genes[self.Xcts[1]._cell_names[1]],
                                                                w12_shape=self.Xcts[1].w12_shape,
                                                                dist_metric=dist_metric,
                                                                rank=rank,
                                                                verbose=self.verbose)
        df_nn_all = pd.concat([self._aligned_result_A, self._aligned_result_B.drop(['ligand', 'receptor'], axis=1)], axis=1) 
        # print('adding column \'diff2\'...')
        df_nn_all['diff2'] = np.square(df_nn_all['dist'].iloc[:, 0] - df_nn_all['dist'].iloc[:, 1]) #there are two 'dist' cols
        if rank:
            # print('adding column \'diff2_rank\'...')
            df_nn_all = df_nn_all.sort_values(by=['diff2'], ascending=False)
            df_nn_all['diff2_rank'] = np.arange(len(df_nn_all)) + 1
        self._aligned_diff_result = df_nn_all
        if self.verbose:
            print("merged pair-wise distances")
        # return df_nn_all

    def chi2_diff_test(self,
                  dof=1,
                  pval=0.05,
                  cal_FDR=True,
                  plot_result=False,
                  ):
        return chi2_diff_test(df_nn=self._aligned_diff_result, 
                        df=dof,
                        pval=pval,
                        FDR=cal_FDR,
                        candidates=self._merge_candidates,
                        plot=plot_result)
    
    @property
    def aligned_diff_result(self):
        if self._aligned_diff_result is None:
            raise AttributeError("No aligned_diff_result created yet. "
                                 "Please call train_nn() to train the neural network to get embeddings first.")
        return self._aligned_diff_result

    # def plot_merge_pcNet_graph(self, 
    #                         gene_names: list[str], 
    #                         sample: int = 0, 
    #                         view: str ="sender", 
    #                         **kwargs):
    #     if view not in ["sender", "receiver"]:
    #         raise ValueError("view needs to be sender or receiver")

    #     g = plot_pcNet_method(self.Xcts[sample]._net_A if view == "sender" else self.Xcts[sample]._net_B,
    #                       gene_names=gene_names,
    #                       tf_names=self.Xcts[sample]._TFs["TF_symbol"].to_list(),
    #                       **kwargs)
    #     return g

