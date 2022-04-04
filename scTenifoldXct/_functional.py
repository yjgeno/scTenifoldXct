import numpy as np
import pandas as pd


def filtered_nn_aligned_dist(df_nn, candidates):
    '''filter and rank only L-R pairs'''
    if 'diff2' in df_nn.columns:
        df_nn_filtered = df_nn.loc[candidates].sort_values(by=['diff2'])  # dist difference^2 ranked L-R candidates
    else:
        df_nn_filtered = df_nn.loc[candidates].sort_values(by=['dist'])  # dist ranked L-R candidates
    print('manifold aligned # of L-R pairs:', len(df_nn_filtered))
    df_nn_filtered['rank_filtered'] = np.arange(len(df_nn_filtered)) + 1

    return df_nn_filtered


def KnK(Xct_object, ko_genelist, copy=True):
    '''set rows and cols of ko_genelist to zeros in correspondence w'''
    from copy import deepcopy
    # ko_gene index in cell A and B
    ko_idx1 = [Xct_object.genes_names[0].index(g) for g in ko_genelist]
    ko_idx2 = [Xct_object.genes_names[1].index(g) + len(Xct_object.genes_names[0]) for g in ko_genelist]

    if copy:
        Xct_object = deepcopy(Xct_object)
    for idx in [ko_idx1, ko_idx2]:
        Xct_object._w[idx, :] = 0
        Xct_object._w[:, idx] = 0

    return Xct_object


def nn_aligned_dist_diff(df_nn1, df_nn2, rank=False):
    '''pair-wise difference of aligned distance'''
    print(f'computing pair-wise distance differences...')
    df_nn_all = pd.concat([df_nn1, df_nn2.drop(['ligand', 'receptor'], axis=1)], axis=1)
    print('adding column \'diff2\'...')
    df_nn_all['diff2'] = np.square(
        df_nn_all['dist'].iloc[:, 0] - df_nn_all['dist'].iloc[:, 1])  # there are two 'dist' cols
    if rank:
        print('adding column \'diff2_rank\'...')
        df_nn_all = df_nn_all.sort_values(by=['diff2'], ascending=False)
        df_nn_all['diff2_rank'] = np.arange(len(df_nn_all)) + 1

    return df_nn_all


def get_genelist(df_enriched, saveas=None):
    '''get a list of single genes from enriched pairs'''
    targets = np.ravel([n.split('_') for n in df_enriched.index])  # .tolist()
    targets = list(set(targets))
    if saveas is not None:
        with open(f'{saveas}.txt', 'w') as file:
            for g in targets:
                file.write(g + '\n')
    return targets