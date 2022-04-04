import numpy as np
from numpy.linalg import svd
import scipy


def pcNet(X, nComp = 3, scale = True, symmetric = False, q = 0): # X: cell * gene, q: 0-100
    X = X.toarray() if scipy.sparse.issparse(X) else X
    if not isinstance(X, np.ndarray):
        raise TypeError('input a numpy array with cells as rows and genes as columns')
    elif nComp < 2 or nComp >= X.shape[1]:
        raise ValueError('nComp should be greater or equal than 2 and lower than the total number of genes')
    
    else:
        n = X.shape[1] # genes
        def pcCoefficients(K):
            y = X[:, K] 
            Xi = np.delete(X, K, 1)
            U, s, VT = svd(Xi, full_matrices=False) 
            #print ('U:', U.shape, 's:', s.shape, 'VT:', VT.shape)
            V = VT[:nComp, :].T
            #print('V:', V.shape)

            score = Xi@V
            t = np.sqrt(np.sum(score**2, axis=0))
            score_lsq = ((score.T / (t**2)[:, None])).T
            beta = np.sum(y[:, None]*score_lsq, axis=0)
            beta = V@beta

            return list(beta)
        
        B = []
        for k in range(n):
            B.append(pcCoefficients(k))
        B = np.array(B)
        
        A = np.ones((n, n), dtype=float)
        np.fill_diagonal(A, 0)
        for i in range(n):
            A[i, A[i, :]==1] = B[i, :]
                 
        if scale:
            absA = abs(A)
            A = A / np.max(absA)
        if q > 0:
            A[absA < np.percentile(absA, q)] = 0
        if symmetric: # place in the end
            A = (A + A.T)/2
        #diag(A) <- 0
        
        return A
      
      
