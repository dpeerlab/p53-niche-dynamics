from scipy.sparse import csr_matrix, find, issparse
from scipy.sparse.linalg import eigs
import matplotlib
import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

def get_diffusion_operator_from_knn(indices, distances):
    print("Computing diffusion kernel using CPU")
    N = indices.shape[0]
    k = np.min(np.ravel(indices.sum(axis=1)))
    adaptive_k = int(np.floor(k / 3))
    adaptive_std = np.zeros(N)
    
    for i in np.arange(len(adaptive_std)):
        if(scipy.sparse.issparse(distances)):
            current_distances = np.ravel(distances[i,:].todense())
            current_distances = current_distances[current_distances > 0]
        else:
            current_distances = distances[i,:]
        adaptive_std[i] = np.sort(current_distances)[adaptive_k - 1]

    # Kernel
    if(scipy.sparse.issparse(indices)):
        sources, targets = indices.nonzero()
        dists = np.ravel(distances[sources, targets])
    else:
        targets = np.ravel(indices)
        sources = np.ravel([[i]*k for i in range(N)])
        dists = np.ravel(distances)

    # X, y specific stds
    dists = dists / adaptive_std[sources]
    W = scipy.sparse.csr_matrix((np.exp(-dists), (sources, targets)), shape=[N, N])
    
    # Diffusion components
    kernel = W + W.T
    
    print("Computing diffusion operator using CPU")
    # Markov
    D = np.ravel(kernel.sum(axis=1))
    D[D != 0] = 1 / D[D != 0]
    D = scipy.sparse.csr_matrix((D, (range(N), range(N))), shape=[N, N])
    T = D.dot(kernel)
    
    return(T, kernel)

def getDiffusionComponents(T, kernel, n_components = 10):
    D, V = eigs(T, n_components, tol=1e-4, maxiter=1000)
    D = np.real(D)
    V = np.real(V)
    inds = np.argsort(D)[::-1]
    D = D[inds]
    V = V[:, inds]
    # Normalize
    for i in range(V.shape[1]):
        V[:, i] = V[:, i] / np.linalg.norm(V[:, i])
    # Create are results dictionary
    res = {"T": T, "EigenVectors": V, "EigenValues": D, 'kernel': kernel}
    res["EigenVectors"] = pd.DataFrame(res["EigenVectors"])
    res["EigenValues"] = pd.Series(res["EigenValues"])
    return(res)

def imputeData(adata, T, n_steps = 1):
    if(issparse(adata.X)):
        my_matrix = adata.X.todense().copy()
    else:
        my_matrix = adata.X.copy()

    data_df = pd.DataFrame(my_matrix, index = adata.obs_names, columns = adata.var_names)
    T_steps = T ** n_steps
    
    imputed_data = pd.DataFrame(
        np.dot(T_steps.todense(), data_df), index=data_df.index, columns=data_df.columns
    )
    return(imputed_data)

def get_imputed_adata_subset(adata, subcells, k, t_steps = 3):
    adata_sub = adata[subcells,:].copy()
    print("Getting knn graph")
    sc.pp.neighbors(adata_sub, n_neighbors=k)
    distances = adata_sub.obsp['distances']
    indices = (distances > 0).astype(np.uint8)
    print("Getting diffusion operator")
    diffusion_operator, kernel = get_diffusion_operator_from_knn(indices, distances)
    print("Imputing data")
    T = diffusion_operator.todense() ** t_steps
    #imputed_data = (diffusion_operator @ adata_sub.X).todense()
    imputed_data = T @ adata_sub.X.todense()
    adata_imputed = sc.AnnData(imputed_data)
    adata_imputed.obs_names = adata_sub.obs_names.copy()
    adata_imputed.var_names = adata_sub.var_names.copy()
    return(adata_imputed)

def plot_diffusion_components(X_embed, DCs, vmax=None):
    """ Plots the diffusion components on obsm embedding
    :return: fig, ax
    """

    # Please run tSNE before plotting diffusion components. #
    # Please run diffusion maps using run_diffusion_map before plotting #

    # Plot
    fig = FigureGrid(DCs.shape[1], 5)

    for i, ax in enumerate(fig):
        ax.scatter(
            X_embed[:,0],
            X_embed[:,1],
            c=DCs[:, i],
            cmap=matplotlib.cm.Spectral_r,
            edgecolors="none",
            s=3,
        )
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.set_aspect("equal")
        ax.set_title("Component %d" % i, fontsize=10)
        ax.set_axis_off()

class FigureGrid:
    """
    Generates a grid of axes for plotting
    axes can be iterated over or selected by number. e.g.:
    >>> # iterate over axes and plot some nonsense
    >>> fig = FigureGrid(4, max_cols=2)
    >>> for i, ax in enumerate(fig):
    >>>     plt.plot(np.arange(10) * i)
    >>> # select axis using indexing
    >>> ax3 = fig[3]
    >>> ax3.set_title("I'm axis 3")
    """

    # Figure Grid is favorable for displaying multiple graphs side by side.

    def __init__(self, n: int, max_cols=3, scale=3):
        """
        :param n: number of axes to generate
        :param max_cols: maximum number of axes in a given row
        """

        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * scale, self.nrows * scale

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]