import numpy as np
import pandas as pd
import scipy
import scanpy as sc

def marker_genes_for_cell_group(adata, query_cells, reference_cells = None, num_genes = 20):
    if reference_cells is None:
        reference_cells = np.array([True] * adata.n_obs)
    adata.obs['group_temp'] = query_cells
    valid_cells = np.logical_or(query_cells, reference_cells)
    adata_sub = adata[valid_cells,:].copy()
    adata_sub.obs['group_temp'] = ['query' if x else 'reference' for x in adata_sub.obs['group_temp']]
    sc.tl.rank_genes_groups(adata_sub, groupby='group_temp', n_genes = num_genes)
    sc.tl.dendrogram(adata_sub, groupby='group_temp')
    sc.pl.rank_genes_groups_dotplot(adata_sub, groupby='group_temp', n_genes=50)

# This function removes genes with 0 variance
def zscore_with_subset(adata, selected_cells):
    mean_val = np.ravel(adata[selected_cells,:].X.mean(axis=0))
    std_val = np.ravel(np.std(adata[selected_cells,:].X.todense(), axis=0))
    adata = adata[:,std_val > 0].copy()
    mean_val = np.ravel(adata[selected_cells,:].X.mean(axis=0))
    std_val = np.ravel(np.std(adata[selected_cells,:].X.todense(), axis=0))
    X = adata.X.todense()
    zscore = np.divide(np.subtract(X, mean_val), std_val)
    adata.layers['zscore_subcells'] = zscore.copy()
    return(adata)
    
def zscore_with_subset(adata, selected_cells, selected_layer = 'X'):
    if selected_layer == 'X':
        mat = adata.X
    else:
        mat = adata.layers[selected_layer]
    if scipy.sparse.issparse(mat):
        mat = mat.todense()
    mean_val = np.ravel(mat[selected_cells,:].mean(axis=0))
    std_val = np.ravel(np.std(mat[selected_cells,:], axis=0))
    mean_val = np.ravel(mat[selected_cells,:].mean(axis=0))
    std_val = np.ravel(np.std(mat[selected_cells,:], axis=0))
    zscore = np.divide(np.subtract(mat, mean_val), std_val)
    zscore[np.isinf(zscore)] = np.nan
    adata.layers['zscore_subcells'] = zscore.copy()
    return(adata)

def obsm_to_df(adata, obsm):
    df = pd.DataFrame(adata_sub.obsm[obsm])
    df.index = adata.obs_names
    return(df)
    
def get_gene(adata, gene_name, layer = None):
    if(layer == None):
        selected_matrix = adata.X
    else:
        selected_matrix = adata.layers[layer]
    if scipy.sparse.issparse(selected_matrix):
        v = np.array(np.ravel(selected_matrix[:,adata.var_names == gene_name].todense()))
    else:
        v = np.array(np.ravel(selected_matrix[:,adata.var_names == gene_name]))
    return(v)


def get_genes(adata, gene_names, layer = None):
    if(layer == None):
        selected_matrix = adata[:,gene_names].X
    else:
        selected_matrix = adata[:,gene_names].layers[layer]
    if scipy.sparse.issparse(selected_matrix):
        v = np.ravel(selected_matrix.todense())
    else:
        v = np.ravel(selected_matrix)
    return(v)

def get_features(adata, features):
    features = np.array(features)
    result = []
    for current_feature in features:
        if np.sum(ismember([current_feature], adata.var_names)[1]) > 0:
            result.append(list(get_gene(adata, current_feature)))
        else:
            result.append(adata.obs[current_feature].to_list())
    return(np.array(result).T)

def annotate_groups(adata_sub, grouping):
    adata_sub.obs['group_id'] = adata_sub.obs[grouping].astype(str).apply(lambda x: '__'.join(x), axis=1)
    group_counts = adata_sub.obs['group_id'].value_counts()
    group_dictionary = {x: x + "_(n = " + str(group_counts.loc[x]) + ")" for x in np.unique(adata_sub.obs['group_id'])}
    adata_sub.obs['group_id'] = [group_dictionary[x] for x in adata_sub.obs['group_id']]

def filter_coordinates(coord, xlim, ylim):
    if isinstance(coord, pd.DataFrame):
        coord = coord.to_numpy()
    valid_points = (coord[:,0] > xlim[0]) & (coord[:,0] < xlim[1]) & (coord[:,1] > ylim[0]) & (coord[:,1] < ylim[1])
    return valid_points

def query_region(adata, coord, xlim, ylim):
    valid_points = filter_coordinates(adata.obsm[coord], xlim, ylim)
    adata_sub = adata[valid_points,:].copy()
    return(adata_sub)

def random_adata_subset(adata, scaling = 0.1):
    import random
    cell_indexes = list(range(adata.n_obs))
    random_indexes = random.sample(cell_indexes, int(adata.n_obs * scaling))
    adata_sub = adata[random_indexes,:].copy()
    return(adata_sub)
