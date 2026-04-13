import scipy
import numpy as np
import pandas as pd
import cupy
import cuml
import cudf
import cugraph
import os

import importlib.util
import sys

def lazy_import(module_name, path_to_file):
    spec = importlib.util.spec_from_file_location(module_name,path_to_file)
    foo = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = foo
    spec.loader.exec_module(foo)
    return foo

utils_dir = "../../utils"
general_utils =  lazy_import("general_utils",os.path.join(utils_dir, "general_utils.py"))
anndata_utils = lazy_import("anndata_utils",os.path.join(utils_dir, "anndata_utils.py"))
plotting_functions = lazy_import("plotting_functions",os.path.join(utils_dir, "plotting_functions.py"))

from general_utils import ismember, grep, array_to_csr_matrix, array_to_edges

def create_diagonal_matrix(v, total_cells):
    cell_indexes = np.array(list(range(total_cells)))
    diagonal_matrix = scipy.sparse.csr_matrix((v,(cell_indexes, cell_indexes)), shape = [total_cells, total_cells])
    return(diagonal_matrix)

def get_average_neighborhood_expression(expression_matrix, spatial_neighborhood, anchor_cells = None, neighbor_cells = None):
    if anchor_cells is None:
        anchor_cells = np.array([True] * expression_matrix.shape[0])
    if neighbor_cells is None:
        neighbor_cells = np.array([True] * expression_matrix.shape[0])
    X_neighbors = expression_matrix[neighbor_cells,:].copy()
    spatial_neighborhood = spatial_neighborhood[anchor_cells,:].copy()[:,neighbor_cells].copy()
    num_neighbors_per_anchor = np.ravel(spatial_neighborhood.sum(axis=1))
    n_anchors = spatial_neighborhood.shape[0]
    normalization_matrix = create_diagonal_matrix(1/num_neighbors_per_anchor, n_anchors)
    average_expression_matrix = normalization_matrix @ spatial_neighborhood @ X_neighbors
    return(average_expression_matrix)

def get_total_neighborhood_expression(expression_matrix, spatial_neighborhood, anchor_cells = None, neighbor_cells = None):
    if anchor_cells is None:
        anchor_cells = np.array([True] * expression_matrix.shape[0])
    if neighbor_cells is None:
        neighbor_cells = np.array([True] * expression_matrix.shape[0])
    X_neighbors = expression_matrix[neighbor_cells,:].copy()
    spatial_neighborhood = spatial_neighborhood[anchor_cells,:].copy()[:,neighbor_cells].copy()
    total_expression_matrix = spatial_neighborhood @ X_neighbors
    return(total_expression_matrix)    

def get_neighborhood_fraction(spatial_neighborhood, boolean_variable):
    n_cells = spatial_neighborhood.shape[0]
    cell_indexes = np.array(list(range(n_cells)))
    masked_matrix = mask_neighborhood_columns(spatial_neighborhood, boolean_variable)
    boolean_fraction = np.ravel(masked_matrix.sum(axis=1)) / np.ravel(spatial_neighborhood.sum(axis=1))
    return(boolean_fraction)

def get_neighborhood_count(spatial_neighborhood, boolean_variable):
    n_cells = spatial_neighborhood.shape[0]
    cell_indexes = np.array(list(range(n_cells)))
    masked_matrix = mask_neighborhood_columns(spatial_neighborhood, boolean_variable)
    boolean_count = np.ravel(masked_matrix.sum(axis=1)) 
    total_count = np.ravel(spatial_neighborhood.sum(axis=1))
    return(boolean_count, total_count)

def mask_neighborhood_columns(spatial_neighborhood, mask):
    mask_matrix = create_diagonal_matrix(mask, spatial_neighborhood.shape[0])
    masked_neighborhood = spatial_neighborhood @ mask_matrix
    return(masked_neighborhood)

def mask_neighborhood_rows(spatial_neighborhood, mask):
    mask_matrix = create_diagonal_matrix(mask, spatial_neighborhood.shape[0])
    masked_neighborhood = mask_matrix @ spatial_neighborhood
    return(masked_neighborhood)

def row_normalize_sparse_matrix(mat):
    normalization_vector = np.ravel(mat.sum(axis=1))
    normalization_factor = create_diagonal_matrix(1/normalization_vector, mat.shape[0])
    normalized_matrix = normalization_factor @ mat * np.median(normalization_vector)
    return(normalized_matrix)

def neighbors_gpu(mat, k=30):
    X_mat = cupy.array(mat)
    model = cuml.neighbors.NearestNeighbors(n_neighbors=k)
    model.fit(X_mat)
    distances, indices = model.kneighbors(X_mat)
    return(distances, indices)

def get_cugraph_from_csr_matrix(mat):
    sources, targets = mat.nonzero()
    edges = cudf.concat([
        cudf.Series(sources, name='source'), # sources
        cudf.Series(targets, name='destination'), # targets
    ], axis=1)
    G = cugraph.from_cudf_edgelist(edges)
    return(G)

def get_cugraph_from_array(mat):
    targets = np.ravel(mat)
    n_neighbors_per_cell = [len(x) for x in mat]
    sources = np.ravel([[i]*k for i,k in zip(range(mat.shape[0]), n_neighbors_per_cell)])
    edges = cudf.concat([
        cudf.Series(sources, name='source'), # sources
        cudf.Series(targets, name='destination'), # targets
    ], axis=1)
    G = cugraph.from_cudf_edgelist(edges)
    return(G)

def get_spatial_knn(adata, k = 30):
    # Compute spatial knn
    spatial_knn = neighbors_gpu(adata.obsm['X_spatial'])[1]
    targets = np.ravel(spatial_knn).get()
    n_neighbors_per_cell = [len(x) for x in spatial_knn]
    sources = np.ravel([[i]*k for i,k in zip(range(spatial_knn.shape[0]), n_neighbors_per_cell)])
    df = pd.DataFrame({'sources': sources, 'targets': targets})
    return(df)

def compute_sample_neighborhoods(adata, slide_key, cell_type_key, selected_cell_types = [], k = 30):
    
    # Pre-compute cell_type_filter
    if(len(selected_cell_types) == 0):
        selected_cell_types = np.unique(adata.obs[cell_type_key])
    cell_type_filter = ismember(adata.obs[cell_type_key], selected_cell_types)[1]

    # Compute spatial neighborhoods per sample
    sources = np.array([], dtype=np.uint32)
    targets = np.array([], dtype=np.uint32)

    unique_slides = np.unique(adata.obs[slide_key])
    for my_slide in tqdm(unique_slides):
        slide_filter = adata.obs[slide_key] == my_slide
        selected_cell_indexes_slide = np.where(slide_filter)[0]
        adata_sub = adata[slide_filter,:].copy()
        
        selected_cells = np.logical_and(slide_filter, cell_type_filter)
        selected_cell_indexes_query = np.where(selected_cells)[0]
        
        result = get_spatial_knn(adata_sub, k=k)        
        sources = np.hstack((sources, selected_cell_indexes_slide[result['sources'].to_numpy().astype(np.uint32)]))
        targets = np.hstack((targets, selected_cell_indexes_slide[result['targets'].to_numpy().astype(np.uint32)]))
        
    spatial_neighborhood = scipy.sparse.csr_matrix(([1] * len(sources), (sources, targets)), shape = [adata.n_obs] * 2).astype(np.uint8)
    return(spatial_neighborhood)

def compute_distances_within_neighborhood(adata, obsm_field, neighbor_graph):
    """Computes euclidian distances within neighborhoods using arbitrary feature matrix

    Args:
        - adata: base adata object
        - obsm_field: feature matrix within the adata object (e.g. "X_spatial")
        - neighbor_graph: scipy.sparse.csr_matrix containing binary neighbor assignments. The number of rows is equal to adata.n_obs

    Returns
    -------
        - scipy.sparse.csr_matrix with distances between source and query cells, as specified in neighbor_graph.
    """
    assert (
        adata.n_obs == neighbor_graph.shape[0]
    ), "The number of cells in adata is different from number of cells in neighbor_graph"
    sources, targets = neighbor_graph.nonzero()
    source_vectors = adata.obsm[obsm_field][sources, :]
    target_vectors = adata.obsm[obsm_field][targets, :]
    dist = np.sqrt(np.power(source_vectors - target_vectors, 2).sum(axis=1))
    distance_neighborhood = scipy.sparse.csr_matrix((dist, (sources, targets)), shape=[adata.n_obs] * 2)
    return distance_neighborhood

