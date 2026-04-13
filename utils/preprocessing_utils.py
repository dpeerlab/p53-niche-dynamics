import scipy
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
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
from general_utils import ismember, grep

def norm_log_transform(adata_local, size_factors, cell_filter):
    if(len(cell_filter) > 0):
        adata_local = adata_local[cell_filter,:].copy()
        size_factors = np.array(size_factors)[cell_filter]
    else:
        adata_local = adata_local.copy()
    print("Normalizing data")
    if(scipy.sparse.issparse(adata_local.X)):
        c = scipy.sparse.diags(1/size_factors)
        adata_local.X = c * adata_local.X * np.median(size_factors)
    else:
        adata_local.X = adata_local.X / size_factors[:,None] * np.median(size_factors)
    
    print("Log transform with pseudocount 1")
    sc.pp.log1p(adata_local)
    return(adata_local)

def HVG_selection(adata_local, num_HVG = None, exclude_genes = [], defined_HVGs = [], mode='lowess', layer='cellbender', n_bin=40):
    if(len(exclude_genes) != 0):
        print("Highly variable gene selection with custom gene exclusion")
        adata_temp = adata_local[:,np.logical_not(exclude_genes)].copy()
    else:
        print("Highly variable gene selection")
        adata_temp = adata_local.copy()
    if(mode=='defined'):
        print("Using predined list of " + str(len(defined_HVGs)) + " highly variable genes")
        adata_temp.var['highly_variable'] = ismember(adata_temp.var_names, defined_HVGs)[1]
    elif(mode=='scanpy'):
        if(num_HVG != None):
            print("Selecting " + str(num_HVG) + " highly variable genes using Scanpy Seurat v3")
            adata_temp.X = adata_temp.layers[layer].copy()
            sc.pp.highly_variable_genes(adata_temp, n_top_genes=num_HVG, flavor='seurat_v3')
        else:
            print("Selecting " + str(num_HVG) + " highly variable genes using Scanpy Seurat v1")
            sc.pp.highly_variable_genes(adata_temp, flavor='seurat')
    else:
        print("Selecting highly variable genes using lowess regression")
        if(num_HVG == None):
            num_HVG = 100
        adata_temp = lowess_HVG(adata_temp, plot_show=True, n_bin=n_bin, min_mean=0.001, n_genes_per_bin = num_HVG, exclude_genes=adata_temp.var['excluded'])
    adata_local.var['highly_variable'] = ismember(adata_local.var_names, adata_temp.var_names[adata_temp.var['highly_variable']])[1]
    actual_HVG = np.sum(adata_local.var['highly_variable'])
    print("Using " + str(actual_HVG) + " highly variable genes")
    return(adata_local)

def kneepoint(vals, xcoords=None):

    # Convert array of values to 2D vectors by position in array
    vals = np.array(vals)
    descending = vals[-1] < vals[0]
    if descending:
        vals = -vals  # Should go from 0 to 1
    n_pts = vals.shape[0]
    if xcoords is None:
        xcoords = range(n_pts)
    pts = np.vstack([xcoords, vals]).T

    last_vec = pts[-1] - pts[0]  # Get vec b/n first and last points
    last_vec_norm = last_vec / np.linalg.norm(last_vec)  # L2 normalized vector
    all_vecs = pts - pts[0]  # Move to origin at Point 1

    scalars = np.dot(all_vecs, last_vec_norm)  # Scalar of projection to last vector
    projections = np.outer(
        scalars, last_vec_norm
    )  # Projection of original vectors onto last vector

    # Get point w/ max difference b/n projection to diagonal and the original vector
    # This occurs at the knee point (similar to ROC graphs)
    vecs = all_vecs - projections
    dists = np.linalg.norm(vecs, axis=1)
    idx = np.argmax(dists)

    return idx

def process_pca_var_explained(cms, var_to_explain):
    try:
        pc_threshold = np.min(np.where(cms >= var_to_explain)[0])
    except:
        print("Not enough PCs to explain %f of variance. Default to 50 pcs." % (var_to_explain))
        pc_threshold = 50
    return(pc_threshold)
    
def process_pca_inflection(cms, epsilon = 0.00001):
    average_window = 10
    print("Detecting inflection point with average smoothing window of %d and epsilon %f" % (average_window, epsilon))
    d1 = np.diff(pd.Series(cms).rolling(average_window).mean()[average_window:])
    d2 = np.diff(pd.Series(d1).rolling(average_window).mean()[average_window:])
    try:
        inflection_pt = np.min(np.where(np.abs(d2) < epsilon))
    except:
        print("Not enough PCs to reach inflection point with %f epsilon. Default to 50 percent explained variance." % (epsilon))
        inflection_pt = process_pca_var_explained(cms, 0.5)
    return(inflection_pt)

def process_pca_kneepoint(cms):
    pc_kneepoint = kneepoint(cms)
    return(pc_kneepoint)

def plot_pca_results(cms, selected_pcs, title):
    f, ax = plt.subplots(1,1,figsize=[3,3])
    
    plt.plot(range(len(cms)), list(cms), color='black')
    
    plt.plot([selected_pcs, selected_pcs], [0, cms[selected_pcs]], 'o--', color='red')
    plt.plot([0, selected_pcs], [cms[selected_pcs], cms[selected_pcs]], 'o--', color='red')
    plt.text(selected_pcs + 20, 0, str(selected_pcs), color='red')
    plt.text(selected_pcs + 20, cms[selected_pcs], "%.2f" % (cms[selected_pcs]), color='red')
    plt.xlabel("PC")
    plt.ylabel("Cumulative variance explained")
    plt.title(title)

def process_pca_results_from_cumulative_variance(cms, pca_mode = 'explainedvar50', min_pcs = 50, epsilon = 0.00001):
    if pca_mode.startswith('explainedvar'):
        var_to_explain = float(re.sub("explainedvar", "", pca_mode)) / 100
        selected_pcs = process_pca_var_explained(cms, var_to_explain)
    if pca_mode.startswith('fixed'):
        selected_pcs = int(re.sub("fixed", "", pca_mode))
    if pca_mode == 'inflection':
        selected_pcs = process_pca_inflection(cms, epsilon = epsilon)
    if pca_mode == 'kneepoint':
        selected_pcs = process_pca_kneepoint(cms)
    if selected_pcs < min_pcs:
        print("Found %d PCs that explain %f variance" % (selected_pcs, np.sum(cms[selected_pcs-1])))
        print("Fewer than minimum number of PCs. Default to %d that explain %f variance" % (min_pcs, np.sum(cms[min_pcs-1])))
        selected_pcs = min_pcs
    var_explained = np.sum(cms[selected_pcs-1])
    plot_pca_results(cms, selected_pcs, pca_mode)
    return(selected_pcs, var_explained)

def process_pca_results(adata_sub, pca_mode = 'explainedvar50', min_pcs = 50, epsilon = 0.00001):
    if 'pca' in adata_sub.uns:
        pca_explained_variance = adata_sub.uns['pca']['variance_ratio']
    if 'pca_explained_var' in adata_sub.uns:
        pca_explained_variance = adata_sub.uns['pca_explained_var']
    cms = np.cumsum(pca_explained_variance)
    selected_pcs, var_explained = process_pca_results_from_cumulative_variance(cms, pca_mode, min_pcs, epsilon = epsilon)
    return(selected_pcs, var_explained)

def consolidate_var_names(adata, return_sparse = False):
    df = pd.DataFrame(adata.X.todense(), index = adata.obs_names, columns = adata.var_names)
    singleton = df.columns.value_counts() == 1
    df_single = df.loc[:,singleton]
    df_multiple = df.loc[:,np.logical_not(singleton)]
    df_multiple = df_multiple.T.groupby(df_multiple.columns).agg('sum').T
    np.sum(df_multiple.index == df_single.index) == df_single.shape[0]
    df = pd.concat((df_single, df_multiple), axis=1)
    obs = adata.obs.copy()
    adata = sc.AnnData(df)
    if(return_sparse):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.obs = obs.copy()
    return(adata)

def add_obs_names(adata, obs_names):
    first_check = ismember(adata.obs_names, obs_names)[1]
    if(np.sum(first_check) / adata.n_obs != 1):
        print("Some obs_names in adata are absent in master obs_names list. Aborting!")
        return

    present_filter = ismember(obs_names, adata.obs_names)[1]    
    absent_filter = np.logical_not(present_filter)
    cellbender_0 = adata.obs_names[np.ravel(adata.X.sum(axis=1) == 0)]
    
    print(str(np.sum(present_filter)) + '\t' + str(np.sum(absent_filter)))
    
    if(np.sum(absent_filter) > 0):
        added_rows = np.zeros([sum(absent_filter), adata.n_vars])
        finalList = adata.obs_names.tolist() + list(obs_names[absent_filter])
        newMatrix = pd.DataFrame(np.concatenate((adata.X.todense(), added_rows), axis=0), index = finalList, columns=adata.var_names)
    else:
        finalList = obs_names[present_filter]
        newMatrix = pd.DataFrame(adata.X.copy(), index = finalList, columns=adata.var_names)
    
    newAdata = sc.AnnData(scipy.sparse.csr_matrix(newMatrix.to_numpy()))
    newAdata.obs_names = finalList
    newAdata.var_names = adata.var_names
    
    newAdata = newAdata[obs_names,:]
    
    valid_obs_names = obs_names[np.logical_not(ismember(obs_names, cellbender_0)[1])]
    oldIndexes = np.sort(ismember(adata.obs_names[ismember(adata.obs_names, valid_obs_names)[1]], obs_names)[0])
    rowSums = newAdata.X.sum(axis=1)
    newIndexes = np.where(rowSums > 0)[0]
    newIndexes2 = np.sort(ismember(obs_names[np.logical_and(present_filter, np.logical_not(ismember(obs_names, cellbender_0)[1]))], obs_names)[0])

    display(np.sum([x == y for x,y in zip(oldIndexes, newIndexes)])/len(oldIndexes))
    display(np.sum([x == y for x,y in zip(oldIndexes, newIndexes2)])/len(oldIndexes))
    
    return(newAdata)
