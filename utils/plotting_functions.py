import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import matplotlib
import scipy
from matplotlib.colors import Colormap
from pkg_resources import parse_version
import os
from shapely import Point, Polygon
from geopandas import GeoDataFrame
import shapely.geometry

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
anndata_utils =  lazy_import("anndata_utils",os.path.join(utils_dir, "anndata_utils.py"))
from general_utils import ismember, grep
from anndata_utils import query_region

cmap0 = matplotlib.colors.LinearSegmentedColormap.from_list(
        'blue2red', ['navy', 'red'])
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list(
        'unevently divided', [(0, 'b'), (.3, 'gray'), (1, 'green')])
cmap0 = matplotlib.colors.LinearSegmentedColormap.from_list(
    'blue2orange', ['royalblue', 'orange'])

def plot_individual_category(
    adata, 
    coord_name, 
    obs_name, 
    selected_values, 
    condition_colors, 
    alpha = 1, 
    s=0.1, 
    frameon=True, 
    figsize_max = 10, 
    ax = None,
    rasterized = True, 
    background_color = 'lightgray', 
    use_outlines = True, 
    axes_color = 'white',
    scale_1 = 5,
    scale_2 = 10,
    scale_3 = 20
):
    x = adata.obsm[coord_name][:,0]
    y = adata.obsm[coord_name][:,1]
    max_x = np.max(x)
    max_y = np.max(y)
    if(max_x > max_y):
        figsize = [figsize_max, figsize_max*max_y/max_x]
    else:
        figsize = [figsize_max*max_x/max_y, figsize_max]
    figsize = [figsize_max, figsize_max]
    conditions = adata.obs[obs_name]
    colors = np.array([condition_colors[x] if x != None else 'lightgray' for x in ismember(conditions, selected_values)[0]])    
    subcells = ismember(conditions, selected_values)[1]
    point_order = np.random.permutation(np.sum(subcells))
    if ax is None:
        f, ax = plt.subplots(1,1,figsize=figsize, dpi=150)
        if not frameon:
            ax.set_xticks([])
            ax.set_yticks([])
        ax.set_facecolor(axes_color)
        if(use_outlines):
            ax.scatter(x,y, facecolor='black', edgecolor='none', marker = 'o', linewidth=0, s=s*scale_3, alpha = 1, rasterized = rasterized)
            ax.scatter(x,y, facecolor='white', edgecolor='none', marker = 'o', linewidth=0, s=s*scale_2, alpha = 1, rasterized = rasterized)
        ax.scatter(x,y, facecolor=background_color, edgecolor='none', marker = 'o', linewidth=0, s=s*scale_1, alpha = 1, rasterized = rasterized)
    ax.scatter(x[subcells][point_order],y[subcells][point_order], facecolor=colors[subcells][point_order], edgecolor=None, linewidth=0, s=s*scale_1, alpha = alpha, rasterized = rasterized)
    ax.grid(False)
    return(ax)

def plot_scatter(adata, obsm_name, colormap, exclude_cells = [], grayout_cells = [], figsize=[2,2], s=1, grayout_color='lightgray', cmap='viridis', vmin=None, vmax=None, use_outlines = True):
    point_order = random.sample(list(range(adata.n_obs)), k=adata.n_obs)
    xvalue = np.ravel(adata.obsm[obsm_name][point_order,0])
    yvalue = np.ravel(adata.obsm[obsm_name][point_order,1])
    if len(exclude_cells) == 0:
        exclude_cells = [False] * adata.n_obs
    if len(grayout_cells) == 0:
        grayout_cells = [False] * adata.n_obs
    exclude_cells = np.array(exclude_cells)[point_order]
    grayout_cells = np.array(grayout_cells)[point_order]
    colormap = colormap[point_order]
    if vmin==None:
        vmin = np.quantile(colormap, 0.00)
    if vmax==None:
        vmax = np.quantile(colormap, 0.99)
    f, ax = plt.subplots(1,1,figsize=figsize)
    if(use_outlines):
        plt.scatter(xvalue, yvalue, color='black', s=s*5, rasterized=True)
        plt.scatter(xvalue, yvalue, color='white', s=s*2, rasterized=True)
    subcells = np.logical_not(exclude_cells)
    plt.scatter(xvalue[subcells], yvalue[subcells], color=grayout_color, s=s, rasterized=True)
    subcells = np.logical_not(np.logical_or(exclude_cells, grayout_cells))
    if(use_outlines):
        plt.scatter(xvalue[subcells], yvalue[subcells], color='black', s=s*5, rasterized=True)
        plt.scatter(xvalue[subcells], yvalue[subcells], color='white', s=s*2, rasterized=True)
    im = plt.scatter(xvalue[subcells], yvalue[subcells], c = colormap[subcells], s=s, vmin=vmin, vmax=vmax, rasterized=True, cmap=cmap)
    plt.grid(False)
    plt.axis(False)
    cb_ax = f.add_axes([0.45,0.85,0.25,0.02])
    f.colorbar(im,orientation='horizontal',cax=cb_ax, fraction=0.01, label='Average z-score')
    return(ax, cb_ax)

def plot_signature(adata, obsm_name, signature_name, exclude_cells = [], grayout_cells = [], figsize=[2,2], s=1, grayout_color='lightgray', cmap='viridis', vmin=None, vmax=None, use_outlines = True):
    colormap = adata.obs[signature_name].to_numpy()
    ax, cb_ax = plot_scatter(adata, obsm_name, colormap, exclude_cells = exclude_cells, grayout_cells = grayout_cells, figsize = figsize, s=s, grayout_color=grayout_color, cmap=cmap, vmin=vmin, vmax=vmax, use_outlines = use_outlines)
    return(ax, cb_ax)

def plot_gene(adata, obsm_name, gene_name, exclude_cells = [], grayout_cells = [], figsize=[2,2], s=1, grayout_color='lightgray', cmap='viridis', vmin=None, vmax=None, use_outlines = True):
    colormap = np.ravel(adata.X[:,adata.var_names == gene_name].todense())
    ax, cb_ax = plot_scatter(adata, obsm_name, colormap, exclude_cells = exclude_cells, grayout_cells = grayout_cells, figsize = figsize, s=s, grayout_color=grayout_color, cmap=cmap, vmin=vmin, vmax=vmax, use_outlines = use_outlines)
    return(ax, cb_ax)

def get_pseudo_colors(adata, feature_names, rgb_cols, quantile_threshold = [0, 0.95]):
    colors_all = get_features(adata, feature_names)
    clip_low = np.quantile(colors_all, quantile_threshold[0], axis = 0)
    colors_all = colors_all-clip_low
    clip_high = np.quantile(colors_all, quantile_threshold[1], axis = 0)
    colors_all = colors_all/(clip_high + 0.01)
    Pseudo_colors = colors_all.dot(rgb_cols)
    Pseudo_colors[Pseudo_colors > 1] = 1
    Pseudo_colors[Pseudo_colors < 0] = 0
    return np.array(Pseudo_colors)

def get_features(adata, features):
    features = np.array(features)
    result = []
    for current_feature in features:
        if np.sum(ismember([current_feature], adata.var_names)[1]) > 0:
            result.append(list(get_gene(adata, current_feature)))
        else:
            result.append(adata.obs[current_feature].to_list())
    return(np.array(result).T)

def get_quant(x):
  " for each value in x, return which quantile it corresponds to "
  return np.array([np.count_nonzero(x<x_i)/(len(x)-1) for x_i in x])


# https://codeberg.org/tfardet/NNGT/commit/8d0e5b6c27796b534f9567a28e1c2b90440cf72c
def get_cmap(colormap, n=None):
    '''
    Get a colormap.

    Parameters
    ----------
    colormap : str or colormap
        Colormap to return.
    n : int, optional
        Take `n` samples from the colormap.
    '''
    if not isinstance(colormap, Colormap):
        colormap = matplotlib.colormaps[colormap]

    if n is None:
        return colormap

    # check version for call to resampled
    # @TODO require matplotlib > 3.6.0 in 2024 or something
    mpl_version = parse_version(matplotlib.__version__)
    min_version = parse_version("3.6.0")

    if mpl_version < min_version:
        return colormap._resample(n)

    return colormap.resampled(n)

def get_custom_legend_patches(condition_colors, selected_values):
    import matplotlib.patches as mpatches
    patchList = []
    for color, value in zip(condition_colors, selected_values):
        data_key = mpatches.Patch(facecolor=color, label=value, edgecolor='k')
        patchList.append(data_key)
    return(patchList)

def get_custom_legend_circles(selected_radius, color = 'navy'):
    patchList = []
    for radius in selected_radius:
        my_circle = Line2D(
            [0], [0], 
            marker='o', 
            color='w', 
            label='Circle',
            markerfacecolor='r', 
            markersize=radius)
        patchList.append(my_circle)    
    return(patchList)

def pseudo_color(
    adata_original, 
    coord, 
    feature_names, 
    rgb_cols, 
    selected_cells = [], 
    figsize = (10,10), 
    s = 1, 
    bounding_box = [], 
    use_outlines = False,
    quantile_threshold = [0, 0.95]
):
    if(len(bounding_box) > 0):
        adata = query_region(adata_original, coord, bounding_box)
    else:
        adata = adata_original
    Pseudo_colors = get_pseudo_colors(adata_original, feature_names, rgb_cols, quantile_threshold = quantile_threshold)
        
    fig,ax=plt.subplots(1,1,figsize=figsize)
    
    if(len(selected_cells) > 0):
        selected_cells = ismember(adata.obs_names, adata_original.obs_names[selected_cells])[1]
    else:
        selected_cells = [True] * adata.n_obs
    
    point_order = np.array(random.sample(list(range(adata.n_obs)), k=adata.n_obs))
    if use_outlines:
        ax.scatter(adata.obsm[coord][:,0][point_order],
                   adata.obsm[coord][:,1][point_order], 
                   color='black',s=s*5, rasterized=True)
        ax.scatter(adata.obsm[coord][:,0][point_order],
                   adata.obsm[coord][:,1][point_order], 
                   color='white',s=s*2, rasterized=True)
    point_order = np.array(random.sample(list(range(adata[selected_cells].n_obs)), k=np.sum(selected_cells)))
    Pseudo_colors = Pseudo_colors[selected_cells,:]
    ax.scatter(adata.obsm[coord][selected_cells,0][point_order],
               adata.obsm[coord][selected_cells,1][point_order], 
               c = Pseudo_colors[point_order,:],s=s, rasterized=True)
    #ax.set_axis_off()

def make_density_map_simple(
    xval, yval, subcells, 
    selected_contours, foreground_color, 
    background_color = None,
    log = False, total_bins = [100,100], gaussian_sigma = 1
):
    if log:
        xval[xval == 0] = np.min(xval[xval > 0])
        yval[yval == 0] = np.min(yval[yval > 0])
        xval = np.log(xval)
        yval = np.log(yval)
    
    import matplotlib.ticker as ticker

    h, binx, biny = np.histogram2d(xval, yval, bins=total_bins)
    
    diffx = np.diff(binx)[0]
    diffy = np.diff(biny)[0]
    
    extension = int(total_bins / 10)
    extend_x = np.array([binx[0] -diffx * (x+1) for x in range(extension)] + [binx[-1] + diffx * (x+1) for x in range(extension)])
    extend_y = np.array([biny[0] -diffy * (x+1) for x in range(extension)] + [biny[-1] + diffy * (x+1) for x in range(extension)])
    binx = np.sort(np.append(binx, extend_x))
    biny = np.sort(np.append(biny, extend_y))

    h_contrast, binx, biny = np.histogram2d(xval[subcells], yval[subcells], bins=[binx, biny])
    h_all, binx, biny = np.histogram2d(xval, yval, bins=[binx, biny])

    h_contrast_smoothed = scipy.ndimage.gaussian_filter(h_contrast, gaussian_sigma)
    h_all_smoothed = scipy.ndimage.gaussian_filter(h_all, gaussian_sigma)
    
    if(not background_color is None):
        if(type(background_color) == 'str'):
            background_color = matplotlib.colors.to_rgb(background_color)
        x_, y_ = np.meshgrid(binx[1:], biny[1:])
        z_grid = h_all_smoothed.T
        z_linear = np.ravel(z_grid)
        contour_all = plt.contour(
            x_,y_,z_grid, 
            levels=np.linspace(np.min(z_linear), np.max(z_linear), selected_contours), 
            linewidths=[1] * (selected_contours+1), 
            colors = np.array([background_color] * (selected_contours+1)), 
            locator = ticker.MaxNLocator(prune = 'lower')
        )    
    x_, y_ = np.meshgrid(binx[1:], biny[1:])
    z_grid = h_contrast_smoothed
    z_linear = np.ravel(z_grid)
    if(type(foreground_color) == 'str'):
        foreground_color = matplotlib.colors.to_rgb(background_color)
    contour_contrast = plt.contour(
        x_,y_,z_grid.T, 
        levels=np.linspace(np.min(z_linear), np.max(z_linear), selected_contours), 
        linewidths=[1] * (selected_contours+1),
        colors = np.array([foreground_color] * (selected_contours+1)), 
        locator = ticker.MaxNLocator(prune = 'lower')
    )
    
    plt.grid(False)
    plt.xticks([])
    plt.yticks([])
    return(x_, y_, z_grid)

def make_density_map(
    adata, 
    subcells, 
    obsm, 
    selected_contours, 
    foreground_color, 
    background_color = None,
    gaussian_sigma = 1
):
    import matplotlib.ticker as ticker
    total_bins = 100

    h, binx, biny = np.histogram2d(adata.obsm[obsm][:,0], adata.obsm[obsm][:,1], bins=total_bins)
    
    diffx = np.diff(binx)[0]
    diffy = np.diff(biny)[0]
    
    extension = int(total_bins / 10)
    extend_x = np.array([binx[0] -diffx * (x+1) for x in range(extension)] + [binx[-1] + diffx * (x+1) for x in range(extension)])
    extend_y = np.array([biny[0] -diffy * (x+1) for x in range(extension)] + [biny[-1] + diffy * (x+1) for x in range(extension)])
    binx = np.sort(np.append(binx, extend_x))
    biny = np.sort(np.append(biny, extend_y))

    h_contrast, binx, biny = np.histogram2d(adata.obsm[obsm][subcells,0], adata.obsm[obsm][subcells,1], bins=[binx, biny])
    h_all, binx, biny = np.histogram2d(adata.obsm[obsm][:,0], adata.obsm[obsm][:,1], bins=[binx, biny])
    
    h_contrast_smoothed = scipy.ndimage.gaussian_filter(h_contrast, gaussian_sigma)
    h_all_smoothed = scipy.ndimage.gaussian_filter(h_all, gaussian_sigma)
    
    if(not background_color is None):
        if(type(background_color) == 'str'):
            background_color = matplotlib.colors.to_rgb(background_color)
        x_, y_ = np.meshgrid(binx[1:], biny[1:])
        z_grid = h_all_smoothed.T
        z_linear = np.ravel(z_grid)
        contour_all = plt.contour(
            x_,y_,z_grid, 
            levels=np.linspace(np.min(z_linear), np.max(z_linear), selected_contours), 
            linewidths=[0] + [1] * (selected_contours), 
            colors = np.array([background_color] * (selected_contours+1)), 
            locator = ticker.MaxNLocator(prune = 'lower')
        )    
    x_, y_ = np.meshgrid(binx[1:], biny[1:])
    z_grid = h_contrast_smoothed
    z_linear = np.ravel(z_grid)
    if(type(foreground_color) == 'str'):
        foreground_color = matplotlib.colors.to_rgb(background_color)
    contour_contrast = plt.contour(
        x_,y_,z_grid.T, 
        levels=np.linspace(np.min(z_linear), np.max(z_linear), selected_contours), 
        linewidths=[0] + [1] * (selected_contours), 
        colors = np.array([foreground_color] * (selected_contours+1)), 
        locator = ticker.MaxNLocator(prune = 'lower')
    )
    
    plt.grid(False)
    plt.xticks([])
    plt.yticks([])
    return(x_, y_, z_grid)

def simple_beeswarm(y, nbins=None):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """
    y = np.asarray(y)
    if nbins is None:
        nbins = len(y) // 6

    # Get upper bounds of bins
    x = np.zeros(len(y))
    ylo = np.min(y)
    yhi = np.max(y)
    dy = (yhi - ylo) / nbins
    ybins = np.linspace(ylo + dy, yhi - dy, nbins - 1)

    # Divide indices into bins
    i = np.arange(len(y))
    ibs = [0] * nbins
    ybs = [0] * nbins
    nmax = 0
    for j, ybin in enumerate(ybins):
        f = y <= ybin
        ibs[j], ybs[j] = i[f], y[f]
        nmax = max(nmax, len(ibs[j]))
        f = ~f
        i, y = i[f], y[f]
    ibs[-1], ybs[-1] = i, y
    nmax = max(nmax, len(ibs[-1]))

    # Assign x indices
    dx = 1 / (nmax // 2)
    for i, y in zip(ibs, ybs):
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(y)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x

def get_shapely_polygon(df):
    coords = [Point(x,y) for x,y in zip(df['vertex_x'], df['vertex_y'])]
    result = Polygon(coords)
    return(result)

def read_10x_nucleus_boundaries(tenx_path, prefix, selected_cells = None, mode = 'nucleus'):
    if mode == 'nucleus':
        nucleus_boundaries = pd.read_parquet(os.path.join(tenx_path, 'nucleus_boundaries.parquet'))
    else:
        nucleus_boundaries = pd.read_parquet(os.path.join(tenx_path, 'cell_boundaries.parquet'))
    if not selected_cells is None:
        nucleus_boundaries = nucleus_boundaries.loc[ismember(nucleus_boundaries['cell_id'], selected_cells)[1],:].copy()
    grouped_nuclei = nucleus_boundaries.groupby('cell_id')
    unique_cells = np.unique(nucleus_boundaries['cell_id'])
    print("Extracting polygons")
    all_nuclei = [get_shapely_polygon(grouped_nuclei.get_group(current_cell)) for current_cell in unique_cells]
    nuclei = GeoDataFrame(geometry = all_nuclei, index = unique_cells)
    nuclei.index = [prefix + "__" + x for x in nuclei.index]
    return(nuclei)

def plot_nuclei_simple(nuclei, ax, color_value = None, subcells = None, cmap='viridis', edgecolor='lightgray', frameon=True, alpha = 0.5, linewidth = 1, background_color = 'white', vmin=None, vmax = None):
    if subcells is None:
        subcells = np.array([True] * nuclei.shape[0])
    current_cells = nuclei.loc[subcells,:].copy()
    current_cells.plot(color=color_value, edgecolor=edgecolor, linewidth = linewidth, ax=ax, alpha = alpha)
    plt.gca().set_facecolor(background_color)
    if(not frameon):
        plt.xticks([])
        plt.yticks([])
    plt.grid(False)
    return(nuclei)

def plot_nuclei_gradient(nuclei, ax, color_value = None, subcells = None, cmap='viridis', edgecolor='lightgray', frameon=True, alpha = 0.5, linewidth = 1, background_color = 'white', vmin=None, vmax = None):
    if subcells is None:
        subcells = np.array([True] * nuclei.shape[0])
    current_cells = nuclei.loc[subcells,:].copy()
    if not color_value is None:
        current_cells['color'] = color_value[subcells]
        current_cells.plot(column = 'color', edgecolor=edgecolor, linewidth = linewidth, ax=ax, alpha = alpha, cmap = cmap, vmin=vmin, vmax = vmax)
    plt.gca().set_facecolor(background_color)
    if(not frameon):
        plt.xticks([])
        plt.yticks([])
    plt.grid(False)
    return(nuclei)

def prepare_adata_window(adata_epi_sub, selected_cell, coordinate_extension = 60):
    selected_slide = adata_epi_sub[selected_cell,:].obs['slide_id'].to_numpy()[0]
    selected_coordinates = adata_epi_sub[selected_cell,:].obsm['X_spatial'][0].copy()
    selected_coordinates = [[x - coordinate_extension, x + coordinate_extension] for x in selected_coordinates]
    
    adata_sub = query_region(adata_epi_sub, "X_spatial", selected_coordinates[0], selected_coordinates[1])
    adata_sub = adata_sub[adata_sub.obs['slide_id'] == selected_slide,:].copy()
    return(adata_sub)

def plot_window(adata_epi_sub, selected_cell, coordinate_extension = 60, obs_key = 'cell_type_1', highlight_categories = [], highlight_colors = [], ax = None):
    adata_sub = prepare_adata_window(adata_epi_sub, selected_cell, coordinate_extension = coordinate_extension)
    ax = plot_individual_category(adata_sub, 'X_spatial', 'cell_type_0', ['Epithelial'], ['gray'], axes_color='white', figsize_max=2, s = 10, use_outlines=False, frameon=False, ax = ax)
    ax = plot_individual_category(adata_sub, 'X_spatial', obs_key, highlight_categories, highlight_colors, axes_color='white', figsize_max=2, s = 5, use_outlines=False, ax = ax)
    #plt.scatter(adata_sub.obsm['X_spatial'][:,0], adata_sub.obsm['X_spatial'][:,1], c = adata_sub.obs['progenitor_DC_original'], cmap = 'magma', s = 20)
    return(adata_sub, ax)

def plot_nuclei_categories_by_cell(adata_epi_sub, selected_cell, coordinate_extension = 60, obs_key = 'cell_type_1', highlight_categories = [], highlight_colors = [], ax = None, figsize = [5,5], edgecolor='k', mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    adata_sub = prepare_adata_window(adata_epi_sub, selected_cell, coordinate_extension = coordinate_extension)
    adata_sub, ax, nuclei = plot_nuclei_categories(adata_sub, obs_key = obs_key, highlight_categories = highlight_categories, highlight_colors = highlight_colors, ax = ax, figsize = figsize, edgecolor=edgecolor, mode = mode, alpha_background = alpha_background, alpha = alpha)
    return adata_sub, ax, nuclei

def plot_nuclei_categories_by_box(adata_epi_sub, selected_slide, bounding_box, obs_key = 'cell_type_1', highlight_categories = [], highlight_colors = [], ax = None, figsize = [5,5], edgecolor='k', mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    adata_sub = query_region(adata_epi_sub, "X_spatial", bounding_box[0], bounding_box[1])
    adata_sub = adata_sub[adata_sub.obs['slide_id'] == selected_slide,:].copy()
    adata_sub, ax, nuclei = plot_nuclei_categories(adata_sub, obs_key = obs_key, highlight_categories = highlight_categories, highlight_colors = highlight_colors, ax = ax, figsize = figsize, edgecolor=edgecolor, mode = mode, alpha_background = alpha_background, alpha = alpha)
    return adata_sub, ax, nuclei

def plot_nuclei_continuum_by_cell(adata_epi_sub, selected_cell, subcells = None, coordinate_extension = 60, obs_key = 'cell_type_1', highlight_categories = [], highlight_colors = [], ax = None, figsize = [5,5], cmap = "magma", vmin = None, vmax = None, edgecolor='none', mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    if subcells is None:
        subcells = np.array([True] * adata_epi_sub.n_obs)
    if vmin is None:
        vmin = np.quantile(adata_epi_sub.obs[obs_key][subcells], 0.05)
    if vmax is None:
        vmax = np.quantile(adata_epi_sub.obs[obs_key][subcells], 0.95)
    
    adata_epi_sub.obs['highlight'] = subcells.copy()
    adata_sub = prepare_adata_window(adata_epi_sub, selected_cell, coordinate_extension = coordinate_extension)
    adata_sub, ax, nuclei = plot_nuclei_continuum(adata_sub, obs_key = obs_key, ax = ax, figsize = figsize, cmap = cmap, vmin = vmin, vmax = vmax, edgecolor=edgecolor, mode = mode, alpha_background = alpha_background, alpha = alpha)
    return adata_sub, ax, nuclei

def plot_nuclei_continuum_by_box(adata_epi_sub, selected_slide, bounding_box, subcells = None, obs_key = 'progenitor_DC_original', cmap = 'magma', vmin = None, vmax = None, ax = None, figsize = [5,5], edgecolor='none', mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    if subcells is None:
        subcells = np.array([True] * adata_epi_sub.n_obs)
    if vmin is None:
        vmin = np.quantile(adata_epi_sub.obs[obs_key][subcells], 0.05)
    if vmax is None:
        vmax = np.quantile(adata_epi_sub.obs[obs_key][subcells], 0.95)
    
    adata_epi_sub.obs['highlight'] = subcells.copy()
    adata_sub = query_region(adata_epi_sub, "X_spatial", bounding_box[0], bounding_box[1])
    adata_sub = adata_sub[adata_sub.obs['slide_id'] == selected_slide,:].copy()
    adata_sub, ax, nuclei = plot_nuclei_continuum(adata_sub, obs_key = obs_key, cmap = cmap, vmin = vmin, vmax = vmax, ax = ax, figsize = figsize, edgecolor=edgecolor, mode = mode, alpha_background = alpha_background, alpha = alpha)
    return adata_sub, ax, nuclei

def plot_nuclei_categories(adata_sub, obs_key = 'cell_type_1', highlight_categories = [], highlight_colors = [], ax = None, figsize = [5,5], edgecolor=None, mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    selected_path = np.unique(adata_sub.obs['slide_path'])[0]
    selected_prefix = np.unique(adata_sub.obs['prefix'])[0]
    nuclei = read_10x_nucleus_boundaries(selected_path, selected_prefix, selected_cells = adata_sub.obs['cell_id'].to_numpy(), mode = mode)
    nuclei = nuclei.loc[adata_sub.obs_names,:].copy()
    subcells = ismember(adata_sub.obs[obs_key], highlight_categories)[1]
    
    if ax is None:
        f, ax = plt.subplots(1,1, figsize = figsize)
    background_cells = np.logical_not(subcells)
    plot_nuclei_simple(nuclei, ax, color_value = 'lightgray', subcells = background_cells, edgecolor='lightgray', frameon=True, alpha = alpha_background, linewidth = 1, background_color = 'white')
    for my_category, my_color in zip(highlight_categories, highlight_colors):
        if edgecolor is None:
            current_edge_color = my_color
        else:
            current_edge_color = edgecolor
        plot_nuclei_simple(nuclei, ax, color_value = my_color, subcells = adata_sub.obs[obs_key] == my_category, edgecolor=current_edge_color, frameon=True, alpha = alpha, linewidth = 1, background_color = 'white')
    plt.xticks([])
    plt.yticks([])
    
    return(adata_sub, ax, nuclei)

def plot_nuclei_continuum(adata_sub, obs_key = 'progenitor_DC_original', cmap = 'magma', vmin = None, vmax = None, ax = None, figsize = [5,5], edgecolor='none', mode = 'nucleus', alpha_background = 0.5, alpha = 1):
    subcells = adata_sub.obs['highlight'].to_numpy()
    
    selected_path = np.unique(adata_sub.obs['slide_path'])[0]
    selected_prefix = np.unique(adata_sub.obs['prefix'])[0]
    nuclei = read_10x_nucleus_boundaries(selected_path, selected_prefix, selected_cells = adata_sub.obs['cell_id'].to_numpy(), mode = mode)
    nuclei = nuclei.loc[adata_sub.obs_names,:].copy()
    nuclei[obs_key] = adata_sub.obs[obs_key]
    
    if ax is None:
        f, ax = plt.subplots(1,1, figsize = figsize)
    background_cells = np.logical_not(subcells)
    plot_nuclei_simple(nuclei, ax, color_value = 'lightgray', subcells = background_cells, edgecolor='lightgray', frameon=True, alpha = alpha_background, linewidth = 1, background_color = 'white')
    plot_nuclei_gradient(nuclei, ax, color_value = nuclei[obs_key].to_numpy(), subcells = subcells, edgecolor=edgecolor, frameon=True, alpha = alpha, linewidth = 1, background_color = 'white', vmin = vmin, vmax = vmax, cmap=cmap)
    plt.xticks([])
    plt.yticks([])
    
    return(adata_sub, ax, nuclei)
