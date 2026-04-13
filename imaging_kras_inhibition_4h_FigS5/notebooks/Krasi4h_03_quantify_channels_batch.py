import os
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
import dask.array as da
import skimage.measure
import spatialdata as sd
import cv2
from multiprocessing import Pool
import itertools
import argparse
from tqdm import tqdm

from skimage.feature import match_template
import skimage.registration
import skimage.restoration
import skimage.util
from skimage.exposure import rescale_intensity
from skimage.segmentation import find_boundaries
from skimage.measure import label, regionprops
from PIL import Image
from skimage import transform, io, exposure, restoration
import scipy.misc

from skimage import exposure, restoration, util, img_as_float
import scipy.ndimage as ndi
from scipy.ndimage._ni_support import _normalize_sequence

def ismember(a, b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    result = np.array([bind.get(itm, None) for itm in a])  # None can be replaced by any other "not in b" value
    result = [result, np.array([not x is None for x in result])]
    return result

def grep(pattern, my_list):
    valid_entries = np.array([pattern in x for x in my_list])
    return np.array(my_list)[valid_entries]
    
def grep_exclude(pattern, my_list):
    valid_entries = np.array([pattern in x for x in my_list])
    return np.array(my_list)[np.logical_not(valid_entries)]

def imnormalize(img, quantile_low = 0.1, quantile_high = 0.99):
    img = img - np.quantile(np.ravel(img), quantile_low)
    img = img / np.quantile(np.ravel(img), quantile_high)
    img = np.clip(img, 0, 1)
    return(img)

def get_image(sdata, key, scale, channel):
    channels = np.array(sd.models.get_channel_names(sdata.images[key])).tolist()
    idx = channels.index(channel)
    return sdata.images[key][scale].image[idx].values

def get_image_slice(sdata, key, scale, channel, row_slice, col_slice):
    channels = np.array(sd.models.get_channel_names(sdata.images[key])).tolist()
    idx = channels.index(channel)
    return sdata.images[key][scale].image[idx][row_slice, col_slice].compute().values    

def erode_labeled_image(labeled_mask, kernel_size = 3):
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
    # To be able to erode per nucleus, I first find the nucleus boundary, then
    # dilate this boundary and us the dilated boundary to punch the original nuclear
    # mask. 
    boundary = find_boundaries(labeled_mask, connectivity=1, mode='inner')
    boundary_mask = (boundary > 0).astype(np.float32)
    boundary_mask = cv2.dilate(boundary_mask, kernel, iterations = 1)
    boundary_mask = -1 * (boundary_mask - 1)
    new_mask = np.multiply(labeled_mask, boundary_mask)
    return(new_mask)

def filter_small_objects(labeled_mask, AREA_THRESHOLD = 20):
    # Exclude nuclei that become too small
    regionprops = skimage.measure.regionprops(label_image = labeled_mask.astype(int))
    labels = np.array([x.label for x in regionprops])
    areas = np.array([x.area for x in regionprops])
    valid_labels = labels[areas >= AREA_THRESHOLD]
    return(valid_labels)

def filter_objects(labeled_image, valid_labels):
    included_mask = np.isin(labeled_image, valid_labels).astype(np.float32)
    labeled_image_filtered = np.multiply(labeled_image, included_mask)
    return(labeled_image_filtered)
    
def rolling_ball_filter(data, ball_radius, spacing=None, top=False, **kwargs):
    """Rolling ball filter implemented with morphology operations

    This implenetation is very similar to that in ImageJ and uses a top hat transform
    with a ball shaped structuring element
    https://en.wikipedia.org/wiki/Top-hat_transform
    #https://stackoverflow.com/questions/29320954/rolling-ball-background-subtraction-algorithm-for-opencv
    
    Parameters
    ----------
    data : ndarray
        image data (assumed to be on a regular grid)
    ball_radius : float
        the radius of the ball to roll
    spacing : int or sequence
        the spacing of the image data
    top : bool
        whether to roll the ball on the top or bottom of the data
    kwargs : key word arguments
        these are passed to the ndimage morphological operations

    Returns
    -------
    data_nb : ndarray
        data with background subtracted
    bg : ndarray
        background that was subtracted from the data
    """
    
    resizeFactor = 1/2
    original_size = data.shape
    data = cv2.resize(data, dsize = [int(x * resizeFactor) for x in original_size])
    radius = int(ball_radius * resizeFactor)
    
    ndim = data.ndim
    if spacing is None:
        spacing = 1

    spacing = _normalize_sequence(spacing, ndim)
    radius = np.asarray(_normalize_sequence(ball_radius, ndim))
    mesh = np.array(np.meshgrid(*[np.arange(-r, r + s, s) for r, s in zip(radius, spacing)], indexing="ij"))
    structure = 2 * np.sqrt(1 - ((mesh / radius.reshape(-1, *((1,) * ndim)))**2).sum(0))
    structure[~np.isfinite(structure)] = 0
    if not top:
        #ndi.white_tophat(data, structure=structure)
        background = ndi.grey_erosion(data, structure=structure, **kwargs)
        background = ndi.grey_dilation(background, structure=structure, **kwargs)
    else:
        #background = ndi.black_tophat(data, structure=structure, output=background)
        background = ndi.grey_dilation(data, structure=structure, **kwargs)
        background = ndi.grey_erosion(background, structure=structure, **kwargs)
    
    background = cv2.resize(background, dsize = original_size[::-1])
    data = cv2.resize(data, dsize = original_size[::-1])
    
    result = data - background
    result[result < 0] = 0
    result[data == 0] = 0
    return result

def top_10percent_px(region_mask, intensity_image):
    region_mask = np.ravel(region_mask)
    intensity_image = np.ravel(intensity_image)
    valid_pixels = intensity_image[region_mask > 0]
    return(np.median(valid_pixels[valid_pixels >= np.quantile(valid_pixels, 0.9)]))

def median_px(region_mask, intensity_image):
    region_mask = np.ravel(region_mask)
    intensity_image = np.ravel(intensity_image)
    valid_pixels = intensity_image[region_mask > 0]
    return(np.median(valid_pixels))

import skimage.measure
def get_intensity_measurement(label_image, intensity_image, channel_name, mask_name):
    regionprops = skimage.measure.regionprops(label_image = label_image.astype(int), intensity_image = intensity_image.astype(int))
    areas = [x.area for x in regionprops]
    solidities = [x.solidity for x in regionprops]
    intensities = [x.intensity_mean for x in regionprops]
    centroids = [x.centroid for x in regionprops]
    centroidx = [x[0] for x in centroids]
    centroidy = [x[1] for x in centroids]
    
    total_intensities = [np.sum(x*y) for x,y in zip(areas, intensities)]
    labels = ["cell_%d" % (x.label) for x in regionprops]
            
    df_intensity = pd.DataFrame({"mean_intensity": intensities, "total_intensity": total_intensities})
    df_intensity.columns = [channel_name + "_" + mask_name + "_" + x for x in df_intensity.columns]
    df_intensity.index = labels

    # Centroids only make sense when dealing with global coordinates
    df_objects = pd.DataFrame({"area": areas, "solidity": solidities, 'centroid_row': centroidx, 'centroid_col': centroidy})
    df_objects = pd.DataFrame({"area": areas, "solidity": solidities})
    df_objects.columns = [mask_name + "_" + x for x in df_objects.columns]
    df_objects.index = labels
    return(df_intensity, df_objects)

def get_foci_measurements(label_image, intensity_image, channel_name, mask_name):
    regionprops = skimage.measure.regionprops(
        label_image = label_image.astype(int), 
        intensity_image = intensity_image.astype(int), 
        extra_properties=[top_10percent_px, median_px]
    )
    top_intensities = np.array([x.top_10percent_px for x in regionprops])
    median_intensities = np.array([x.median_px for x in regionprops])
    mean_intensity = np.array([x.intensity_mean for x in regionprops])
    if(np.sum(median_intensities > 0) > 0):
        pseudocount = np.min(median_intensities[median_intensities > 0])
    else:
        pseudocount = 1
    median_intensities = median_intensities + pseudocount
    top_intensities = top_intensities + pseudocount        
    labels = ["cell_%d" % (x.label) for x in regionprops]
    
    df_intensity = pd.DataFrame({"top_intensity": top_intensities, "median_intensity": median_intensities, 'intensity_ratio': top_intensities/median_intensities})
    df_intensity.columns = [channel_name + "_" + mask_name + "_" + x for x in df_intensity.columns]
    df_intensity.index = labels
    
    return(df_intensity)

def quantify_masks(intensity_image, mask_images, mask_names, channel_name, background_radius = None, quantify_foci = False):
    if not (background_radius is None):
        IM_background = rolling_ball_filter(intensity_image, background_radius)
    all_results_intensities = []
    all_results_objects = []
    for my_mask, my_name in zip(mask_images, mask_names):
        current_result_intensity, current_result_objects = get_intensity_measurement(my_mask, intensity_image, channel_name, my_name)
        if(quantify_foci):
            current_result_foci = get_foci_measurements(my_mask, intensity_image, channel_name, my_name)
            consolidated_table = pd.concat([current_result_foci, current_result_intensity], axis = 1)
        else:
            consolidated_table = current_result_intensity
        all_results_intensities.append(consolidated_table)
        all_results_objects.append(current_result_objects)
    result_final_intensity = pd.concat(all_results_intensities, axis = 1)
    result_final_objects = pd.concat(all_results_objects, axis = 1)
    return([result_final_intensity, result_final_objects])

def process_patch(
    sdata, 
    mask_nucleus_original, 
    mask_cyto_original, 
    df_centroids, 
    row_start, 
    col_start, 
    height=1000, 
    width = 1000, 
    offset = 50
):
    bounding_box = [[row_start, row_start + height - 1], [col_start, col_start + width - 1]]
    bounding_box_offset = [
        [np.max([bounding_box[0][0]-offset, 0]),
         np.min([bounding_box[0][1]+offset, mask_nucleus_original.shape[0]])],
        [np.max([bounding_box[1][0]-offset, 0]),
         np.min([bounding_box[1][1]+offset, mask_nucleus_original.shape[1]])]
    ]
    row_filter = np.logical_and(df_centroids['centroid_row'] >= bounding_box[0][0], df_centroids['centroid_row'] <= bounding_box[0][1])
    col_filter = np.logical_and(df_centroids['centroid_col'] >= bounding_box[1][0], df_centroids['centroid_col'] <= bounding_box[1][1])
    cells_in_box = df_centroids['label'][np.logical_and(row_filter, col_filter)].to_numpy()
    
    row_slice = slice(bounding_box_offset[0][0], bounding_box_offset[0][1])
    col_slice = slice(bounding_box_offset[1][0], bounding_box_offset[1][1])
    
    mask_cyto = mask_cyto_original[row_slice, col_slice].copy()
    mask_nucleus = mask_nucleus_original[row_slice, col_slice].copy()
    
    mask_nucleus_eroded = erode_labeled_image(mask_nucleus, kernel_size = 3)
    valid_labels = filter_small_objects(mask_nucleus_eroded, AREA_THRESHOLD = 20)
    valid_labels = np.intersect1d(valid_labels, cells_in_box) # Filter objects in non-offset box to avoid redundant quantification
    # This last step is important because it ensures that nucleus and cytoplasmic masks are aligned
    mask_nucleus_eroded_filtered = filter_objects(mask_nucleus_eroded, valid_labels)
    mask_cyto_filtered = filter_objects(mask_cyto, valid_labels)
    mask_nucleus_original_filtered = filter_objects(mask_nucleus, valid_labels)

    if(np.sum(np.ravel(mask_nucleus_eroded_filtered) > 0) == 0):
        return(None)

    selected_channels = ['DAPI', 'dpERK', 'MSN', 'HMGA2', 'LGALS4', 'E-cadherin', 'GFP', 'ANXA10', 'VIM', 'CC3', 'CPA1']
    selected_masks = [mask_cyto_filtered, mask_nucleus_eroded_filtered]
    mask_names = ["cyto", "nucleus"]
    
    result = []
    # Quantify all channels
    for current_channel in selected_channels:
        current_intensity_image = get_image_slice(sdata, "COMET", "scale1", current_channel, row_slice, col_slice)
        df_intensity, df_objects = quantify_masks(current_intensity_image, selected_masks, mask_names, current_channel, None)
        result.append(df_intensity)
    # Quantify p19 foci
    current_channel = "p19"
    current_intensity_image = get_image_slice(sdata, "COMET", "scale1", current_channel, row_slice, col_slice)
    df_intensity, df_objects = quantify_masks(current_intensity_image, selected_masks, mask_names, current_channel, None, quantify_foci=True)
    result.append(df_intensity)
    # Consolidate results and add global centroid information
    result_final = pd.concat(result + [df_objects], axis = 1)
    result_final['centroid_row'] = df_centroids.loc[result_final.index, "centroid_row"]
    result_final['centroid_col'] = df_centroids.loc[result_final.index, "centroid_col"]
    return(result_final)

def process_patch_simple(coords):
    result = process_patch(
        sdata, mask_nucleus_original, mask_cyto_original, df_centroids, 
        row_start=coords[0], col_start = coords[1], height = window_size, width = window_size, offset = 50
    )
    return(result)

if __name__ == '__main__':
    # Get slide id from arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--slide')
    args = parser.parse_args()
    SLIDE_ID = args.slide
    # Read sdata object
    main_path = ".."
    sdata_path = os.path.join(main_path, "sdata", SLIDE_ID)
    sdata = sd.read_zarr(sdata_path)
    # Read masks
    print("Reading mask files")
    mask_path = os.path.join(main_path, "cellpose_masks_consolidated")
    possible_files = np.sort(os.listdir(mask_path))
    mask_nucleus_file = grep(SLIDE_ID, grep("nucleus", possible_files))[0]
    mask_cyto_file = grep(SLIDE_ID, grep("cyto", possible_files))[0]
    mask_nucleus_original = np.load(os.path.join(mask_path, mask_nucleus_file))['arr_0']
    mask_cyto_original = np.load(os.path.join(mask_path, mask_cyto_file))['arr_0']
    # Get global coordinate system
    print("Getting global coordinate system")
    regionprops = skimage.measure.regionprops(label_image = mask_nucleus_original.astype(int))
    centroids = [x.centroid for x in regionprops]
    centroidx = [x[0] for x in centroids]
    centroidy = [x[1] for x in centroids]
    labels = [x.label for x in regionprops]
    df_centroids = pd.DataFrame({'centroid_row': centroidx, 'centroid_col': centroidy, 'label': labels})
    df_centroids.index = ["cell_" + str(x) for x in labels]
    # Iterate over image patches
    print("Iterating over windows")
    window_size = 1000
    row_ranges = list(range(0, mask_nucleus_original.shape[0], window_size))
    col_ranges = list(range(0, mask_nucleus_original.shape[1], window_size))
    joint_ranges = itertools.product(row_ranges, col_ranges)
    result = []
    for my_range in tqdm(list(joint_ranges)):
        temp = process_patch(
            sdata, mask_nucleus_original, mask_cyto_original, df_centroids, 
            row_start=my_range[0], col_start = my_range[1], height = window_size, width = window_size, offset = 50
        )
        if not temp is None:
            result.append(temp)
    # Consolidate final results and write to file
    result_final = pd.concat(result, axis=0)
    output_dir = os.path.join(main_path, "measurements")
    os.makedirs(output_dir, exist_ok = True)
    result_final.to_csv(os.path.join(output_dir, SLIDE_ID + "__measurements.csv.gz"))
