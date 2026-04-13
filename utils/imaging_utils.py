import cv2
import numpy as np
from skimage import exposure, restoration, util, img_as_float
import scipy.ndimage as ndi
from scipy.ndimage._ni_support import _normalize_sequence
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

def imnormalize(img, quantile_low = 0.1, quantile_high = 0.99):
    img = img - np.quantile(np.ravel(img), quantile_low)
    img = img / np.quantile(np.ravel(img), quantile_high)
    img = np.clip(img, 0, 1)
    return(img)

def im2rgb(IM, color, quantile_high = 0.99, quantile_low = 0.1):
    norm = imnormalize(IM, quantile_low = quantile_low, quantile_high = quantile_high)
    RGB = np.zeros(list(norm.shape) + [3])
    for i in range(len(color)):
        RGB[...,i] = color[i] * norm
    return(RGB)

def imblend(G, R = [], B = []):
    G = im2rgb(G, [0,1,0])
    if len(R) > 0:
        R = im2rgb(R, [1,0,0])
    else:
        R = np.zeros(R.shape)
    if len(B) > 0:
        B = im2rgb(B, [0,0,1])
    else:
        B = np.zeros(B.shape)
    IM = R + G + B
    return(IM)

def imoverlay(IM, mask, color=[1,1,1], px=0):
    from skimage.segmentation import find_boundaries
    if len(IM.shape) == 2:
        overlay_data = im2rgb(IM, [1,1,1])
    else:
        overlay_data = IM.copy()
    boundary = find_boundaries(mask, connectivity=1, mode='inner')
    boundary_mask = np.zeros(overlay_data.shape[:2])
    boundary_mask[boundary > 0] = 1
    if px > 0:
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(px,px))
        boundary_mask = cv2.dilate(boundary_mask,kernel,iterations = 1)
    for i in range(len(color)):
        overlay_data[boundary_mask > 0,i] = color[i]
    return(overlay_data)

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