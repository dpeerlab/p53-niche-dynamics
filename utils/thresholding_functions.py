import numpy as np
import pandas as pd

# From skimage package
def threshold_otsu(counts, bin_centers):
    """Return threshold value based on Otsu's method.

    Either image or hist must be provided. If hist is provided, the actual
    histogram of the image is ignored.

    Parameters
    ----------
    hist : array, or 2-tuple of arrays, optional
        Histogram from which to determine the threshold, and optionally a
        corresponding array of bin center intensities. If no hist provided,
        this function will compute it from the image.


    Returns
    -------
    threshold : float
        Upper threshold value. All pixels with an intensity higher than
        this value are assumed to be foreground.

    References
    ----------
    .. [1] Wikipedia, https://en.wikipedia.org/wiki/Otsu's_Method

    """
    # class probabilities for all possible thresholds
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(counts * bin_centers) / weight1
    mean2 = (np.cumsum((counts * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of ``weight1``/``mean1`` should pair with zero values in
    # ``weight2``/``mean2``, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[idx]

    return threshold

# From skimage package, motified to handle low or high threshold
def threshold_triangle(counts, bin_centers, upper = True):
    bin_centers = bin_centers
    hist = counts.copy()
    nbins = len(hist)

    # Find peak, lowest and highest gray levels.
    arg_peak_height = np.argmax(hist)
    peak_height = hist[arg_peak_height]
    arg_low_level, arg_high_level = np.flatnonzero(hist)[[0, -1]]

    if(upper):
        flip=True
    else:
        flip=False
    # Flip is True if left tail is shorter.
    #flip = arg_peak_height - arg_low_level < arg_high_level - arg_peak_height
    if flip:
        hist = hist[::-1]
        arg_low_level = nbins - arg_high_level - 1
        arg_peak_height = nbins - arg_peak_height - 1

    # If flip == True, arg_high_level becomes incorrect
    # but we don't need it anymore.
    del(arg_high_level)

    # Set up the coordinate system.
    width = arg_peak_height - arg_low_level
    x1 = np.arange(width)
    y1 = hist[x1 + arg_low_level]

    # Normalize.
    norm = np.sqrt(peak_height**2 + width**2)
    peak_height /= norm
    width /= norm

    # Maximize the length.
    # The ImageJ implementation includes an additional constant when calculating
    # the length, but here we omit it as it does not affect the location of the
    # minimum.
    length = peak_height * x1 - width * y1
    arg_level = np.argmax(length) + arg_low_level

    if flip:
        arg_level = nbins - arg_level - 1

    return bin_centers[arg_level]

def histogram_knee(h, plot=False):
    # This is triangle method for thresholding
    curve =  [1-x for x in h[0]]
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T

    #firstPoint = allCoord[0]
    #Assume that the threshold is higher than half-point the x axis
    min_threshold = int(allCoord.shape[0]/2)
    min_threshold = np.min(np.where(h[1] > 0)[0])
    firstPoint = allCoord[np.argmin(allCoord[:min_threshold,1])]
    #minPoint = allCoord[-1]
    #Draw line to maximum of distribution
    min_point_coordinate = np.argmin(allCoord[:,1])
    minPoint = allCoord[min_point_coordinate]
    lineVec = minPoint - firstPoint

    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))

    vecFromFirst = allCoord - firstPoint

    scalarProduct = np.sum(vecFromFirst * np.tile(lineVecNorm, (nPoints, 1)), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine[min_threshold:min_point_coordinate]) + min_threshold
    if plot:
        f = plt.figure()
        plt.plot(allCoord[:,0], allCoord[:,1])
        plt.plot([firstPoint[0], minPoint[0]], [firstPoint[1], minPoint[1]])
        plt.plot(vecFromFirstParallel[:,0], vecFromFirstParallel[:,1])
        plt.plot(vecFromFirst[:,0], vecFromFirst[:,1])
        plt.plot(vecFromFirst[:,0], distToLine)
        plt.plot([allCoord[idxOfBestPoint,0]]*2, plt.gca().get_ylim(), '--')
    return(idxOfBestPoint)