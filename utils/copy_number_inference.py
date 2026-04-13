import scipy.cluster.hierarchy as shc
import collections
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
from general_utils import ismember

def readGenomeAnnotation_human(annotation_file):
    # Read genome annotation table
    genome_annotation = pd.read_csv(annotation_file, sep='\t', header=0)

    # Sort table by txStart position
    ordering = np.argsort(genome_annotation.loc[:,'txStart'])
    genome_annotation = genome_annotation.iloc[ordering,:]

    # Format names of chromosomes to get a numeric value
    special_region = [len(re.split('_', x)) > 1 for x in genome_annotation.loc[:,'chrom']]
    genome_annotation = genome_annotation.loc[np.logical_not(special_region),:]
    chromosome_names = [re.split('_', x)[0] for x in genome_annotation.loc[:,'chrom']]
    chromosome_names[:] = [s.replace('X', '23') for s in chromosome_names]
    chromosome_names[:] = [s.replace('Y', '24') for s in chromosome_names]
    chromosome_names[:] = [s.replace('M', '25') for s in chromosome_names]
    chromosome_names[:] = [s.replace('Un', '26') for s in chromosome_names]
    chromosome_names[:] = [s.replace('chr', '') for s in chromosome_names]

    # Transform chromosome names to strings and add a new column to the genome annotation table
    chromosome_names[:] = [int(i) for i in chromosome_names]
    genome_annotation.loc[:,'chrNum'] = pd.Series(chromosome_names, index=genome_annotation.index)

    # Sort table on the basis of chromosome number. The mergesort parameter is critical, since 
    # argsort will randomize indexes of duplicated values otherwise
    ordering = np.argsort(genome_annotation.loc[:,'chrNum'], kind='mergesort')
    genome_annotation = genome_annotation.iloc[ordering,:]

    # Transform all gene names to upper case
    genome_annotation.loc[:,'name2'] = [str.upper(s) for s in genome_annotation.loc[:,'name2']]
    
    return(genome_annotation)

def readGenomeAnnotation_mouse(annotation_file):
    # Read genome annotation table
    genome_annotation = pd.read_csv(annotation_file, sep='\t', header=0)

    # Sort table by txStart position
    ordering = np.argsort(genome_annotation.loc[:,'txStart'])
    genome_annotation = genome_annotation.iloc[ordering,:]

    # Format names of chromosomes to get a numeric value
    special_region = [len(re.split('_', x)) > 1 for x in genome_annotation.loc[:,'chrom']]
    genome_annotation = genome_annotation.loc[np.logical_not(special_region),:]
    chromosome_names = [re.split('_', x)[0] for x in genome_annotation.loc[:,'chrom']]
    chromosome_names[:] = [s.replace('X', '20') for s in chromosome_names]
    chromosome_names[:] = [s.replace('Y', '21') for s in chromosome_names]
    chromosome_names[:] = [s.replace('M', '22') for s in chromosome_names]
    chromosome_names[:] = [s.replace('Un', '23') for s in chromosome_names]
    chromosome_names[:] = [s.replace('chr', '') for s in chromosome_names]

    # Transform chromosome names to strings and add a new column to the genome annotation table
    chromosome_names[:] = [int(i) for i in chromosome_names]
    genome_annotation.loc[:,'chrNum'] = pd.Series(chromosome_names, index=genome_annotation.index)

    # Sort table on the basis of chromosome number. The mergesort parameter is critical, since 
    # argsort will randomize indexes of duplicated values otherwise
    ordering = np.argsort(genome_annotation.loc[:,'chrNum'], kind='mergesort')
    genome_annotation = genome_annotation.iloc[ordering,:]

    # Transform all gene names to upper case
    genome_annotation.loc[:,'name2'] = [str.upper(s) for s in genome_annotation.loc[:,'name2']]
    
    return(genome_annotation)

def slideWindowAverage(sourceMat, window, separator):
    targetMat = sourceMat.T.groupby(separator).rolling(window).mean().T.to_numpy()
    return(targetMat)

def getCopyNumberMatrix(dataMatrix, annotationMatrix, gene_ids, sample_ids, reference_cells, min_threshold = 0.1, window_size = 100, blacklisted=[]):
    """
    Params
    - dataMatrix: counts matrix in pandas data.frame format. Assumes linear counts (i.e. not log or z-scored)
    - annotationMatrix: data frame of genome annotation. Important fields are 'name2' (gene symbol) and 'chrNum'. Genes should be pre-sorted by chromosome, and position within chromosome (see readGenomeAnnotation function)
    - gene_ids: vector of gene symbols corresponding to columns in dataMatrix
    - sample_ids: vector of sample ids corresponding to rows of dataMatrix
    - reference_cells: cells as diploid reference.
    - min_threshold: threshold for minimum average counts in diploid reference to include gene in copy number inference
    - window_size: number of genes to consider when averaging sliding windows.
    - blacklisted: vector of gene symbols to exclude from CNV inference (optional)
    
    Output
    - smoothed_counts: Matrix of smoothed counts (log2)
    - annotationMatrix: subset of genome annotation matrix (rows correspond to columns of smoothed counts)
    """
    gene_ids = [str.upper(s) for s in gene_ids]
    dataMatrix.columns = gene_ids

    # Filter genes in the matrix on the basis of expression, presence in annotation table and uniqueness
    expression_threshold = min_threshold
    diploid_means = (dataMatrix.loc[reference_cells,:]).mean(axis=0) # Obtain the mean expression of each gene in the reference
    expression_filter = diploid_means > expression_threshold
    
    annotationMatrix = annotationMatrix.loc[~annotationMatrix.loc[:,'name2'].duplicated(),:]
    annotated_filter = ismember(dataMatrix.columns, annotationMatrix.loc[:,'name2'])[1]
    
    id_count = collections.Counter(dataMatrix.columns)
    unique_genes = [x == 1 for x in id_count.values()]
    
    validGenes = unique_genes
    validGenes = np.logical_and(validGenes, expression_filter)
    validGenes = np.logical_and(validGenes, annotated_filter)
    
    if(len(blacklisted) > 0):
        excluded_genes = ismember(gene_ids, blacklisted)[1]
        validGenes = np.logical_and(validGenes, np.logical_not(excluded_genes))

    dataMatrix = dataMatrix.loc[:,validGenes]
    diploid_means = diploid_means[validGenes]
    
    annotationMatrix = annotationMatrix.loc[ismember(annotationMatrix.loc[:,'name2'], dataMatrix.columns)[1],:]
    
    # Reorder matrix so to match the indexes in annotationMatrix 
    # (assumes that this has already been sorted by chromosome position and transcription start)
    dataMatrix = dataMatrix.loc[:,annotationMatrix.loc[:,'name2']]

    source_counts = np.array(dataMatrix)

    # Obtain the mean expression of each gene in the reference
    diploid_means = (source_counts[reference_cells,:]).mean(axis=0)

    # For each gene in each cell, calculate the ratio of counts relative to the average count in the reference sample
    # Compute log2 of this ratio
    pseudocount = 0.1
    normalized_counts = (source_counts + pseudocount) / (diploid_means[:,None].T + pseudocount)
    normalized_counts = np.log(normalized_counts) / np.log(2)

    # Set the minimum and maximum log2 ratio to be -3 and 3, respectively.
    normalized_counts = normalized_counts.clip(-3,3)
    normalized_counts = pd.DataFrame(normalized_counts, columns = annotationMatrix.loc[:,'name2'])
    smoothed_counts = slideWindowAverage(normalized_counts, window_size, np.array(annotationMatrix['chrNum']))
    validColumns = np.isnan(smoothed_counts).sum(axis=0) == 0;
    smoothed_counts = smoothed_counts[:,validColumns]
    annotationMatrix = annotationMatrix.loc[validColumns,:]
    print('Recentering medians')
    
    # Recenter smoothed counts to a median of 0
    median_counts = np.nanmedian(smoothed_counts,axis=1)
    smoothed_counts = smoothed_counts - median_counts[:,None]
    print('done')
    return smoothed_counts, annotationMatrix

def fastClusterSmoothedCounts(smoothed_counts, sample_ids, chrNum, selected_chromosomes):
    unique_samples = np.unique(sample_ids)
    indexes = range(len(sample_ids))
    indexes = np.array(indexes)
    cluster_assignment = np.array([0.] * len(sample_ids))
    
    unique_chromosomes = np.unique(chrNum)
    chromosome_filter = ismember(unique_chromosomes, selected_chromosomes)[1]

    smoothed_counts_clustered = smoothed_counts.copy()
    
    chrom_average = []
    for chrom in unique_chromosomes:
        chrom_average.append(np.ravel(np.median(smoothed_counts[:,chrNum == chrom], axis=1)))
    chrom_average = np.array(pd.DataFrame(chrom_average).T)
    
    for i in range(len(unique_samples)):
        try:
            s = unique_samples[i]
            display(s)
            validRows = sample_ids == s
            currentIndexes = indexes[validRows]
            temp = chrom_average[validRows,:]
            linkage = shc.linkage(temp[:,np.where(chromosome_filter)[0]], method='ward')
            dend = shc.dendrogram(linkage, no_plot=True)
            smoothed_counts_clustered[validRows,:] = smoothed_counts[currentIndexes[dend['leaves']],:]
            indexes[validRows] = currentIndexes[dend['leaves']]
        except ValueError:
            print("Did not cluster " + s)

    return smoothed_counts_clustered, indexes

def clusterSmoothedCounts(smoothed_counts, sample_ids, chrNum, selected_chromosomes):
    unique_samples = list(set(sample_ids))
    indexes = range(len(sample_ids))
    indexes = np.array(indexes)
    cluster_assignment = np.array([0.] * len(sample_ids))
    
    chromosome_filter = ismember(chrNum, selected_chromosomes)[1]

    smoothed_counts_clustered = smoothed_counts.copy()

    for i in range(len(unique_samples)):
        try:
            s = unique_samples[i]
            display(s)
            validRows = sample_ids == s
            currentIndexes = indexes[validRows]
            temp = smoothed_counts[validRows,:]
            linkage = shc.linkage(temp[:,np.where(chromosome_filter)[0]], method='ward')
            dend = shc.dendrogram(linkage, no_plot=True)
            smoothed_counts_clustered[validRows,:] = temp[dend['leaves'],:]
            indexes[validRows] = currentIndexes[dend['leaves']]
            subclusters = shc.cut_tree(linkage, height = 100)
            cluster_assignment[validRows] = [format(x, '2.2f') for x in np.ravel(np.float32(subclusters))[dend['leaves']] / 100 + i]
        except ValueError:
            print("Did not cluster " + s)

    return smoothed_counts_clustered, indexes    

def plot_CNV_matrix(data_matrix, row_variable, col_variable, selected_cells, ax):
    selected_cells = np.array(selected_cells)
    row_variable = np.array(row_variable)[selected_cells]
    col_variable = np.array(col_variable)
    
    col_change = np.where(col_variable[:-1] != col_variable[1:])[0]
    
    row_num = ismember(row_variable, np.unique(row_variable))[0]
    row_change = np.where([a!=b for a,b in zip(row_num[:-1], row_num[1:])])[0]

    plt.rcParams["axes.grid"] = True
    
    im1 = ax.imshow(data_matrix[selected_cells,:], vmin=-1, vmax=1, aspect = 'auto', interpolation = None)
    ax.set_xticks(col_change)
    ax.set_yticks(row_change)
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
    im1.set_cmap("RdBu_r")
    ax.set_xlabel('Genomic position')
    cbar = plt.colorbar(im1, ax=ax, orientation='horizontal', shrink=0.6, pad=0.1)
    cbar.ax.set_xlabel('Smoothed counts (log2)')
