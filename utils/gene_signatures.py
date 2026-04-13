from tqdm import tqdm
import os
import pandas as pd
import numpy as np
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

def toupper(my_list):
    return(np.array([str(x).upper() for x in my_list]))

def convert_to_mouse(gene_list):
    try:
        h2m
    except NameError:
        from gseapy import Biomart
        bm = Biomart()
        h2m = bm.query(dataset='hsapiens_gene_ensembl', attributes=['external_gene_name','mmusculus_homolog_associated_gene_name'])
        h2m = h2m.astype(str)
        h2m = h2m.loc[np.logical_not(h2m['mmusculus_homolog_associated_gene_name'] == 'nan'),:].copy()
    converted_list = toupper(np.unique(h2m['mmusculus_homolog_associated_gene_name'][ismember(h2m['external_gene_name'], gene_list)[1]]))
    return(converted_list)

def attach_selected_signatures(adata, signature_path, selected_signatures, layer = 'zscore', human_symbols = False):
    possible_signatures = np.array(os.listdir(signature_path))
    for current_signature in tqdm(selected_signatures):
        current_file = grep(current_signature, possible_signatures)[0]
        current_genes = toupper(pd.read_csv(os.path.join(signature_path, current_file))['gene_id'].to_numpy())
        if human_symbols:
            current_genes = convert_to_mouse(current_genes)
        attach_custom_signature(adata, current_genes, current_signature, layer)

def attach_custom_signature(adata, signature_genes, signature_name, layer):
    valid_genes = ismember(adata.var_names,signature_genes)[1]
    adata.var[str(signature_name)] = valid_genes
    if layer == 'X':
        adata.obs[str(signature_name)] = np.array(np.nanmean(adata.X[:,valid_genes], axis=1))
    else:
        adata.obs[str(signature_name)] = np.array(np.nanmean(adata.layers[layer][:,valid_genes], axis=1))

def attach_signatures_low_memory(adata, signature_path, selected_signatures):
    data_genes = np.array([str(x).upper() for x in adata.var_names])
    my_signatures = [os.path.join(signature_path, x) for x in os.listdir(signature_path)]
    signature_names = []
    signature_dict = {}
    for current_signature in tqdm(my_signatures):
        if os.path.isdir(current_signature):
            continue
        current_file = grep(current_signature, possible_signatures)[0]
        gene_set = pd.read_csv(os.path.join(signature_path, current_file), index_col = None, sep = "\t")
        gene_set = np.array(gene_set.iloc[:,0])
        gene_set = np.array([str(x).upper() for x in gene_set])
        adata_temp = adata[:,ismember(data_genes, gene_set)[1]].copy()
        sc.pp.scale(adata_temp)
        adata.obs[current_name] = np.ravel(adata_temp.X.mean(axis=1))
        adata.var[current_name] = ismember(data_genes, gene_set)[1]
        signature_names.append(current_name)
        signature_dict[current_name] = gene_set.copy()
    return(signature_names, signature_dict)
