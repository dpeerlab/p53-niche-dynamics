from itertools import chain
import os
import numpy as np
import scipy

def ismember(a, b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    result = [bind.get(itm, None) for itm in a]  # None can be replaced by any other "not in b" value
    result = [result, [not x is None for x in result]]
    return result

def precompute_for_ismember(b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return(bind)

def precompute_for_ismember_comprehensive(b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = [i]
        else:
            bind[elt] += [i]
    return(bind)

def ismember_precomputed(a,bind):
    result = np.array([bind.get(itm, None) for itm in a])  # None can be replaced by any other "not in b" value
    result = [result, np.array([not x is None for x in result])]
    return result

def ismember_precomputed_comprehensive(a,bind):
    result = [bind.get(itm, None) for itm in a]  # None can be replaced by any other "not in b" value
    result = [result, [not x is None for x in result]]
    return result
    
def flatten(my_list):
    return(list(itertools.chain(*my_list)))
    
def grep(pattern, my_list):
    valid_entries = np.array([pattern in x for x in my_list])
    return np.array(my_list)[valid_entries]

def grep_indexes(my_pattern, my_list):
    my_list = np.array(my_list)
    my_indexes = np.array([my_pattern in x for x in my_list])
    return(my_indexes)

def grep_exclude(pattern, my_list):
    valid_entries = np.array([pattern in x for x in my_list])
    return np.array(my_list)[np.logical_not(valid_entries)]

def write_object(myObject, file):
    import dill
    with open(file, 'wb') as f:
        f.write(dill.dumps(myObject))

def read_object(file):
    import dill
    with open(file, 'rb') as f:
        myObject = dill.loads(f.read())
    return myObject

# http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def flatten_chain(matrix):
    return list(chain.from_iterable(matrix))

def array_to_csr_matrix(mat):
    """Get csr_matrix from list of list graph representation

    Params:
        - mat: list of lists in np.array format, each row in the array corresponds to the source node, values correspond to indexes the connected nodes. It assumes that the row indexes are in the same space as the values in the array.

    Returns
    -------
        csr_matrix representation of graph
    """
    sources, targets = array_to_edges(mat)
    result = scipy.sparse.csr_matrix(([1] * len(sources), (sources, targets)), shape=[mat.shape[0]] * 2)
    return result

def array_to_edges(mat):
    """Get csr_matrix from list of list graph representation

    Params:
        - mat: list of lists in np.array format, each row in the array corresponds to the source node, values correspond to indexes the connected nodes. It assumes that the row indexes are in the same space as the values in the array.

    Returns
    -------
        sources and targets of the graph
    """
    targets = np.ravel(mat)
    n_neighbors_per_cell = [len(x) for x in mat]
    sources = np.ravel([[i] * k for i, k in zip(range(mat.shape[0]), n_neighbors_per_cell)])
    return sources, targets

