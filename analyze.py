"""TODO: docstring

"""
from scipy import mat, zeros, arange, array
from data import SeqList, Seq
from math import log

def calc_all(index_function, **kwargs):
    seq_list = kwargs['seq_list']
    for seq in seq_list:
        seq.set_attr(index_function.__name__, index_function(seq, **kwargs))
        
def calc_dist_matrix(seq_list):
    dim = len(seq_list)
    out_matrix = mat(zeros((dim, dim)))
    seq_str_list = []
    for i in range(dim):
        seq_str_list += [str(seq_list[i])]
        for j in range(dim):
            out_matrix[i,j] = seq_list[i].dist(seq_list[j])
    return out_matrix

def get_dists_for(seq, seq_list, dist_matrix):
    seq_i = seq_list.index(seq)
    dists = dist_matrix[...,seq_i]
    return dists

def get_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1):
    output_seq_list = SeqList()
    dists = get_dists_for(seq, seq_list, dist_matrix)
    bool_array = array(dists == dist)[...,0]
    dist_x_indices = arange(len(seq_list))[bool_array]
    for i in dist_x_indices:
        output_seq_list += seq_list[i]
    return output_seq_list

def get_num_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1):
    dists = get_dists_for(seq, seq_list, dist_matrix)
    bool_array = array(dists == dist)[...,0]
    dist_x_indices = arange(len(seq_list))[bool_array]
    num_dist_x_neighbors = len(dist_x_indices)
    return num_dist_x_neighbors

def get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1):
    output_abund = 0
    dists = get_dists_for(seq, seq_list, dist_matrix)
    bool_array = array(dists == dist)[...,0]
    dist_x_indices = arange(len(seq_list))[bool_array]
    for i in dist_x_indices:
        output_abund += float(seq_list[i])
    return output_abund

def centrality_index(seq, **kwargs):
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    dists = get_dists_for(seq, seq_list, dist_matrix)
    some_index = 0
    for dist in dists:
        some_index += 1.0 / (dist + 1)
    return float(some_index)