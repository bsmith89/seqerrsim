"""

"""
from scipy import mat, zeros
from data import SeqList, Seq
from math import log

def dist_matrix(seq_list):
    dim = len(seq_list)
    out_matrix = mat(zeros((dim, dim)))
    seq_str_list = []
    for i in range(dim):
        seq_str_list += [str(seq_list[i])]
        for j in range(dim):
            out_matrix[i,j] = seq_list[i].dist(seq_list[j])
    return out_matrix

def neighbors(seq, seq_list, dist_matrix, dist = 1):
    """Returns a SeqList object of all sequences neighboring seq."""
    assert seq in seq_list
    neighbors = SeqList()
    seq_i = seq_list.index(seq)
    # The distance of every sequence compared to seq
    comparisons = dist_matrix[...,seq_i]
    assert len(comparisons) == len(seq_list)
    for j in range(len(comparisons)):
        if comparisons[j] == dist:
            neighbors += seq_list[j]
    return neighbors
    

def num_neighbors(seq, seq_list, dist_matrix, dist = 1):
    return len(neighbors(seq, seq_list, dist_matrix, dist))

def abund_neighbors(seq, seq_list, dist_matrix, dist = 1):
    return neighbors(seq, seq_list, dist_matrix, dist).total_abund()
    

def abund_index(seq, seq_list):
    """Returns an index value for seq.
    
    The abund_index of a sequence is its relative abundance.

    """
    return seq_list.norm[seq].abund

def log_abund_index(seq, seq_list):
    """Returns an index value for seq.
    
    The log_abund_index of a sequence is the log of its normalized abundance.

    """
    return log(seq_list.norm[seq].abund)

def neighbor_index(seq, seq_list, dist_matrix):
    """Returns an index value for seq.
    
    The neighbor_index of a sequence is the number of neighboring
    (only 1 nucleotides off) sequences present in the sample.

    """
    return num_neighbors(seq, seq_list, dist_matrix, dist = 1)

def weighted_neighbor_index(seq, seq_list, dist_matrix):
    """Returns an index value for seq.
    
    The weighted_neighbor_index of a sequence is the total abundance
    of neighboring sequences.

    """
    return abund_neighbors(seq, seq_list, dist_matrix, dist = 1)

def inverse_weighted_neighbor_index(seq, seq_list, dist_matrix):
    """Returns an index value for seq.
    
    The inverse_weighted_neighbor_index of a sequence is the
    abundance of the sequence divided by the total abundance of
    neighboring sequences.

    """
    return seq_list[seq].abund / abund_neighbors(seq, seq_list, dist_matrix, dist = 1)

def inverse_neighbor_index(seq, seq_list, dist_matrix):
    """Returns an index value for seq.
    
    The inverse_neighbor_index of a sequence is its abundance
    divided by the number of neighbors.
    
    """
    return seq_list[seq].abund / num_neighbors(seq, seq_list, dist_matrix, dist = 1)