"""

"""
from scipy import mat, zeros
from data import SeqList, Seq
from math import log

def print_index_tuples(index_tuples):
    for item in index_tuples:
        print("%f\t%s" % (item[0], str(item[1])))
        
def calc_dist_matrix(seq_list):
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

def dists_for(seq, seq_list, dist_matrix):
    return dist_matrix[...,seq_list.index(seq)]


def use_index(function, seq_list, **kwargs):
    """Returns a list of sequences in seq_list sorted by index value.
    
    """
    try:
        dists = kwargs['dist_matrix']
    except KeyError:
        dists = calc_dist_matrix(seq_list)
    index_list = []
    kwargs = {}
    kwargs['dist_matrix'] = dists
    kwargs['seq_list'] = seq_list
    for seq in seq_list:
        index_list.append(function(seq, **kwargs))
    return sorted(zip(index_list, seq_list))

def abund_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The abund_index of a sequence is its relative abundance.

    """
    seq_list = kwargs['seq_list']
    return seq_list.norm()[seq].abund

def log_abund_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The log_abund_index of a sequence is the log of its normalized abundance.

    """
    seq_list = kwargs['seq_list']
    return log(float(seq_list.norm()[seq].abund))

def neighbor_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The neighbor_index of a sequence is the number of neighboring
    (only 1 nucleotides off) sequences present in the sample.

    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    return num_neighbors(seq, seq_list, dist_matrix, dist = 1)

def weighted_neighbor_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The weighted_neighbor_index of a sequence is the total abundance
    of neighboring sequences.

    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    return abund_neighbors(seq, seq_list, dist_matrix, dist = 1)

def inverse_weighted_neighbor_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The inverse_weighted_neighbor_index of a sequence is the
    abundance of the sequence divided by the total abundance of
    neighboring sequences.

    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    return float(seq_list[seq].abund) / abund_neighbors(seq, seq_list, dist_matrix, dist = 1)

def inverse_neighbor_index(seq, **kwargs):
    """Returns an index value for seq.
    
    The inverse_neighbor_index of a sequence is its abundance
    divided by the number of neighbors.
    
    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    return float(seq_list[seq].abund) / num_neighbors(seq, seq_list, dist_matrix, dist = 1)

def abund_times_neighbors_index(seq, **kwargs):
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    return float(seq_list[seq].abund) * num_neighbors(seq, seq_list, dist_matrix, dist = 1)

def index0(seq, **kwargs):
    """Returns the index0 of the sequence.
    
    The index0 is the abundance of the sequence minus the abundance
    of all neighbors divided by their distance from the
    sequence plus 1.
    
    """
    seq_list = kwargs['seq_list'].norm()
    dist_matrix = kwargs['dist_matrix']
    dists = dists_for(seq, seq_list, dist_matrix)
    abund = float(seq_list[seq].abund)
    index0_curr = abund
    for i in range(len(seq_list)):
        other_seq = seq_list[i]
        if other_seq is not seq:
            index0_curr -= other_seq.abund / (dists[i] + 1.0)
    return index0_curr
    