"""

"""
from scipy import mat, zeros

def dist_matrix(seq_list):
    dim = len(seq_list)
    out_matrix = mat(zeros((dim, dim)))
    for i in range(dim):
        for j in range(dim):
            out_matrix[i,j] = seq_list[i].dist(seq_list[j])
    return out_matrix

def abund_index(seq_list):
    """Returns a list of tuples, (Seq, index_value).
    
    The abund_index of a sequence is its abundance.

    """
    pass

def log_abund_index(seq_list):
    """Returns a list of tuples, (Seq, index_value).
    
    The log_abund_index of a sequence is the log of its abundance.

    """
    
