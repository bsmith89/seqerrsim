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
