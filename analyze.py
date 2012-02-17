"""TODO: docstring

"""
from scipy import mat, zeros, arange, array
from data import SeqList

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

def local_maximum(seq, **kwargs):
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    dist_x_neighbors = get_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1)
    for neighbor in dist_x_neighbors:
        if neighbor > seq:
            return False
    return True

def dist_closest_greater(seq, **kwargs):
    """Returns closest distance of a sequence with higher abundance.
    
    TODO: Fix this.  Plus it's *really slow.
    
    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    dist = 1
    while True:
        dist_x_neighbors = get_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = dist)
        for neighbor in dist_x_neighbors:
            if neighbor > seq:
                return dist
        dist += 1
    return None

def adjusted_abundance3(seq, **kwargs):
    """Abundance of the sequence minus the abundance of 1-neighbors.
    
    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    abund_dist_1_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1)
    abund_dist_2_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 2)
    abund_dist_3_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 3)
    return float(seq) - abund_dist_1_neighbors / 2.0**1 - abund_dist_2_neighbors / 2.0**2 - abund_dist_3_neighbors / 2.0**3

def adjusted_abundance2(seq, **kwargs):
    """Abundance of the sequence minus the abundance of 1-neighbors.
    
    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    abund_dist_1_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1)
    abund_dist_2_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 2)
    return float(seq) - abund_dist_1_neighbors / 2.0**1 - abund_dist_2_neighbors / 2.0**2

def adjusted_abundance1(seq, **kwargs):
    """Abundance of the sequence minus the abundance of 1-neighbors.
    
    """
    seq_list = kwargs['seq_list']
    dist_matrix = kwargs['dist_matrix']
    abund_dist_1_neighbors = get_abund_dist_x_neighbors_for(seq, seq_list, dist_matrix, dist = 1)
    return float(seq) - abund_dist_1_neighbors / 2.0**1


class SeqCompMatrix(object):
    def __init__(self, seq_list):
        self._matrix_dict = {}
        for seq1 in seq_list:
            self._matrix_dict[seq1] = {}
            for seq2 in seq_list:
                self._matrix_dict[seq1][seq2] = {}
                
    def contains(self, seq):
        return seq in self._matrix_dict
        
    def _set_comp(self, seq1, seq2, key, value):
        self._matrix_dict[seq1][seq2][key] = value
        
    def _get_comp(self, seq1, seq2, key):
        return self._matrix_dict[seq1][seq2][key]
    
    def _del_comp(self, seq1, seq2, key):
        del self._matrix_dict[seq1][seq2][key]
        
    def setc(self, seq1, seq2, key, value):
        assert self.contains(seq1)
        assert self.contains(seq2)
        self._set_comp(seq1, seq2, key, value)
        
    def getc(self, seq1, seq2, key):
        assert self.contains(seq1)
        assert self.contains(seq2)
        try:
            value = self._get_comp(seq1, seq2, key)
        except KeyError: # try to calculate the comparison from a seq1 method
            try:
                value = seq1.__getattr__(key)(seq2)
            except AttributeError:
                value = seq2.__getattr__(key)(seq1)
            self.set(seq1, seq2, key, value)
            self.set(seq2, seq1, key, value)
        return self._get_comp(seq1, seq2, key)
        
    def delc(self, seq1, seq2, key):
        assert self.contains(seq1)
        assert self.contains(seq2)
        self._del_comp(seq1, seq2, key)
        
        
class DistMatrix(SeqCompMatrix):
    
    def get(self, seq1, seq2):
        return self.getc(seq1, seq2, key = "dist")
        
        
        
        
        
        
        
        
        
        
        