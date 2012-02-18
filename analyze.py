"""TODO: docstring

"""
from scipy import mat, zeros, arange, array
from data import SeqList


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
        except KeyError:
            # try to calculate the comparison from a seq1 method
            try:
                value = seq1.__getattribute__(key)(seq2)
            except AttributeError: 
                # but if that fails try a seq2 method instead
                value = seq2.__getattribute__(key)(seq1)
            self._set_comp(seq1, seq2, key, value)
            self._set_comp(seq2, seq1, key, value)
        return self._get_comp(seq1, seq2, key)
        
    def delc(self, seq1, seq2, key):
        assert self.contains(seq1)
        assert self.contains(seq2)
        self._del_comp(seq1, seq2, key)
        
    def seq_comp_array(self, seq1, key):
        """
        
        TODO: Make this less of a bottleneck
        
        """
        out_dict = {}
        for seq2 in self._matrix_dict[seq1]:
            out_dict[seq2] = self.getc(seq1, seq2, key)
        return out_dict
    
    def seqs_for_which_true(self, seq1, key, test, **testargs):
        out_list = []
        seq_comp_array = self.seq_comp_array(seq1, key = key)
        for seq2 in seq_comp_array:
            if test(seq_comp_array[seq2], **testargs):
                out_list += [seq2]
        return out_list
        
        
class DistMatrix(SeqCompMatrix):
    
    def get(self, seq1, seq2):
        return self.getc(seq1, seq2, key = "dist")
        
        
def calc_all(index_function, **kwargs):
    seq_list = kwargs['seq_list']
    for seq in seq_list:
        seq.__setattr__(index_function.__name__, index_function(seq, **kwargs))
        
def binomial_corrected_abundance(seq, **kwargs):
    working_abund = seq.abund
    dist_mat = kwargs['dist_mat']
    e = kwargs['error_rate']
    l = len(seq)
    dist_array = dist_mat.seq_comp_array(seq, key = 'dist')
    for other_seq in dist_array:
        a = other_seq.abund
        d = dist_array[other_seq]
        working_abund -= a * ((e / 3) ** d) * ((1 - e) ** (l - d))
    return working_abund


        
        
        
        
        
        
        
        