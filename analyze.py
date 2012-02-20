"""Data structures and functions for analyzing a SeqList object.

"""
import scipy
from scipy.misc import comb

class DistMatrix(object):
    """A memory/CPU efficient way to retrieve hamming distances.
    
    Distance values can not be set to 0.0 for any two sequences that
    are not identical, and then they are by definition 0.0.
    
    See Also:
    data.Seq.dist
    
    """
    def __init__(self, seq_list):
        n = len(seq_list)
        self._matrix = scipy.zeros([n, n], dtype = float)
        self._seq_order = list(seq_list)
        self.seq_list = seq_list
        
    def update(self, seq_list):
        """Update self to work for a new seq list.
        
        If possible conserves the internal sequence order so that it
        doesn't have to recalculate already calculated values.
        
        """
        n = len(seq_list)
        if len(self._seq_order) is not n:
            self.__init__(seq_list)
            return
        for seq in seq_list:
            if seq not in self:
                self.__init__(seq_list)
                return
        self.seq_list = seq_list
                
    def __contains__(self, seq):
        """
        
        >> for i in len(seq_list):
        >>     assert seq_list[i] in DistMatrix(seq_list) == True
        
        """
        return seq in self._seq_order
        
    def _set(self, seq1, seq2, value):
        """Internal method for setting a distance value.
        
        self._matrix[i, j] == self._matrix[j, i] == value
        
        """
        i = self._seq_order.index(seq1)
        j = self._seq_order.index(seq2)
        assert value is not 0.0
        self._matrix[i, j] = value
        self._matrix[j, i] = value
        
    def _get(self, seq1, seq2):
        """Return the distance between seq1 and seq2 if stored.
        
        Internal method for getting a previously stored distance.
        
        Raises a KeyError if dist(seq1, seq2) has not been set.
        
        """
        i = self._seq_order.index(seq1)
        j = self._seq_order.index(seq2)
        if seq1 == seq2:
            return 0.0
        value = self._matrix[i, j]
        if value == 0.0:
            raise KeyError('That distance has not been set or equals 0.0 (illegal).')
        return self._matrix[i, j]
        
    def get(self, seq1, seq2):
        """Return the distance between seq1 and seq2
        
        If the value has already been stored, retrieves it.
        Otherwise the value is calculated afresh.
        
        """
        try:
            assert seq1 in self
        except AssertionError:
            raise KeyError('seq1 was not in self\nseq1: %s' % str(seq1))
        try:
            assert seq2 in self
        except AssertionError:
            raise KeyError('seq2 was not in self\nseq2: %s' % str(seq2))
        try:
            value = self._get(seq1, seq2)
        except KeyError:
            # try to calculate the distance from seq1
            try:
                value = seq1.dist(seq2)
            except AttributeError: 
                # but if that fails try seq2
                # if this raises a AttributeError let it pass
                value = seq2.dist(seq1)
            self._set(seq1, seq2, value)
        return self._get(seq1, seq2)

    def __getitem__(self, seq1):
        """Return a DistArray for seq1.
        
        for i in len(seq_list):
            for j in len(seq_list):
                DistMatrix(seq_list)[seq_list[i]][seq_list[j]] == \
                seq_list[i].dist(seq_list[j])
        
        """
        return DistArray(seq1, self)
    
    def __iter__(self):
        """Iterate over DistArrays for every seq.
        
        """
        for seq in self.seq_list:
            yield self[seq]
            
        
class DistArray(object):
    """A reference to a particular row of DistMatrix
    
    """
    def __init__(self, seq, dist_matrix):
        self._matrix = dist_matrix
        self._seq = seq
        
    def __getitem__(self, other_seq):
        """Return the distance between self.seq and other_seq
        
        See Also:
        DistMatrix.get(self._seq, other_seq)
        
        """
        return self._matrix.get(self._seq, other_seq)
    
    def __iter__(self):
        """Iterate over the sequences being compared.
        
        iter(dist_mat) should be identical to
        iter(dist_mat[i]) for all i in [0, len(dist_mat))
        
        """
        return iter(self._matrix.seq_list)
        
        
def calc_all(index_function, **kwargs):
    """Calculate an index_function For every sequence in seq_list.
    
    Stores the result of calling index_function(seq, **kwargs) for seq
    as an attribute of seq with name index_function.__name__
    
    seq_list and any additional arguments required by index_function
    must be provided as key-word arguments.
    
    """
    seq_list = kwargs['seq_list']
    for seq in seq_list:
        seq.__setattr__(index_function.__name__, index_function(seq, 
                                                                **kwargs))
        
def binomial_corrected_abundance(seq, **kwargs):
    """Calculate the abundance of seq less the amount expected by error
    
    Based on the expected missequencing of other_seq as seq (A_ij)
    by the rule:
    A_ij = A_j * comb(L, d) * (e / 3)^d * (1 - e)^(L-d)
    
    where:
        A_j = the true abundance of other_seq
        L = length of seq = length other_seq
        d = nucleotide differences between seq and other_seq
        e = per nucleotide error-rate
        
    Originally meant to be applied iteratively to determine which
    sequences were at too high an abundance to be artifacts of
    sequencing error.
    
    """
    working_abund = seq.abund
    dist_mat = kwargs['dist_mat']
    e = kwargs['error_rate']
    l = len(seq)
    dist_array = dist_mat[seq]
    for other_seq in dist_array:
        a = other_seq.abund
        d = dist_array[other_seq] * l
        working_abund -= a * (1 / comb(l, d)) * ((e / 3.0) ** d) * \
                             ((1.0 - e) ** (l - d))
    return working_abund
        