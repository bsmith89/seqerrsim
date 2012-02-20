"""Simulation functions for sampling and sequencing errors.

"""
import random
from data import Seq, SeqList, ALPHABET

def sample(seq_list, n):
    """Return a random sample from the provided SeqList.
    
    Randomly selects n sequences from seq_list with a uniform
    distribution weighted by relative seq.abund values.
    
    """
    sample = SeqList()
    for i in range(n):
        sample += seq_list.random_seq()
    return sample
        
def sample_with_errors(seq_list, n, error_rate, alph = None):
    """Return a random sample with errors from seq_list.
    
    Returned values have random errors at a per-nucleotide error rate
    of error_rate.
    
    alph is a list object of characters that can be found or
    erroneously sequenced.  Default value can be found in data.ALPHABET
    but is most likely set to ['A', 'C', 'G', 'T']
    (data.DEFAULT_ALPHABET).
    
    See Also:
    simulate.sample
    
    """
    if alph is None:
        alph = ALPHABET
    sample = SeqList()
    for i in range(n):
        pick = seq_list.random_seq()
        pick_str_with_errors = with_errors(str(pick), error_rate, alph)
        sample += Seq(pick_str_with_errors, abund = 1.0)
    return sample
        
            
def with_errors(seq_str, error_rate, alph = ALPHABET):
    """Returns a sequence based on seq_str with errors.
    
    See Also:
    simulate.sample_with_errors
    
    """
    out_str = ""
    for index in range(len(seq_str)):
        nucl = seq_str[index]
        seqd_nucl = nucl
        if random.random() < error_rate:
            while seqd_nucl == nucl:
                seqd_nucl = random.choice(alph)
        out_str += seqd_nucl
    return out_str
    
