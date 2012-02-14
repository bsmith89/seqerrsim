"""Simulates sampling and sequencing errors

"""
import random
from data import Seq, SeqList, ALPHABET, parse_file

def sample(seqs, n):
    sample = SeqList()
    for i in range(n):
        sample += seqs.random_seq()
    return sample
        
def sample_with_errors(seq_list, n, error_rate, alph = None):
    if alph is None:
        alph = ALPHABET
    sample = SeqList()
    for i in range(n):
        pick = seq_list.random_seq()
        pick_str_with_errors = with_errors(str(pick), error_rate, alph)
        sample += Seq(pick_str_with_errors, 1.0)
    return sample
        
            
def with_errors(seq_str, error_rate, alph = ALPHABET):
    out_str = ""
    for index in range(len(seq_str)):
        nucl = seq_str[index]
        seqd_nucl = nucl
        if random.random() < error_rate:
            while seqd_nucl == nucl:
                seqd_nucl = random.choice(alph)
        out_str += seqd_nucl
    return out_str
    
