# !usr/bin/env python
"""Simple sequence data simulation

This module contains methods for parsing fasta files of the correct
format and outputting a simulated sequence data set.

"""
import sys
import optparse
import random

DEFAULT_ALPHABET = ['A','C','G','T']

def parse_fa(fa_file, abund_data=False):
    """Take a FASTA file and return a list of seqs.

    If abund_data == True (default = false) The list contains each
    seq proportionally to its abundance.  This abundance data is read
    from the FASTA file as the last number on the title line.  If the
    last word on the title line is not a integer abundance is assumed
    to be 1.

    fa_file must be an open file object.

    """
    out_seqs = []
    curr_seq = ""
    curr_abund = 0
    for line in fa_file:
        # If this is the title line, then add current seq to out_seqs
        # and reset curr_seq and curr_abund.
        if line.startswith('>'):
            if curr_seq == "": # if this is the very first title line
                continue
            seqs += [curr_seq]*curr_abund
            curr_seq = ""
            # Sequence abundances should be an integer at the end of
            # the title line.
            if len(line.split()) > 1:
                try:
                    curr_abund = int(line.split()[-1])
                except ValueError:
                    curr_abund = 1                
        else:
            # Add the current line of sequence, removing white space
            # and setting all letters uppercase.
            curr_seq.join(line.strip().upper())
    # Add the final sequence and return output
    if curr_seq != "": # in case the file was empty
        seqs += [curr_seq]*curr_abund
    return out_seqs

def sample_from_seqs(n, seqs_list, rand_seed = None):
    """Choose with replacement n sequences from seqs_list.

    Return a list of n sequences chosen randomly from seqs_list.  If
    rand_seed is set, it is used as the seed value for sequence
    picking.
    
    """
    out_seqs = []
    for i in range(n):
        seq = random.choice(seqs_list)
        out_seqs += [seq]
    return out_seqs

def simulate_errors(e, seqs_list, alphabet = DEFAULT_ALPHABET):
    """Return a list of seqs where errors have been introduced.

    Errors are introduced at a per nucleotide error rate of e.  All
    erroneous nucleotides are picked from the character list alphabet
    with the true nucleotide removed.
    
    """
    # The maximum integer to be randomly generated should be no
    # higher than error-rate^(-1) so that an error rate of 0.01
    # means that integers should be randomly picked up to 100.
    range_len = int(1.0/e)
    out_seqs = []
    for true_seq in seqs:
        observed_seq = ""
        for true_nucl in true_seq:
            observed_nucl = None
            # This assertion shouldn't really be necessary
            assert true_nucl in alphabet
            # if the random integer (from 0 to range_len-1) is 0
            if random.randrange(range_len) == 0:
                obs_nucl = random.choice(alphabet)
                # The erronious nucleotide should never be the true
                while obs_nucl == nucl:
                    obs_nucl = random.choice(alphabet)
            else:
                obs_nucl = true_nucl
            obs_seq.join(obs_nucl)
        out_seqs += [obs_seq]
    return out_seqs

def quick_data(file_path, n=10000, e=0.01):
    return simulate_errors(e, sample_from_seqs(n, parse_fa(open(file_path), abund_data = True)))

if __name__ == '__main__':
    p=optparse.OptionParser()
    p.add_option('--fasta', '-f')
    p.add_option('--output', '-o')
    p.add_option('--sample-size', '-n', default = '10000')
    p.add_option('--error-rate', '-e', default = '0.01')
    opts, args = p.parse_args()
    if opts.fasta:
        infile = open(opts.fasta)
    else:
        infile = sys.stdin
    if opts.output:
        outfile = open(opts.output)
    else:
        outfile = sys.stdout
    error_rate = float(opts.error_rate)
    sample_size = int(opts.sample_size)
    seq_data = simulate_errors(error_rate, sample_from_seqs(sample_size, parse_fa(infile)))
    infile.close()
    index = 0
    for seq in seq_data:
        seq_name = "seq%d" % index
        outfile.write(">%s\n%s\n" % (seq_name, seq))
        index += 1
    outfile.close()
