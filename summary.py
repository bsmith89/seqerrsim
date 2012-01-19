import sys
import optparse

def list_seqs(seqs):
    """Returns a sorted list of all sequences in seqs.

    Somewhat computationally expensive, so if you'll ever need
    abundance data, you should probably use count(seqs) and cache
    the resulting sequence list with its abundance data.

    """
    return sorted(count(seqs).keys())

def calc_diffs_matrix(seqs):
    """Calculate the distances between all seqs.
    
    Returns a dictionary of dictionaries (or matrix) where
    dict[seq1][seq2] is the number of nucleotide differences b/t
    seq1 and seq2.
    
    """
    out_matrix = {}
    num_done = 0
    for seq1 in seqs:
        seq1_vers = {}
        for seq2 in seqs:   
            if seq2 in out_matrix and seq1 in out_matrix[seq2]:
                seq1_vers[seq2] = out_matrix[seq2][seq1]
            else:
                seq1_vers[seq2] = nucl_diffs(seq1, seq2)
        out_matrix[seq1] = seq1_vers
        num_done += 1
    return out_matrix

def nucl_diffs(seq1, seq2):
    """Return the number of loci with different nucleotides in seqs.

    seq1 and seq2 must be the same length.
    """
    count = 0
    assert len(seq1) == len(seq2)
    if seq1 == seq2:
        return 0
    for char1, char2 in zip(seq1, seq2):
        if char1 != char2:
            count += 1
    return count

def hamming_dist(seq1, seq2):
    return nucl_diffs(seq1, seq2)

def abund(seqs):
    """Returns the abundance of each seq in a list of seqs."""
    out_dict = {}
    for seq in seqs:
        if seq in out_dict:
            out_dict[seq] += 1
        else:
            out_dict[seq] = 1
    return out_dict

def calc_index0(seqs, diffs_matrix = None, abund = None):
    """Calculate the index0 for each seq in seqs.

    This is just a prototype function, where the result is a dict
    where every key points to a zero.

    """
    if diffs_matrix is None:
        diffs_matrix = calc_diffs_matrix(seqs)
    if abund is None:
        abund = abund(seqs)
    out_dict = {}
    for seq1 in diffs_matrix.keys():
        out_dict[seq1] = 0
    return out_dict

def calc_index1(seqs, diffs_matrix = None, abund = None):
    """Calculate the index1 for each seq in seqs.

    Index1 = The total number of sequences within dist=1 of seq.
    """
    if diffs_matrix is None:
        diffs_matrix = calc_diffs_matrix(seqs)
    if abund is None:
        abund = abund(seqs)
    out_dict = {}
    for seq1 in diffs_matrix.keys():
        for seq2 in diffs_matrix[seq1].keys():
            if diffs_matrix[seq1][seq2] == 1:
                try:
                    out_dict[seq1] += abund[seq2]
                except KeyError:
                    out_dict[seq1] = abund[seq2]
    return out_dict

def calc_index2(seqs, diffs_matrix = None, abund = None):
    """Calculate the index2 for each seq in seqs.

    Index2 = The total number of sequences within dist<1 of seq
    """
    if diffs_matrix is None:
        diffs_matrix = calc_diffs_matrix(seqs)
    if abund is None:
        abund = abund(seqs)
    out_dict = {}
    for seq1 in diffs_matrix.keys():
        for seq2 in diffs_matrix[seq1].keys():
            if diffs_matrix[seq1][seq2] <= 1:
                out_dict[seq1] += abund[seq2]
    return out_dict

def calc_index3(seqs, diffs_matrix = None, abund = None):
    """Calculate the index3 for each seq in seqs.

    Index3 = sum(sequences within dist=n of seq divided by n+1)
    """
    if diffs_matrix is None:
        diffs_matrix = calc_diffs_matrix(seqs)
    if abund is None:
        abund = abund(seqs)
    out_dict = {}
    for seq1 in diffs_matrix.keys():
        for seq2 in diffs_matrix[seq1].keys():
            dist = nucl_diffs(seq1, seq2)
            out_dict[seq1] += 1.0 / dist + 1
    return out_dict

if __name__ == '__main__':
    p=optparse.OptionParser()
    p.add_option('--fasta', '-f')
    p.add_option('--output', '-o')
    opts, args = p.parse_args()
    if opts.fasta:
        infile = open(opts.fasta)
    else:
        infile = sys.stdin
    if opts.out:
        outfile = open(opts.output)
    else:
        outfile = sys.stdout

    infile.close()
    outfile.close()
