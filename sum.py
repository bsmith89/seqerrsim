import sys
import optparse

def list_seqs(seqs):
    return count(seqs).keys()

def diffs_matrix(seqs, max_dist=None):
    out = {}
    seqs = count(seqs).keys()
    num_done = 0
    for seq1 in seqs:
        counts = {}
        for seq2 in seqs:   
            dist = None
            if seq2 in out and seq1 in out[seq2]:
                dist = out[seq2][seq1]
            else:
                dist = count_diffs(seq1, seq2)
            if max_dist != None and dist > max_dist:
                continue
            else:
                counts[seq2] = dist
        out[seq1] = counts
        num_done += 1
        print("%s done %d" % (seq1, num_done))
    return out

def calc_centr_index(diffs_matrix, counts):
    """Calculate the centrality index for each seq in diffs_matrix

    The centrality index of seq1 is the total sum of seqs found
    weighted by 1 over their distance from seq1.

    Unfortunately this doesn't seem to be the solution for a sample
    of size 10000.  I found that a 1% sequence had a centrality index
    of only ~70, compared to ~1400 for a 99% abundance sequence.
    To me, this means that the index is not independent of abundance,
    but I would still expect it to become increasingly independent as
    sample size -> inf.  This index should really be re-labeled,
    because it is not really an index of centrality, although I wanted
    it to be one.

    A centrality index should be independent of sample size, and
    reflect how closely "neighboring" sequences are distributed to
    the perfectly binomial distribution predicted.

    """
    out = {}
    for seq1 in diffs_matrix.keys():
        centr_index = 0
        for seq2 in diffs_matrix[seq1]:
            if seq2 == seq1:
                continue
            else:
                centr_index += (1.0 / diffs_matrix[seq1][seq2])
        out[seq1] = centr_index
    return out

def count_diffs(seq1, seq2):
    count = 0
    assert len(seq1) == len(seq2)
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count

def count(seqs):

    out = {}
    for seq in seqs:
        if seq in out:
            out[seq] += 1
        else:
            out[seq] = 1
    return out

def rabund(seqs, size = None):
    out = {}
    seq_counts = count(seqs)
    if size == None:
        size = sum(seq_counts.values())
    for seq in count(seq_counts):
        out[seq] = float(seq_counts[seq])/size
    return out

if __name__ == '__main__':
    p=optparse.OptionParser()
    p.add_option('--fasta', '-f')
    p.add_option('--output', '-o')
    p.add_option('--sample-size', '-n', default = '1000')
    p.add_option('--error-rate', '-e', default = '0.01')
    opts, args = p.parse_args()
    if opts.fasta:
        infile = open(opts.fasta)
    else:
        infile = sys.stdin
    if opts.out:
        outfile = open(opts.output)
    else:
        outfile = sys.stdout
    error_rate = float(opts.error_rate)
    sample_size = int(sample_size)


