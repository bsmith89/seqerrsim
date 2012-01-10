import sys
import optparse
import random

DEFAULT_ALPHABET = ['A','C','G','T']

def pars_fasta(fa_file, abund_data=False):
    seqs = []
    curr_seq = ""
    curr_abund = 0
    for line in fa_file:
        if line.startswith('>'):
            if curr_seq != "":
                seqs += [curr_seq]*curr_abund
            curr_seq = ""
            try:
                curr_abund = int(line.split()[-1]) # Seq abundances
            except:
                curr_abund = 1                 # should be an int
                                               # at the
                                               # end of the title
                                               # line.
        else:
            curr_seq += line.strip()
    if curr_seq != "":
        seqs += [curr_seq]*curr_abund
    return seqs

def sample(n, seqs):
    out = []
    for i in range(n):
        out += [random.choice(seqs).upper()]
    return out

def add_errors(e, seqs, alphabet = DEFAULT_ALPHABET):
    max = 1.0/e
    out = []
    for seq in seqs:
        obs_seq = ""
        for nucl in seq:
            obs_nucl = None
            assert nucl in alphabet
            if random.randrange(max) == 0:
                obs_nucl = random.choice(alphabet)
                while obs_nucl == nucl:
                    obs_nucl = random.choice(alphabet)
            else:
                obs_nucl = nucl
            obs_seq += obs_nucl
        out += [obs_seq]
    return out

def quick_data(filename, n=10000, e=0.01):
    return add_errors(e, sample(n, pars_fasta(open(filename))))

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


