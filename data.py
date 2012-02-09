"""Package which defines the Seq, Seqs classes as well as sequence I/O

"""
import random
from scipy import array
from scipy.spatial.distance import hamming


DEFAULT_ALPHABET = ['A', 'C', 'G', 'T']
DEFAULT_ABUND = 1

ALPHABET = DEFAULT_ALPHABET

class Seq():
    """Sequence and its abundance.
    
    Also defines methods for comparing sequences
    
    """
    def __init__(self, seq_str, abund):
        self.seq = array(list(seq_str.strip().upper()))
        self.abund = float(abund)

    def __str__(self):
        seq_str = ""
        for nucl in self.seq:
            seq_str += str(nucl)
        return seq_str

    def __iter__(self):
        return iter(str(self))

    def __len__(self):
        return len(self.seq)

    def __lt__(self, other):
        return str(self) < str(other)

    def dist(self, other):
        """Returns the hamming distance of self vs. other)"""
        assert len(self) == len(other)
        return hamming(self.seq, other.seq) * len(self)
    
    def __repr__(self):
        return ("Seq(seq='%s', abund=%f)" % (str(self), float(self)))
        
    def __eq__(self, other):
        if scipy.all(self.seq == other.seq):
            return True
        else:
            return False
        
    def __iadd__(self, value):
        self.abund += float(value)
        return self

    def __add__(self, value):
        return Seq(str(self), self.abund + float(value))
        
    def __float__(self):
        return float(self.abund)

class SeqList():
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seqs = None):
        # the dict attribute is a dictionary with the sequence
        # string as the key and the Seq object as the value.
        self.seq_dict = {}
        if seqs is not None:
            for seq in seqs:
                # add the sequence to the dictionary
                if seq in self:
                    self[seq] += seq.abund
                else:
                    self[seq] = seq
            
    def __iter__(self):
        return iter(self.seq_dict.values())    
    
    def seq_strs(self):
        return self.seq_dict.keys()

    def __repr__(self):
        return "SeqList(%s)" % repr(self.seq_dict.values())
    
    def __contains__(self, seq):
        """Check if the sequence or Seq object are in SeqList."""
        if str(seq) in self.seq_strs():
            return True
        else:
            return False
        
    def __getitem__(self, key):
        """Return the Seq object which has the same sequence as key
        
        """
        return self.seq_dict[str(key)]
    
    def __len__(self):
        return len(self.seq_dict)
    
    def __setitem__(self, key, value):
        """Add/Replace the item at str(key) with value."""
        self.seq_dict[str(key)] = value
        return self
        
    def __iadd__(self, seq_object):
        """Add a single Seq object to self."""
        if seq_object in self:
            self[seq_object] = Seq(str(seq_object),
                                   self[seq_object].abund +
                                       float(seq_object))
        else:
            self[seq_object] = seq_object
        return self
    
    def append(self, seq_objects):
        """Add a list of seqs to self."""
        for seq in seq_objects:
            self += seq
        return self

    def total_abund(self):
        running_total = 0.0
        for seq in self:
            running_total += float(seq)
        return running_total
        
    def norm(self):
        """Return a SeqList with normalized abundances, Total = 1.0"""
        out_seq_list = SeqList()
        abund_total = self.total_abund()
        for seq in self:
            out_seq_list += Seq(seq_str = str(seq), abund = float(seq) / abund_total)
        return out_seq_list
            
    def random_seq(self):
        """Return a random Seq object from self.seq_dict.
        
        This algorithm could probably be sped up a good deal but
        for now I don't need to.  Sampling doesn't seem to be the
        bottleneck.
        
        """
        abund_total = self.total_abund()
        index = random.random()
        adjusted_index = index * abund_total
        tally = 0.0
        for seq in self:
            tally += float(seq)
            if tally > adjusted_index:
                return Seq(str(seq), abund = 1)

def parse_file(fasta):
    """Takes a open file object and returns a Seqs object.
    
    """
    seqs = SeqList()
    curr_seq = ""
    curr_abund = 0.0
    for line in fasta:
        if line[0] == ">":
            # add the current sequence to seqs
            if curr_seq is not "":
                seqs += Seq(curr_seq, curr_abund)
            # reset curr_*
            curr_seq = ""
            curr_abund = 0.0
            # now, get the current abundance of this new sequence
            try:
                curr_abund = float(line.split()[-1])
            except (ValueError, IndexError):
                # There is no abundance on the title line
                curr_abund = DEFAULT_ABUND
        else:
            # append the sequence to curr_seq
            curr_seq += line.strip().upper()
    seqs += Seq(curr_seq, curr_abund)
    return seqs
