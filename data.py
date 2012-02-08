"""Package which defines the Seq, Seqs classes as well as sequence I/O

"""
import random

DEFAULT_ALPHABET = ['A', 'C', 'G', 'T']
DEFAULT_ABUND = 1

ALPHABET = DEFAULT_ALPHABET

class Seq():
    """Sequence and its abundance.
    
    Also defines methods for comparing sequences
    
    """
    def __init__(self, seq, abund):
        self.str = seq.strip().upper()
        self.abund = float(abund)
        
    def __str__(self):
        return self.str
    
    def __repr__(self):
        return ("Seq(seq='%s', abund=%f)" % (str(self), float(self)))
        
    def __eq__(self, other):
        if self.str == str(other):
            return True
        else:
            return False
        
    def __iadd__(self, value):
        self.abund += float(value)
        return self
        
    def __float__(self):
        return float(self.abund)

class SeqList():
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seqs = None):
        # the dict attribute is a dictionary with the sequence
        # string as the key and the Seq object as the value.
        self.dict = {}
        if seqs is not None:
            for seq in seqs:
                # add the sequence to the dictionary
                if seq in self:
                    self[seq] += seq.abund
                else:
                    self[seq] = seq
            
    def __iter__(self):
        return iter(self.dict.values())    
    
    def __repr__(self):
        return "SeqList(%s)" % repr(self.dict.values())
    
    def __contains__(self, seq):
        """Check if the sequence or Seq object are in SeqList."""
        if str(seq) in self.dict.keys():
            return True
        elif seq in self.dict.values():
            return True
        else:
            return False
        
    def __getitem__(self, key):
        """Return the Seq object which has the same sequence as key
        
        """
        return self.dict[str(key)]
    
    def __len__(self):
        return len(self.dict.keys())
    
    def __setitem__(self, key, value):
        """Add/Replace the item at str(key) with value."""
        self.dict[str(key)] = value
        return self
        
    def __iadd__(self, value_or_list):
        """Add the Seq object to SeqList[str(seq)]"""
        try:
            items = list(value_or_list)
        except TypeError:
            items = [value_or_list]
        for item in items:
            if item in self:
                self[item] += item
            else:
                self[item] = item
        return self
    
    def append(self, seq):
        """See __iadd__"""
        self += seq
        return self
        
    def norm_abunds(self):
        """Normalize the abundances to a total of 1.0"""
        abund_total = self.get_abund_total()
        for seq in self:
            seq.abund = seq.abund/abund_total
            
    def get_abund_total(self):
        running_total = 0.0
        for seq in self:
            running_total += float(seq)
        return running_total
            
#    def dereplicate(self):
#        """Remove duplicate sequences and combine their abundances.
#        
#        """
#        unique = SeqList()
#        for seq in self:
#            if seq in unique:
#                unique[seq] += seq.abund
#            else:
#                unique[seq] = seq
#        self.dict = unique.dict
        
    def random_seq(self):
        abund_total = self.get_abund_total()
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
                seqs.append(Seq(curr_seq, curr_abund))
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
    seqs.append(Seq(curr_seq, curr_abund))
    return seqs