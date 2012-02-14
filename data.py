"""Package which defines the Seq, Seqs classes as well as sequence I/O

"""
import random
from scipy import array
from scipy.spatial.distance import hamming


DEFAULT_ALPHABET = ['A', 'C', 'G', 'T']
DEFAULT_ABUND = 1.0

ALPHABET = DEFAULT_ALPHABET

class Seq():
    """Sequence and its abundance.
    
    Also defines methods for comparing sequences
    
    """
    def __init__(self, seq_str, abund, attrs = None):
        self._seq = array(list(seq_str.strip().upper()))
        self._abund = float(abund)
        if attrs is None:
            self._attrs = {}
                    
    def set_attr(self, attr, value):
        self._attrs[attr] = value
        
    def get_attr(self, attr):
        return self._attrs[attr]
    
    def get_attrs(self):
        return self._attrs
    
    def del_attr(self, attr):
        del self._attrs[attr]
        
    def get_seq(self):
        return self._seq

    def __str__(self):
        seq_str = ""
        for nucl in self._seq:
            seq_str += str(nucl)
        return seq_str

    def __len__(self):
        return len(self._seq)

    def __lt__(self, other):
        return float(self) < float(other)

    def dist(self, other):
        """Returns the hamming distance of self vs. other)"""
        assert len(self) == len(other)
        return hamming(self._seq, other._seq) * len(self)
    
    def __repr__(self):
        return ("Seq(seq='%s', abund=%f, attrs=%s)" % (str(self), float(self), str(self._attrs)))
        
    def __eq__(self, other):
        return str(self) == str(other)
        
    def __iadd__(self, value):
        self._abund += float(value)
        return self

    def __add__(self, value):
        return Seq(str(self), self._abund + float(value))
        
    def __float__(self):
        return float(self._abund)


class SeqList():
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seqs = None):
        # the dict attribute is a dictionary with the sequence
        # string as the key and the Seq object as the value.
        self._seq_dict = {}
        self._seq_list = []
        if seqs is not None:
            for seq in seqs:
                # add the sequence to the dictionary
                if seq in self:
                    self[seq] += seq._abund
                else:
                    try:
                        abund = float(seq)
                    except TypeError:
                        abund = DEFAULT_ABUND
                    try:
                        attrs = seq.get_attrs()
                    except TypeError:
                        attrs = None
                    self[seq] = Seq(str(seq), abund, attrs)
            
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]
            
    def __reversed__(self):
        for i in range(len(self)-1, -1, -1):
            yield self[i]
            
    def sort_by(self, attr, reverse = False):
        working_seq_list = []
        working_attr_list = []
        sorted_seq_list = []
        for seq in self:
            working_seq_list += [str(seq)]
            if attr == "abund":
                working_attr_list += [float(seq)]
            elif attr == "seq":
                working_attr_list += [str(seq)]
            else:
                try:
                    working_attr_list += [seq.get_attr(attr)]
                except KeyError:
                    working_attr_list += None
        sorted_zipper = sorted(zip(working_attr_list, working_seq_list))
        for tup in sorted_zipper:
            sorted_seq_list += [str(tup[1])]
        if reverse is not True:
            self._seq_list = sorted_seq_list
        else:
            self._seq_list = list(reversed(sorted_seq_list))
        
    
    def seq_strs(self):
        return self._seq_list

    def __repr__(self):
        seqs_str = ""
        for seq in self:
            seqs_str += "%s\n" % repr(seq)
        return "SeqList(%s)" % seqs_str
    
    def __str__(self):
        seqs_str = ""
        for seq in self:
            seqs_str += "%s\t%f\n" % (str(seq), float(seq))
        return "SeqList\n(\n%s)" % seqs_str
    
    def __contains__(self, seq):
        """Check if the sequence or Seq object are in SeqList."""
        if str(seq) in self._seq_list:
            return True
        else:
            return False
        
    def __getitem__(self, key):
        """Return the Seq object which has the same sequence as key
        
        """
        try:
            return self._seq_dict[self._seq_list[key]]
        except TypeError:
            return self._seq_dict[str(key)]
    
    def __len__(self):
        return len(self._seq_dict)
    
    def __setitem__(self, key, value):
        """Add/Replace the item at str(key) with value."""
        self._seq_dict[str(key)] = value
        if not value in self._seq_list:
            self._seq_list.append(str(value))
        return self
        
    def __iadd__(self, seq_object):
        """Add a single Seq object to self."""
        if seq_object in self:
            self[seq_object] = Seq(str(seq_object),
                                   self[seq_object]._abund +
                                       float(seq_object))
        else:
            self[seq_object] = seq_object
        return self
    
    def index(self, value):
        assert value in self
        return self._seq_list.index(str(value))
    
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
        """Return a random Seq object from self.
        
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
