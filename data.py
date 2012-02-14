"""Package which defines the Seq, Seqs classes as well as sequence I/O

"""
import random
import scipy as sp
from scipy.spatial.distance import hamming


DEFAULT_ALPHABET = ['A', 'C', 'G', 'T']
DEFAULT_ABUND = 1.0

ALPHABET = DEFAULT_ALPHABET

class Seq():
    """Sequence and its abundance.
    
    Also defines methods for comparing sequences
    
    """
    def __init__(self, sequence, **kwargs):
        self._sequence = sp.array(list(sequence))
        for attr in kwargs.keys():
            self.__setattr__(attr, kwargs[attr])
    
    def __iter__(self):
        for nucl in self._sequence:
            yield nucl
    
    def __len__(self):
        return len(self._sequence)
    
    def __eq__(self, other):
        try:
            other_seq_array = other.get_seq_array()
        except AttributeError:
            other_seq_array = sp.array(list(other))
        return all(self.get_seq_array() == other_seq_array)
    
    def get_attrs(self):
        attr_dict = dict(self.__dict__)
        del attr_dict['_sequence']
        return attr_dict
            
    def get_seq_array(self):
        return self._sequence
    
    def __lt__(self, other):
        return str(self) < str(other)
    
    def __dist__(self, other):
        try:
            other_seq_array = other.get_seq_array()
        except AttributeError:
            other_seq_array = sp.array(list(other))
        return hamming(self.get_seq_array(), other_seq_array)
        
    def __str__(self):
        return ''.join(self.get_seq_array())
    
    def __repr__(self):
        other_attrs = self.get_attrs()
        out_str = "Seq('%s', **%s)" % (str(self), str(other_attrs))
        return out_str
    
    def get_abund(self):
        return self.abund
    
    def set_abund(self, value):
        self.abund = float(value)


class SeqList():
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seq_list = []):
        self._seq_list = []
        self.append(seq_list)
    
    def __iter__(self):
        return iter(self._seq_list)
    
    def __reversed__(self):
        return reversed(self._seq_list)
    
    def __repr__(self):
        seqs_str = ""
        for seq in self:
            seqs_str += "%s\n" % repr(seq)
        repr_str = "SeqList(\n%s)" % seqs_str
        return repr_str
    
    def __str__(self):
        header_str = ""
        body_str = ""
        curr_index = 0
        attr_lists = {}
        for seq in self:
            attrs_for_this_seq = seq.get_attrs()
            for attr in attrs_for_this_seq:
                if attr not in attr_lists:
                    attr_lists[attr] += ["NA"] * curr_index
            for attr in attr_lists:
                try:
                    value_for_this_seq = attrs_for_this_seq[attr]
                except KeyError:
                    value_for_this_seq = "NA"
                finally:
                    attr_lists[attr] += [value_for_this_seq]
            curr_index += 1
        for attr in attr_lists:
            header_str += "%s\t" % attr
        header_str += "sequence"
        attr_zip = zip(*attr_lists.values())
        for line in attr_zip:
            body_str += "\n%s" % "\t".join(line)
        return header_str + body_str
                
    def fasta(self):
        fasta_str = ""
        seq_count = 1
        for seq in self:
            title_line = ">seq%d" % seq_count
            attrs = seq.get_attrs()
            seq_str = str(seq)
            for attr in attrs:
                title_line += "|%s=%s" % (attr, str(attrs[attr]))
            fasta_str += "%s\n%s\n" % (title_line, seq_str)
        return fasta_str
    
    def __contains__(self, seq):
        if seq in self._seq_list:
            return True
        return False
    
    def __getitem__(self, key_or_index):
        try:
            index = self.index(key_or_index)
        except ValueError:
            try:
                assert (index % 1.0 == 0)
            except AssertionError:
                raise KeyError("No sequence %s found in this SeqList" % str(key_or_index))
            finally:
                index = key_or_index
        finally:
            return self._seq_list[index]
    
    def __setitem__(self, key_or_index, seq):
        if self.contains(key_or_index):
            self._seq_list[self.index(key_or_index)] = seq
        elif key_or_index % 1.0 == 0.0:
            self._seq_list[key_or_index] = seq
        else:
            raise KeyError("No ")
    
    def __len__(self):
        return len(self._seq_list)
    
    def _add_one_seq(self, seq):
        pass
    
    def __iadd__(self, seq_or_seqlist):
        # Pseudocode:
        # if seq_or_seqlist is a seqlist:
        #     for seq in seqlist:
        #         self._add_one_seq(seq)
        # else:
        #     self._add_one_seq(seq)
        pass
    
    def index(self):
        pass
    
    def append(self, seq_or_seqlist):
        self += seq_or_seqlist
    
    def abund_total(self):
        pass
    
    def normalize(self):
        pass
    
    def random_seq(self):
        pass
    
    def sort_by(self, *args):
        pass
    
    def copy(self):
        pass
    
class SeqMatrix():
    """Encodes all-by-all relationships for a list of sequences.
    
    """
    def __init__(self):
        pass

def parse_file(fasta):
    """Takes a open file object and returns a Seqs object.
    
    """
    pass
