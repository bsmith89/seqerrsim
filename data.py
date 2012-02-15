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
        try: # see if sequence already has attributes
            seq_attrs = sequence.get_attrs
        except AttributeError: # if not, move on.
            pass
        else:
            for attr in seq_attrs:
                self.__setattr__(attr, seq_attrs[attr])
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
        
    def is_seq_object(self):
        return True


class SeqList():
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seq_list = []):
        self._seq_list = []
        self.append_copy(seq_list)
    
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
        else:
            return False
        
    def __len__(self):
        return len(self._seq_list)
        
    def __index__(self, value):
        if not value in self:
            raise ValueError("%s is not in the SeqList" % value)
        else:
            return self._seq_list.index(value)
    
    def __getitem__(self, key_or_index):
        try: # act like it's a key
            index = self.seq_list.index(key_or_index)
            return self._seq_list[index]
        except ValueError:
            try: # act like it's an index
                index = key_or_index
                return self._seq_list[index]
            except TypeError:
                raise KeyError("%s is not in the SeqList" % key_or_index)
    
    def __setitem__(self, key_or_index, value):
        try: # act like it's a key
            index = self.seq_list.index(key_or_index)
            self._seq_list[index] = value
        except ValueError:
            try: # act like it's an index
                index = key_or_index
                self._seq_list[index] = value
            except TypeError:
                key = key_or_index
                raise KeyError("%s is not in the SeqList" % key)
        
    def __delitem__(self, key_or_index):
        try: # act like it's a key
            index = self.seq_list.index(key_or_index)
            del self._seq_list[index]
        except ValueError:
            try: # act like it's an index
                index = key_or_index
                del self._seq_list[index]
            except TypeError:
                key = key_or_index
                raise KeyError("%s is not in the SeqList" % key)
    
    def append(self, seq):
        """Add a single Seq object to SeqList
        
        """
        if seq in self:
            self[seq] = Seq(seq, abund = self[seq] + seq.abund)
        else:
            self._seq_list.append(seq)
        return self
        
    def extend(self, seqlist):
        """Add a list of Seq objects to SeqList.
        
        TODO: this
        
        """
        for seq in seqlist:
            self.append(seq)
        return self
            
    def __iadd__(self, seq_or_list):
        """Append the Seq object(s) to the SeqList.
        
        SeqList() += Seq() <==> SeqList().append(Seq())
        SeqList() += SeqList() <==> SeqList().extend(SeqList())
        
        """
        if seq_or_list.is_seq_object():
            self.append(seq_or_list)
        else:
            for seq in seq_or_list:
                self.append(seq)
        return self
    
    def abund_total(self):
        abund_total = 0.0
        for seq in self:
            abund_total += seq.abund
        return abund_total
    
    def normalize(self):
        abund_total = self.abund_total()
        for seq in self:
            self[seq] = Seq(seq, abund = seq.abund / abund_total)
    
    def random_seq(self):
        abund_total = self.total_abund()
        index = random.random()
        adjusted_index = index * abund_total
        tally = 0.0
        for seq in self:
            tally += seq.abund
            if tally > adjusted_index:
                return Seq(str(seq), abund = 1)
    
    def sort_by(self, reverse = [False], *args):
        """Sorts the SeqList by the arguments in *args.
        
        Forward or reverse is listed in reverse as booleans
        reverse[i] == True: args[i] should be sorted in reverse.
        
        TODO: implement
        
        """
        pass
    
    def _sort_by_one_attr(self, reverse = False, attr):
        """Hidden method for a single attribute sort    
        
        """
        pass
    
    def _attr_zip(self, attr):
        attr_zip = []
        for seq in self:
            attr_zip.append(seq.__getattr__[attr], )
        
    
    def copy(self):
        return SeqList(self)

def parse_file(fasta):
    """Takes a open file object and returns a Seqs object.
    
    """
    pass
