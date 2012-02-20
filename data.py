"""Data structures and functions for retrieving/storing sequence data

"""
import random
import scipy as sp
from scipy.spatial.distance import hamming


DEFAULT_ALPHABET = ['A', 'C', 'G', 'T']
DEFAULT_ABUND = 1.0

ALPHABET = DEFAULT_ALPHABET

class Seq(object):
    """Sequence and its abundance.
    
    Also defines methods for comparing sequences.
    
    """
    def __init__(self, sequence, **kwargs):
        self._sequence = sp.array(list(sequence))
        try: # see if sequence already has attributes
            seq_attrs = sequence.get_attrs()
        except AttributeError: # if not, move on.
            pass
        else:
            for attr in seq_attrs:
                self.__setattr__(attr, seq_attrs[attr])
        for attr in kwargs.keys():
            self.__setattr__(attr, kwargs[attr])
    
    def __iter__(self):
        """Iterate over the nucleotides in Seq
        
        """
        for nucl in self._sequence:
            yield nucl
    
    def __len__(self):
        """Return the length of Seq
        
        """
        return len(self._sequence)
    
    def __eq__(self, other):
        """Return if other is same sequence as self.
        
        Supposed to work on any sequence of nucleotides,
        whether a Seq object or not.
        TODO: Test this
        
        Does not consider anything except for nucleotide order.
        
        """
        try:
            assert len(self) == len(other)
        except (AssertionError, TypeError):
            return False
        try:
            other_seq_array = other.get_seq_array()
        except AttributeError:
            other_seq_array = sp.array(list(other))
        return all(self.get_seq_array() == other_seq_array)
    
    def get_attrs(self):
        """Returns a list of all the attributes except the Seq array.
        
        """
        attr_dict = dict(self.__dict__)
        del attr_dict['_sequence']
        return attr_dict
            
    def get_seq_array(self):
        """Returns a copy of the numpy array of Seq's sequence
        
        """
        return sp.array(self._sequence)
    
    def __lt__(self, other):
        """Return if the sequence is alphabetically less than other
        
        """
        return str(self) < str(other)
        
    def __str__(self):
        """Return a str of sequence
        
        """
        return ''.join(self.get_seq_array())
    
    def fasta(self, name = None):
        """Return a string in FASTA format of self
        
        Also includes all attribute values on self.
        """
        title_str = ">"
        if name is None:
            try:
                name = self.name
            except AttributeError:
                pass
        title_str += str(name)
        attr_list = self.get_attrs()
        for attr in attr_list:
            title_str += "|%s=%s" % (str(attr), str(attr_list[attr]))
        title_str += "\n"
        return title_str + str(self) + "\n"
    
    def __repr__(self):
        """Return a string representation, evaluates identical to self
        
        """
        other_attrs = self.get_attrs()
        out_str = "Seq('%s', **%s)" % (str(self), str(other_attrs))
        return out_str
    
    def get_abund(self):
        """Return the abundance of self
        
        Assumes that self.abund has been set, although that is not
        necessarily the case.  Raises a AttributeError if self.abund
        has not been set.
        
        """
        return self.abund
    
    def set_abund(self, value):
        """Sets self.abund to value.
        
        """
        self.abund = float(value)
        
    def dist(self, other):
        """Return the hamming distance between self and other
        
        Attempts to have other provide its own numpy array, but able
        to produce one if necessary.
        
        """
        try:
            other_seq_array = other.get_seq_array()
        except AttributeError:
            other_seq_array = sp.array(list(other))
        return hamming(self.get_seq_array(), other_seq_array)
        
    def is_seq_object(self):
        """Allows external functions to make sure a seq is a Seq object
        
        """
        return True


class SeqList(object):
    """List of Seq objects.
    
    Also defines methods for analyzing the sequence sample
    
    """
    def __init__(self, seq_list = []):
        self._seq_list = []
        self.extend_with_copy(seq_list)
    
    def __iter__(self):
        """Iterate over the sequences in SeqList
        
        """
        return iter(self._seq_list)
    
    def __reversed__(self):
        """Iterate over the reversed sequences in SeqList
        
        """
        return reversed(self._seq_list)
    
    def __repr__(self):
        """Return a string representation, evaluates identical to self
        
        """
        seqs_str = ""
        for seq in self:
            seqs_str += "%s\n" % repr(seq)
        repr_str = "SeqList(\n%s)" % seqs_str
        return repr_str
    
    def __str__(self):
        """Return a tabular view of all seqs in self
        
        Returns a tab delimited table of the Seq's in self along with
        their attributes.
        
        TODO: Break this up into smaller processes.
        
        """
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
        """Return a string of self which follows FASTA format.
        
        """
        fasta_str = ""
        seq_count = 1
        for seq in self:
            try:
                name = seq.ident
            except AttributeError:
                name = "seq%d" % seq_count
            seq_fasta = seq.fasta(name = name)
            fasta_str += seq_fasta
            seq_count += 1
        return fasta_str
    
    def __contains__(self, seq):
        """Return True if sequence seq is found in self
        
        Uses Seq.__eq__ which does not consider anything except for
        the nucleotide sequence in seq.
        
        """
        if seq in self._seq_list:
            return True
        else:
            return False
        
    def __len__(self):
        """Return the number of Seq objects in self
        
        """
        return len(self._seq_list)
        
    def __index__(self, value):
        """Return the index of self which contains value
        
        Dependent on the Seq.__eq__() function, which is supposed to
        consider any two sequences 
        
        """
        if not value in self:
            raise ValueError("%s is not in the SeqList" % value)
        else:
            return self._seq_list.index(value)
    
    def __getitem__(self, key_or_index):
        """Set the item with the provided index or sequence
        
        SeqList['AAAA'] <==> SeqList[index_of_AAAA]
        
        Evaluates the key_or_index first as a key
        TODO: account for negative indices
        
        """
        
        try: # act like it's a key
            index = self._seq_list.index(key_or_index)
        except ValueError: # oops, it wasn't a key.
            pass
        else:
            return self._seq_list[index]
        try: # act like it's an index
            index = key_or_index
            return self._seq_list[index]
        except TypeError: # doesn't look like a valid index either
            raise KeyError("%s is not in the SeqList" % key_or_index)
    
    def __setitem__(self, key_or_index, value):
        """Set the item with the provided index or sequence to value
        
        If the key provided is different from str(value) the index now
        points to value, but with a new key (i.e. str(value)).
        TODO: account for negative indices

        
        """
        try: # act like it's a key
            index = self._seq_list.index(key_or_index)
        except ValueError: # oops, look like it wasn't a key
            pass
        else:
            self._seq_list[index] = value
            return
        try: # act like it's an index
            self._seq_list[key_or_index] = value
        except TypeError: # not a valid index either
            raise KeyError("%s is not in the SeqList" % str(key_or_index))
        
    def __delitem__(self, key_or_index):
        """Delete the item with provided index or sequence
        
        TODO: account for negative indices
        
        """
        try: # act like it's a key
            index = self._seq_list.index(key_or_index)
        except ValueError: # oops, look like it wasn't a key
            pass
        else:
            del self._seq_list[index]
            return
        try: # act like it's an index
            del self._seq_list[key_or_index]
        except TypeError: # not a valid index either
            raise KeyError("%s is not in the SeqList" % str(key_or_index))
    
    def append(self, seq):
        """Add a single Seq object to SeqList
        
        """
        if seq in self:
            self[seq] = Seq(seq, abund = self[seq].abund + seq.abund)
        else:
            self._seq_list.append(seq)
        
    def extend(self, seqlist):
        """Add a list of Seq objects to SeqList.

        """
        for seq in seqlist:
            self.append(seq)
            
    def __iadd__(self, seq_or_list):
        """Append the Seq object(s) to the SeqList.
        
        SeqList() += Seq() <==> SeqList().append(Seq())
        SeqList() += SeqList() <==> SeqList().extend(SeqList())
        
        TODO: test this
        TODO: write a doc-test
        
        """
        try: # treat it like a Seq object
            assert seq_or_list.is_seq_object()
        except AttributeError:
            pass
        else:
            self.append(seq_or_list)
            return self
        # okay, so it must be a Seq list
        for seq in seq_or_list:
            self.append(seq)
        return self
    
    def abund_total(self):
        """Returns the sum of all sequence abundances
        
        """
        abund_total = 0.0
        for seq in self:
            abund_total += seq.abund
        return abund_total
    
    def normalize(self):
        """Repopulate the SeqList with Seq's with normalized abunds
        
        """
        abund_total = self.abund_total()
        for seq in self:
            self[seq] = Seq(seq, abund = seq.abund / abund_total)
    
    def random_seq(self):
        """Return a random Seq object from SeqList
        
        """
        abund_total = self.abund_total()
        index = random.random()
        adjusted_index = index * abund_total
        tally = 0.0
        for seq in self:
            tally += seq.abund
            if tally > adjusted_index:
                return Seq(str(seq), abund = 1)
    
    def sort_by(self, *attrs, **kwargs):
        """Sorts the SeqList by the arguments in *args.
        
        Forward or reverse is listed in reverse as booleans
        reverse[i] == True: args[i] should be sorted in reverse.
        
        TODO: implement
        
        """
        try:
            reverse = kwargs['reverse']
        except KeyError:
            reverse = None
        if reverse is None:
            reverse = [False]*len(attrs)
        elif reverse == True:
            reverse = [True]*len(attrs)
        elif reverse == False:
            reverse = [False]*len(attrs)
        for attr, if_reverse in reversed(zip(attrs, reverse)):
            self._sort_by_one_attr(attr, reverse = if_reverse)
    
    def _sort_by_one_attr(self, attr, reverse = False):
        """Hidden method for a single attribute sort    
        
        """
        if reverse is not True or False:
            raise ValueError('Reverse must be a boolean, not "%s"' %
                             str(reverse))
        sorted_zip = sorted(self._attr_zip(attr), reverse = reverse)
        new_seq_list = []
        for pair in sorted_zip:
            new_seq_list.append(pair[-1])
        self._seq_list = new_seq_list
            
    def _attr_zip(self, attr):
        """Return a list of tuples, [(attr_value, seq_object),...].
        
        """
        attr_zip = []
        for seq in self:
            if attr == 'sequence':
                value = str(seq)
            else:
                value = seq.__getattribute__(attr)
            attr_zip.append((value, seq))
        return attr_zip
        
    
    def copy(self):
        """Return a copy of self.
        
        """
        return SeqList(self)

def parse_file(fasta):
    """Takes a open file object and returns a SeqList object.
    
    """
    out_list = SeqList()
    curr_seq = ""
    curr_attrs = {}
    for line in fasta:
        if line.startswith(';'): # ignore comment lines
            continue
        elif line.startswith(">"): # is it a title line?
            if curr_seq is not "": # unless this is the first Seq
                # add the previously stored information as a Seq object
                out_list.append(Seq(curr_seq, **curr_attrs))
                # and reset the curr_ values
                curr_seq = ""
                curr_attrs = {}
            # collect new information
            attr_strs = line[1:].split('|') # everything's split by '|'
            # The identifier (or name) is found right after the
            # '>'
            identifier = attr_strs[0]
            curr_attrs['ident'] = identifier
            attr_index = 1 # just for naming unnamed attributes
            for attr_str in attr_strs[1:]:
                try:
                    # assume that attributes are named, and connected
                    # by a '='
                    attr_name, value = attr_str.split('=')
                except ValueError: # but if not...
                    value = attr_str
                    # name it 'attr(01/02/.../10/...)'
                    attr_name = "attr%02d" % attr_index
                attr_name = attr_name.strip()
                value = value.strip()
                try:
                    # treat value like a float first
                    value = float(value)
                except ValueError:
                    # otherwise it's just gonna stay a string
                    pass
                curr_attrs[attr_name] = value
                attr_index += 1
        else:
            curr_seq += line.strip()
    if curr_seq is not "":
        out_list.append(Seq(curr_seq, **curr_attrs))
    return out_list
                    
