from __future__ import division
import collections, itertools, re, pysam
from . import parse_sam

DEFAULT_ALPHABET = 'ACGT' # expected characters in UMI sequences
DEFAULT_SEPARATOR = '+' # what separates the two UMIs in paired-end read names
RE_EXCLUSION = re.compile('[^%s]' % (DEFAULT_ALPHABET + DEFAULT_SEPARATOR)) # match any unexpected character (like N)

def umi_is_good (umi):
	return (RE_EXCLUSION.search(umi) is None)

def get_separator_position (umi):
	try:
		return umi.index(DEFAULT_SEPARATOR)
	except ValueError:
		return None

def make_umi_list (length, separator_position = None, alphabet = DEFAULT_ALPHABET):
	for sequence in itertools.product(alphabet, repeat = length):
		umi = ''.join(sequence)
		if separator_position is not None: umi = umi[:separator_position] + DEFAULT_SEPARATOR + umi[separator_position:]
		yield umi

class UmiValues:
	def __init__ (self,
		initial_data =       None, # should be list of (key, value) pairs
		length =             None,
		separator_position = None,
		alphabet =           DEFAULT_ALPHABET
	):
		self.alphabet = alphabet
		if initial_data is not None:
			assert length is None and separator_position is None
			initial_data = list(initial_data)
			example_umi = initial_data[0][0]
			assert umi_is_good(example_umi)
			self.length = len(example_umi) - example_umi.count(DEFAULT_SEPARATOR)
			self.separator_position = get_separator_position(example_umi)
			self.data = collections.Counter()
			for pair in initial_data:
				assert(self.is_valid(pair[0])) # verify valid entries
				if pair[1] != 0: self.data[pair[0]] = pair[1]
		else:
			assert length is not None
			self.length =             length
			self.separator_position = separator_position
			self.data =               collections.Counter()

	def __len__ (self):
		return len(self.alphabet) ** self.length

	def is_valid (self, key):
		return (len(key) == self.length + (self.separator_position is not None) and get_separator_position(key) == self.separator_position)

	def __getitem__ (self, key):
		if not self.is_valid(key): raise KeyError(key)
		return self.data[key]

	def __setitem__ (self, key, value):
		if not self.is_valid(key): raise KeyError(key)
		if value == 0: # delete zeroes; conveniently __setitem__ also covers incrementing/decrementing
			try:
				del self.data[key]
			except KeyError:
				pass
		else:
			self.data[key] = value

	# standard dict functions
	def keys (self):
		return make_umi_list(self.length, self.separator_position, self.alphabet)
	def values (self):
		return (self.data[key] for key in self.keys())
	def items (self):
		return ((key, self.data[key]) for key in self.keys())
	
	# special dict functions for nonzero values only
	def nonzero_keys (self):
		return (key for key in sorted(self.data.keys()) if self.data[key]) # double check in case the Counter contains zeroes
	def nonzero_values (self):
		return (self.data[key] for key in self.nonzero_keys())
	def nonzero_items (self):
		return ((key, self.data[key]) for key in self.nonzero_keys())
	
	# convenience function
	def n_nonzero (self): # return the number of nonzero counts
		return len(self.data)

def parse_umi (read_name, truncate = None):
	for label in read_name.split(' ')[:2]: # to allow NCBI format or regular Illumina
		if label.count(':') in (5, 7): # Casava pre-1.8: should be 5 (4 + the UMI hack); Casava 1.8+ / bcl2fastq 2.17+: should be 7 (with optional UMI field)
			umi = label.partition('#')[0].partition('/')[0].rpartition(':')[2] # don't include the space or # and the stuff after it, if present
			return (umi if truncate is None else umi[:truncate + umi.count(DEFAULT_SEPARATOR)]) # don't count the pair separator when truncating
	# only get here if nothing was found
	raise RuntimeError('read name %s does not contain UMI in expected Casava/bcl2fastq format' % label)

def set_umi (alignment, umi = None, truncate = None): # set an alignment's MI tag to a specific sequence, or parse it from the read name if no specific sequence provided; then return the updated alignment
	if umi is None: umi = parse_umi(alignment.query_name, truncate)
	alignment.set_tag('MI', umi)
	return alignment

def get_umi (alignment):
	return alignment.get_tag('MI')

def read_umi_counts_from_table (in_file, truncate = None):
	result = None
	for line in in_file:
		split_line = line.split()
		try:
			umi = split_line[0]
			if truncate is not None: umi = umi[:truncate]
			try:
				result[umi] = int(split_line[1])
			except TypeError: # no result yet
				result = UmiValues([(umi, 1)])
		except IndexError: # insufficient data in this line
			pass
	if not result: raise RuntimeError('bad format in UMI table')
	return result

def read_umi_counts_from_reads (in_file, truncate = None): # in_file should be a pysam.Samfile or a Bio.SeqIO.parse in 'fastq' format, or at least contain an Illumina-formatted name in either 'query_name' or 'id'
	umi_totals = None
	for read in in_file:
		if hasattr(read, 'is_paired') and read.is_paired and not parse_sam.alignment_is_properly_paired(read): continue
		try:
			read_name = read.query_name
		except AttributeError:
			read_name = read.id # EAFP; if this isn't found either, AttributeError is still raised
		umi = parse_umi(read_name, truncate)
		if not umi_is_good(umi): continue
		try:
			umi_totals[umi] += 1
		except TypeError: # no data yet
			umi_totals = UmiValues([(umi, 1)])
	return umi_totals

def mark_duplicates (reads, n):
	'''
	mark 'n' reads from 'reads' as duplicates
	reads to mark as the duplicates are chosen by lowest base quality
	'''
	assert len(reads) >= n
	if n > 0:
		sorted_reads = sorted(reads, key = parse_sam.get_quality)
		for i in range(n): sorted_reads[i].is_duplicate = True
	return reads

