from __future__ import division
import collections, itertools, re, pysam, parse_sam

alphabet = 'ACGT' # expected characters in UMI sequences
default_pair_separator = '+' # what separates the two UMIs in paired-end read names
re_exclusion = re.compile('[^%s]' % (alphabet + default_pair_separator)) # match any unexpected character (like N)

def umi_is_good (umi):
	return (re_exclusion.search(umi) is None)

def make_umi_list (length, separator_position = None, alphabet = alphabet):
	for sequence in itertools.product(alphabet, repeat = length):
		umi = ''.join(sequence)
		if separator_position is not None: umi = umi[:separator_position] + default_pair_separator + umi[separator_position:]
		yield umi

def get_umi (read_name, truncate = None):
	for label in read_name.split(' ')[:2]: # to allow NCBI format or regular Illumina
		if label.count(':') in (5, 7): # Casava pre-1.8: should be 5 (4 + the UMI hack); Casava 1.8+ / bcl2fastq 2.17+: should be 7 (with optional UMI field)
			umi = label.partition('#')[0].partition('/')[0].rpartition(':')[2] # don't include the space or # and the stuff after it, if present
			return (umi if truncate is None else umi[:truncate + umi.count(default_pair_separator)]) # don't count the pair separator when truncating
	# only get here if nothing was found
	raise RuntimeError('read name %s does not contain UMI in expected Casava/bcl2fastq format' % label)

def get_pair_separator_position (umi):
	try:
		return umi.index(default_pair_separator)
	except ValueError:
		return None

def read_umi_counts_from_table (in_file, truncate = None):
	result = collections.OrderedDict()
	for line in in_file:
		split_line = line.split()
		try:
			umi = split_line[0]
			if truncate is not None: umi = umi[:truncate]
			try:
				result[umi] = int(split_line[1])
			except IndexError: # no count given
				result[umi] = 0
		except IndexError: # empty line
			pass
	if not result: raise RuntimeError('bad format in UMI table')
	return result

def read_umi_counts_from_reads (in_file, truncate = None): # in_file should be a pysam.Samfile or a Bio.SeqIO.parse in 'fastq' format, or at least contain an Illumina-formatted name in either 'query_name' or 'id'
	umi_totals = collections.Counter()
	umi_length = None
	for read in in_file:
		if hasattr(read, 'is_paired') and read.is_paired and not parse_sam.alignment_is_properly_paired(read): continue
		try:
			read_name = read.query_name
		except AttributeError:
			read_name = read.id # EAFP; if this isn't found either, AttributeError is still raised
		umi = get_umi(read_name, truncate)
		if len(umi) - umi.count(default_pair_separator) != umi_length:
			if umi_length is None:
				umi_length = len(umi) - umi.count(default_pair_separator)
				umi_totals = make_umi_counts(make_umi_list(umi_length, get_pair_separator_position(umi)))
			else:
				raise RuntimeError('different UMI length in read ' + read_name)
		umi_totals[umi] += 1
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

