from __future__ import division
import collections, itertools, re, pysam, parse_sam

alphabet = 'ACGT'
re_exclusion = re.compile('[^%s]' % alphabet)

def make_umi_list (length, alphabet = alphabet):
	return (''.join(umi) for umi in itertools.product(alphabet, repeat = length))

def make_umi_counts (umi_list, counts = None):
	try:
		return collections.OrderedDict((umi, count) for umi, count in zip(umi_list, counts))
	except TypeError:
		return collections.OrderedDict((umi, 0) for umi in umi_list)

def get_umi (read_name): # Illumina-specific
	if read_name.count(':') != 7: raise RuntimeError('read name %s does not contain UMI in expected Casava 1.8+ / bcl2fastq 2.17+ format' % read_name)
	return read_name.partition(' ')[0].rpartition(':')[2]

def umi_is_good (umi):
	return (len(umi) > 0 and re_exclusion.search(umi) is None)

def read_umi_counts_from_table (in_file):
	result = collections.OrderedDict()
	for line in in_file:
		split_line = line.split()
		if len(split_line) >= 2: result[split_line[0]] = int(split_line[1])
	if len(result) == 0: raise RuntimeError('bad format in UMI table')
	return result

def read_umi_counts_from_reads (in_file): # in_file should be a pysam.Samfile or a Bio.SeqIO.parse in 'fastq' format, or at least contain an Illumina-formatted name in either 'query_name' or 'id'
	umi_totals, umi_length = None, None
	for read in in_file:
		try:
			read_name = read.query_name
		except AttributeError:
			read_name = read.id # EAFP; if this isn't found either, AttributeError is still raised
		umi = get_umi(read_name)
		if umi_length is None:
			umi_length = len(umi)
			umi_totals = make_umi_counts(make_umi_list(umi_length))
		elif len(umi) != umi_length:
			raise RuntimeError('different UMI length in read ' + read_name)
		try:
			umi_totals[umi] += 1
		except KeyError:
			pass # exclude bad UMIs
	if umi_totals is None: raise RuntimeError('no valid reads detected')
	return umi_totals

def mark_duplicates (reads, target_umi_counts):
	'''
	assumes 'reads' is a dict where each key is a UMI and each value is a list of pysam.AlignedSegment objects (pointers), all marked as is_duplicate = False
	assumes 'target_umi_counts' is a dict where each key is a UMI and each value is the number of non-duplicate reads
	so if len(reads[umi]) == 5 and target_umi_counts[umi] == 3, two reads from this UMI will be marked as duplicates
	reads to mark as the duplicates are chosen by lowest base quality
	'''
	for umi in reads.keys():
		assert len(reads[umi]) >= target_umi_counts[umi]
		if len(reads[umi]) > target_umi_counts[umi]:
			sorted_reads = sorted(reads[umi], key = parse_sam.get_quality)
			for i in range(len(reads[umi]) - target_umi_counts[umi]):
				sorted_reads[i].is_duplicate = True

