from __future__ import division
import collections, itertools, parse_sam

def make_umi_list (length, alphabet = 'ACGT'):
	return (''.join(umi) for umi in itertools.product(alphabet, repeat = length))

def make_umi_counts (umi_list, counts = None):
	if counts:
		return collections.OrderedDict((umi, count) for umi, count in zip(umi_list, counts))
	else:
		return collections.OrderedDict((umi, 0) for umi in umi_list)

def read_umi_counts_from_table (infile):
	result = collections.OrderedDict()
	for line in infile:
		split_line = line.split()
		if len(split_line) >= 2: result[split_line[0]] = int(split_line[1])
	if len(result) == 0: raise RuntimeError('bad format in UMI table')
	return result

def read_umi_counts_from_sam (infile):
	umi_totals = None
	umi_length = 0
	file_pos = infile.tell() # save position so we can return there afterward
	infile.reset()
	for read in infile:
		umi = parse_sam.get_umi(read)
		if not umi_length:
			umi_length = len(umi)
			umi_totals = make_umi_counts(make_umi_list(umi_length))
		elif len(umi) != umi_length:
			raise RuntimeError('different UMI length in read ' + read.query_name)
		if umi in umi_totals: umi_totals[umi] += 1 # exclude bad UMIs (containing N)
	infile.seek(file_pos)
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

