from __future__ import division
import collections, itertools, parse_sam

def make_umi_list (length, alphabet = 'ACGT'):
	return (''.join(umi) for umi in itertools.product(alphabet, repeat = length))

def make_umi_counts (umi_list, counts = None):
	if counts:
		return collections.OrderedDict((umi, count) for umi, count in zip(umi_list, counts))
	else:
		return collections.OrderedDict((umi, 0) for umi in umi_list)

def deduplicate (reads, p): # assumes 'reads' is a dict where each key is a UMI and each value is a list of pysam.AlignedSegment objects, all marked as is_duplicate = False
	n_hit = sum(len(umi_hits) for umi_hits in reads.values())
	assert p >= len(reads) / n_hit
	
	n_dup = round(n_hit * (1 - p))
	
	# this is stupid; an ideal solution would try to keep the number of retained reads from each UMI roughly proportional to how many hits it had, except always at least 1 for every UMI that had a nonzero count, and prioritize the best reads (longest with highest total score)
	for umi_reads in reads.values():
		try:
			for read in umi_reads:
				if n_dup == 0: raise StopIteration
				read.is_duplicate = True
				n_dup -= 1
		except StopIteration:
			break

