from __future__ import division
from . import umi_data

def estimate_p (umi_counts):
	n_hit, n_umi = 0
	for count in umi_counts.nonzero_values():
		n_hit += count
		n_umi += 1
	return n_umi / n_hit

def deduplicate_counts (umi_counts):
	return umi_data.UmiValues([(key, 1) for key in umi_counts.nonzero_keys()])

