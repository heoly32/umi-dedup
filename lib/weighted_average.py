from __future__ import division
from . import umi_data, apportion_counts

# this approach was proposed by Fumiaki Katagiri

def estimate_count (umi_counts):
	c0 = len(umi_counts) - umi_counts.n_nonzero()
	c1 = total = 0
	for count in umi_counts.nonzero_values():
		total += count
		if count == 1: c1 += 1
	return round( (c0 * umi_counts.n_nonzero() + c1 * total) / (c0 + c1) )

def deduplicate_counts (umi_counts):
	return apportion_counts.apportion_umi_values(umi_counts, estimate_count(umi_counts))

