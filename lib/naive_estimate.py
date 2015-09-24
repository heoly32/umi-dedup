from __future__ import division

def estimate_p (umi_counts):
	n_hit = sum(umi_counts.values())
	n_umi = sum(count > 0 for count in umi_counts.values())
	return n_umi / n_hit

