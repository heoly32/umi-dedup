from __future__ import division
import collections
from . import umi_data, apportion_counts

# extended from Fumiaki Katagiri's approach

def estimate_count (umi_counts):
	count_frequencies = collections.Counter(umi_counts.nonzero_values())
	count_frequencies[0] = len(umi_counts) - sum(count_frequencies.values()) # infer the frequency of zero counts
	weighted_sum = 0
	for count_i, frequency_i in count_frequencies.items():
		truncated_sum = 0	
		for count_j, frequency_j in count_frequencies.items():
			truncated_sum += frequency_j * min(count_j, count_i + 1)
		weighted_sum += frequency_i * truncated_sum
	return round(weighted_sum / len(umi_counts))

def deduplicate_counts (umi_counts):
	return apportion_counts.apportion_umi_values(umi_counts, estimate_count(umi_counts))

