from __future__ import division
import collections
from . import umi_data, apportion_counts

# extended from Fumiaki Katagiri's approach

def estimate_count (umi_counts):
	count_frequencies = collections.Counter(umi_counts.nonzero_values())
	n_umi = len(umi_counts)
	count_frequencies[0] = n_umi - sum(count_frequencies.values()) # infer the frequency of zero counts
	sorted_counts = sorted(count_frequencies.keys())
	
	weighted_sum = 0
	for count_i, frequency_i in count_frequencies.items():
		n_counted = truncated_sum = 0	
		for count_j in sorted_counts:
			if count_j > count_i: break # skip the counts that will be truncated
			frequency_j = count_frequencies[count_j]
			n_counted += frequency_j
			truncated_sum += count_j * frequency_j
		truncated_sum += (count_i + 1) * (n_umi - n_counted) # infer the sum of all the truncated counts
		weighted_sum += frequency_i * truncated_sum
	return round(weighted_sum / n_umi)

def deduplicate_counts (umi_counts):
	return apportion_counts.apportion_umi_values(umi_counts, estimate_count(umi_counts))

