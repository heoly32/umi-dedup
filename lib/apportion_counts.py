from __future__ import division
from . import umi_data

def apportion_counts (counts, target_sum):
	counts = list(counts)
	result = [int(count > 0) for count in counts]
	divisor = float(target_sum) / sum(counts)
	quotients = (count / divisor for count in counts)
	residuals = [quotient - new_count for quotient, new_count in zip(quotients, result)]
	remaining_counts = target_sum - sum(result)
	while remaining_counts > 0:
		which_to_increment = residuals.index(max(residuals))
		result[which_to_increment] += 1
		residuals[which_to_increment] -= 1
		remaining_counts -= 1
	return result

