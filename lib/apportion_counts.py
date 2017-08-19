from __future__ import division
import numpy
from . import umi_data

def safe_round (x):
	assert(x > 0)
	if x < 1:
		return 1
	else:
		return round(x)

def apportion_counts (counts, target_sum):
	counts = list(counts)
	original_sum = sum(counts)
	assert(target_sum <= original_sum)
	assert(target_sum >= sum(int(count > 0) for count in counts))
	
	# initialize: start with simple proportionality but ensure all nonzero counts remain at least 1
	divisor = float(original_sum) / target_sum
	perfect_targets = [count / divisor for count in counts]
	results = list(map(safe_round, perfect_targets))
	residuals = [result - target for result, target in zip(results, perfect_targets)]
	remaining_counts = target_sum - sum(results)
	
	# adjust: add or subtract counts according to the most extreme residuals
	while remaining_counts > 0:
		which_to_increment = numpy.argmin(residuals)
		results[which_to_increment] += 1
		residuals[which_to_increment] += 1
		remaining_counts -= 1
	
	if remaining_counts < 0:
		for i in range(len(counts)):
			if results[i] <= 1: residuals[i] = -numpy.inf # set residuals of counts 1 or 0 to negative infinity so they're never reduced further
		while remaining_counts < 0:
			which_to_decrement = numpy.argmax(residuals)
			results[which_to_decrement] -= 1
			if results[which_to_decrement] == 1:
				residuals[which_to_decrement] = -numpy.inf
			else:
				residuals[which_to_decrement] -= 1
			remaining_counts += 1
	assert(sum(results) == target_sum)
	return results

def apportion_umi_values (counts, target_sum):
	return umi_data.UmiValues(zip(counts.nonzero_keys(), apportion_counts(counts.nonzero_values(), target_sum)))

