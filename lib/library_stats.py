from __future__ import division
import math

EPSILON = 1e-9 # max error allowed for Newton's method

def mean(x):
	total = n = 0
	for i in x:
		total += i
		n += 1
	return total / n

def entropy (x):
	x = list(x)
	x_sum = sum(x)
	p = (x_i / x_sum for x_i in x)
	return - sum(p_i * math.log(p_i) for p_i in p if p_i != 0)

def w0(x): # Lambert W function using Newton's method, http://code.activestate.com/recipes/577729-lambert-w-function/
	w = x
	while True:
		ew = math.exp(w)
		wNew = w - (w * ew - x) / (w * ew + ew)
		if abs(w - wNew) <= EPSILON: break
		w = wNew
	return w

def estimate_library_size(distinct_reads, total_reads): # from Lander-Waterman equation, http://sourceforge.net/p/samtools/mailman/samtools-help/thread/DUB405-EAS154589A1ACEF2BE4C573D4592180%40phx.gbl/#msg31513744
	assert distinct_reads <= total_reads
	return int(round((distinct_reads * total_reads) / (distinct_reads * w0(- math.exp(- total_reads / distinct_reads) * total_reads / distinct_reads) + total_reads)))

