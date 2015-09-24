import collections
from . import parse_sam

class ReadBuffer(collections.deque()):
	'''
	Queue of SAM reads grouped by start (5') position. Expects to have reads added in order of leftmost position (standard SAM sorting) and produces them in the same order. This allows collecting reads by start position without breaking the sort order.
	'''
	def popleft():
		
