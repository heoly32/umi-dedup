from __future__ import division
import time, sys

DEFAULT_INTERVAL = 1

class ProgressTrackerByPosition:
	'''
	given pysam.AlignedSegment objects, report progress along a pysam.AlignmentFile according to its header
	'''
	
	def __init__ (self,
		alignment_file,
		interval = DEFAULT_INTERVAL
	):
		self.lengths =					alignment_file.lengths
		self.total =						sum(self.lengths)
		self.interval =					interval
		self.creation_time =		time.time()
		self.start_time =				self.creation_time
		self.last_update_time =	self.start_time
	
	def update (self, alignment):
		current_time = time.time()
		if current_time - self.last_update_time >= self.interval:
			progress = (sum(self.lengths[:alignment.reference_id]) + (alignment.reference_start + 1)) / self.total # alignment.reference_start + 1 because it's 0-based, so this will create a division by 0 if you hit the very first position of the very first reference
			elapsed_time = current_time - self.start_time
			remaining_time = elapsed_time / progress - elapsed_time
			sys.stderr.write('\r' + ' ' * 50 + '\r%i%% (%i s remaining)' % (100 * progress, remaining_time))
			self.last_update_time = current_time
	
	def reset (self):
		self.start_time = time.time()

	def __del__ (self):
		sys.stderr.write('\r' + ' ' * 50 + '\rcompleted in %i s\n\n' % (time.time() - self.creation_time))

