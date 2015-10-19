import collections, parse_sam, umi_data, optical_duplicates, naive_estimate, bayes_estimate

class DuplicateMarker:
	'''
	a generator of duplicate-marked alignments (like pysam.AlignedSegment objects), given an iterable of alignments (like a pysam.AlignmentFile) as input
	assumes the input iterable is sorted by leftmost coordinate (SAMtools style) and yields alignments in the same order
	does not yield unusable alignments
	central concept: cheat a little, detecting duplicate reads by the fact that they align to the same (start) position rather than by their sequences
	implementation: as you traverse the coordinate-sorted input reads, add each read to both a FIFO buffer (so they can be output in the same order) and a dictionary that groups reads by start position and strand (which is what you need for deduplication); this is not inefficient because both data structures contain pointers to the same alignment objects
	then, after each input read, check the left end of the buffer to tell whether the oldest read is in a position that will never accumulate any more hits (because the input is sorted by coordinate); if so, estimate the duplication at that position, mark all the reads there accordingly, then output the read with the appropriate marking
	'''
	def __init__(self, 
		alignments,
		umi_frequency = None,
		algorithm = 'bayes',
		optical_dist = optical_duplicates.DEFAULT_DIST,
		truncate_umi = None,
		nsamp = bayes_estimate.DEFAULT_NSAMP,
		nthin = bayes_estimate.DEFAULT_NTHIN,
		nburn = bayes_estimate.DEFAULT_NBURN
	):
		self.alignments = alignments
		self.umi_frequency = umi_frequency
		self.optical_dist = optical_dist
		self.truncate_umi = truncate_umi
		if algorithm == 'naive':
			self.umi_dup_function = naive_estimate.deduplicate_counts
		elif algorithm in ('bayes', 'uniform-bayes'):
			self.umi_dup_function = lambda counts: bayes_estimate.deduplicate_counts(umi_counts = counts, nsamp = nsamp, nthin = nthin, nburn = nburn, uniform = (algorithm == 'uniform-bayes'))
		else:
			raise NotImplementedError
		self.alignment_buffer = collections.deque()
		self.pos_tracker = ({}, {}) # data structure containing alignments by position rather than sort order; top level is by strand (0 = forward, 1 = reverse), then next level is by 5' read start position (dict since these will be sparse and are only looked up by identity), and that contains a variety of data; there is no level for reference ID because there is no reason to store more than one chromosome at a time
		self.counts = collections.Counter()
		self.output_generator = self.get_marked_alignment()
	
	def __iter__(self):
		return self
	
	def __next__(self):
		return next(self.output_generator)
	
	def iter(self):
		return self
	
	def next(self):
		return next(self.output_generator)
	
	def pop_buffer(self):
		alignment = self.alignment_buffer.popleft()
		start_pos = parse_sam.get_start_pos(alignment)
		pos_data = self.pos_tracker[alignment.is_reverse][start_pos]
		
		# deduplicate reads
		if not pos_data['deduplicated']:
			alignments = pos_data['alignments']
			
			# first pass: mark optical duplicates
			if self.optical_dist != 0:
				for opt_dups in optical_duplicates.get_optical_duplicates(alignments, self.optical_dist):
					for alignment in umi_data.mark_duplicates(opt_dups, len(opt_dups) - 1):
						# remove duplicate reads from the tracker so they won't be considered later (they're still in the read buffer)
						if alignment.is_duplicate: alignments.remove(alignment)
					self.counts['optical duplicate'] += len(opt_dups) - 1
			
			# second pass: mark PCR duplicates
			alignments_by_umi = collections.defaultdict(list)
			for this_alignment in alignments: alignments_by_umi[umi_data.get_umi(this_alignment.query_name, self.truncate_umi)] += [this_alignment]
			dedup_counts = self.umi_dup_function(umi_data.make_umi_counts(alignments_by_umi.keys(), map(len, alignments_by_umi.values())))
			for umi, alignments, dedup_count in zip(alignments_by_umi.keys(), alignments_by_umi.values(), dedup_counts.values()):
				assert alignments
				n_dup = len(alignments) - dedup_count
				umi_data.mark_duplicates(alignments, n_dup)
				self.counts['PCR duplicate'] += n_dup
				self.counts['UMI rescued'] += 1
				self.counts['algorithm rescued'] += dedup_count - 1
			self.counts['UMI rescued'] -= 1 # count the first read at this position as distinct
			self.counts['distinct'] += 1
			
			pos_data['deduplicated'] = True
		
		# garbage collection
		if alignment is pos_data['last alignment']: del self.pos_tracker[alignment.is_reverse][start_pos]
		
		return alignment
	
	def get_marked_alignment(self):
		for alignment in self.alignments:
			self.counts['alignment'] += 1
			if not ( # verify sorting
				(not self.alignment_buffer) or (
					alignment.reference_id > self.alignment_buffer[-1].reference_id or (
						alignment.reference_id == self.alignment_buffer[-1].reference_id and
						alignment.reference_start >= self.alignment_buffer[-1].reference_start
					)
				)
			): raise RuntimeError('alignment %s out of order: verify sorting' % alignment.query_name)
			if not parse_sam.alignment_is_good(alignment): continue
			umi = umi_data.get_umi(alignment.query_name, self.truncate_umi)
			if not umi_data.umi_is_good(umi): continue
			alignment.is_duplicate = False # not sure how to handle alignments that have already been deduplicated somehow, so just ignore previous annotations
			start_pos = parse_sam.get_start_pos(alignment)
			self.counts['usable alignment'] += 1
			
			# advance the buffer
			while self.alignment_buffer and (
				self.alignment_buffer[0].reference_id < alignment.reference_id or # new chromosome
				parse_sam.get_start_pos(self.alignment_buffer[0]) < alignment.reference_start # oldest buffer member is now guaranteed not to get any more hits at its position
			): yield self.pop_buffer()
			
			# add the alignment to the tracking data structures
			self.alignment_buffer.extend([alignment])
			try:
				self.pos_tracker[alignment.is_reverse][start_pos]['alignments'] += [alignment]
			except KeyError: # first time we've seen this position+strand
				self.pos_tracker[alignment.is_reverse][start_pos] = {'alignments': [alignment], 'deduplicated': False}
			self.pos_tracker[alignment.is_reverse][start_pos]['last alignment'] = alignment
		
		# flush the buffer
		while self.alignment_buffer: yield self.pop_buffer()

