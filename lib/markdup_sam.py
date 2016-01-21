import collections, copy, parse_sam, umi_data, optical_duplicates, naive_estimate, bayes_estimate

DUP_CATEGORIES = ['optical duplicate', 'PCR duplicate']

class DuplicateMarker:
	'''
	a generator of duplicate-marked alignments (like pysam.AlignedSegment objects), given an iterable of alignments (like a pysam.AlignmentFile) as input
	assumes the input iterable is sorted by leftmost coordinate (SAMtools style) and yields alignments in the same order
	does not yield unusable alignments
	central concept: cheat a little, detecting duplicate reads by the fact that they align to the same (start) position rather than by their sequences
	implementation: as you traverse the coordinate-sorted input reads, add each read to both a FIFO buffer (so they can be output in the same order) and a dictionary that groups reads by strand, start position, and mate start position (None if single-end); this is not inefficient because both data structures contain pointers to the same alignment objects
	at the strand+position level (before mate start position), track metadata: whether this strand+position has been deduplicated and which alignment was added to it last (when this alignment is processed, the tracker can be deleted safely)
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
			self.umi_dup_function = lambda counts: bayes_estimate.deduplicate_counts(umi_counts = counts, nsamp = nsamp, nthin = nthin, nburn = nburn, uniform = (algorithm == 'uniform-bayes'),  total_counts = self.umi_frequency)
		else:
			raise NotImplementedError
		self.alignment_buffer = collections.deque()
		self.pos_tracker = ({}, {}) # data structure containing alignments by position rather than sort order; top level is by strand (0 = forward, 1 = reverse), then next level is by 5' read start position (dict since these will be sparse and are only looked up by identity), then next level is by 5' start position of mate read, and each element of that contains a variety of data; there is no level for reference ID because there is no reason to store more than one chromosome at a time
		self.counts = collections.Counter()
		self.output_generator = self.get_marked_alignment()
		self.current_reference_id = 0
		self.most_recent_left_pos = 0

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
			for mate_start_pos, alignments_with_this_mate in pos_data['alignments by mate start'].iteritems(): # iterate over mate start positions
				alignments_to_dedup = copy.copy(alignments_with_this_mate)
				
				# if mate is already deduplicated, use those results
				if mate_start_pos in pos_data['already processed']:
					umi_nondup_counts = collections.Counter()
					for already_processed_alignment in alignments_with_this_mate:
						umi = umi_data.get_umi(already_processed_alignment.query_name, self.truncate_umi)
						try:
							category = pos_data['already processed'][mate_start_pos][umi][already_processed_alignment.query_name]
							if category in DUP_CATEGORIES:
								already_processed_alignment.is_duplicate = True
								self.counts[category] += 1
							else:
								already_processed_alignment.is_duplicate = False
								umi_nondup_counts[umi] += 1
							alignments_to_dedup.remove(already_processed_alignment)
						except KeyError: # this alignment is missing from the list, so its mate must have been missing from the data, so it still needs deduplication
							pass
					for umi_nondup_count in umi_nondup_counts.values():
						self.counts['UMI rescued'] += 1
						self.counts['algorithm rescued'] += umi_nondup_count - 1
					self.counts['UMI rescued'] -= bool(umi_nondup_counts) # count the first read at this position as distinct
					self.counts['distinct'] += bool(umi_nondup_counts)
				
				if alignments_to_dedup: # false if they had all been deduplicated already
					alignment_categories = {} # key = alignment query_name, value = which category it is (from DUP_CATEGORIES or otherwise)
					
					# first pass: mark optical duplicates
					if self.optical_dist != 0:
						for opt_dups in optical_duplicates.get_optical_duplicates(alignments_to_dedup, self.optical_dist):
							for dup_alignment in umi_data.mark_duplicates(opt_dups, len(opt_dups) - 1):
								# remove duplicate reads from the tracker so they won't be considered later (they're still in the read buffer)
								if dup_alignment.is_duplicate:
									alignments_to_dedup.remove(dup_alignment)
									alignment_categories[dup_alignment.query_name] = 'optical duplicate'
									self.counts['optical duplicate'] += 1
					
					# second pass: mark PCR duplicates
					alignments_by_umi = collections.defaultdict(list)
					for this_alignment in alignments_to_dedup: alignments_by_umi[umi_data.get_umi(this_alignment.query_name, self.truncate_umi)] += [this_alignment]
					dedup_counts = self.umi_dup_function(umi_data.make_umi_counts(alignments_by_umi.keys(), map(len, alignments_by_umi.values())))
					for umi, alignments_with_this_umi, dedup_count in zip(alignments_by_umi.keys(), alignments_by_umi.values(), dedup_counts.values()):
						assert alignments_with_this_umi
						n_dup = len(alignments_with_this_umi) - dedup_count
						for marked_alignment in umi_data.mark_duplicates(alignments_with_this_umi, n_dup):
							alignment_categories[marked_alignment.query_name] = ('PCR duplicate' if marked_alignment.is_duplicate else 'nonduplicate')
						self.counts['PCR duplicate'] += n_dup
						self.counts['UMI rescued'] += 1
						self.counts['algorithm rescued'] += dedup_count - 1
					self.counts['UMI rescued'] -= bool(alignments_by_umi) # count the first read at this position as distinct
					self.counts['distinct'] += bool(alignments_by_umi)
					
					# pass duplicate marking to mates
					if mate_start_pos is not None:
						for categorized_alignment, category in alignment_categories.iteritems():
							categorized_alignment_umi = umi_data.get_umi(categorized_alignment)
							try:
								self.pos_tracker[not alignment.is_reverse][mate_start_pos]['already processed'][start_pos][categorized_alignment_umi][categorized_alignment] = category
							except KeyError: # first time we've seen this mate strand+position+mate+UMI
								try:
									self.pos_tracker[not alignment.is_reverse][mate_start_pos]['already processed'][start_pos][categorized_alignment_umi] = {categorized_alignment: category}
								except KeyError: # first time we've seen this mate strand+position+mate
									try:
										self.pos_tracker[not alignment.is_reverse][mate_start_pos]['already processed'][start_pos] = {categorized_alignment_umi: {categorized_alignment: category}}
									except KeyError: # first time we've seen this mate strand+position
										self.pos_tracker[not alignment.is_reverse][mate_start_pos] = {
											'alignments by mate start': collections.defaultdict(list),
											'already processed': {start_pos: {categorized_alignment_umi: {categorized_alignment: category}}},
											'deduplicated': False,
											'last alignment': None
										}
			
			pos_data['deduplicated'] = True
		
		# garbage collection
		if pos_data['last alignment'] is alignment:	del self.pos_tracker[alignment.is_reverse][start_pos]

		return alignment

	def tracker_is_empty(self): # verify that any remaining positions are only "already processed" reads that never appeared
		for tracker in self.pos_tracker:
			for pos in tracker.values():
				if pos['alignments by mate start'] or pos['last alignment']: return False
		return True

	def get_marked_alignment(self):
		for alignment in self.alignments:
			self.counts['alignment'] += 1
			if not ( # verify sorting
				(not self.alignment_buffer) or (
					alignment.reference_id > self.current_reference_id or (
						alignment.reference_id == self.current_reference_id and
						alignment.reference_start >= self.most_recent_left_pos
					)
				)
			): raise RuntimeError('alignment %s out of order: verify sorting' % alignment.query_name)
			if not parse_sam.alignment_is_good(alignment): continue
			if alignment.is_paired and not parse_sam.alignment_is_properly_paired(alignment): continue
			umi = umi_data.get_umi(alignment.query_name, self.truncate_umi)
			if not umi_data.umi_is_good(umi): continue
			alignment.is_duplicate = False # not sure how to handle alignments that have already been deduplicated somehow, so just ignore previous annotations
			start_pos = parse_sam.get_start_pos(alignment)
			mate_start_pos = parse_sam.get_mate_start_pos(alignment)
			self.counts['usable alignment'] += 1

			# advance the buffer
			while self.alignment_buffer and (
				self.current_reference_id < alignment.reference_id or # new chromosome
				parse_sam.get_start_pos(self.alignment_buffer[0]) < alignment.reference_start # oldest buffer member is now guaranteed not to get any more hits at its position
			): yield self.pop_buffer()
			if alignment.reference_id > self.current_reference_id: # clear the tracker
				assert self.tracker_is_empty()
				self.pos_tracker = ({}, {})
				self.current_reference_id = alignment.reference_id

			# add the alignment to the tracking data structures
			self.alignment_buffer.extend([alignment])
			try:
				self.pos_tracker[alignment.is_reverse][start_pos]['alignments by mate start'][mate_start_pos] += [alignment]
			except KeyError: # first time we've seen this strand+position
				self.pos_tracker[alignment.is_reverse][start_pos] = {
					'alignments by mate start': collections.defaultdict(list, [(mate_start_pos, [alignment])]),
					'already processed': collections.defaultdict(dict),
					'deduplicated': False
				}
			self.pos_tracker[alignment.is_reverse][start_pos]['last alignment'] = alignment
			self.most_recent_left_pos = alignment.reference_start

		# flush the buffer
		while self.alignment_buffer: yield self.pop_buffer()
		
		assert self.tracker_is_empty()
		assert self.counts['usable alignment'] == sum(self.counts[x] for x in ['distinct', 'optical duplicate', 'PCR duplicate', 'UMI rescued', 'algorithm rescued'])

