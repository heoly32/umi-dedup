import collections, copy, multiprocessing
from . import parse_sam, umi_data, optical_duplicates, naive_estimate, bayes_estimate, sequence_error, library_stats

# Initiate sequence correction functor
sequence_correcter = sequence_error.ClusterAndReducer()

PROCESSES = 4
DUP_CATEGORIES = ['optical duplicate', 'PCR duplicate']

class PosTracker:
	'''
	data structure for tracking alignments at a certain start position
	meant to be used internally in DuplicateMarker
	'''

	__slots__ = ['alignments_by_mate', 'alignments_already_processed', 'last_alignment', 'deduplicated']

	def __init__(self,
		alignments_by_mate = None,
		alignments_already_processed = None,
		last_alignment = None,
		deduplicated = False
	):
		self.alignments_by_mate = alignments_by_mate if alignments_by_mate is not None else collections.defaultdict(list)
		self.alignments_already_processed = alignments_already_processed if alignments_already_processed is not None else collections.defaultdict(lambda: collections.defaultdict(dict))
		self.last_alignment = last_alignment
		self.deduplicated = deduplicated

	def __repr__(self):
		return '%s(alignments_by_mate = %s, alignments_already_processed = %s, last_alignment = %s, deduplicated = %s)' % (self.__class__, self.alignments_by_mate, self.alignments_already_processed, self.last_alignment, self.deduplicated)

def dedup_pos(pos_data, umi_dup_function, optical_dist = 0, sequence_correction = None):
	# trackers for summary statistics
	category_counts = collections.Counter()
	pos_counts = {'before': [], 'after': []}
	for mate_start_pos, alignments_with_this_mate in pos_data.alignments_by_mate.iteritems(): # iterate over mate start positions
		alignments_to_dedup = copy.copy(alignments_with_this_mate)

		# if mate is already deduplicated, use those results
		if mate_start_pos in pos_data.alignments_already_processed:
			umi_counts = collections.Counter()
			umi_nondup_counts = collections.Counter()
			for already_processed_alignment in alignments_with_this_mate:
				umi_counts[already_processed_alignment.umi] += 1
				try:
					category = pos_data.alignments_already_processed[mate_start_pos][already_processed_alignment.umi][already_processed_alignment.name]
					if category in DUP_CATEGORIES:
						already_processed_alignment.is_duplicate = True
						category_counts[category] += 1
					else:
						already_processed_alignment.is_duplicate = False
						umi_nondup_counts[already_processed_alignment.umi] += 1
					alignments_to_dedup.remove(already_processed_alignment)
				except KeyError: # this alignment is missing from the list, so its mate must have been missing from the data, so it still needs deduplication
					pass
			for umi_nondup_count in umi_nondup_counts.values():
				category_counts['UMI rescued'] += 1
				category_counts['algorithm rescued'] += umi_nondup_count - 1

			category_counts['UMI rescued'] -= bool(umi_nondup_counts) # count the first read at this position as distinct
			category_counts['distinct'] += bool(umi_nondup_counts)
			pos_counts['before'].append(umi_counts.values())
			pos_counts['after'].append(umi_nondup_counts.values())

		if alignments_to_dedup: # false if they had all been deduplicated already
			alignment_categories = {} # key = alignment query_name, value = which category it is (from DUP_CATEGORIES or otherwise)
			alignments_by_umi = collections.defaultdict(list)
			for this_alignment in alignments_to_dedup: alignments_by_umi[this_alignment.umi] += [this_alignment]
			count_by_umi = umi_data.UmiValues([(umi, len(hits)) for umi, hits in alignments_by_umi.iteritems()])
			pos_counts['before'].append(count_by_umi.nonzero_values())

			# first pass: mark optical duplicates
			if optical_dist != 0:
				for opt_dups in optical_duplicates.get_optical_duplicates(alignments_to_dedup, optical_dist):
					for dup_alignment in umi_data.mark_duplicates(opt_dups, len(opt_dups) - 1):
						# remove duplicate reads from the tracker so they won't be considered later (they're still in the read buffer)
						if dup_alignment.is_duplicate:
							alignments_by_umi[dup_alignment.umi].remove(dup_alignment)
							if len(alignments_by_umi[dup_alignment.umi]) == 0: del alignments_by_umi[dup_alignment.umi]
							count_by_umi[dup_alignment.umi] -= 1
							alignment_categories[dup_alignment.name] = 'optical duplicate'
							category_counts['optical duplicate'] += 1

			# second pass: mark PCR duplicates
			dedup_counts = umi_dup_function(count_by_umi)
			pos_counts['after'].append(dedup_counts.nonzero_values())
			if sequence_correction is not None:
				pre_correction_dict = {umi: len(hits) for umi, hits in alignments_by_umi.iteritems()}
				pre_correction_count = sum(map(len, alignments_by_umi.values()))
				alignments_with_new_umi, first_clusters, second_clusters = sequence_correcter(alignments_by_umi)
				obsolete_umis = set()
				for alignment, umi in alignments_with_new_umi:
					obsolete_umis.add(alignment.umi)
					alignments_by_umi[umi].append(alignment)
					category_counts['sequence correction'] += 1
				for umi in obsolete_umis:
					del alignments_by_umi[umi]
				post_correction_count = sum(map(len, alignments_by_umi.values()))
				assert pre_correction_count == post_correction_count
			for umi, alignments_with_this_umi in alignments_by_umi.iteritems():
				dedup_count = dedup_counts[umi]
				assert alignments_with_this_umi and dedup_count
				n_dup = len(alignments_with_this_umi) - dedup_count
				for marked_alignment in umi_data.mark_duplicates(alignments_with_this_umi, n_dup):
					alignment_categories[marked_alignment.name] = ('PCR duplicate' if marked_alignment.is_duplicate else 'nonduplicate')
				category_counts['PCR duplicate'] += n_dup
				category_counts['UMI rescued'] += 1
				category_counts['algorithm rescued'] += dedup_count - 1
			category_counts['UMI rescued'] -= bool(alignments_by_umi) # count the first read at this position as distinct
			category_counts['distinct'] += bool(alignments_by_umi)

#				# pass duplicate marking to mates DANGER DANGER THIS IS TEMPORARILY DISABLED
#				if mate_start_pos is not None:
#					for categorized_alignment, category in alignment_categories.iteritems():
#						self.pos_tracker[not alignment.is_reverse][mate_start_pos].alignments_already_processed[start_pos][categorized_alignment.umi][categorized_alignment] = category

	pos_data.deduplicated = True
	return (pos_data, category_counts)


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
		nburn = bayes_estimate.DEFAULT_NBURN,
		alpha2 = bayes_estimate.DEFAULT_ALPHA2,
		prior = None,
		filter_counts = True,
		sequence_correction = None
	):
		self.alignment_source = alignments
		self.umi_frequency = umi_frequency
		self.optical_dist = optical_dist
		self.truncate_umi = truncate_umi
		self.sequence_correction = sequence_correction
		if algorithm == 'naive':
			self.umi_dup_function = naive_estimate.deduplicate_counts
		elif algorithm in ('bayes', 'uniform-bayes'):
			self.umi_dup_function = lambda counts: bayes_estimate.deduplicate_counts(umi_counts = counts, nsamp = nsamp, nthin = nthin, nburn = nburn, uniform = (algorithm == 'uniform-bayes'),  alpha2 = alpha2, total_counts = self.umi_frequency, prior = prior, filter_counts = filter_counts)
		else:
			raise NotImplementedError
		self.raw_alignments = {}
		self.alignment_buffer = collections.deque()
		self.pos_tracker = (collections.defaultdict(PosTracker), collections.defaultdict(PosTracker)) # data structure containing alignments by position rather than sort order; top level is by strand (0 = forward, 1 = reverse), then next level is by 5' read start position (dict since these will be sparse and are only looked up by identity), then next level is by 5' start position of mate read, and each element of that contains a variety of data; there is no level for reference ID because there is no reason to store more than one chromosome at a time
		self.output_generator = self.get_marked_alignment()
		self.current_reference_id = 0
		self.most_recent_left_pos = 0
		self.pool = multiprocessing.Pool(processes = PROCESSES)

		# Initiate sequence correction functor
		if sequence_correction is not None:
			sequence_correcter = sequence_error.ClusterAndReducer(sequence_correction)

		# trackers for summary statistics
		self.category_counts = collections.Counter()
		self.pos_counts = {'before': [], 'after': []} # position hit counts before or afterdeduplication


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
		pos_data = self.pos_tracker[alignment.is_reverse][alignment.start_pos]

		# deduplicate reads
		if not pos_data.deduplicated:
			pos_data, category_counts = dedup_pos(pos_data, self.umi_dup_function, self.optical_dist, self.sequence_correction)
			self.pos_tracker[alignment.is_reverse][alignment.start_pos] = pos_data
			self.category_counts.update(category_counts)

		# garbage collection
		if pos_data.last_alignment is alignment:	del self.pos_tracker[alignment.is_reverse][alignment.start_pos]
		self.raw_alignments[alignment.name].is_duplicate = alignment.is_duplicate
		result = self.raw_alignments[alignment.name]
		del self.raw_alignments[alignment.name]
		
		return result

	def tracker_is_empty(self): # verify that any remaining positions are only "already processed" reads that never appeared
		for tracker in self.pos_tracker:
			for pos in tracker.values():
				if pos.alignments_by_mate or pos.last_alignment: return False
		return True

	def get_marked_alignment(self):
		for raw_alignment in self.alignment_source:
			alignment = parse_sam.ParsedAlignment(raw_alignment, self.truncate_umi)
			self.category_counts['alignment'] += 1
			if not ( # verify sorting
				(not self.alignment_buffer) or (
					alignment.reference_id > self.current_reference_id or (
						alignment.reference_id == self.current_reference_id and
						alignment.left_pos >= self.most_recent_left_pos
					)
				)
			): raise RuntimeError('alignment %s out of order: verify sorting' % alignment.query_name)
			if not alignment.is_good: continue
			if alignment.is_paired and not alignment.is_properly_paired: continue
			if not umi_data.umi_is_good(alignment.umi): continue
			alignment.is_duplicate = False # not sure how to handle alignments that have already been deduplicated somehow, so just ignore previous annotations
			self.category_counts['usable alignment'] += 1
			
			# advance the buffer
			while self.alignment_buffer and (
				self.current_reference_id < alignment.reference_id or # new chromosome
				alignment.start_pos < alignment.left_pos # oldest buffer member is now guaranteed not to get any more hits at its position
			): yield self.pop_buffer()
			if alignment.reference_id > self.current_reference_id: # clear the tracker
				assert self.tracker_is_empty()
				self.pos_tracker = (collections.defaultdict(PosTracker), collections.defaultdict(PosTracker))
				self.current_reference_id = alignment.reference_id

			# add the alignment to the tracking data structures
			assert not alignment.name in self.raw_alignments
			self.raw_alignments[alignment.name] = raw_alignment
			self.alignment_buffer.extend([alignment])
			self.pos_tracker[alignment.is_reverse][alignment.start_pos].alignments_by_mate[alignment.mate_start_pos] += [alignment]
			self.pos_tracker[alignment.is_reverse][alignment.start_pos].last_alignment = alignment
			self.most_recent_left_pos = alignment.left_pos

		# flush the buffer
		while self.alignment_buffer: yield self.pop_buffer()

		assert self.tracker_is_empty()
		assert self.category_counts['usable alignment'] == sum(self.category_counts[x] for x in ['distinct', 'optical duplicate', 'PCR duplicate', 'UMI rescued', 'algorithm rescued'])

	def get_mean_pos_entropy(self, which): # which: 'before' or 'after' deduplication
		return library_stats.mean(library_stats.entropy(pos) for pos in self.pos_counts[which])

	def get_library_entropy(self, which): # which: 'before' or 'after' deduplication
		return library_stats.entropy(map(sum, self.pos_counts[which]))

	def estimate_library_size(self):
		return library_stats.estimate_library_size(self.category_counts['distinct'] + self.category_counts['UMI rescued'] + self.category_counts['algorithm rescued'], self.category_counts['usable alignment'])

