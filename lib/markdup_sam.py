import collections, copy, multiprocessing
import time
from . import parse_sam, umi_data, optical_duplicates, naive_estimate, bayes_estimate, sequence_error, library_stats

test = False # test


DUP_CATEGORIES = ['optical duplicate', 'PCR duplicate']

class PosTracker:
	'''
	data structure for tracking alignments at a certain start position
	meant to be used internally in DuplicateMarker
	'''

	__slots__ = ['alignments_by_mate', 'alignments_already_processed', 'last_alignment_name']

	def __init__(self,
		alignments_by_mate =            collections.defaultdict(list),
		alignments_already_processed =  {},
		last_alignment_name =           None,
	):
#		self.alignments_already_processed = alignments_already_processed if alignments_already_processed is not None else collections.defaultdict(lambda: collections.defaultdict(dict))
		self.alignments_by_mate =            alignments_by_mate
		self.alignments_already_processed =  alignments_already_processed
		self.last_alignment_name =           last_alignment_name
	
	def __repr__(self):
		return '%s(alignments_by_mate = %s, alignments_already_processed = %s, last_alignment_name = %s)' % (self.__class__, self.alignments_by_mate, self.alignments_already_processed, self.last_alignment_name)

def dedup_counts(counts, algorithm, *args, **kwargs):
	if algorithm == 'naive':
		return naive_estimate.deduplicate_counts(counts)
	elif algorithm in ('bayes', 'uniform-bayes'):
		return bayes_estimate.deduplicate_counts(umi_counts = counts, uniform = (algorithm == 'uniform-bayes'), *args, **kwargs)
	else:
		raise NotImplementedError

def dedup_pos(pos_data, sequence_corrector = None, optical_dist = 0, *dedup_args, **dedup_kwargs):
	# trackers for summary statistics
	category_counts = collections.Counter()
	pos_counts = {'before': [], 'after': []}
	for mate_start_pos, alignments_with_this_mate in pos_data.alignments_by_mate.items(): # iterate over mate start positions
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

			category_counts['UMI rescued'] -= bool(umi_nondup_counts) # count the first alignment at this position as distinct
			category_counts['distinct'] += bool(umi_nondup_counts)
			pos_counts['before'].append(umi_counts.values())
			pos_counts['after'].append(umi_nondup_counts.values())

		if alignments_to_dedup: # false if they had all been deduplicated already
			alignment_categories = {} # key = alignment query_name, value = which category it is (from DUP_CATEGORIES or otherwise)
			alignments_by_umi = collections.defaultdict(list)
			for this_alignment in alignments_to_dedup: alignments_by_umi[this_alignment.umi] += [this_alignment]
			count_by_umi = umi_data.UmiValues([(umi, len(hits)) for umi, hits in alignments_by_umi.items()])
			pos_counts['before'].append(count_by_umi.nonzero_values())

			# first pass: mark optical duplicates
			if optical_dist != 0:
				for opt_dups in optical_duplicates.get_optical_duplicates(alignments_to_dedup, optical_dist):
					for dup_alignment in umi_data.mark_duplicates(opt_dups, len(opt_dups) - 1):
						# remove duplicate alignments from the tracker so they won't be considered later (they're still in the alignment buffer)
						if dup_alignment.is_duplicate:
							alignments_by_umi[dup_alignment.umi].remove(dup_alignment)
							if len(alignments_by_umi[dup_alignment.umi]) == 0: del alignments_by_umi[dup_alignment.umi]
							count_by_umi[dup_alignment.umi] -= 1
							alignment_categories[dup_alignment.name] = 'optical duplicate'
							category_counts['optical duplicate'] += 1

			# second pass: mark PCR duplicates
			dedupped_counts = dedup_counts(count_by_umi, *dedup_args, **dedup_kwargs)
			pos_counts['after'].append(dedupped_counts.nonzero_values())
			if sequence_corrector is not None:
				pre_correction_count = sum(map(len, alignments_by_umi.values()))
				alignments_with_new_umi = sequence_corrector(alignments_by_umi)
				try:
					for alignment, umi in alignments_with_new_umi:
						alignments_by_umi[umi].append(alignment)
						category_counts['sequence correction'] += 1
						try:
							del alignments_by_umi[alignment.umi]
						except KeyError:
							pass
				except TypeError:
					pass
				post_correction_count = sum(map(len, alignments_by_umi.values()))
				assert pre_correction_count == post_correction_count
				post_correction_count = sum(map(len, alignments_by_umi.values()))
				assert pre_correction_count == post_correction_count
			for umi, alignments_with_this_umi in alignments_by_umi.items():
				dedup_count = dedupped_counts[umi]
				assert alignments_with_this_umi and dedup_count
				n_dup = len(alignments_with_this_umi) - dedup_count
				for marked_alignment in umi_data.mark_duplicates(alignments_with_this_umi, n_dup):
					alignment_categories[marked_alignment.name] = ('PCR duplicate' if marked_alignment.is_duplicate else 'nonduplicate')
				category_counts['PCR duplicate'] += n_dup
				category_counts['UMI rescued'] += 1
				category_counts['algorithm rescued'] += dedup_count - 1
			category_counts['UMI rescued'] -= bool(alignments_by_umi) # count the first alignment at this position as distinct
			category_counts['distinct'] += bool(alignments_by_umi)

#				# pass duplicate marking to mates DANGER DANGER THIS IS TEMPORARILY DISABLED
#				if mate_start_pos is not None:
#					for categorized_alignment, category in alignment_categories.iteritems():
#						self.pos_tracker[not alignment.is_reverse][mate_start_pos].alignments_already_processed[start_pos][categorized_alignment.umi][categorized_alignment] = category

	return (pos_data, category_counts)

def dedup_worker(queue_to_dedup, queue_dedupped, sequence_correction, *args, **kwargs):
	sequence_corrector = (None if sequence_correction is None else sequence_error.ClusterAndReducer(sequence_correction))
	while True:
		is_reverse, pos, pos_data = queue_to_dedup.get()
		if test: print('\t\t\tworking on %i' % pos)
		start_time = time.clock()
		new_pos_data = dedup_pos(pos_data, sequence_corrector, *args, **kwargs)
		queue_dedupped.put((is_reverse, pos, new_pos_data))
		queue_to_dedup.task_done()
		if test: print('\t\t\tfinished %i in %f s' % (pos, time.clock() - start_time))

class DuplicateMarker:
	'''
	a generator of duplicate-marked alignments (like pysam.AlignedSegment objects), given an iterable of alignments (like a pysam.AlignmentFile) as input
	assumes the input iterable is sorted by leftmost coordinate (SAMtools style) and yields alignments in the same order
	does not yield unusable alignments
	central concept: cheat a little, detecting duplicate alignments by the fact that they align to the same (start) position rather than by their sequences
	implementation: as you traverse the coordinate-sorted input alignments, add each alignment to both a FIFO buffer (so they can be output in the same order) and a dictionary that groups alignments by strand, start position, and mate start position (None if single-end); this is not inefficient because both data structures contain pointers to the same alignment objects
	at the strand+position level (before mate start position), track metadata: whether this strand+position has been deduplicated and which alignment was added to it last (when this alignment is processed, the tracker can be deleted safely)
	then, after each input alignment, check the left end of the buffer to tell whether the oldest alignment is in a position that will never accumulate any more hits (because the input is sorted by coordinate); if so, estimate the duplication at that position, mark all the alignments there accordingly, then output the alignment with the appropriate marking
	'''
	def __init__(self,
		alignments,
		truncate_umi = None,
		processes = (1 if test else multiprocessing.cpu_count()), # test
		*dedup_args,
		**dedup_kwargs
	):
		self.truncate_umi =              truncate_umi
		self.dedup_args =                dedup_args
		self.dedup_kwargs =              dedup_kwargs
		self.alignment_source =          alignments
		self.raw_alignments =            {}
		self.alignment_buffer =          collections.deque()
		self.pos_tracker_deduplicated =  (collections.defaultdict(PosTracker), collections.defaultdict(PosTracker)) # data structure containing alignments by position rather than sort order; top level is by strand (0 = forward, 1 = reverse), then next level is by 5' alignment start position (dict since these will be sparse and are only looked up by identity), then next level is by 5' start position of mate alignment, and each element of that contains a variety of data; there is no level for reference ID because there is no reason to store more than one chromosome at a time; this one only contains positions that have completed deduplication
		self.pos_tracker_to_populate =   (collections.OrderedDict(), collections.OrderedDict()) # like self.pos_tracker_deduplicated except it contains positions that are still being populated with new alignments and haven't been deduplicated yet; uses OrderedDict so they can be deduplicated in order of genome position
		self.output_generator =          self.get_marked_alignment()
		self.queue_to_dedup =            multiprocessing.JoinableQueue()
		self.queue_dedupped =            multiprocessing.JoinableQueue()
		for i in range(processes):
			p = multiprocessing.Process(target = dedup_worker, args = [self.queue_to_dedup, self.queue_dedupped] + list(dedup_args), kwargs = dedup_kwargs)
			p.daemon = True
			p.start()
		self.current_reference_id = 0
		self.most_recent_left_pos = 0

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
		pos_data = self.pos_tracker_deduplicated[alignment.is_reverse][alignment.start_pos] # KeyError if it hasn't been deduplicated yet

		# garbage collection
		if pos_data.last_alignment_name == alignment.name: del self.pos_tracker_deduplicated[alignment.is_reverse][alignment.start_pos]
		self.raw_alignments[alignment.name].is_duplicate = alignment.is_duplicate
		result = self.raw_alignments[alignment.name]
		del self.raw_alignments[alignment.name]
		
		if test: print('\toutput %i:left%i start:%i %s' % (alignment.reference_id, alignment.left_pos, alignment.start_pos, alignment.name)) # test
		return alignment.unparse(result)
	
	def enqueue_positions(self, new_reference_id, new_pos):
		# enqueue position data eligible for deduplication
		if test: print('\tadvance to ref%i:left%i' % (new_reference_id, new_pos)) # test
		for is_reverse in (0, 1):
			while self.pos_tracker_to_populate[is_reverse] and (				
				self.current_reference_id < new_reference_id or # entire chromosome is done
				list(self.pos_tracker_to_populate[is_reverse].keys())[0] < new_pos # position has been passed and is now guaranteed not to get any more hits
			):
				pos = list(self.pos_tracker_to_populate[is_reverse].keys())[0]
				if max(map(len, self.pos_tracker_to_populate[is_reverse][pos].alignments_by_mate.values())) == 1: # singletons skip the queue
					self.pos_tracker_deduplicated[is_reverse][pos] = self.pos_tracker_to_populate[is_reverse][pos]
					self.category_counts['distinct'] += len(self.pos_tracker_to_populate[is_reverse][pos].alignments_by_mate)
					self.pos_counts['before'] += [[1]] * len(self.pos_tracker_to_populate[is_reverse][pos].alignments_by_mate)
					self.pos_counts['after'] += [[1]] * len(self.pos_tracker_to_populate[is_reverse][pos].alignments_by_mate)
				else:
					self.queue_to_dedup.put((is_reverse, pos, self.pos_tracker_to_populate[is_reverse][pos]))
					if test: print('\t\tenqueue start%i' % pos) # test
				del self.pos_tracker_to_populate[is_reverse][pos]
	
	def dedup_until(self, new_reference_id, new_pos):
		self.enqueue_positions(new_reference_id, new_pos)
		
		# update the tracking data with any positions that have finished deduplication
		while (
			self.queue_dedupped._unfinished_tasks.get_value() > 0
			) or(
				self.current_reference_id < new_reference_id and self.queue_to_dedup._unfinished_tasks.get_value() > 0 # if this is the end of the chromosome, continue until the entire queue is processed
			):
			if test: print('\t\t%i tasks in work queue, %i in finished queue' % (self.queue_to_dedup._unfinished_tasks.get_value(), self.queue_dedupped._unfinished_tasks.get_value()))
			is_reverse, pos, (pos_data, category_counts) = self.queue_dedupped.get() # this will wait until self.queue_dedupped is no longer empty, if it was empty in this cycle but we got here anyway at the end of the chromosome
			self.pos_tracker_deduplicated[is_reverse][pos] = pos_data
			self.category_counts.update(category_counts)
			self.queue_dedupped.task_done()
			if test: print('\t\tfinished start%i' % pos) # test
		if test: print('\t\t%i tasks in work queue, %i in finished queue' % (self.queue_to_dedup._unfinished_tasks.get_value(), self.queue_dedupped._unfinished_tasks.get_value()))
	
	def tracker_is_empty(self): # verify that any remaining positions are only "already processed" alignments that never appeared
		for tracker in (self.pos_tracker_to_populate, self.pos_tracker_deduplicated):
			for strand_tracker in tracker:
				for pos_data in strand_tracker.values():
					if pos_data.alignments_by_mate or pos_data.last_alignment_name:
						return False
		return True
	
	def add_alignment(self, alignment, raw_alignment):
		# add the alignment to the tracking data structures
		assert not alignment.name in self.raw_alignments
		self.raw_alignments[alignment.name] = raw_alignment
		self.alignment_buffer.extend([alignment])
		try: # this position is already in the tracker
			self.pos_tracker_to_populate[alignment.is_reverse][alignment.start_pos].alignments_by_mate[alignment.mate_start_pos].append(alignment)
			self.pos_tracker_to_populate[alignment.is_reverse][alignment.start_pos].last_alignment_name = alignment.name
		except KeyError: # first time seeing this position
			# first pop higher positions from the tracker, to keep it ordered if there are any
			n_to_pop = 0
			for pos in reversed(self.pos_tracker_to_populate[alignment.is_reverse]):
				if pos > alignment.start_pos:
					assert alignment.is_reverse # this should never happen on the forward strand, where start position is left position so they're already in order
					n_to_pop += 1
				else:
					break
			popped = []
			for i in range(n_to_pop): popped.append(self.pos_tracker_to_populate[alignment.is_reverse].popitem())
			
			# now add the new position
			self.pos_tracker_to_populate[alignment.is_reverse][alignment.start_pos] = PosTracker(
				alignments_by_mate = {alignment.mate_start_pos: [alignment]},
				last_alignment_name = alignment.name
			)
			
			# finally, put the popped positions back on
			for pos, pos_data in popped: self.pos_tracker_to_populate[alignment.is_reverse][pos] = pos_data
		
		self.most_recent_left_pos = alignment.left_pos
		if test: print('\tadd %i:start%i' % (alignment.reference_id, alignment.start_pos)) # test
	
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
			if test: print('read %i:left%i start%i %s' % (alignment.reference_id, alignment.left_pos, alignment.start_pos, alignment.name)) # test
			
			# deduplicate eligible positions
			self.dedup_until(alignment.reference_id, alignment.left_pos)
			
			# advance the alignment buffer as far as possible
			while self.alignment_buffer and self.alignment_buffer[0].start_pos in self.pos_tracker_deduplicated[self.alignment_buffer[0].is_reverse]: yield self.pop_buffer()
			
			# reset the tracker if starting a new chromosome
			if alignment.reference_id > self.current_reference_id:
				assert self.tracker_is_empty()
				self.pos_tracker_deduplicated =  (collections.defaultdict(PosTracker), collections.defaultdict(PosTracker))
				self.pos_tracker_to_populate =   (collections.OrderedDict(), collections.OrderedDict())
				self.current_reference_id =      alignment.reference_id

			self.add_alignment(alignment, raw_alignment)

		# flush the buffer
		self.dedup_until(self.current_reference_id + 1, 0)
		while self.alignment_buffer: yield self.pop_buffer()

		assert self.tracker_is_empty()
		assert self.category_counts['usable alignment'] == sum(self.category_counts[x] for x in ['distinct', 'optical duplicate', 'PCR duplicate', 'UMI rescued', 'algorithm rescued'])

	def get_mean_pos_entropy(self, which): # which: 'before' or 'after' deduplication
		return library_stats.mean(library_stats.entropy(pos) for pos in self.pos_counts[which])

	def get_library_entropy(self, which): # which: 'before' or 'after' deduplication
		return library_stats.entropy(map(sum, self.pos_counts[which]))

	def estimate_library_size(self):
		return library_stats.estimate_library_size(self.category_counts['distinct'] + self.category_counts['UMI rescued'] + self.category_counts['algorithm rescued'], self.category_counts['usable alignment'])
