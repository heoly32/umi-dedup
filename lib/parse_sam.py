import collections, warnings, copy, pysam

# much of this assumes Illumina-formatted read names


def alignment_is_good (alignment):
	return not (alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary) # only count primary alignments, and discard unmapped reads since it's too expensive (and useless?) to find their duplicates

def alignment_is_properly_paired (alignment): # return True if pair is in expected form, False if unexpected pairing or just unpaired
	return (
		alignment.template_length != 0 and
		alignment.reference_id == alignment.next_reference_id and
		(
			(not alignment.is_reverse) and
			alignment.mate_is_reverse and
			alignment.next_reference_start >= alignment.reference_start
		) or (
			alignment.is_reverse and
			(not alignment.mate_is_reverse) and
			alignment.next_reference_start <= alignment.reference_start and
			alignment.next_reference_start == alignment.reference_end + alignment.template_length # template_length should be negative
		)
	)

def get_start_pos (alignment): # alignment start position in 0-based counting
	return alignment.reference_end - 1 if alignment.is_reverse else alignment.reference_start

def get_mate_start_pos (alignment): # paired alignment's start position in 0-based counting (returns None if unpaired)
	if not alignment_is_properly_paired(alignment): return None
	# streamlined by the assumption of proper pairing
	return (alignment.reference_start + alignment.template_length - 1) if not alignment.is_reverse else alignment.next_reference_start

def is_first_mate (alignment): # check if this is the first or second in the pair; returns False if unpaired
	return alignment.is_paired and get_mate_start_pos(alignment) > get_start_pos(alignment)	

def get_quality_sum (alignment):
	return sum(alignment.query_qualities)

Coords = collections.namedtuple('Coords', ['tile', 'x', 'y']) # "tile" actually also contains machine name, flow cell name, lane, etc. - enough to specify this tile all in one unique string

def get_cluster_coords (alignment):
	tile, x_string, y_string = alignment.query_name.rsplit(':', 3)[:3]
	return Coords(tile, int(x_string), int(y_string))

def get_umi (read_name, truncate = None):
	for label in read_name.split(' ')[:2]: # to allow NCBI format or regular Illumina
		if label.count(':') in (5, 7): # Casava pre-1.8: should be 5 (4 + the UMI hack); Casava 1.8+ / bcl2fastq 2.17+: should be 7 (with optional UMI field)
			umi = label.partition('#')[0].partition('/')[0].rpartition(':')[2] # don't include the space or # and the stuff after it, if present
			return (umi if truncate is None else umi[:truncate + umi.count(DEFAULT_SEPARATOR)]) # don't count the pair separator when truncating
	# only get here if nothing was found
	raise RuntimeError('read name %s does not contain UMI in expected Casava/bcl2fastq format' % label)

class ParsedAlignment:
	'''
	container for useful information parsed out of a pysam.AlignedSegment
	'''
	def __init__ (self, alignment, truncate = None):
		self.is_good = 							alignment_is_good(alignment)
		self.is_properly_paired = 	alignment_is_properly_paired(alignment)
		self.start_pos = 						get_start_pos(alignment)
		self.mate_start_pos = 			get_mate_start_pos(alignment)
		self.is_first_mate = 				is_first_mate(alignment)
		self.quality_sum = 					get_quality_sum(alignment)
		self.cluster_coords = 			get_cluster_coords(alignment)
		self.umi = 									get_umi(alignment.query_name, truncate)
		self.is_reverse = 					alignment.is_reverse
		self.is_paired = 						alignment.is_paired
		self.is_duplicate = 				alignment.is_duplicate
		self.left_pos = 						alignment.reference_start
		self.reference_id = 				alignment.reference_id
		self.name = 								alignment.query_name
	
	def unparse (self, alignment):
		alignment.set_tag('MI', self.umi)
		alignment.is_duplicate = self.is_duplicate
		return alignment

