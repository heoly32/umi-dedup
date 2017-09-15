import collections, warnings, copy, pysam

MAX_READ_LENGTH = 200 # maximum expected read length, to know when you're beyond soft-clipping range (set this high to be cautious)

def alignment_is_good (alignment):
	return not (alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary) # only count primary alignments, and discard unmapped reads since it's too expensive (and useless?) to find their duplicates

def alignment_is_properly_paired (alignment): # return True if pair is in expected form, False if unexpected pairing or just unpaired
	if not alignment.is_paired: return False
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

def get_start_pos (alignment): # alignment start position in 0-based counting; do some arithmetic to ignore soft-clipping
	return alignment.reference_end + (alignment.query_length - alignment.query_alignment_end) - 1 if alignment.is_reverse else alignment.reference_start - alignment.query_alignment_start

def get_mate_start_pos (alignment): # paired alignment's start position in 0-based counting (returns None if unpaired)
	if not alignment_is_properly_paired(alignment): return None
	# streamlined by the assumption of proper pairing
	return (alignment.reference_start + alignment.template_length - 1) if not alignment.is_reverse else alignment.next_reference_start

def is_first_mate (alignment): # check if this is the first or second in the pair; returns False if unpaired
	return alignment.is_paired and get_mate_start_pos(alignment) > get_start_pos(alignment)	

def get_quality (alignment):
	return sum(alignment.query_qualities)

Coords = collections.namedtuple('Coords', ['tile', 'x', 'y']) # "tile" actually also contains machine name, flow cell name, lane, etc. - enough to specify this tile all in one unique string

def get_coords (alignment):
	tile, x_string, y_string = alignment.query_name.rsplit(':', 3)[:3]
	return Coords(tile, int(x_string), int(y_string))

