import collections

def alignment_is_good (alignment):
	if alignment.is_paired or alignment.is_read2: raise RuntimeError('paired-end reads are currently unsupported')
	return not (alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary) # only count primary alignments, and discard unmapped reads since it's too expensive (and useless?) to find their duplicates

def get_start_pos (alignment): # alignment start position in 0-based counting
	return alignment.reference_end - 1 if alignment.is_reverse else alignment.reference_start

def get_quality (alignment):
	return sum(alignment.query_qualities)

Coords = collections.namedtuple('Coords', ['tile', 'x', 'y']) # "tile" actually also contains machine name, flow cell name, lane, etc. - enough to specify this tile all in one unique string

def get_coords (alignment):
	tile, x_string, y_string = alignment.query_name.rsplit(':', 3)[:3]
	return Coords(tile, int(x_string), int(y_string))

