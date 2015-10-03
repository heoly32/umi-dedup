import collections

def read_is_good (read):
	if read.is_paired or read.is_read2: raise RuntimeError('paired-end reads are currently unsupported')
	return not (read.is_unmapped or read.is_secondary or read.is_supplementary) # only count primary alignments, and discard unmapped reads since it's too expensive (and useless?) to find their duplicates

def get_start_pos (read): # read start position in 0-based counting
	return read.reference_end - 1 if read.is_reverse else read.reference_start

def get_quality (read):
	return sum(read.query_qualities)

Coords = collections.namedtuple('Coords', ['tile', 'x', 'y']) # "tile" actually also contains machine name, flow cell name, lane, etc. - enough to specify this tile all in one unique string

def get_coords (read):
	tile, x_string, y_string = read.query_name.rsplit(':', 3)[:3]
	return Coords(tile, int(x_string), int(y_string))

