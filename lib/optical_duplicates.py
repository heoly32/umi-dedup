import collections, parse_sam

DEFAULT_DIST = 100

def are_optical_duplicates (coords1, coords2, max_dist):
	return (coords1.x - coords2.x) ** 2 + (coords1.y - coords2.y) ** 2 <= max_dist

def get_optical_duplicates (reads, max_dist): # return a list where each element is a list of the reads that are optical duplicates of each other (each read is within max_dist of at least one other read in the set)
	coords_by_tile = {} # group reads by which tile they came from; each value is another dictionary, in which the key is the coordinates and the value is the read itself
	for read in reads:
		coords = parse_sam.get_coords(read)
		try:
			coords_by_tile[coords.tile] += [(read, coords)]
		except KeyError:
			coords_by_tile[coords.tile] = [(read, coords)]
	
	result = [] # each element will be a list of reads that are optical duplicates of each other (within max_dist pixels)
	for tile_reads in coords_by_tile.values():
		if len(tile_reads) > 1:
			which_group = list(range(len(tile_reads))) # directory of which group each read is in; initially each read is in its own group
			groups = [[read[0]] for read in tile_reads] # list of all the groups
			for i in range(len(tile_reads) - 1):
				for j in range(i + 1, len(tile_reads)):
					if are_optical_duplicates(tile_reads[i][1], tile_reads[j][1], max_dist):
						groups[which_group[i]] += groups[which_group[j]] # move the second read's entire group into the first read's group
						groups[which_group[j]] = [] # delete the second read's group so it won't be duplicated
						which_group[j] = which_group[i] # point the second read to the first read's group
			for group in groups:
				if len(group) > 1: result += [group] # return only the groups with multiple elements

	return result

