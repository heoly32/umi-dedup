import collections, parse_sam
	
def are_optical_duplicates (coords1, coords2, max_dist):
	return (coords1.tile == coords2.tile and (coords1.x - coords2.x) ** 2 + (coords1.y - coords2.y) ** 2 <= max_dist)

def get_optical_duplicates (reads, max_dist): # return a list where each element is a list of the reads that are optical duplicates of each other (each read is within max_dist of at least one other read in the set)
	coords = map(parse_sam.get_coords, reads)
	which_group = list(range(len(reads))) # directory of which group each read is in; initially each read is in its own group
	groups = [[read] for read in reads] # list of all the groups
	for i in range(len(reads) - 1):
		for j in range(i + 1, len(reads)):
			if are_optical_duplicates(coords[i], coords[j], max_dist):
				groups[which_group[i]] += groups[which_group[j]] # move the second read's entire group into the first read's group
				groups[which_group[j]] = [] # delete the second read's group so it won't be duplicated
				which_group[j] = which_group[i] # point the second read to the first read's group
	return [group for group in groups if len(group) > 1] # return only the groups with multiple elements

