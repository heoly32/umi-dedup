import collections
from . import parse_sam
import numpy as np
from scipy.spatial.distance import pdist, squareform

DEFAULT_DIST = 100

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
			# Calculate all distance of points in tile reads
			# Create coordinates matrix in numpy array
			coord_mat = np.array([[tile_reads[i][1].x, tile_reads[i][1].y] for i in range(len(tile_reads))])
			# Calculate all distance of points in matrix using the squared Euclidean distance
			dist = squareform(pdist(coord_mat, 'sqeuclidean')) # (x1 - x2) ** 2 + (y1 - y2) ** 2
			# Set values at (0, 0), (1,1), (2,2) .. which have zero to max_dist + 1
			for i in range(len(tile_reads)):
				dist[i][i] = max_dist + 1
			# Compare if there are values which equal to or smaller than max_dist
			dist_logi = dist <= max_dist

			which_group = list(range(len(tile_reads))) # directory of which group each read is in; initially each read is in its own group
			groups = [[read[0]] for read in tile_reads] # list of all the groups
			for i in range(len(tile_reads) - 1):
				if not dist_logi[i].any(): # If there is no True, skip i
					continue

				for j in range(i + 1, len(tile_reads)):
					if dist_logi[i][j]: # See if distance between ith and jth is equal to or smaller than max_dist
						groups[which_group[i]] += groups[which_group[j]] # move the second read's entire group into the first read's group
						groups[which_group[j]] = [] # delete the second read's group so it won't be duplicated
						which_group[j] = which_group[i] # point the second read to the first read's group
			for group in groups:
				if len(group) > 1: result += [group] # return only the groups with multiple elements

	return result

