import collections
	
def are_optical_duplicates (coords1, coords2, max_dist):
	return (coords1.tile == coords2.tile and (coords1.x - coords2.x) ** 2 + (coords1.y - coords2.y) ** 2 <= max_dist)

def which_optical_duplicates (coord_list, max_dist):
	is_optical_duplicate = [False] * len(coord_list)
	for i in range(len(coord_list) - 1):
		for j in range(i + 1, len(coord_list)):
			if are_optical_duplicates(coord_list[i], coord_list[j], max_dist):
				is_optical_duplicate[i] = is_optical_duplicate[j] = True
	return is_optical_duplicate

