import itertools, Bio.SeqIO
from . import umi_data

def get_umi (seq, length, before = 0, mask_pos = []):
	'''
	extract UMI from a read sequence
	optionally specify offset (number of bases left of UMI to discard)
	optionally specify mask positions (0-based indices of bases within the UMI to discard)
	'''
	if length == 0: return ''
	
	umi = str(seq[before:(before + length)])
	if mask_pos is not None:
		for pos, next_pos in itertools.izip_longest(mask_pos[::-1], mask_pos[-2::-1], fillvalue = None):
			if next_pos >= pos: raise RuntimeError('mask_pos must be in ascending order')
			umi = umi[:pos] + umi[pos + 1:]
		
	return umi

def add_umi_to_read (read, umi, trim_length = 0):
	'''
	add UMI to Illumina-format FASTQ read name
	read is expected to be a Bio.SeqIO.SeqRecord object
	optionally also remove a specified number of bases from the beginning of the read (both sequence and qualities)
	'''
	
	labels = read.description.split(' ')
	assert labels[0] == read.id # this will be important later
	which_read_name = None
	for i in range(2):
		if labels[i].count(':') in (4, 6): # Illumina format
			which_read_name = i
			break
	if which_read_name is None: raise RuntimeError('unrecognized format in read ' + read.description)
	# protect against Casava version < 1.8 without breaking 1.8+
	part1 = labels[which_read_name].partition('#')
	part2 = part1[0].partition('/')
	labels[which_read_name] = part2[0] + ':' + umi + ''.join(part2[1:] + part1[1:])
	read.description = ' '.join(labels)
	read.id = labels[0] # if this is the part of the name that changed, read.id must also change or the Bio.SeqIO.parse object will give bad output
	
	qualities = read.letter_annotations['phred_quality']
	read.letter_annotations = {} # letter annotations must be emptied before changing sequence
	read.seq = read.seq[trim_length:]
	read.letter_annotations['phred_quality'] = qualities[trim_length:]
	
	return read

def get_read_umis (
	in_file,
	umi_length,
	before = 0,
	after = 0,
	mask_pos = [],
	relabel = True
):
	'''
	returns a generator of pairs of Bio.SeqIO.SeqRecord objects and UMIs extracted from them
	moves UMI from read sequence to read name unless disabled
	defensively sorts mask_pos just in case you're naughty
	'''
	trim_length = before + umi_length + after
	mask_pos = sorted(mask_pos)
	if len(mask_pos) > 0 and mask_pos[-1] > umi_length: raise RuntimeError('UMI is only %i bases; can\'t mask position %i' % (umi_length, mask_pos[-1]))
	for read in Bio.SeqIO.parse(in_file, 'fastq'):
		umi = get_umi(read.seq, umi_length, before, mask_pos)
		yield (
			(add_umi_to_read(read, umi, trim_length) if relabel else read),
			umi
		)

def get_umi_labeled_reads (*args, **kwargs):
	'''
	like get_read_umis but only yield the relabeled reads, not the UMIs too
	'''
	return (read for read, umi in get_read_umis(*args, relabel = True, **kwargs))

def get_read_pair_umis (
	in_file1,
	in_file2,
	umi_length1,
	umi_length2,
	before1 = 0,
	before2 = 0,
	after1 = 0,
	after2 = 0,
	mask_pos1 = [],
	mask_pos2 = [],
	pair_separator = umi_data.DEFAULT_SEPARATOR,
	relabel = True
):
	'''
	returns a generator of pairs of pairs of Bio.SeqIO.SeqRecord objects and their UMIs
	moves UMIs from read sequences to read names unless disabled
	gives both members of the pair the same UMI field, a combination of both UMIs with a separator between them
	or if only one read has the UMI, use 0 as the other one's length, and the field will just contain the single UMI with no separator
	'''
	trim_length1, trim_length2 = before1 + umi_length1 + after1, before2 + umi_length2 + after2
	if umi_length1 == 0 or umi_length2 == 0: pair_separator = ''
	
	for read_and_umi1, read_and_umi2 in zip(get_read_umis(in_file1, umi_length1, before1, after1, mask_pos1, False), get_read_umis(in_file2, umi_length2, before2, after2, mask_pos2, False)):
		read1, umi1 = read_and_umi1
		read2, umi2 = read_and_umi2
		if read1.name != read2.name: raise RuntimeError('mismatched reads\n%s (%s)\n%s (%s)\n' % (read1.name, in_file1.name, read2.name, in_file2.name))

		if relabel:
			combined_umi = umi1 + pair_separator + umi2
			yield ((add_umi_to_read(read1, combined_umi, trim_length1), umi1), (add_umi_to_read(read2, combined_umi, trim_length2), umi2))
		else:
			yield ((read1, umi1), (read2, umi2))

def get_umi_labeled_read_pairs (*args, **kwargs):
	'''
	like get_read_pair_umis but only yield a pair of relabeled reads, not the UMIs too
	'''
	return ((pair1[0], pair2[0]) for pair1, pair2 in get_read_pair_umis(*args, relabel = True, **kwargs))

