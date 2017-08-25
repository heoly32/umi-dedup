import collections, itertools
from . import umi_data

Read = collections.namedtuple('Read', ['name', 'seq', 'qualities'])

def readfq(fp): # this is a generator function
	'''
	from Readfq by Heng Li, https://github.com/lh3/readfq
	yields reads as tuples of name, sequence, qualities (all strings)
	'''
	last = None # this is a buffer keeping the last unprocessed line
	while True: # mimic closure; is it a bad idea?
		if not last: # the first record or a record following a fastq
			for l in fp: # search for the start of the next record
				if l[0] in '>@': # fasta/q header line
					last = l[:-1] # save this line
					break
		if not last: break
		name, seqs, last = last[1:].partition(" ")[0], [], None
		for l in fp: # read the sequence
			if l[0] in '@+>':
				last = l[:-1]
				break
			seqs.append(l[:-1])
		if not last or last[0] != '+': # this is a fasta record
			yield Read(name, ''.join(seqs), None) # yield a fasta record
			if not last: break
		else: # this is a fastq record
			seq, leng, seqs = ''.join(seqs), 0, []
			for l in fp: # read the quality
				seqs.append(l[:-1])
				leng += len(l) - 1
				if leng >= len(seq): # have read enough quality
					last = None
					yield Read(name, seq, ''.join(seqs)); # yield a fastq record
					break
			if last: # reach EOF before reading enough quality
				yield Read(name, seq, None) # yield a fasta record instead
				break

def writefq (read):
	return '@%s\n%s\n+\n%s\n' % (read.name, read.seq, read.qualities)

def get_umi (seq, length, before = 0, mask_pos = []):
	'''
	extract UMI from a read sequence
	optionally specify offset (number of bases left of UMI to discard)
	optionally specify mask positions (0-based indices of bases within the UMI to discard)
	'''
	if length == 0: return ''
	
	umi = str(seq[before:(before + length)])
	if mask_pos is not None:
		for pos, next_pos in itertools.zip_longest(mask_pos[::-1], mask_pos[-2::-1], fillvalue = None):
			if next_pos >= pos: raise RuntimeError('mask_pos must be in ascending order')
			umi = umi[:pos] + umi[pos + 1:]
		
	return umi

def add_umi_to_read (read, umi, trim_length = 0):
	'''
	add UMI to Illumina-format FASTQ read name
	optionally also remove a specified number of bases from the beginning of the read (both sequence and qualities)
	'''
	labels = read.name.split(' ')
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
	name = ' '.join(labels)
	
	seq = read.seq[trim_length:]
	qualities = read.qualities[trim_length:]
	
	return Read(name, seq, qualities)

def get_read_umis (
	in_file,
	umi_length,
	before =   0,
	after =    0,
	mask_pos = [],
	relabel =  True
):
	'''
	returns a generator of pairs of Read instances and UMIs extracted from them
	moves UMI from read sequence to read name unless disabled
	defensively sorts mask_pos just in case you're naughty
	'''
	trim_length = before + umi_length + after
	mask_pos = sorted(mask_pos)
	if len(mask_pos) > 0 and mask_pos[-1] > umi_length: raise RuntimeError('UMI is only %i bases; can\'t mask position %i' % (umi_length, mask_pos[-1]))
	for read in readfq(in_file):
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
	before1 =        0,
	before2 =        0,
	after1 =         0,
	after2 =         0,
	mask_pos1 =      [],
	mask_pos2 =      [],
	pair_separator = umi_data.DEFAULT_SEPARATOR,
	relabel =        True
):
	'''
	returns a generator of pairs of pairs of Read objects and their UMIs
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

