#!/usr/bin/env python

import argparse, sys, collections, Bio.SeqIO

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a FASTQ file containing UMIs prepended to alignable sequences, and output the same file with the UMIs moved into the read names. Report the base frequencies.')
parser_data = parser.add_argument_group('data files')
parser_umi = parser.add_argument_group('UMI format')
parser_umi.add_argument('umi_length', action = 'store', type = int, help = 'length of UMI sequence (these bases become the UMI label)')
parser_umi.add_argument('-b', '--before', action = 'store', type = int, default = 0, help = 'number of bases before UMI to discard')
parser_umi.add_argument('-a', '--after', action = 'store', type = int, default = 0, help = 'number of bases after UMI to discard')
parser_umi.add_argument('-m', '--mask', action = 'append', type = int, help = 'position within UMI to discard (can be used multiple times)')
parser_data.add_argument('in_file', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin, help = 'input FASTQ')
parser_data.add_argument('out_file', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout, help = 'output FASTQ')
args = parser.parse_args()

read_counter = 0
base_counter = [collections.Counter() for i in range(args.umi_length)]

trim_length = args.before + args.umi_length + args.after
try:
	mask_pos = sorted(args.mask, reverse = True) # sorted in descending order so they can be done from the right end, i.e. each mask doesn't change the relative position of the next
	if mask_pos[0] > args.umi_length: raise RuntimeError('UMI is only %i bases; can\'t mask position %i' % (args.umi_length, mask_pos[0]))
except TypeError:
	mask_pos = []

for read in Bio.SeqIO.parse(args.in_file, 'fastq'):

	if len(read.seq) <= trim_length: raise RuntimeError('%i bases to be trimmed but read %s is only %i bases' % (trim_length, read.id, len(read.seq)))
	
	# parse UMI
	umi = str(read.seq[args.before:(args.before + args.umi_length)])
	for i in range(args.umi_length): base_counter[i][umi[i]] += 1
	
	# mask UMI
	for pos in mask_pos: umi = umi[:(pos - 1)] + umi[pos:]
	
	# modify read
	read.id += ':' + umi
	if ' ' in read.description: # problems with spaces
		read.description = read.description[:read.description.index(' ')] + ':' + umi + read.description[read.description.index(' '):]
	else:
		read.description += ':' + umi
	qualities = read.letter_annotations['phred_quality']
	read.letter_annotations = {} # letter annotations must be emptied before changing sequence
	
	read.seq = read.seq[trim_length:]
	read.letter_annotations['phred_quality'] = qualities[trim_length:]
	
	# output read
	args.out_file.write(read.format('fastq'))
	read_counter += 1


# generate summary statistics
sys.stderr.write('%i reads processed\n\n' % read_counter)

alphabet = sorted(set.union(*(set(i) for i in base_counter))) # detect from data since it may include N
sys.stderr.write('UMI base frequency by position\n')
sys.stderr.write('\t'.join(alphabet) + '\n')
for pos in base_counter:	sys.stderr.write('\t'.join(str(pos[base]) for base in alphabet) + '\n')
sys.stderr.write('\n')

