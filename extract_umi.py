#!/usr/bin/env python

import argparse, sys, collections
from lib import parse_fastq

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a FASTQ file containing UMIs prepended to alignable sequences, and output the same file with the UMIs moved into the read names. Report the base frequencies.')
parser_data = parser.add_argument_group('data files')
parser_umi = parser.add_argument_group('UMI format')
parser_umi.add_argument('umi_length', action = 'store', type = int, help = 'length of UMI sequence (these bases become the UMI label)')
parser_umi.add_argument('-b', '--before', action = 'store', type = int, default = 0, help = 'number of bases before UMI to discard')
parser_umi.add_argument('-a', '--after', action = 'store', type = int, default = 0, help = 'number of bases after UMI to discard')
parser_umi.add_argument('-m', '--mask', action = 'append', type = int, default = [], help = 'position within UMI to discard (can be used multiple times)')
parser_data.add_argument('in_file', action = 'store', nargs = '?', type = argparse.FileType('r'), default = sys.stdin, help = 'input FASTQ')
parser_data.add_argument('out_file', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout, help = 'output FASTQ')
args = parser.parse_args()

read_counter = 0
base_counter = [collections.Counter() for i in range(args.umi_length)]

for read, umi in parse_fastq.get_read_umis(args.in_file, args.umi_length, args.before, args.after, args.mask):
	args.out_file.write(read.format('fastq'))
	read_counter += 1
	for i in range(args.umi_length): base_counter[i][umi[i]] += 1
args.in_file.close()
args.out_file.close()

# generate summary statistics
sys.stderr.write('%i reads processed\n\n' % read_counter)

alphabet = sorted(set.union(*(map(set, base_counter)))) # detect from data since it may include N
sys.stderr.write('UMI base frequency by position\n')
sys.stderr.write('\t'.join(alphabet) + '\n')
for pos in base_counter: sys.stderr.write('\t'.join(str(pos[base]) for base in alphabet) + '\n')
sys.stderr.write('\n')

