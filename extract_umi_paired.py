#!/usr/bin/env python

import argparse, sys, collections
from lib import parse_fastq

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a pair of FASTQ files containing UMIs prepended to alignable sequences, and output the same files with the UMIs moved into the read names. Give both reads in a pair the same UMI label, combined from both UMIs if there are two. Report the base frequencies.')
parser_data = parser.add_argument_group('data files')
parser_umi = parser.add_argument_group('UMI format')
parser_umi.add_argument('umi_length1', action = 'store', type = int, help = 'length of UMI sequence in read 1 (these bases become the UMI label)')
parser_umi.add_argument('umi_length2', action = 'store', type = int)
parser_umi.add_argument('-b1', '--before1', action = 'store', type = int, default = 0, help = 'number of bases before UMI to discard from read 1')
parser_umi.add_argument('-b2', '--before2', action = 'store', type = int, default = 0)
parser_umi.add_argument('-a1', '--after1', action = 'store', type = int, default = 0, help = 'number of bases after UMI to discard from read 1')
parser_umi.add_argument('-a2', '--after2', action = 'store', type = int, default = 0)
parser_umi.add_argument('-m1', '--mask1', action = 'append', type = int, default = [], help = 'position within UMI to discard from read 1 (can be used multiple times)')
parser_umi.add_argument('-m2', '--mask2', action = 'append', type = int, default = [])
parser_data.add_argument('in_file1', action = 'store', type = argparse.FileType('r'), help = 'input FASTQ, read 1')
parser_data.add_argument('in_file2', action = 'store', type = argparse.FileType('r'))
parser_data.add_argument('out_file1', action = 'store', type = argparse.FileType('w'), help = 'output FASTQ, read 1')
parser_data.add_argument('out_file2', action = 'store', type = argparse.FileType('w'))
args = parser.parse_args()

read_pair_counter = 0
base_counter1 = [collections.Counter() for i in range(args.umi_length1)]
base_counter2 = [collections.Counter() for i in range(args.umi_length2)]

for pair1, pair2 in parse_fastq.get_read_pair_umis(
	args.in_file1,
	args.in_file2,
	args.umi_length1,
	args.umi_length2,
	args.before1,
	args.before2,
	args.after1,
	args.after2,
	args.mask1,
	args.mask2
):
	args.out_file1.write(pair1[0].format('fastq'))
	args.out_file2.write(pair2[0].format('fastq'))
	read_pair_counter += 1
	for i in range(args.umi_length1): base_counter1[i][pair1[1][i]] += 1
	for i in range(args.umi_length2): base_counter2[i][pair2[1][i]] += 1
map(lambda x: x.close(), (args.in_file1, args.in_file2, args.out_file1, args.out_file2))

# generate summary statistics
sys.stderr.write('%i read pairs processed\n\n' % read_pair_counter)

alphabet = sorted(set.union(*(map(set, base_counter1 + base_counter2)))) # detect from data since it may include N
sys.stderr.write('UMI base frequency by position\n')
sys.stderr.write('\t'.join(alphabet) + '\n')
for pos in base_counter1: sys.stderr.write('\t'.join(str(pos[base]) for base in alphabet) + '\n')
if args.umi_length1 > 0 and args.umi_length2 > 0: sys.stderr.write('\t'.join(['+'] * len(alphabet)) + '\n')
for pos in base_counter2: sys.stderr.write('\t'.join(str(pos[base]) for base in alphabet) + '\n')
sys.stderr.write('\n')

