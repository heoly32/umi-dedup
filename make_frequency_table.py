#!/usr/bin/env python

import collections, itertools, argparse, sys, Bio.SeqIO, pysam
from lib import umi_data

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a BAM or FASTQ file and generate a table containing the number of times each UMI was observed.')
parser_data = parser.add_argument_group('data files')
parser_format = parser.add_argument_group('format')
parser_perf = parser.add_argument_group('performance testing')
parser_format.add_argument('-f', '--fastq', action = 'store_true', help = 'input file is FASTQ format rather than BAM')
parser_perf.add_argument('--truncate_umi', action = 'store', type = int, default = None, help = 'truncate UMI sequences to this length')
parser_data.add_argument('in_file', action = 'store', nargs = '?', default = sys.stdin, help = 'input BAM or FASTQ')
parser_data.add_argument('out_file', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout, help = 'output table')
args = parser.parse_args()

umi_totals = umi_data.read_umi_counts_from_reads((Bio.SeqIO.parse(args.in_file, 'fastq') if args.fastq else pysam.Samfile('-' if args.in_file is sys.stdin else args.in_file, 'rb')), args.truncate_umi)
try: args.in_file.close()
except AttributeError: pass # if it's a string it won't close, but that's okay because it's already garbage-collected

# generate summary statistics
sys.stderr.write('%i UMIs read\n' % sum(umi_totals.values()))
for umi, count in umi_totals.items():
	args.out_file.write('%s\t%i\n' % (umi, count))
args.out_file.close()

