#!/usr/bin/env python

import collections, itertools, argparse, sys, Bio.SeqIO, pysam
from lib import umi_data

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a BAM or FASTQ file and generate a table containing the number of times each UMI was observed.')
parser.add_argument('-f', '--fastq', action = 'store_true', help = 'input file is FASTQ format rather than BAM')
parser.add_argument('in_file', action = 'store', nargs = '?', default = sys.stdin)
parser.add_argument('out_file', action = 'store', nargs = '?', type = argparse.FileType('w'), default = sys.stdout)
args = parser.parse_args()

umi_totals = umi_data.read_umi_counts_from_reads(Bio.SeqIO.parse(args.in_file, 'fastq') if args.fastq else pysam.Samfile('-' if args.in_file is sys.stdin else args.in_file, 'rb'))

# generate summary statistics
sys.stderr.write('%i UMIs read\n' % sum(umi_totals.values()))
for umi, count in umi_totals.items():
	args.out_file.write('%s\t%i\n' % (umi, count))

