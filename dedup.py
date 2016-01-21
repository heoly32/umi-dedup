#!/usr/bin/env python

import copy, collections, argparse, pysam, sys
from lib import parse_sam, umi_data, optical_duplicates, naive_estimate, bayes_estimate, markdup_sam

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a coordinate-sorted BAM file with labeled UMIs and mark or remove duplicates due to PCR or optical cloning, but not duplicates present in the original library. When PCR/optical duplicates are detected, the reads with the highest total base qualities are marked as non-duplicate - note we do not discriminate on MAPQ, or other alignment features, because this would bias against polymorphisms.')
parser_data = parser.add_argument_group('data files')
parser_format = parser.add_argument_group('format')
parser_alg = parser.add_argument_group('algorithm')
parser_perf = parser.add_argument_group('performance testing')
parser_format.add_argument('-r', '--remove', action = 'store_true', help = 'remove PCR/optical duplicates instead of marking them')
parser_alg.add_argument('-d', '--dist', action = 'store', type = int, default = optical_duplicates.DEFAULT_DIST, help = 'maximum pixel distance for optical duplicates (Euclidean); set to 0 to skip optical duplicate detection')
parser_alg.add_argument('-a', '--algorithm', action = 'store', default = 'naive', choices = ['naive', 'bayes', 'uniform-bayes'], help = 'algorithm for duplicate identification')
parser_alg.add_argument('--nsamp', action = 'store', type = int, default = bayes_estimate.DEFAULT_NSAMP)
parser_alg.add_argument('--nthin', action = 'store', type = int, default = bayes_estimate.DEFAULT_NTHIN)
parser_alg.add_argument('--nburn', action = 'store', type = int, default = bayes_estimate.DEFAULT_NBURN)
parser_perf.add_argument('--truncate_umi', action = 'store', type = int, default = None, help = 'truncate UMI sequences to this length')
parser_data.add_argument('in_file', action = 'store', nargs = '?', default = '-', help = 'input BAM')
parser_data.add_argument('out_file', action = 'store', nargs = '?', default = '-', help = 'output BAM')
parser_data.add_argument('-u', '--umi_table', action = 'store', type = argparse.FileType('r'), help = 'table of UMI sequences and (optional) prior frequencies')
args = parser.parse_args()
if args.algorithm == 'bayes' and args.umi_table is None and args.in_file == '-':
	raise RuntimeError('for the Bayesian algorithm, you must provide a UMI table filename, a BAM filename, or both')


in_bam = pysam.Samfile(args.in_file, 'rb')
if in_bam.header['HD'].get('SO') != 'coordinate': raise RuntimeError('input file must be sorted by coordinate')
out_bam = pysam.Samfile(args.out_file, 'wb', template = in_bam) # should add a line to the header indicating it was processed

# first pass through the input: get total UMI counts (or use table instead, if provided)
try:
	umi_totals = umi_data.read_umi_counts_from_table(args.umi_table, args.truncate_umi)
	args.umi_table.close()
except TypeError:
	if args.algorithm == 'bayes':
		umi_totals = umi_data.read_umi_counts_from_reads(in_bam, args.truncate_umi)
		sys.stderr.write('%i\tusable alignments read\n\n' % sum(umi_totals.values()))
		in_bam.reset()
	else:
		umi_totals = None

# second pass: mark duplicates
dup_marker = markdup_sam.DuplicateMarker(
	alignments = in_bam,
	umi_frequency = umi_totals,
	algorithm = args.algorithm,
	optical_dist = args.dist,
	truncate_umi = args.truncate_umi,
	nsamp = args.nsamp,
	nthin = args.nthin,
	nburn = args.nburn
)
for alignment in dup_marker:
	if not (args.remove and alignment.is_duplicate): out_bam.write(alignment)
in_bam.close()
out_bam.close()

# report summary statistics
if args.algorithm == 'bayes' and args.umi_table is None: # would already have reported alignments read
	assert sum(umi_totals.values()) == dup_marker.counts['usable alignment']
sys.stderr.write(
	('%i\talignments read\n%i\tusable alignments read\n\n' % (dup_marker.counts['alignment'], dup_marker.counts['usable alignment']) if not (args.algorithm == 'bayes' and args.umi_table is None) else '') +
	'%i\tdistinct alignments\n' % dup_marker.counts['distinct'] +
	('%i\toptical duplicates\n' % dup_marker.counts['optical duplicate'] if args.dist != 0 else '') +
	'%i\tPCR duplicates\n%i\tpre-PCR duplicates rescued by UMIs\n%i\tpre-PCR duplicates rescued by algorithm\n' % tuple(dup_marker.counts[x] for x in ['PCR duplicate', 'UMI rescued', 'algorithm rescued'])
)


