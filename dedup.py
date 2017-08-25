#!/usr/bin/env python3

import copy, collections, argparse, pysam, sys
from lib import parse_sam, umi_data, optical_duplicates, naive_estimate, markdup_sam, pysam_progress

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a coordinate-sorted BAM file with labeled UMIs and mark or remove duplicates due to PCR or optical cloning, but not duplicates present in the original library. When PCR/optical duplicates are detected, the reads with the highest total base qualities are marked as non-duplicate - note we do not discriminate on MAPQ, or other alignment features, because this would bias against polymorphisms.')
parser_data = parser.add_argument_group('data files')
parser_format = parser.add_argument_group('format')
parser_alg = parser.add_argument_group('algorithm')
parser_perf = parser.add_argument_group('performance testing')
parser_reporting = parser.add_argument_group('reporting')
parser_format.add_argument('-r', '--remove', action = 'store_true', help = 'remove PCR/optical duplicates instead of marking them')
parser_alg.add_argument('-c', '--sequence_correction', action = 'store', choices = ['directional', 'kmeans'], help = 'correct UMI sequences before deduplication')
parser_alg.add_argument('-d', '--dist', action = 'store', type = int, default = optical_duplicates.DEFAULT_DIST, help = 'maximum pixel distance for optical duplicates (Euclidean); set to 0 to skip optical duplicate detection')
parser_alg.add_argument('-a', '--algorithm', action = 'store', default = 'naive', choices = ['naive'], help = 'algorithm for duplicate identification')
parser_perf.add_argument('--truncate_umi', action = 'store', type = int, default = None, help = 'truncate UMI sequences to this length')
parser_data.add_argument('in_file', action = 'store', nargs = '?', default = '-', help = 'input BAM')
parser_data.add_argument('out_file', action = 'store', nargs = '?', default = '-', help = 'output BAM')
parser_data.add_argument('-u', '--umi_table', action = 'store', type = argparse.FileType('r'), help = 'table of UMI sequences and (optional) prior frequencies')
parser_reporting.add_argument('-s', '--stats', action = 'store_true', help = 'compute additional library stats')
parser_reporting.add_argument('-q', '--quiet', action = 'store_true', help = 'don\'t show progress updates')
args = parser.parse_args()
# if args.algorithm == 'bayes' and args.umi_table is None and args.in_file == '-':
# 	raise RuntimeError('for the Bayesian algorithm, you must provide a UMI table filename, a BAM filename, or both')


in_bam = pysam.Samfile(args.in_file, 'rb')
if not args.quiet: progress = pysam_progress.ProgressTrackerByPosition(in_bam)
if in_bam.header['HD'].get('SO') != 'coordinate': raise RuntimeError('input file must be sorted by coordinate')
# create output file with modified header
bam_header = in_bam.header
bam_header['PG'].append({
	'ID': 'umi-bayes',
	'PN': 'umi-bayes',
	'VN': 'devel',
	'CL': ' '.join(sys.argv)
})
out_bam = pysam.Samfile(args.out_file, 'wb', header = bam_header)

# first pass through the input: get total UMI counts (or use table instead, if provided)
try:
	umi_totals = umi_data.read_umi_counts_from_table(args.umi_table, args.truncate_umi)
	args.umi_table.close()
except TypeError:
	umi_totals = None

# second pass: mark duplicates
dup_marker = markdup_sam.DuplicateMarker(
	alignments =          in_bam,
	umi_frequency =       umi_totals,
	algorithm =           args.algorithm,
	optical_dist =        args.dist,
	truncate_umi =        args.truncate_umi,
	sequence_correction = args.sequence_correction
)
if not args.quiet: progress.reset()
for alignment in dup_marker:
	if not (args.remove and alignment.is_duplicate): out_bam.write(alignment)
	if not args.quiet: progress.update(alignment)
in_bam.close()
out_bam.close()
if not args.quiet: del progress

# report summary statistics
# alignment counts
# if args.algorithm == 'bayes' and args.umi_table is None: # would already have reported alignments read
# 	assert sum(umi_totals.nonzero_values()) == dup_marker.category_counts['usable alignment']
sys.stderr.write(
	'%i\talignments read\n%i\tusable alignments read\n\n' %                                                       (dup_marker.category_counts['alignment'], dup_marker.category_counts['usable alignment']) +
	'%i\tUMI sequence corrections\n' %                                                                            dup_marker.category_counts['sequence correction'] +
	'\n' +
	'%i\tdistinct alignments\n' %                                                                                 dup_marker.category_counts['distinct'] +
	('%i\toptical duplicates\n' %                                                                                 dup_marker.category_counts['optical duplicate'] if args.dist != 0 else '') +
	'%i\tPCR duplicates\n%i\tpre-PCR duplicates rescued by UMIs\n%i\tpre-PCR duplicates rescued by algorithm\n' % tuple(dup_marker.category_counts[x] for x in ['PCR duplicate', 'UMI rescued', 'algorithm rescued'])
)
# library stats
if args.stats:
	sys.stderr.write(
		'\n' +
		'%.3f\tmean position entropy before deduplication\n' % dup_marker.get_mean_pos_entropy('before') +
		'%.3f\tmean position entropy after deduplication\n' %  dup_marker.get_mean_pos_entropy('after') +
		'%.3f\tlibrary entropy before deduplication\n' %       dup_marker.get_library_entropy('before') +
		'%.3f\tlibrary entropy after deduplication\n' %        dup_marker.get_library_entropy('after') +
		'%i\testimated library size\n' %                       dup_marker.estimate_library_size()
	)

