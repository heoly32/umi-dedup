#!/usr/bin/env python

# get ERCC hit counts, separated by UMI, from a BAM file and print a nice matrix
# usage: ercc_umi_hits.py < infile.bam > outfile.tsv
import collections, itertools, sys, pysam

bases = 'ACGT'
umi_length = 5
umi_list = [''.join(umi) for umi in itertools.product(bases, repeat = umi_length)]

counts = {}

sam = pysam.Samfile('/dev/stdin', 'rb')
which_ercc = [tid for tid in range(0, sam.nreferences) if sam.getrname(tid).startswith('ERCC')]

for ercc_tid in which_ercc:
	ercc_name = sam.getrname(ercc_tid)
	if not ercc_name in counts:
		counts[ercc_name] = {}

n_reads = n_ercc_hits = n_ercc_forward_hits = 0
for read in sam:
	if read.is_secondary or read.is_unmapped or read.is_read2:
		continue
	n_reads += 1
	if read.reference_id in which_ercc:
		n_ercc_hits += 1
		if not read.is_reverse:
			n_ercc_forward_hits += 1
			transcript = sam.getrname(read.reference_id)
			pos = read.reference_start
			umi = read.query_name[read.query_name.rindex(':') + 1:]
			if not pos in counts[transcript]: counts[transcript][pos] = collections.Counter()
			counts[transcript][pos][umi] += 1

sys.stderr.write('%i aligned reads, %i ERCC hits, %i forward-oriented\n' % (n_reads, n_ercc_hits, n_ercc_forward_hits))

print('\t'.join(['position'] + umi_list))
for transcript in sorted(counts.keys()):
	for pos in sorted(counts[transcript].keys()):
		print('\t'.join(['%s:%i' % (transcript, pos + 1)] + [str(counts[transcript][pos][umi]) for umi in umi_list]))

