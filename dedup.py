#!/usr/bin/env python

import collections, argparse, pysam, sys
from lib import parse_sam, umi_data, naive_estimate

# parse arguments
parser = argparse.ArgumentParser(description = 'Read a coordinate-sorted SAM file with labeled UMIs and mark or remove duplicates due to PCR or optical cloning, but not duplicates present in the original library. When PCR/optical duplicates are detected, the reads with the highest total base qualities are marked as non-duplicate - note we do not discriminate on MAPQ, or other alignment features, because this would bias against polymorphisms.')
parser.add_argument('-r', '--remove', action = 'store_true', help = 'remove PCR/optical duplicates instead of marking them')
#parser.add_argument('-d', '--dist', action = 'store', help = 'maximum pixel distance for optical duplicates (Euclidean); set to 0 to skip optical duplicate detection', type = int, default = 100)
parser.add_argument('-u', '--umi_table', action = 'store', help = 'table of UMI sequences and prior frequencies')
parser.add_argument('infile', action = 'store', nargs = '?', default = '-') # actually this filename should be required since we can't use a stream
parser.add_argument('outfile', action = 'store', nargs = '?', default = '-')
args = parser.parse_args()
# maybe should have some options about the output format
if args.umi_table is None and args.infile == '-':
	raise RuntimeError('cannot read from stream unless UMI table is provided')

in_sam = pysam.Samfile(args.infile, 'rb')
if in_sam.header['HD'].get('SO') != 'coordinate': raise RuntimeError('input file must be sorted by coordinate')
out_sam = pysam.Samfile(args.outfile, 'wb', template = in_sam) # should add a line to the header indicating it was processed
read_counter = collections.Counter()

'''
central concept: cheat a little, detecting duplicate reads by the fact that they align to the same (start) position rather than by their sequences
implementation: as you traverse the coordinate-sorted input reads, add each read to both a FIFO buffer (so they can be output in the same order) and a dictionary that groups reads by start position and strand (which is what you need for deduplication); this is not inefficient because both data structures contain pointers to the same pysam.AlignedSegment objects
then, after each input read, check the left end of the buffer to tell whether the oldest read is in a position that will never accumulate any more hits (because the input is sorted by coordinate); if so, estimate the duplication at that position, mark all the reads there accordingly, then output the read with the appropriate marking
'''


read_buffer = collections.deque()
pos_tracker = ({}, {}) # data structure containing observed UMIs and corresponding tracking information; top level is by strand (0 = forward, 1 = reverse), then next level is by 5' read start position (dict since these will be sparse and are only looked up by identity), then at each position the next level is by UMI (dict), and that contains a variety of data; there is no level for reference ID because there is no reason to store more than one chromosome at a time
def pop_buffer(): # pop the oldest read off the buffer (into the output), but first make sure its position has been deduplicated
	read = read_buffer.popleft()
	start_pos, umi = parse_sam.get_start_pos(read), parse_sam.get_umi(read)
	this_pos = pos_tracker[read.is_reverse][start_pos]
	
	# deduplicate reads at this position
	if not this_pos['deduplicated']:
		
		umi_counts = umi_data.make_umi_counts(umi_totals.keys(), (len(this_pos['reads'][umi]) if umi in this_pos['reads'] else 0 for umi in umi_totals.keys()))
		
		# P ESTIMATION / DEDUPLICATION GOES HERE
		dedup_counts = naive_estimate.deduplicate_counts(umi_counts)
		
		umi_data.mark_duplicates(this_pos['reads'], dedup_counts)
		this_pos['deduplicated'] = True
	
	# output read
	read_counter['duplicate' if read.is_duplicate else 'nonduplicate'] += 1
	if not (args.remove and read.is_duplicate): out_sam.write(read)
	
	# prune the tracker
	if read is this_pos['last read']: del pos_tracker[read.is_reverse][start_pos]


# first pass through the input: get total UMI counts (or use table instead, if provided)
umi_totals = umi_data.read_umi_counts_from_sam(in_sam) if args.umi_table is None else umi_data.read_umi_counts_from_table(open(args.umi_table))
if args.umi_table is None: sys.stderr.write('%i\tusable alignments read\n' % sum(umi_totals.values()))


# second pass through the input
for read in in_sam:
	if not parse_sam.is_good(read):	continue
	read_counter['read'] += 1
	
	start_pos, umi = parse_sam.get_start_pos(read), parse_sam.get_umi(read)
	
	# advance the buffer
	while read_buffer and (read_buffer[0].reference_id < read.reference_id or parse_sam.get_start_pos(read_buffer[0]) < read.reference_start): # if the top of the buffer is at a position that's definitely not going to get any more hits
		pop_buffer()
	
	# add read to buffer and tracking data structure	
	read_buffer.extend([read])
	if start_pos not in pos_tracker[read.is_reverse]: # first time we've seen this position+strand
		pos_tracker[read.is_reverse][start_pos] = {'reads': {umi: [read]}, 'deduplicated': False}
	elif umi not in pos_tracker[read.is_reverse][start_pos]['reads']: # first time we've seen this UMI there
		pos_tracker[read.is_reverse][start_pos]['reads'][umi] = [read]
	else:
		pos_tracker[read.is_reverse][start_pos]['reads'][umi] += [read]
	pos_tracker[read.is_reverse][start_pos]['last read'] = read

# flush the buffer
while read_buffer: pop_buffer()


# generate summary statistics
if args.umi_table is None: assert sum(umi_totals.values()) == read_counter['read']
if args.umi_table is not None: sys.stderr.write('%i\tusable alignments read\n' % read_counter['read'])
sys.stderr.write('%i\tunduplicated\n%i\tduplicates\n' % (read_counter['nonduplicate'], read_counter['duplicate']))

