#!/usr/bin/env python3

import pysam
import argparse
import os
from wig_files_writer import parse_wig2


#############################################
# Converts counts from a wig file into 		#
#	a bam file (according to iCLIP results)	#
#	Each positions will results in the 		#
#	number of bam records according to the 	#
# 	count number 							#
#############################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Wig2 file used as input')
parser.add_argument('--output', '-o', help='prefix of the output files')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

#######################
#######################
###    variables    ###
#######################
#######################


########################
########################
###    parameters    ###
########################
########################


##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main():
	chromosome_lengths = []
	parsed_wig = parse_wig2(args.input)
	for chrom in parsed_wig:
		chromosome_lengths.append({'LN': max([int(k) for k in parsed_wig[chrom].keys()]),
									'SN': chrom})

	header = { 'HD': {'VN': '1.0'},
				'SQ': chromosome_lengths }

	index_chrom_by_name = build_dict(chromosome_lengths, key = 'SN')
	print(header)


	with pysam.AlignmentFile(args.output, "wb", header=header) as outf:
		for chrom in parsed_wig:
			for pos in parsed_wig[chrom].keys():
				for strand in ['forward', 'reverse']:
					for count in list(range(1, int( parsed_wig[chrom][pos][strand] ) + 1)):
						a = pysam.AlignedSegment()
						a.query_name = chrom + "_" + pos + "_" + strand + "_" + str(count)
						a.query_sequence = "N"
						a.flag = 0 if strand == 'forward' else 16
						a.reference_id = index_chrom_by_name.get(chrom)["index"]
						a.reference_start = int(pos) -1
						a.mapping_quality = 20
						a.cigar = [(0,1)]
						#a.next_reference_id = 0
						#a.next_reference_start=199
						a.template_length=0
						a.query_qualities = pysam.qualitystring_to_array("<")
						a.tags = [("NM", 0)]
						outf.write(a)

	

#######################
#######################
###    functions    ###
#######################
#######################

def build_dict(seq, key):
    return dict((d[key], dict(d, index=index)) for (index, d) in enumerate(seq))

##########################
### starts main script ###
##########################
main()

