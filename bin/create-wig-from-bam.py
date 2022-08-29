#!/usr/bin/env python

from __future__ import division
import argparse
import os
import re
import pybam #needs python2
from wig_files_writer import *

#############################################
# Creates a wig-file out of a bam-file		#
# Forward and reverse strand will be split	#
# Python 2 is needed due to pybam			#
#############################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='Bam file used as input')
parser.add_argument('--mapq', '-q', help='Minimum mapping quality')
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

if args.mapq:
	minQuality = int(args.mapq)
else:
	minQuality = 2

#outName =  os.path.splitext( os.path.basename( args.input ) )[0]
output = args.output
cwd = os.getcwd()

##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main():
	peaks = {}
	# Both variables will be used to calculate how many reads are thrown out due to the quality filtering
	countLowQualReads = 0
	countTotalReads = 0
	# go through every entry in the bam file
	for alignment in pybam.read(args.input):
		chromosome = alignment.sam_rname
		flag = alignment.sam_flag
		position = alignment.sam_pos0 + 1
		mapq = alignment.sam_mapq
		cigar = alignment.sam_cigar_string
		countTotalReads += 1
		# only mappings with a quality score above or equal to the minimum quality will be taken.
		# The rest is discarded
		if int(mapq) >= minQuality:
			# formats the flag into a binary value -> easier to parse
			binaryFlag = "{0:016b}".format(flag)
			# the 12th position says if the read mapped to the forward (0) or the reverse (1) strand
			if int(binaryFlag[11]) == 0:
				position = position - 1 # the peak is one nucleotide in front of the read
				# There is no position 0 in wig files. Needs to be removed or else it will trigger an error
				if position <= 0:
					continue
				# This whole block is only for the forward strand
				# Check if the chromosome was already initialized. If not it will now be initialized
				if chromosome in peaks:
					# Check if the position inside of the chromosome was already initialized.
					#  If yes, the value of the position is increased by one. If not it will be initialized with a value of 1
					if position in peaks[chromosome]:
						if "forward" in peaks[chromosome][position].keys():
							peaks[chromosome][position]["forward"] += 1
						else:
							peaks[chromosome][position]["forward"] =  1
					else:
						peaks[chromosome][position] = {"forward": 1}
				else:
					peaks[chromosome] = {position: {"forward": 1}}
			elif int(binaryFlag[11]) == 1:
				# This whole block is only for the reverse strand.
				# It basically functions like the above block with the difference that the values are negative (due to the strand being reverse) and that the positions are taken -1
				position = position + 1 # The peak for reverse reads the nucleotide behind the alignment not the last nucleotide of the alignment
				#//TODO: read in chromosome lengths and check if the allocated position is behind the end of the chromosome
				reversePosition = calculateStart(position, cigar) # - 1 ; need in case a cl-site is detected off the chromosome 
				if chromosome in peaks:
					if reversePosition in peaks[chromosome]:
						if "reverse" in peaks[chromosome][reversePosition].keys():
							peaks[chromosome][reversePosition]["reverse"] -= 1
						else:
							peaks[chromosome][reversePosition]["reverse"] = -1
					else: 
						peaks[chromosome][reversePosition] = {"reverse": -1}
				else:
					peaks[chromosome] = {reversePosition: {"reverse": -1}}
			else:
				print("Erroneous flag: %s" % (binaryFlag))
		else:
			countLowQualReads += 1

	# print how many of the reads are thrown due to low quality
	print("%i (%.2f%%) of the %i reads were discarded due to their mapq-score (below %i)" % (
		countLowQualReads,
		100 * (countLowQualReads/countTotalReads),
		countTotalReads,
		minQuality
	))

	peaks = fill_entries(peaks)

	write_wig2(peaks, output)

#######################
#######################
###    functions    ###
#######################
#######################

#####################################
# Calculates the actual start 		#
#  of the read out of its end 		#
#  position and the cigar string 	#
#####################################

def calculateStart( endPos, cigar ):
	matches = re.findall(r"\d+\D", cigar)
	currentPos = endPos
	for cigarPart in matches:
		number = re.findall(r"\d+", cigarPart)[0]
		operator = re.findall(r"\D", cigarPart)[0]
		if operator == 'D' or operator == 'M' or operator == 'N' or operator == 'X' or operator == '=':
			currentPos += int(number)

	return currentPos

#################################
# Fills missing entries with 0 	#
#################################

def fill_entries(wig):

	for chromosome in wig:
		for position in wig[chromosome]:
			if not "forward" in wig[chromosome][position].keys():
				wig[chromosome][position]["forward"] = 0
			elif not "reverse" in wig[chromosome][position].keys():
				wig[chromosome][position]["reverse"] = 0

	return wig


##########################
### starts main script ###
##########################
main()
