#!/usr/bin/env python3

from __future__ import division
import argparse
import os
from wig_files_writer import *




#####################################################
# splits a wig2 file into 2 separate wig files		#
#####################################################

parser = argparse.ArgumentParser()

parser.add_argument('--input', '-i', help='Name of input wig2 file')
parser.add_argument('--output', '-o', help='Prefix of output files')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

#######################
#######################
###    variables    ###
#######################
#######################


########################
########################
###    arguments     ###
########################
########################

outputPrefix = args.output


##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main():
	wig = parse_wig2(args.input)

	forwardPeaks = {}
	reversePeaks = {}

	for chromosome in wig:
		for position in wig[chromosome]:
			if not float(wig[chromosome][position]["forward"]) == 0:
				if not chromosome in forwardPeaks.keys():
					forwardPeaks[chromosome] = {}
				forwardPeaks[chromosome][position] = wig[chromosome][position]["forward"]
			if not float(wig[chromosome][position]["reverse"]) == 0:
				if not chromosome in reversePeaks.keys():
					reversePeaks[chromosome] = {}
				reversePeaks[chromosome][position] = wig[chromosome][position]["reverse"]  


	if not len(forwardPeaks) == 0:
		write_wig(forwardPeaks, outputPrefix + "_forward.wig")
	if not len(reversePeaks) == 0:
		write_wig(reversePeaks, outputPrefix + "_reverse.wig")

	


#######################
#######################
###    functions    ###
#######################
#######################




#######################
#######################
###      tests      ###
#######################
#######################




##########################
### starts main script ###
##########################
main()