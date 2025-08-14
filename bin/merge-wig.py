#!/usr/bin/env python3

from __future__ import division
import argparse
import os
import math
from wig_files_writer import *




#####################################################
# Merges several replicates (as wig files) into a	#
#  representative version							#
#####################################################

parser = argparse.ArgumentParser()

parser.add_argument('--wig', '-w', nargs='+', help='Wig files to be merged into a single representative file')
parser.add_argument('--min_number_signal', '-m', type=int, default=None, help='Minimum amount of replicates with signal (at least strength 1) necessary to show peak in merged output. If no value is provided at least 51 percent of the samples need to show a signal.')
parser.add_argument('--output', '-o', help='Name of the output file')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

#######################
#######################
###    variables    ###
#######################
#######################

default_percent_of_samples_with_peak = 0.51


########################
########################
###    arguments     ###
########################
########################

replicates = args.wig

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(min_rep):
	print("\nStarting merging wig files")
	print( str(  len(replicates) ) + " replicates used for merging\n" )

	parsedReplicates = []
	

	# Parses all wig files used as input. Per cycle one replicate consisting of a forward and a reverse file will be parsed.
	# Additionally, the total amount of hits per replicate is calculated

	for replicate in replicates:
		parsedReplicates.append(parse_wig2(replicate))

	# get the names of all chromosomes
	chromosomeNames = []
	for k in parsedReplicates:
		chromosomeNames += k.keys()
	chromosomeNames = list(set(chromosomeNames))
	mergedChromosomes = {}
	numberReplicates = len(parsedReplicates)
	#TODO: Check min amount of samples with 
	# Check minimum amount of replicates with signal are necessary in order to add peak to merged output
	min_sample_with_peak = 0
	## Default mode: at least 51% of samples need signal
	if min_rep is None:
		min_sample_with_peak = math.ceil(numberReplicates * default_percent_of_samples_with_peak)
	## If chosen amount of replicates with peak is higher tahn the actual number of replicates then all replicates need to have a peak present
	elif min_rep >= numberReplicates:
		min_sample_with_peak = numberReplicates
	else:
		min_sample_with_peak = min_rep


	#go through every chromosome and get a unique list off all positions of all involved data sets. 
	for chromosome in chromosomeNames:
		mergedChromosomes[chromosome] = {}
		allPositions = []
		for wig in parsedReplicates:
			if chromosome in wig.keys():
				allPositions += wig[chromosome].keys()
		allPositions = list(set(allPositions))
		# create the mean hitCount for every position
		for position in allPositions:
			collectCountsForward = []
			collectCountsReverse = []
			for wig in parsedReplicates:
				if chromosome in wig.keys():
					if position in wig[chromosome].keys():
						collectCountsForward.append(wig[chromosome][position]["forward"])
						collectCountsReverse.append(wig[chromosome][position]["reverse"])
			peaksForward = sum(1 for x in collectCountsForward if x != 0)
			peaksReverse = sum(1 for x in collectCountsReverse if x != 0)
			if peaksForward >= min_sample_with_peak or peaksReverse >= min_sample_with_peak:
				mergedChromosomes[chromosome][position] = {}
				if peaksForward >= min_sample_with_peak:
					mergedChromosomes[chromosome][position]["forward"] = sum(collectCountsForward) / numberReplicates
				if peaksReverse >= min_sample_with_peak:
					mergedChromosomes[chromosome][position]["reverse"] = sum(collectCountsReverse) / numberReplicates

	write_wig2(mergedChromosomes, args.output)
	
	print("\nmerging complete\n")

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
main(min_rep = args.min_number_signal)
