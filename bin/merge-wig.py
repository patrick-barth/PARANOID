#!/usr/bin/env python3

from __future__ import division
import argparse
import os
from wig_files_writer import *




#####################################################
# Merges several replicates (as wig files) into a	#
#  representative version							#
#####################################################

parser = argparse.ArgumentParser()

parser.add_argument('--wig', '-w', nargs='+', help='Wig files to be merged into a single representative file')
parser.add_argument('--output', '-o', help='Name of the output file')
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

replicates = args.wig

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main():
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
			addedCountsForward = 0
			addedCountsReverse = 0
			for wig in parsedReplicates:
				if chromosome in wig.keys():
					if position in wig[chromosome].keys():
						addedCountsForward += wig[chromosome][position]["forward"]
						addedCountsReverse += wig[chromosome][position]["reverse"]
			mergedChromosomes[chromosome][position] = {}
			mergedChromosomes[chromosome][position]["forward"] = addedCountsForward / numberReplicates
			mergedChromosomes[chromosome][position]["reverse"] = addedCountsReverse / numberReplicates

	write_wig2(mergedChromosomes, args.output)
	
	print("\nmerging complete\n")

#######################
#######################
###    functions    ###
#######################
#######################


#####################################
# Merges all given wig files.		#
# Returns a representative wig file #
#####################################

def merge_wig(wigFiles):
	
	return mergedChromosomes





#######################
#######################
###      tests      ###
#######################
#######################




##########################
### starts main script ###
##########################
main()
