#!/usr/bin/env python3

import argparse
import os
import numpy
from wig_files_writer import *


import sys
#####################################################
# Calculates the distances between peaks from 
# 	iCLIP-experiments
#####################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', nargs='+', help='WIG files used for distance calculations')
parser.add_argument('--output', '-o', help='Name of the output file')
parser.add_argument('--percentile', '-p', type=float, help='Percentile to use for cutoff')
parser.add_argument('--distance', '-d', type=int, help='Maximum distance to be considered')
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

percentile = args.percentile
print("Percentile used for cutoff determination: " + str(percentile))
maxDistance = args.distance
outFile = args.output


##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main():
	wigFiles = []
	for file in args.input:
		wigFiles.append( parse_wig2( file ) )



	allCounts = bundle_counts( wigFiles ) 

	numpy.set_printoptions(threshold = sys.maxsize)

	cutoff = numpy.percentile(allCounts, percentile)
	print("The value cutoff is set to " + str(cutoff) )

	# remove all elements with a value below the cutoff
	for file in list(wigFiles):
		for chrom in dict(file):
			for position in dict(file[chrom]):
				for strand, value in dict(file[chrom][position]).items():
			# a copy of the dictionary has to be created as the size changes during the loop
					if abs(value) < cutoff:
						del file[chrom][position][strand]
				if not file[chrom][position]:
					del file[chrom][position]
			if not file[chrom]:
				del file[chrom]
	wigFiles = list(filter(None,wigFiles))

	distances = [0] * maxDistance

	splitFiles = splitForwardReverse(wigFiles)

	for file in splitFiles:
		for chrom in file:
			allPositions = sorted( [ int(i) for i in file[chrom].keys() ] )
			for position in allPositions:
				# Ends the loop if it reached the last
				if position == allPositions[-1]:
					break
				indexToTest = allPositions.index(position) + 1

				while indexToTest < len(allPositions) and allPositions[indexToTest] - position <= maxDistance :

					currentDistance = allPositions[indexToTest] - position - 1
					distances[currentDistance] = distances[currentDistance] + 1

					indexToTest = indexToTest + 1

	write_output(distances, outFile)



	
#######################
#######################
###    functions    ###
#######################
#######################


def bundle_counts(wigFiles):
	allCounts = []
	floatCounts = numpy.array(allCounts, dtype = numpy.float32)
	for file in wigFiles:
		for chromosome in file:
			for pos in file[chromosome]:
				for strand in file[chromosome][pos]:
					if not file[chromosome][pos][strand] == 0:
						floatCounts = numpy.append(floatCounts, abs(file[chromosome][pos][strand]))
	return floatCounts

def write_output(distances, fileName):

	outFile = open(os.path.realpath(fileName), 'w')

	outFile.write( "%s\t%s\n" % (
		"distance",
		"value"
	) )

	count = 1
	for entry in distances:
		outFile.write( "%s\t%s\n" % (
			str(count),
			str(entry)
		) )
		count += 1

#######################
#######################
###      tests      ###
#######################
#######################




##########################
### starts main script ###
##########################
main()
