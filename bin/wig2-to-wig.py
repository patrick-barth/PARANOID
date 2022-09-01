#!/usr/bin/env python3

from __future__ import division
import argparse
import os




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
		forwardPeaks[chromosome] = {}
		reversePeaks[chromosome] = {}
		for position in wig[chromosome]:
			if not float(wig[chromosome][position]["forward"]) == 0:
				forwardPeaks[chromosome][position] = wig[chromosome][position]["forward"]
			if not float(wig[chromosome][position]["reverse"]) == 0:
				reversePeaks[chromosome][position] = wig[chromosome][position]["reverse"]  


	write_wig(forwardPeaks, outputPrefix + "_forward.wig")
	write_wig(reversePeaks, outputPrefix + "_reverse.wig")


#######################
#######################
###    functions    ###
#######################
#######################

#####################################
# parses a wig file. 				#
# Returns a dictionary containing	#
#	a dictionary:					#
# 	key -> chromsome				#
#	value -> 						#
#		key -> position 			#
#		value -> count 				#
#####################################

def parse_wig2( wigFile ):
	chromosomes = {}
	currentChromosome = ""
	with open ( wigFile, 'r' ) as wigFile:
		for line in wigFile:
			line = line.rstrip()
			# a file might contain several chromosomes. New chromosomes start with a line like "variableStep chrom=DQ380154.1 span=1". Due to that they are divided by this line. 
			# Also the chromosome name is parsed from that line
			if line.startswith("variableStep"):
				currentChromosome = line.split(" ")[1].split("=")[1]
				chromosomes[currentChromosome] = {}
			else:
				position, forward, reverse = line.split("\t")
				chromosomes[currentChromosome][position] = {"forward": forward, "reverse": reverse}
	return chromosomes


#################################
# Writes wig file to target 	#
#  directory with given name	#
#################################

def write_wig(wig, fileName):
	print("Writing file:\t" + fileName)

	outFile = open(os.path.realpath( fileName ), 'w')

	for chromosome in wig:
		# write header for file
		outFile. write( "%s\n" % (
			"variableStep chrom=" + chromosome + " span=1"
		) )
		for key, value in sorted(wig[chromosome].items(), key=lambda x: int(x[0])):
			outFile.write( "%s %s\n" % (
				str(key),
				str(value)
				) )




#######################
#######################
###      tests      ###
#######################
#######################




##########################
### starts main script ###
##########################
main()