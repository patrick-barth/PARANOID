#!/usr/bin/env python3

import argparse
import os
import re
import sys
import numpy
from Bio import SeqIO
from wig_files_writer import *

#################################################
# Counts the nucleotides occuring at cross-link	#
#  sites. Takes wig file as input and writes	#
#  results to stdout.							#
#  Allowes to define a cutoff for values that	#
#  are considered. 								#
#################################################


parser = argparse.ArgumentParser()
parser.add_argument('--input', 			'-i', type=str,	nargs='+', 					help='Input file(s). Supported formats: WIG, BED')
parser.add_argument('--reference', 		'-r', type=str,								help='Reference in fasta format from which the nucleotides are read')
parser.add_argument('--length', 		'-l', type=int,				default=10,		help='Number of nucleotides to each side of the cross-link')
parser.add_argument('--percentile', 	'-p', type=float,			default=90,		help='Percentile used to calculate the cutoff')
parser.add_argument('--output', 		'-o', type=str,								help='Output file for extracted sequences')
parser.add_argument('--outfmt_fasta', 	'-f', action='store_true',	default=False,  help='If true, sequences will be put out in fasta format')
parser.add_argument('--omit_cl', 		'-u', action='store_true',	default=False,  help='BOOLEAN: If true nucleotide at cross/link site will be omitted')
parser.add_argument('--generate_bed', 	'-b', 										help='If used a BED file is generated')
args = parser.parse_args()

#######################
#######################
###    variables    ###
#######################
#######################

countOutputs = 0 
	
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

def main(input,reference,length,percentile,output,outfmt_fasta,omit_cl):
	countOutOfBounds = 0
	referenceSequences = parse_fasta(reference)

	# calculate percentile cutoff
	wigFiles = []
	for file in input:
		wigFiles.append( parse_wig( file ) )

	allCounts = bundle_counts( wigFiles ) 

	numpy.set_printoptions(threshold = sys.maxsize)
	cutoff = numpy.percentile(allCounts, percentile)
	print("The value cutoff is set to " + str(cutoff) )

	if( generate_bed ):
		bed_file = open(generate_bed, 'w', encoding="utf-8")


	for file in input:
		parsedFile = parse_wig(file)
		for chromosome in sorted(parsedFile):
			# check if the current chromosome of the wig file is found in the reference. If not it is skipped and a warning is printed
			if chromosome not in referenceSequences.keys():
				print("Warning: " + chromosome + " from file " + file + " not found in Reference file")
				print("Reference contains following chromosomes: " + str(referenceSequences.keys()))
				continue
			for entry in parsedFile[chromosome]:
				if( abs(parsedFile[chromosome][entry]) >= cutoff ):
					#check if cross-link is not too close to the end or the beginning of the reference
					if int(entry) - length - 1 >= 0 and int(entry) + length <= len(referenceSequences[chromosome]) : 
						extractionStart = int(entry) - length - 1
						extractionEnd 	= int(entry) + length # no '-1' since the last number is exclusive when getting a substring
						strand 			= "+" if parsedFile[chromosome][entry] > 0 else "-"

						if(generate_bed):
							bed_file.write(chromosome + "\t" + str(extractionStart) + "\t" + str(extractionEnd) + "\tname\t0\t" + strand + "\n")

						if(omit_cl):
							extractedSequence = referenceSequences[chromosome][extractionStart:int(entry)-1] + 'n' + referenceSequences[chromosome][int(entry):extractionEnd]
						else:
							extractedSequence = referenceSequences[chromosome][extractionStart:extractionEnd]

						if(strand == "-"):
							extractedSequence = extractedSequence.reverse_complement()

						write_sequence(extractedSequence, output)
					else:
						countOutOfBounds += 1
	
	if( generate_bed ):
		bed_file.close()

	print(str(countOutOfBounds) + " peaks were too close to one of the chromosome ends and were thus ignored")



					

			
		

#######################
#######################
###    functions    ###
#######################
#######################

def parse_fasta( fastaFile ):
	sequences = {}
	with open (fastaFile, "r") as fastaFile:
		for record in SeqIO.parse(fastaFile, "fasta"):
			sequences[record.id] = record.seq 
	return sequences

def bundle_counts(wigFiles):
	allCounts = []
	floatCounts = numpy.array(allCounts, dtype = numpy.float32)
	for file in wigFiles:
		for chromosome in file:
			for pos in file[chromosome]:
				if not file[chromosome][pos] == 0:
					floatCounts = numpy.append(floatCounts, abs(file[chromosome][pos]))
	return floatCounts

def write_sequence(sequence, outputFile):
	global countOutputs
	with open(outputFile, 'a') as file:
		if args.outfmt_fasta:
			file.write('>sequence' + str(countOutputs) + os.linesep )
			countOutputs += 1 
		file.write(str(sequence) + os.linesep)
		

##########################
### starts main script ###
##########################
main(input=args.input,
	reference=args.reference,
	length=args.length,
	percentile=args.percentile,
	output=args.output,
	outfmt_fasta=args.outfmt_fasta,
	omit_cl=args.omit_cl)

