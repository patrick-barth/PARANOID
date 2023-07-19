#!/usr/bin/env python3

import argparse
import os
import sys
from typing import Dict, List
import numpy
from Bio import SeqIO
from wig_files_writer import *

#################################################
# Counts nucleotides occurring at cross-link	#
#  sites. Takes wig file as input and writes	#
#  results to stdout.							#
#  Allows to define a cutoff for values that	#
#  are considered. 								#
#################################################


parser = argparse.ArgumentParser()
parser.add_argument('--input', 			'-i', type=str,	nargs='+', 					help='Input file(s). Supported formats: WIG, BED')
parser.add_argument('--reference', 		'-r', type=str,								help='Reference in fasta format from which the nucleotides are read')
parser.add_argument('--length', 		'-l', type=int,				default=20,		help='Number of nucleotides to each side of the cross-link')
parser.add_argument('--percentile', 	'-p', type=float,			default=90,		help='Percentile used to calculate the cutoff')
parser.add_argument('--output', 		'-o', type=str,								help='Output file for extracted sequences')
parser.add_argument('--omit_cl', 		'-u', action='store_true',	default=False,  help='BOOLEAN: If true nucleotide at cross/link site will be omitted')
parser.add_argument('--omit_width',		'-w', type=int,				default=0,  	help='INT: Allows to omit nucleotides besides the CL site. Only works when --omit_cl is stated.')
parser.add_argument('--generate_bed', 	'-b', 										help='If used a BED file is generated')
args = parser.parse_args()

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(input,reference,length,percentile,output,omit_cl,omit_width):
	referenceSequences: Dict 	= parse_fasta(reference)
	extension:			str 	= os.path.splitext(input[0])[1][1:]
	# Check if all files have the same ending
	if not all([x.endswith(extension) for x in input]):
		errx("Input files contain different extensions")
	# Check if omitted nucleotides don't exceed the total extracted length
	if args.omit_width >= args.length:
		errx("Length of omitted nucleotides exceeds the actual length of extracted nucleotides")

	chromosomes: Dict[str,List[Dict[str,str]]] = {}
	for k in sorted(referenceSequences):
		chromosomes[k]: List[tuple[int,str]] = []

	if 		extension in ['wig']:
		# Calculate cutoff (only for wig files)
		wigFiles: List[Dict[Dict[str,float]]] = []
		for file in input:
			wigFiles.append(parse_wig(file))
			collectedPeakValues: List[float] = bundle_counts(wigFiles=wigFiles)
		numpy.set_printoptions(threshold = sys.maxsize)
		cutoff = numpy.percentile(collectedPeakValues, percentile)
		print("The value cutoff is set to " + str(cutoff) )

		for file in wigFiles:
			firstValue: int = list(list(file.values())[0].values())[0]
			strand: 	str = '+' if firstValue > 0 else '-'
			for chromosome in file:
				if not chromosome in chromosomes:
					errx('Chromosome %s not found in reference' % (chromosome))
				for position,value in file[chromosome].items():
					if abs(value) >= cutoff:
						chromosomes[chromosome].append((int(position),strand))

	elif 	extension in ['bed']:
		for file in input:
			with open (file, 'r') as file:
				for line in file:
					(chromosome,pos1,pos2,_,_,strand,_) = line.split("\t")
					if not chromosome in chromosomes:
						errx('Chromosome %s not found in reference' % (chromosome))
					chromosomes[chromosome].append((int(pos1),strand))

	# Extract sequences around determined peaks
	FASTA_output: 		str = ''
	sequenceCounter: 	int = 1
	for chromosome,entries in chromosomes.items():
		reference_current: 	str = referenceSequences[chromosome]
		reference_length: 	int = len(reference_current)
		for position,strand in entries:
			if position - length - 1 < 0 or position + length > reference_length:
				continue
			sequenceStart: 	int = position - length - 1
			sequenceEnd: 	int = position + length
			sequence:		str = reference_current[sequenceStart:sequenceEnd]
			if omit_cl == True:
				sequence = sequence[:length - omit_width] + 'n' + ('n' * omit_width * 2) + sequence[length + omit_width + 1:]
			if 	strand == '-':
				sequence = sequence.reverse_complement()
			FASTA_output += '>sequence%s chromosome=%s position=%s strand=%s\n' % (sequenceCounter,chromosome,position,strand)
			FASTA_output += '%s\n' % (sequence)
			sequenceCounter += 1

	out=open(output,'w')
	out.write(FASTA_output)
	out.close()

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

def errx(message):
    print(message)
    exit(1)

##########################
### starts main script ###
##########################
main(input=args.input,
	reference=args.reference,
	length=args.length,
	percentile=args.percentile,
	output=args.output,
	omit_cl=args.omit_cl,
	omit_width=args.omit_width)