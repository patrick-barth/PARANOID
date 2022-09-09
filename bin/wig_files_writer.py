#!/usr/bin/env python3

import os

def parse_wig( wigFile ):
	chromosomes = {}
	currentChromosome = ""
	with open ( wigFile, 'r', encoding='utf-8' ) as wigFile:
		for line in wigFile:
			# a file might contain several chromosomes. New chromosomes start with a line like "variableStep chrom=DQ380154.1 span=1". Due to that they are divided by this line. 
			# Also the chromosome name is parsed from that line
			if line.startswith("variableStep"):
				currentChromosome = line.split(" ")[1].split("=")[1]
				chromosomes[currentChromosome] = {}
			else:
				position, count = line.split(" ")
				chromosomes[currentChromosome][position] = float(count)
	return chromosomes

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
				chromosomes[currentChromosome][position] = {"forward": float(forward), "reverse": float(reverse)}
	return chromosomes

def write_wig2(wig, fileName):
	print("Writing file:\t" + fileName)

	outFile = open(os.path.realpath( fileName ), 'w')

	for chromosome in wig:
		# write header for file
		outFile. write( "%s\n" % (
			"variableStep chrom=" + chromosome + " span=1"
		) )
		for key, value in sorted(wig[chromosome].items(), key=lambda x: int(x[0])):
			outFile.write( "%s\t%s\t%s\n" % (
				str(key),
				str(wig[chromosome][key]["forward"]),
				str(wig[chromosome][key]["reverse"])
			) )

def splitForwardReverse(wig2):
	
	splitWig = []

	for wig in wig2:
		forward = {}
		reverse = {}
		for chromosome in wig:
			for position in wig[chromosome]:
				for strand in wig[chromosome][position]:
					if strand == "forward":
						if not forward:
							forward[chromosome] = {position: wig[chromosome][position]["forward"]}
						elif not chromosome in forward.keys():
							forward[chromosome] = {position: wig[chromosome][position]["forward"]}
						else:
							forward[chromosome][position] = wig[chromosome][position]["forward"]
					elif strand == "reverse":
						if not reverse:
							reverse[chromosome] = {position: wig[chromosome][position]["reverse"]}
						elif not chromosome in reverse.keys():
							reverse[chromosome] = {position: wig[chromosome][position]["reverse"]}
						else:
							reverse[chromosome][position] = wig[chromosome][position]["reverse"]
		splitWig.append(forward)
		splitWig.append(reverse)
	return splitWig