#!/usr/bin/env python3

import argparse
import os
from typing import Dict,Any,List

#####################################################
# Calculates distribution of RNA-{regions,species}	#
#####################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', nargs='+', help='Input file')
parser.add_argument('--rna_species','-r',nargs='+',help='All RNA species to be included')
parser.add_argument('--output', '-o', help='Output file')
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

rna_species = args.rna_species


##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main():
	samples_collected = {}

	# iterate through every file given as input
	for file in args.input:
		[name_of_current_file, rna_of_current_file, _] = file.split(".")

		# Check if the current file contains reads assigned to a feature that is wanted in the analysis
		if not rna_of_current_file in rna_species:
			print("Warning: " + file + " has species not included in this analysis (not stated in --rna_species)")
			print("RNA species stated in --rna_species are: " + " ".join(rna_species))
			continue

		# Go thtough every line and split it
		with open(file,"r") as parse_file:
			for line in parse_file:
				line: str = line.strip()
				[read_name, aligned, number_alignments, transcript] = line.split("\t")

				if not read_name in samples_collected:
					samples_collected[read_name]: Dict[str,bool] = {}
					for rna in rna_species:
						samples_collected[read_name][rna]: bool = False
				if aligned == "Assigned":
					samples_collected[read_name][rna_of_current_file]: bool = True

	# initiate dictionary to count how many reads are assigned to which species
	rna_species_counts: Dict[str,int] = {}
	for species in rna_species:
		rna_species_counts[species]: int = 0
	read_count: int = 0
	not_assigned: int = 0
	ambiguous: int = 0

	for read_name, assignments in samples_collected.items():
		read_count += 1
		number_of_assignments = get_amount_of_assigned_rna_species(assignments)
		if number_of_assignments == 0:
			not_assigned += 1
		elif number_of_assignments == 1:
			assigned_species = [k for k,v in assignments.items() if v][0]
			rna_species_counts[assigned_species] += 1
		elif number_of_assignments >= 2:
			# IF it will be necessary to adapt the script to split the amounts of multimapped reads: here would probably the best part to implement it
			ambiguous += 1

	rna_species_counts["not_assigned"]: int = not_assigned
	rna_species_counts["ambiguous"]: int = ambiguous

	if not test_read_amounts(rna_species_counts, read_count):
		print("Warning: assignments and total read counts do not match:")
		print("Total count: " + str(read_count))
		print("RNA-species assigments: ")
		print(rna_species_counts)


	out_file = open(os.path.realpath(args.output), 'w')
	out_file.write("%s\t%s\t%s" % (
		"RNA_species",
		"number_assignments",
		"percentage"
	))
	for i in rna_species_counts:
		out_file.write("\n%s\t%s\t%s" % (
			i,
			str(rna_species_counts[i]),
			str("{:.2f}".format(get_percentage_amount(rna_species_counts[i], read_count)))
		))
	out_file.write("\n%s\t%s\t%s" % (
		"total",
		str(read_count),
		str("{:.2f}".format(get_percentage_amount(read_count,read_count)))
	))


#######################
#######################
###    functions    ###
#######################
#######################

def get_amount_of_assigned_rna_species( assignments ):
	count = 0
	for rna_species in assignments:
		if assignments[rna_species]:
			count = count + 1
	return count

def test_read_amounts( assignments, total_count ):
	count_species: int = 0
	for i in assignments:
		count_species += assignments[i]

	if count_species == total_count:
		return True
	else:
		return False

#TODO: finish method
def get_percentage_amount(div,total):
	return (div * 100)/total

##########################
### starts main script ###
##########################
main()

