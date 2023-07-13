#!/usr/bin/env python3

from __future__ import division
import argparse
import os
import re
from wig_files_writer import *

parser = argparse.ArgumentParser()

parser.add_argument('--wig',    '-w', nargs=2, metavar=('forward','reverse'),   help='input files in wig format; First the forward and then the reverse file. Both should be from the same sample')
parser.add_argument('--output', '-o',                                           help='prefix for output files')
args = parser.parse_args()

def main():
	replicate = args.wig

	parsedReplicate = { "forward": parse_wig( replicate[0] ),
					"reverse": parse_wig( replicate[1] ) }
	#print(parsedReplicate)
	peaks = {}

	for chromosome in parsedReplicate['forward']:
		print(chromosome)
		if not chromosome in peaks:
			peaks[chromosome] = {}
		for position in parsedReplicate['forward'][chromosome]:
			if not position in peaks[chromosome]:
				peaks[chromosome][position] = {	'forward': parsedReplicate['forward'][chromosome][position],
												'reverse': 0}
			else:
				peaks[chromosome][position]['forward'] = parsedReplicate['forward'][chromosome][position]
	
	for chromosome in parsedReplicate['reverse']:
		if not chromosome in peaks:
			peaks[chromosome] = {}
		for position in parsedReplicate['reverse'][chromosome]:
			if not position in peaks[chromosome]:
				peaks[chromosome][position] = {	'forward': 0,
												'reverse': parsedReplicate['reverse'][chromosome][position]}
			else:
				peaks[chromosome][position]['reverse'] = parsedReplicate['reverse'][chromosome][position]
	
	if not len(peaks) == 0:
		write_wig2(peaks, args.output)

main()