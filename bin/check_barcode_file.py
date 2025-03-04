#!/usr/bin/env python3

import sys
import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', 	'-i', default=None,	         								help='Input file')
parser.add_argument('--output', '-o', default=None,											help='Output file')
args = parser.parse_args()

def ensure_can_read(fname):
    if not (os.path.exists(fname) and os.path.isfile(fname) and os.access(fname, os.R_OK)):
        print(f"Unable to read file {fname}")
        exit(1)

def errx(message):
    print(message)
    exit(1)

def main():
    file_input = args.input
    file_output = args.output

    print(file_input)
    
    if file_input is None or file_output is None:
        errx("Usage: check_barcode_file.py --input barcodes.tsv --output checked_barcodes.tsv")

    nucl_pattern = re.compile(r"^[ATGCatgc]+$")
    repl_pattern = re.compile(r"^(.*)(_rep_\d+)$")

    ensure_can_read(file_input)

    exp_names = list()
    barcodes = list()
    replicates = dict()

    validated = dict()
    
    with open(file_input) as tsv:
        for line in tsv:

            # skip empty lines
            if len(line.lstrip()) == 0:
                continue

            if line != line.lstrip():
                errx("Invalid whitespace at beginning of line: " + line)

            # silently ignore trailing whitespace
            line = line.rstrip()

            splitted = line.split("\t")
            if (len(splitted) != 2):
                errx("Wrong number of columns in line: " + line)
            
            expname = splitted[0].replace('-','_')
            barcode = splitted[1]

            # check for duplicate experiment names
            if expname in exp_names:
                errx("Duplicate experiment name: " + expname)
            exp_names.append(expname)

            # count replicates for experiments
            repl_match = re.search(repl_pattern, expname)
            if repl_match is not None:
                exp = repl_match.group(1)
                if exp not in replicates.keys():
                    replicates[exp] = 1
                else:
                    replicates[exp] = replicates[exp] + 1

            # check barcode uniqueness
            if barcode in barcodes:
                errx("Duplicate barcode: " + barcode)
            barcodes.append(barcode)

            # only nucleotides allowed in barcode
            nucl_match = re.search(nucl_pattern, barcode)
            if nucl_match is None:
                errx("Invalid characters in barcode \"" + barcode + "\" for experiment " + expname)


            # assume ok
            validated[expname] = barcode

    # if any experiment occurs with only one replicate print warning to stdout and continue
    for exp in replicates.keys():
        if replicates[exp] == 1:
            print("Warning: Experiment " + exp + " has only one replicate.")


    # output
    with open(file_output, 'w') as outfile:
        for exp, bc in validated.items():
            outfile.write(f"{exp}\t{bc}\n")


if __name__== "__main__":
    main()
