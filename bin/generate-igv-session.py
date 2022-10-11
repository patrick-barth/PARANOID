#!/usr/bin/env python3

#####################################
# Writes IGV-session file (XML)		#
#####################################

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--reference', 		'-r', type=str, 				help='Reference file')
parser.add_argument('--input_path',		'-i', type=str,			 		help='Path in which depicted files will be stored')
parser.add_argument('--tracks', 		'-t', type=str, nargs='+', 		help='Files displayed as tracks')
parser.add_argument('--output', 		'-o', type=str, 				help='Name of the output file')
parser.add_argument('--panel_height', 	'-p', type=int, default=938, 	help='Height of each track')
parser.add_argument('--panel_width', 	'-w', type=int, default=1901, 	help='Width of each track')
parser.add_argument('--font_size', 		'-f', type=int, default=10, 	help='Font size of track name')
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


##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main(reference,input_path,tracks,output,panel_height,panel_width,font_size):

	XML_string = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
	XML_string += '<Session genome="%s" hasGeneTrack="false" hasSequenceTrack="true" locus="All" version="8">\n' % (reference)

	# Write Resources block
	XML_string += '\t<Resources>\n'
	for i in tracks:
		file_name = os.path.basename(i)
		extension = os.path.splitext(i)[1][1:]
		XML_string += '\t\t<Resource path="./%s/%s" type="%s"/>\n' % (input_path,file_name,extension)
	XML_string += '\t</Resources>\n'

	# Writes Panels block -> contains information about tracks
	XML_string += '\t<Panel height="%s" name="DataPanel" width="%s">\n' % (panel_height,panel_width)
	XML_string += '\t\t<Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="%s" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="POSITIVE" shouldShowTranslation="false" visible="true"/>\n' % (font_size)

	for i in tracks:
		file_name 					= os.path.basename(i)
		file_name_without_extension = os.path.splitext(file_name)[0]
		extension 					= os.path.splitext(i)[1][1:]
		
		if extension in ['wig']:
			strand	= '+' if '_forward' in file_name else '-'
			maximum = max(get_values_wig(i)) 	if strand == '+' else 0
			minimum = 0 						if strand == '+' else min(get_values_wig(i))
			XML_string += '\t\t<Track attributeKey="%s" autoScale="true" clazz="org.broad.igv.track.DataSourceTrack" fontSize="%s" id="./%s/%s" name="%s" renderer="BAR_CHART" visible="true" windowFunction="mean">\n' % (file_name,font_size,input_path,file_name,file_name_without_extension)
			XML_string += '\t\t\t<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="%s" minimum="%s" type="LINEAR"/>\n' % (maximum,minimum)
			XML_string += '\t\t</Track>\n'
		else:
			print("not supposed to happen")

	XML_string += '\t</Panel>\n'

	XML_string += '</session>'

	out = open(output, 'w')
	out.write(XML_string)
	out.close

#######################
#######################
###    functions    ###
#######################
#######################

# returns all values as a sorted list of floats
def get_values_wig(wig_file):
	collect_values = []
	with open (wig_file) as wig_file:
		for line in wig_file:
			line = line.rstrip()
			if line.startswith("variableStep"):
				continue
			else:
				(position,value) = line.split(" ")
				collect_values.append(float(value))
	collect_values = sorted(collect_values)
	return(collect_values)



##########################
### starts main script ###
##########################
main(reference=args.reference,
	input_path=args.input_path,
	tracks=args.tracks,
	output=args.output,
	panel_height=args.panel_height,
	panel_width=args.panel_width,
	font_size=args.font_size)

