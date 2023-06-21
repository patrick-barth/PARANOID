/*
 * Prepares Annotation file in order to include it into the IGV
 * Input: [GTF]|[GFF3] Annotation file
 * Output: name_annotation 			-> [STR] Name of the resulting zipped annotation file]
 *		files_annotation_for_igv 	-> [MISC] All files generated by tabix to prepare the annotation file for the IGV
 */
process prepare_annotation_for_igv {
	publishDir "${params.output}", mode: 'copy', pattern: "${annotation.baseName}.sorted.gff.gz*"

	input:
	path(annotation)

	output:
	val("${annotation.baseName}.sorted.gff.gz"), emit: name_annotation
	path("${annotation.baseName}.sorted.gff.gz*"), emit: files_annotation_for_igv

	"""
	~/software/gff3sort/gff3sort.pl ${annotation} > ${annotation.baseName}.sorted.gff
	bgzip ${annotation.baseName}.sorted.gff
	tabix ${annotation.baseName}.sorted.gff.gz
	"""
}

/*
 * Generates an IGV-session that can be directly imported into the IGV
 * Input: 	[bigWIG] Peak files in bigWig format
 *			[BAM] Aligned sequences
 *			[STR] String pointing to the directory which contains the peak files 
 *			[FASTA] Reference sequence 
 *			[STR] Name of annotation file - if none is used 'NO_NAME' should be given here
 * Output: [XML] IGV-session to import into the IGV  
 */
process generate_igv_session {
	publishDir "${params.output}", mode: 'copy', pattern: 'igv-session.xml'

	input:
	path(tracks)
	path(bam)
	val(track_path)
	path(ref)
	val(annotation)

	output:
	path('igv-session.xml')

	script:
	if(params.annotation == 'NO_FILE')
		"""
		generate-igv-session.py --reference ${ref} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
	else
		"""
		generate-igv-session.py --reference ${ref} --annotation ${annotation} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
}