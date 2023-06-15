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