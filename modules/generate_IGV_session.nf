/*
 * Prepares Annotation file in order to include it into the IGV
 */
process prepare_annotation_for_igv {
	publishDir "${params.output}", mode: 'copy', pattern: "${annotation.baseName}.sorted.gff.gz*"

	input:
	path(annotation)

	output:
	val("${annotation.baseName}.sorted.gff.gz"), 	emit: name_annotation
	path("${annotation.baseName}.sorted.gff.gz*"), 	emit: files_annotation_for_igv
	path("${task.process}.version.txt"), 			emit: version

	"""
	gff3sort.pl ${annotation} > ${annotation.baseName}.sorted.gff
	bgzip ${annotation.baseName}.sorted.gff
	tabix ${annotation.baseName}.sorted.gff.gz

	echo -e "${task.process}\tgff3sort.pl\tno_version_available" > ${task.process}.version.txt
    echo -e "${task.process}\tbgzip\t\$(bgzip --help 2>&1 | head -2 | tail -1 | cut -d' ' -f2)" >> ${task.process}.version.txt
	echo -e "${task.process}\ttabix\t\$(tabix --help 2>&1 | head -2 | tail -1 | cut -d' ' -f2)" >> ${task.process}.version.txt
	"""
}

/*
 * Generates an IGV-session that can be directly imported into the IGV
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
	path('igv-session.xml'), 				emit: xml_to_output_dir
	path("${task.process}.version.txt"), 	emit: version

	script:
	def local_annotation = !params.annotation == 'NO_FILE' ? '--annotation ' + gtf : ''

	"""
	generate-igv-session.py \
		--reference ${ref} \
		${local_annotation} \
		--input_path ${track_path} \
		--tracks ${tracks} \
		${bam} \
		--output igv-session.xml

	echo -e "${task.process}\tgenerate-igv-session.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}
