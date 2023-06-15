/*
 * Splits a BAM files into several smaller ones to reduce the overall computation time
 * Input: [BAM] Alignments 
 * Output: [BAM] Several smaller alignment files 
 */
process split_bam_by_chromosome{
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("*.bam")

	"""
	bamtools split -in ${query} -reference
	"""
}

/*
 * Sorts alignments according to the coordinates and generates index files
 * Input: [BAM] Alignments 
 * Output: bam_sorted_alignments -> [BAM] Alignments sorted by coordinates
 * 	bai_index_files -> [BAI] Index file 
 */
process sort_and_index_alignment{
	tag {query.simpleName}
	publishDir "${params.output}/alignments", mode: 'copy', pattern: "${query.simpleName}.sorted.bam*"

	input:
	path(query)

	output:
	path("${query.simpleName}.sorted.bam"), emit: bam_sorted_alignments
    path("${query.simpleName}.sorted.bam.bai"), emit: bai_index_files 

	"""
	samtools sort ${query} > ${query.simpleName}.sorted.bam
	samtools index ${query.simpleName}.sorted.bam
	"""
}

/*
 * Calculates distances between peaks
 * Params: params.percentile -> Percentile of peaks being omited in order to filter out background noise. Omits peaks according to their height
 *		params.distance -> Maximum distance between peaks allowed
 * Input: [WIG2] Peak files
 * Output: [TSV] Tab separated files showing the amounts of found distances
 */
process output_reference {
	publishDir "${params.output}", mode: 'copy', pattern: "${query}"

	input:
	path(query)

	output:
	path(query)

	"""
	"""
}

/*
 * Calculates distances between peaks
 * Params: params.percentile -> Percentile of peaks being omited in order to filter out background noise. Omits peaks according to their height
 *		params.distance -> Maximum distance between peaks allowed
 * Input: [WIG2] Peak files
 * Output: [TSV] Tab separated files showing the amounts of found distances
 */
process multiqc{
	publishDir "${params.output}/statistics", mode: 'move'

	input:
	path(adapter)
	path(qual)
	path(qc1)
	path(qc2)
	path(split)
	path(mapping)
	path(deduplication)

	output:
	path("multiqc_*")

	"""
	multiqc .
	"""
}

/*
 * Calculates distances between peaks
 * Params: params.percentile -> Percentile of peaks being omited in order to filter out background noise. Omits peaks according to their height
 *		params.distance -> Maximum distance between peaks allowed
 * Input: [WIG2] Peak files
 * Output: [TSV] Tab separated files showing the amounts of found distances
 */
process collect_workflow_metrics{
	publishDir "${params.output}/execution_metrics", mode: "move"

	input:
	val(command_line)

	output:
	path("execution_information.txt"), emit: execution_information
	path("container_information.txt"), emit: container_information
	path("parameter_information.txt"), emit: parameter_information

	"""
	echo "Execution command: ${command_line}" > execution_information.txt
	echo "Project directory: ${workflow.projectDir}" >> execution_information.txt
	echo "Path to executed script: ${workflow.scriptFile}" >> execution_information.txt
	echo "Configuration file: ${workflow.configFiles}" >> execution_information.txt
	echo "Configuration profile: ${workflow.profile}" >> execution_information.txt
	echo "Nextflow version: ${nextflow.version}" >> execution_information.txt
	echo "Workflow version: ${workflow.manifest.version}" >> execution_information.txt
	echo "Execution directory: ${workflow.launchDir}" >> execution_information.txt
	echo "Unique session ID: ${workflow.sessionId} >> execution_information.txt"

	echo "Container engine used: ${workflow.containerEngine}" > container_information.txt
	echo "Containers used: ${workflow.container}" >> container_information.txt

	echo "${params}" >> parameter_information.txt
	"""
}