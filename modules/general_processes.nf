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
 * Outputs the reference to the output directory
 * Input: [FASTA] Reference file 
 * Output: [FASTA] Reference file 
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
 * Collects statistics of a variety of tools and runs multiqc on them to summarize them
 * Input: 	[TXT] Report Adapter trimming from Trim Galore
 *			[TXT] Report Quality filtering from cutadapt
 *			[HTML]&[ZIP] General read report from FastQC (before preprocessing)
 *			[HTML]&[ZIP] General read report from FastQC (after preprocessing)	
 *			[LOG] Report experimental barcode splitting from fastx_barcode_splitter.pl
 *			[TXT] Report alignment of reads from Bowtie2 or STAR
 *			[TSV] Report PCR deduplication from umi_tools dedup
 * Output: [HTML] Report overview and [DIR] Directory containing information about the reports
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
 * Collects a variety of different metrics such as the command line used to execute the workflow and the parameters used
 * Params: params
 * Input: [VAL] Peak files
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