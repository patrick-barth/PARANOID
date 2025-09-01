/*
 * Splits a BAM files into several smaller ones to reduce the overall computation time
 */
process split_bam_by_chromosome{
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("*.bam"), emit: bam_split

	"""
	bamtools split \
		-in ${query} \
		-reference
	"""
}

/*
 * Sorts alignments according to the coordinates and generates index files
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
	path(query), emit: ref_to_output_dir

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
	path(mapping)
	path(deduplication)

	output:
	path("multiqc_*"), emit: collected_stats_to_output_dir

	"""
	multiqc .
	"""
}

/*
 * Collects a variety of different metrics such as the command line used to execute the workflow and the parameters used
 */
process collect_workflow_metrics{
	publishDir "${params.output}/metadata", mode: "move"

	output:
	path("workflow_metrics.txt"), emit: output

	"""
	cat <<EOF > workflow_metrics.txt
    Author: ${params.manifest.author}
    Pipeline version: ${params.manifest.version}
    Nextflow version: ${nextflow.version}
    Working directory: ${workflow.workDir}
    Project directory: ${workflow.projectDir}
    User name: ${workflow.userName}
    Execution directory: ${workflow.launchDir}
    Executed script: ${workflow.scriptFile}
    Configuration file: ${workflow.configFiles}
    Configuration profile: ${workflow.profile}
    Command line: ${workflow.commandLine}
    Parameters: ${params}
    Unique session ID: ${workflow.sessionId}

    Container engine: ${workflow.containerEngine}    
    Containers used: ${workflow.container}

    Git repository: ${workflow.repository}
    Repository revision: ${workflow.revision}
    EOF
	"""
}

/*
 * Calculates the md5sum of all files given via input
 */
process get_md5sum {
    publishDir "${params.output}/metadata", mode: 'copy', pattern: "md5sums.txt"

    input:
    path(query)

    output:
    path("md5sums.txt"), emit: txt_to_output_dir

    script:
    """
    for i in ${query}
	do
		if test -f \$i; then
			md5sum \$i >> md5sums.txt
		fi
	done

    #echo -e "${task.process}\tmd5sum\t\$(md5sum --version | head -1 | rev | cut -f 1 -d' ' | rev)" > ${task.process}.version.txt
    """
}

/*
 * Collects all version outputs and merges them into a single file
 */
process collect_versions {
    publishDir "${params.output}/metadata", mode: 'copy', pattern: "tool_versions.txt"

    input:
    path(query)

    output:
    path('tool_versions.txt'), emit: txt_to_output_dir

    script:
    """
    echo -e "Process\ttool\tversion" > tool_versions.txt
    for i in ${query}
	do
		cat \$i >> tool_versions.txt
	done
    """
}
