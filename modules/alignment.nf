/*
 * Prepares index for bowtie2
 */
process build_index_bowtie {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*")

	"""
	bowtie2-build ${ref} ${ref}
	"""
}

/*
 * Alignes reads to a reference via bowtie2, filters out unaligned sequeces
 * and converts output to BAM. Reference needs to be indexed.
 */
process mapping_bowtie{
	tag {query.simpleName}

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path "${query.baseName}.bam", emit: bam_alignments
	path "${query.simpleName}.statistics.txt", emit: report_alignments

	"""
	bowtie2 --no-unal -q -p ${task.cpus} -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
	"""
}

/*
 * Prepares index for STAR. Includes alignment if one is given
 */
process build_index_STAR {

	input:
	path(referenceGenome)
	path(gtf)

	output:
	path(index)

	script:
	if(params.annotation == 'NO_FILE')
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} 
		"""
	else
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${gtf}
		"""
}

/*
 * Alignes reads to reference via STAR. Reference needs to be index before.
 * Suppresses clipping at 5' end (--alignEndsType Extend5pOfRead1)
 * Output is provided as BAM file that is sorted by coordinates (--outSAMtype BAM SortedByCoordinate)
 */
process mapping_STAR{
	tag {query.simpleName}

	input:
	tuple path(query), path(indexDir)

	output:
	path("${query.baseName}.Aligned.sortedByCoord.out.bam"), emit: bam_alignments
	path("${query.baseName}.Log.*"), emit: report_alignments

	"""
	STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --alignEndsType Extend5pOfRead1 --outSAMtype BAM SortedByCoordinate
	"""
}

/*
 * Filters out BAM files without any valid alignments
 * Names of filtered out files are collected and forwarded 
 */
process filter_empty_bams{
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.filtered.bam"), emit: bam_filtered_empty, optional: true
	path("${query.simpleName}.no_alignments.txt"), emit: report_empty_alignments, optional: true
	"""
	if [[ \$(samtools view ${query} | wc -l) == 0 ]]; then
		echo "${query.simpleName} contains no alignable reads" > ${query.simpleName}.no_alignments.txt
	else
		mv ${query} ${query.baseName}.filtered.bam
	fi
	"""
}

/*
 * Collects all files filtered out via filter_empty_bams and writes
 * their names in a single file.
 */
process collect_experiments_without_alignments {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'experiments-without-alignments.txt'

	input:
	file(query) from log_experiments_without_alignments.flatten().toList()

	output:
	file("experiments-without-alignments.txt") optional true into output_log_exp_without_alignment

	"""
	if [[ ! \$(cat ${query} | wc -l) == 0 ]]; then
		cat ${query} > experiments-without-alignments.txt
	fi
	"""
}