/*
 * Prepares index for bowtie2
 * Input: [FASTA] Reference file 
 * Output: Tuple of [FASTA] reference file and [BT2] index files generated for the reference   
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
 * Input: Tuple of [FASTA] reference file and [BT2] index files generated for the reference
 * 		[FASTA] Read files to be aligned to the reference
 * Output: bam_alignments 	-> [BAM] Aligned sequences
 * 		  report_alignments -> [TXT] Alignment reports
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
 * Input: [FASTA] Reference sequence
 *		[GTF]|[GFF3] Annotation file - if non is given it should say [NO_FILE]
 * Output: [DIR] Directory with index generated for given reference file
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
 * Output is provided as BAM file that is sorted by coordinates
 * Input: Tuple of [FASTQ] Read files to be aligned and [DIR] Directory containing the STAR index 
 * Output: bam_alignments -> [BAM] Aligned reads sorted by coordinates 
 *		report_alignments -> [TXT] Alignment report
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
 * Input: [BAM] Aligned sequences 
 * Output: bam_filtered_empty 	-> [BAM] Aligned sequences without empty files
 *		report_empty_alignments -> [TXT] Name of bam files without any alignments (one file per empty sample)
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
 * Collects all files filtered out via filter_empty_bams and writes their names in a single file.
 * Input: [TXT] Files containing names of samples without alignment 
 * Output: [TXT] Report file containing all names of samples without alignments  
 */
process collect_experiments_without_alignments {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'experiments-without-alignments.txt'

	input:
	file(query)

	output:
	file("experiments-without-alignments.txt"), optional: true

	"""
	if [[ ! \$(cat ${query} | wc -l) == 0 ]]; then
		cat ${query} > experiments-without-alignments.txt
	fi
	"""
}