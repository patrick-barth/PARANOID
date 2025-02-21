process build_index_bowtie {

	input:
	path(ref)

	output:
	tuple path("${ref}"), path("${ref}.*"), emit: index
	path("${task.process}.version.txt"), 	emit: version

	"""
	bowtie2-build ${ref} ${ref}

	echo -e "${task.process}\tbowtie2\t\$(bowtie2-build --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

process mapping_bowtie{
	tag {query.simpleName}

	input:
	tuple path(ref), path(index)
	path(query)

	output:
	path("${query.baseName}.bam"), 				emit: alignments
	path("${query.simpleName}.statistics.txt"), emit: report
	path("${task.process}.version.txt"), 		emit: version

	script:
	def local_all_alignments 	= params.report_all_alignments ? '-a' : ''
	//def some_alignments = params.max_alignments && !params.report_all_alignments ? "-k " + params.max_alignments : '' //TODO: Add when max_alignments is changed

	"""
	bowtie2 \
		--no-unal \
		-q \
		${local_all_alignments} \
		-k ${params.max_alignments} \
		-p ${task.cpus} \
		--seed 0 \
		-U ${query} \
		-x ${ref} \
		2> ${query.simpleName}.statistics.txt | samtools view -bS - \
		> ${query.baseName}.bam

	echo -e "${task.process}\tbowtie2\t\$(bowtie2 --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
		"""
}

process build_index_STAR {

	input:
	path(referenceGenome)
	path(gtf)

	output:
	path(index), 							emit: index
	path("${task.process}.version.txt"), 	emit: version

	script:
	def local_annotation = !params.annotation == 'NO_FILE' ? '--sjdbGTFfile ' + gtf : ''

	"""
	mkdir index
	STAR \
		--runThreadN ${task.cpus} \
		--runMode genomeGenerate \
		--genomeDir ./index \
		--genomeFastaFiles ${referenceGenome} \
		${local_annotation}
	
	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}

process mapping_STAR{
	tag {query.simpleName}

	input:
	tuple path(query), path(indexDir)

	output:
	path("${query.baseName}.Aligned.sortedByCoord.out.bam"), 	emit: alignments
	path("${query.baseName}.Log.*"), 							emit: report
	path("${task.process}.version.txt"), 						emit: version

	script:
	def local_all_alignments = params.report_all_alignments ? '--outSAMmultNmax -1' : ''
	//def some_alignments = params.max_alignments && !params.report_all_alignments ? "--outSAMmultNmax " + params.max_alignments : '' //TODO: Add when max alignment is changed

	"""
	STAR \
		--runThreadN ${task.cpus} \
		--genomeDir ${indexDir} \
		--readFilesIn ${query} \
		--outFileNamePrefix ${query.baseName}. \
		--alignEndsType Extend5pOfRead1 \
		${local_all_alignments} \
		--outSAMtype BAM SortedByCoordinate

	echo -e "${task.process}\tSTAR\t\$(STAR --version)" > ${task.process}.version.txt
	"""
}

process filter_empty_bams{
	tag {query.simpleName}

	input:
	path(query)

	output:
	path("${query.baseName}.filtered.bam"), 		emit: bam_filtered, optional: true
	path("${query.simpleName}.no_alignments.txt"), 	emit: report_empty, optional: true
	"""
	if [[ \$(samtools view ${query} | wc -l) == 0 ]]; then
		echo "${query.simpleName} contains no alignable reads" > ${query.simpleName}.no_alignments.txt
	else
		mv ${query} ${query.baseName}.filtered.bam
	fi
	"""
}

process collect_experiments_without_alignments {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'experiments-without-alignments.txt'

	input:
	path(query)

	output:
	path("experiments-without-alignments.txt"), emit: output, optional: true

	"""
	if [[ ! \$(cat ${query} | wc -l) == 0 ]]; then
		cat ${query} > experiments-without-alignments.txt
	fi
	"""
}
