/*
 * Sorts a BAM file
 */
process sort_bam_before_strand_pref {
	tag {query.baseName}

	input:
	path(query)

	output:
	path("${query.baseName}.sorted.bam"), 	emit: bam_sorted
	path("${task.process}.version.txt"), 	emit: version

	"""
	samtools sort ${query} > ${query.baseName}.sorted.bam

	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Determines the amount of reads aligned to the forward and the reverse strand for each chromosome in the reference
 */
process determine_strand_preference {
	tag {name}
	publishDir "${params.output}/strand-distribution", mode: 'copy', pattern: "${name}.strand_proportion.txt"

	input:
	tuple val(name), path(query), path(reference)

	output:
	path("${name}.strand_proportion.txt"), 	emit: txt_strand_proportions
	path("${task.process}.version.txt"), 	emit: version

	"""
	egrep '^>' ${reference} | cut -f1 -d' ' | cut -c2- > references.txt
	for i in ${query}
	do
		samtools index \$i
	done
	touch ${name}.strand_proportion.txt
	echo -e "chromosome\tforward\treverse" >> ${name}.strand_proportion.txt
	while read r; do
  		echo -e "\$r\t\$(for i in *.sorted.bam; do samtools view -F 20 -q ${params.mapq} \$i \$r | wc -l; done | paste -s -d+ | bc)\t\$(for i in *.sorted.bam; do samtools view -f 16 -q ${params.mapq} \$i \$r | wc -l; done | paste -s -d+ | bc)" >> ${name}.strand_proportion.txt;
	done <references.txt

	echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Visualizes determined strand proportions into bar charts
 */
process visualize_strand_preference {
	publishDir "${params.output}/strand-distribution/visualization", mode: 'copy', pattern: "${strand.simpleName}.png"

	input:
	path(strand)

	output:
	path("${strand.simpleName}.png"), 		emit: png_to_output_dir
	path("${task.process}.version.txt"), 	emit: version

	"""
	visualize_strand_distribution.R \
		--input ${strand} \
		--output ${strand.simpleName} \
		--type png

	echo -e "${task.process}\tvisualize_strand_distribution.R\tcustom_script" >> ${task.process}.version.txt
	"""
}