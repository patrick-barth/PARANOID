process get_chromosome_sizes{
	input:
	path(ref)

	output:
	path("${ref.simpleName}.chromosome_sizes.txt")

	"""
	samtools faidx ${ref}
	cut -f1,2 ${ref}.fai > ${ref.simpleName}.chromosome_sizes.txt
	"""
}

process calculate_crosslink_sites{
	tag {query.simpleName}
	publishDir "${params.output}/cross-link-sites/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

	input:
	tuple path(query), path(chrom_sizes)
	output:
	path("${query.simpleName}.wig2"), emit: wig2_cross_link_sites, optional: true
	tuple val("cross-link-sites"), path("${query.simpleName}_forward.wig"), path("${query.simpleName}_reverse.wig"), emit: wig_cross_link_sites_split, optional: true
	path("${query.simpleName}_{forward,reverse}.wig"), emit: wig_cross_link_sites, optional: true

	"""
	create-wig-from-bam.py --input ${query} --mapq ${params.mapq} --chrom_sizes ${chrom_sizes} --output ${query.simpleName}.wig2
	if [[ -f "${query.simpleName}.wig2" ]]; then
		wig2-to-wig.py --input ${query.simpleName}.wig2 --output ${query.simpleName}
	fi
	"""
}

process merge_wigs{
    tag {name}
    publishDir "${params.output}/cross-link-sites-merged/wig", mode: 'copy', pattern: "${name}_{forward,reverse}.wig"

    input:
    tuple val(name), path(query)

    output:
    path("${name}.wig2"), emit: wig2_merged
    tuple val("cross-link-sites-merged"), path("${name}_forward.wig"), path("${name}_reverse.wig"), emit: wig_merged_cross_link_sites
    path("${name}_{forward,reverse}.wig")

    script:
    if(name !=~ "unmatched*")
        """
        merge-wig.py --wig ${query} --output ${name}.wig2
        wig2-to-wig.py --input ${name}.wig2 --output ${name}
        """
    else
        """
        merge-wig.py --wig ${query} --output unmatched.wig2
        wig2-to-wig.py --input unmatched.wig2 --output unmatched
        """
}

process split_wig2_for_correlation{
    tag {query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}_forward.wig"), emit: wig_split_forward, optional: true
    path("${query.simpleName}_reverse.wig"), emit: wig_split_reverse, optional: true
    path("${query.simpleName}_{forward,reverse}.wig"), emit: wig_split_both_strands, optional: true

    """
    wig2-to-wig.py --input ${query} --output ${query.simpleName}
    """
}

process calc_wig_correlation{

    tag {name}
    publishDir "${params.output}/correlation_of_replicates", mode: 'copy', pattern: "${name}_${strand}_correlation.{png,csv}"

    input:
    tuple val(name),path(query),val(strand),path(chrom_sizes)

    output:
    path("${name}_${strand}_correlation.png"), optional: true
    path("${name}_${strand}_correlation.csv"), optional: true

    script:
    String[] test_size = query
    """
    if [[ ${test_size.size()} > 1 ]]; then
        wig_file_statistics.R --input_path . --chrom_length ${chrom_sizes} --output ${name}_${strand} --type png
    fi
    """
}

process wig_to_bigWig{
	tag {forward.simpleName}
	publishDir "${params.output}/${out_dir}/bigWig", mode: 'copy', pattern: "*.bw"

	input:
	tuple val(out_dir), path(forward), path(reverse), path(chrom_sizes)

	output:
	path("*.bw"), optional: true
	tuple val(out_dir), path("*.bw"), emit: bigWig_both_strands, optional: true
	tuple val(out_dir), path("${reverse.baseName}.bw"), emit: bigWig_reverse, optional: true
	tuple val(out_dir), path("${forward.baseName}.bw"), emit: bigWig_forward, optional: true

	"""
	if [[ \$(cat ${forward} | wc -l) > 1 ]]; then
		wigToBigWig ${forward} ${chrom_sizes} ${forward.baseName}.bw
	fi
	if [[ \$(cat ${reverse} | wc -l) > 1 ]]; then
		wigToBigWig ${reverse} ${chrom_sizes} ${reverse.baseName}.bw
	fi
	"""
}

process bigWig_to_bedgraph{
	tag {bigWig.simpleName}
	publishDir "${params.output}/${out_dir}/bedgraph", mode: 'copy', pattern: "*.bedgraph"

	input:
	tuple val(out_dir), path(bigWig)

	output:
	path("*.bedgraph")

	"""
	bigWigToBedGraph ${bigWig} ${bigWig.baseName}.bedgraph
	"""
}

process index_for_peak_calling {
    tag{query.simpleName}

    input:
    path(query)

    output:
    tuple path("${query.simpleName}.sorted.bam"), path("${query.simpleName}.sorted.bam.bai")

    """
    samtools sort ${query} > ${query.simpleName}.sorted.bam
    samtools index ${query.simpleName}.sorted.bam
    """
}

process pureCLIP {
    tag{bam.simpleName}
    cache false
    errorStrategy 'ignore' //TODO: is supposed to be only temporal. Need to find a solution for: ERROR: Emission probability became 0.0! This might be due to artifacts or outliers.

    publishDir "${params.output}/peak_calling", mode: 'copy', pattern: "${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed"

    input:
    tuple path(bam), path(bai), path(ref)

    output:
    path("${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed")
    path("${bam.simpleName}.pureCLIP_crosslink_sites.bed"), emit: bed_crosslink_sites
    path("${bam.simpleName}.pureCLIP_crosslink_sites.params"), emit: report_pureCLIP

    script:
    if(params.peak_calling_for_high_coverage == true && params.peak_calling_regions == true)
        """
        pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -or ${bam.simpleName}.pureCLIP_crosslink_regions.bed -dm ${params.peak_calling_regions_width} -mtc 5000 -mtc2 5000 -ld
        """
    else if(params.peak_calling_for_high_coverage == true && params.peak_calling_regions == false)
        """
        pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -mtc 5000 -mtc2 5000 -ld
        """
    else if(params.peak_calling_for_high_coverage == false && params.peak_calling_regions == true)
        """
        pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -or ${bam.simpleName}.pureCLIP_crosslink_regions.bed -dm ${params.peak_calling_regions_width}
        """
    else if(params.peak_calling_for_high_coverage == false && params.peak_calling_regions == false)
        """
        ln -s ${ref} ${ref.baseName}_symlink.fasta
        pureclip -i ${bam} -bai ${bai} -g ${ref.baseName}_symlink.fasta -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed
        """
}

process split_wig_2_for_peak_height_hist {
	tag {query.simpleName}

	input:
	path(query)

	output:
	tuple val("${query.simpleName}"), path("${query.simpleName}_forward.wig"), path("${query.simpleName}_reverse.wig"), optional: true

	"""
	wig2-to-wig.py --input ${query} --output ${query.simpleName}
	"""
}

process generate_peak_height_histogram {
	tag {query}
	publishDir "${params.output}/peak_height_distribution", mode: 'copy'

	input:
	tuple val(query), path(forward), path(reverse)

	output:
	path("${query}.png")

	"""
	generate-peak-height-histogram.R --input . --output ${query} --type png --color "${params.color_barplot}" --percentile ${params.percentile}
	"""
}