/*
 * Retrieves sizes of all chromosomes in reference file 
 * Input: [FASTA] Reference file 
 * Output: [TXT] Name and size of all chromosomes 
 */
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

/*
 * Calculates cross-link sites from BAM files and returns them as 2 WIG files - one for the forward and one for the reverse strand.
 * Forward cross-link sites are determined by StartPosition -1
 * Reverse reads are dertermined by StartPosition + CIGAR string
 *  Additionally, alignments are filtered by the mapq score stated via params.mapq
 * Input: Tuple of [BAM] aligned reads and [TXT] file stating chromosome names and sizes 
 * Params: params.mapq -> MAPQ-score to filter out low quality alignments
 * Output:  wig2_cross_link_sites       -> [WIG2] Cross-link sites 
 *          wig_cross_link_sites_split  -> Tuple of [STR] saying "cross-link-sites", [WIG] forward peaks and [WIG] reverse peaks
 *          wig_cross_link_sites        -> [WIG] Cross-link sites with strands being divided
 */
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

/*
 * Performs peak calling via PureCLIP. Can do single peaks or peak regions.
 * Input: Tuple of [BAM] alignment file, [BAI] index file and [FASTA] reference file
 * Params:  params.peak_calling_for_high_coverage   -> Enables special settings I found useful when analysing data with huge amount of peaks all across the reference (-mtc 5000 -mtc2 5000 -ld)
 *          params.peak_calling_regions             -> Can place determined cross-link sites in close proximity into a single cross-link region
 *          params.peak_calling_regions_width       -> Determines the maximum allowed width of cross-link regions
 * Output:  bed_crosslink_sites -> [BED] Cross-link sites determined by PureCLIP
 *          report_pureCLIP     -> [PARAMS] Report of PureCLIP
 */
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

/*
 * Transforms output from pureCLIP to WIG format
 */


process pureCLIP_to_wig{
    tag {query.simpleName}
    publishDir "${params.output}/peak_calling/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

    input:
    path(query)

    output:
    path("${query.simpleName}.wig2"), emit: wig2_peak_called_cl_sites, optional: true
    tuple val("peak_calling"), path("${query.simpleName}_forward.wig"), emit: wig_peak_called_cl_sites_forward, optional: true
    tuple val("peak_calling"), path("${query.simpleName}_reverse.wig"), emit: wig_peak_called_cl_sites_reverse, optional: true
    path('empty_sample,txt'), emit: txt_empty_sample, optional: true

    """
    if [[ "${query.size()}" > 0 ]]; then
        # split into strands
        awk '{if (\$6 == "+") {print \$0 > "${query.simpleName}_forward.bed";} else {print \$0 > "${query.simpleName}_reverse.bed";}}' ${query}
        # Add 1 into the 4th column (since it#s from peak calling the actual size of the peak does not matter anymore)
        #  Then extracts important columns and converts the data into a valif WIG file
        if [[ -e "${query.simpleName}_forward.bed" ]]; then
            awk 'BEGIN {FS = "\\t";OFS = "\\t";}{\$NF = \$NF "\\t1";print \$0;}' ${query.simpleName}_forward.bed > ${query.simpleName}_forward.1.bed
            awk 'BEGIN {FS = "\\t";}{if (prev_header != \$1) {prev_header = \$1; printf("variableStep chrom=%s span=1\\n", \$1);}printf("%s %s\\n", \$3, \$8);}' ${query.simpleName}_forward.1.bed > ${query.simpleName}_forward.wig
        else
            touch ${query.simpleName}_forward.wig
        fi
        if [[ -e "${query.simpleName}_reverse.bed" ]]; then
            awk 'BEGIN {FS = "\\t";OFS = "\\t";}{\$NF = \$NF "\\t-1";print \$0;}' ${query.simpleName}_reverse.bed > ${query.simpleName}_reverse.1.bed
            awk 'BEGIN {FS = "\\t";}{if (prev_header != \$1) {prev_header = \$1; printf("variableStep chrom=%s span=1\\n", \$1);}printf("%s %s\\n", \$3, \$8);}' ${query.simpleName}_reverse.1.bed > ${query.simpleName}_reverse.wig
        else
            touch ${query.simpleName}_reverse.wig
        fi
        wig-to-wig2.py --wig ${query.simpleName}_forward.wig ${query.simpleName}_reverse.wig --output ${query.simpleName}.wig2
    else
        echo "${query.simpleName}" > empty_sample.txt
    fi
    """
}

/*
 * Transforms WIG files to bigWig
 * Input: Tuple of [STR] directory to save results to (merged or not-merged), [WIG] forward cross-link sites, [WIG] reverse cross-link sites, [TXT] chroms names and sizes
 * Output:  bigWig -> Tuple of [STR] output directory, [BW] cross-link-sites
 */
process wig_to_bigWig_peak_called{
	tag {query.simpleName}
	publishDir "${params.output}/${out_dir}/bigWig", mode: 'copy', pattern: "*.bw"

	input:
	tuple val(out_dir), path(query), path(chrom_sizes)

	output:
	path("*.bw"), optional: true
	tuple val(out_dir), path("*.bw"), emit: bigWig, optional: true

	"""
	if [[ \$(cat ${query} | wc -l) > 1 ]]; then
		wigToBigWig ${query} ${chrom_sizes} ${query.baseName}.bw
	fi
	"""
}

/*
 * Transforms bigWig files into the bedGraph format
 * Input: Tuple of [STR] output directory and [BW] cross-link sites
 * Output: [BEDGRAPH] Cross-link sites  
 */
process bigWig_to_bedgraph_peak_called{
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

/*
 * Merges several WIG2 files into a single representative form 
 * Input: Tuple of [STR] Experiment name, [WIG2] several cross-link site files 
 * Output:  wig2_merged                 -> [WIG2] Merged cross-link sites
 *          wig_merged_cross_link_sites -> Tuple of [STR] saying "cross-link-sites-merged", [WIG] merged forward peaks, [WIG] merged reverse peaks
 */
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

/*
 * Splits WIG2 files into 2 WIG files - one for each strand
 * Input: [WIG2] Cross-link sites 
 * Output:  wig_split_forward       -> [WIG] forward cross-link sites
 *          wig_split_reverse       -> [WIG] reverse cross-link sites
 *          wig_split_both_strands  -> [WIG] cross-link sites for both strands
 */
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

/*
 * Calculates correlation for merging of WIG files. Both strands can get their correlation analysed separately or combined
 * Input: Tuple of [STR] name of experiment, [WIG] cross-link sites to be merged, [STR] current strand, [TXT] names of all chromosomes and their sizes
 * Output:  correlation_heatmap -> [PNG] Heatmap showing the correlation of merged WIG files
 *          correlation_matrix  -> [CSV] Matrix showing the correlation of merged WIG files
 */
process calc_wig_correlation{

    tag {name}
    publishDir "${params.output}/correlation_of_replicates", mode: 'copy', pattern: "${name}_${strand}_correlation.{png,csv}"

    input:
    tuple val(name),path(query),val(strand),path(chrom_sizes)

    output:
    path("${name}_${strand}_correlation.png"), emit: correlation_heatmap, optional: true
    path("${name}_${strand}_correlation.csv"), emit: correlation_matrix, optional: true

    script:
    String[] test_size = query
    """
    if [[ ${test_size.size()} > 1 ]]; then
        wig_file_statistics.R --input_path . --chrom_length ${chrom_sizes} --output ${name}_${strand} --type png
    fi
    """
}

/*
 * Transforms WIG files to bigWig
 * Input: Tuple of [STR] directory to save results to (merged or not-merged), [WIG] forward cross-link sites, [WIG] reverse cross-link sites, [TXT] chroms names and sizes
 * Output:  bigWig_both_strands -> Tuple of [STR] output directory, [BW] cross-link-sites for both strands
 *          bigWig_reverse      -> Tuple of [STR] output directory, [BW] reverse cross-link-sites 
 *          bigWig_forward      -> Tuple of [STR] output directory, [BW] forward cross-link-sites
 */
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

/*
 * Transforms bigWig files into the bedGraph format
 * Input: Tuple of [STR] output directory and [BW] cross-link sites
 * Output: [BEDGRAPH] Cross-link sites  
 */
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

/*
 * Sorts and indexes BAM files before peak calling
 * Input: [BAM] Alignment file
 * Output: Tuple of [BAM] Sorted alignment file and [BAI] index file
 */
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

/*
 * Splits WIG2 into 2 WIG files - one for each strand
 * Input: [WIG2] Cross-link sites
 * Output: Tuple of [STR] experiment name, [WIG] forward cross-link sites and [WIG] reverse cross-link sites
 */
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

/*
 * Generates a peak height histogram with a bar showing the percentile cutoff being used by several other analyses
 * Input: Tuple of [STR] experiment name, [WIG] forward cross-link sites and [WIG] reverse cross-link sites
 * Params:  params.color_barplot    -> [STR] Hexcode for color of histogram
 *          params.percentile       -> [FLOAT] Percentile used to show the peak height cutoff employed by other analyses
 * Output: [PNG] Histogram of peak heights  
 */
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