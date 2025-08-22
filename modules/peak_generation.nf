/*
 * Retrieves sizes of all chromosomes in reference file 
 */
process get_chromosome_sizes{
	input:
	path(ref)

	output:
	path("${ref.simpleName}.chromosome_sizes.txt"), emit: chrom_sizes
    path("${task.process}.version.txt"), 	        emit: version

	"""
	samtools faidx ${ref}

	cut \
        -f1,2 \
        ${ref}.fai \
        > ${ref.simpleName}.chromosome_sizes.txt

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}

/*
 * Calculates cross-link sites from BAM files and returns them as 2 WIG files - one for the forward and one for the reverse strand.
 * Forward cross-link sites are determined by StartPosition -1
 * Reverse reads are dertermined by StartPosition + CIGAR string
 *  Additionally, alignments are filtered by the mapq score stated via params.mapq
 */
process calculate_crosslink_sites{
	tag {query.simpleName}
	publishDir "${params.output}/cross_link_sites_raw/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

	input:
	tuple path(query), path(chrom_sizes)

	output:
	path("${query.simpleName}.wig2"), emit: wig2_cross_link_sites, optional: true
	tuple val("cross_link_sites_raw"), path("${query.simpleName}_forward.wig"), emit: wig_cross_link_sites_split_forward, optional: true
    tuple val("cross_link_sites_raw"), path("${query.simpleName}_reverse.wig"), emit: wig_cross_link_sites_split_reverse, optional: true
	path("${query.simpleName}_{forward,reverse}.wig"), emit: wig_cross_link_sites, optional: true
    path("${task.process}.version.txt"), 	emit: version

	"""
	create-wig-from-bam.py \
        --input ${query} \
        --mapq ${params.mapq} \
        --chrom_sizes ${chrom_sizes} \
        --output ${query.simpleName}.wig2

	if [[ -f "${query.simpleName}.wig2" ]]; then
		wig2-to-wig.py \
            --input ${query.simpleName}.wig2 \
            --output ${query.simpleName}
	fi

    echo -e "${task.process}\tcreate-wig-from-bam.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Prepares reference file for peak calling by changing all non ACGTN nucleotides to Ns
 */
process prepare_ref_peak_calling {
    input:
    path(ref)

    output:
    path("${ref.baseName}.change_nuc.fa"), emit: fasta_reference_adapted_for_pureCLIP

    """
    awk '!/^>/ {
        gsub(/[^ACGTNacgtn]/, \"N\");
        print;
        next;
    } 1' ${ref} > ${ref.baseName}.change_nuc.fa
    """
}

/*
 * Performs peak calling via PureCLIP. Can do single peaks or peak regions.
 */
process pureCLIP {
    tag{bam.simpleName}
    cache false
    errorStrategy 'ignore' //TODO: is supposed to be only temporal. Need to find a solution for: ERROR: Emission probability became 0.0! This might be due to artifacts or outliers.

    publishDir "${params.output}/peak_calling", mode: 'copy', pattern: "${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed"
    publishDir "${params.output}/statistics/peak_calling/params", mode: 'copy', pattern: "${bam.simpleName}.pureCLIP_crosslink_sites.params"

    input:
    tuple path(bam), path(bai), path(ref)

    output:
    path("${bam.simpleName}.pureCLIP_crosslink_regions.bed"), optional: true,   emit: bed_crosslink_regions
    path("${bam.simpleName}.pureCLIP_crosslink_sites.bed"),                     emit: bed_crosslink_sites
    path("${bam.simpleName}.pureCLIP_crosslink_sites.params"),                  emit: report_pureCLIP
    path("${task.process}.version.txt"), 	                                    emit: version

    script:
    def local_regions       = params.peak_calling_regions ? '-or ' + ${bam.simpleName} + '.pureCLIP_crosslink_regions.bed' : ''
    def local_region_width  = params.peak_calling_regions ? '-dm ' + ${params.peak_calling_regions_width} : ''
    def local_high_coverage = params.peak_calling_for_high_coverage ? '-mtc 5000 -mtc2 5000 -ld' : ''
    def reference           = ref.endsWith('.fa') || ref.endsWith('.fasta') ? ref : ref.toString()[0..<ref.toString().lastIndexOf('.')] + '.fa'

    """
    pureclip \
        -i ${bam} \
        -bai ${bai} \
        -g ${reference} \
        -nt ${task.cpus} \
        -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed \
        ${local_regions} \
        ${local_region_width} \
        ${local_high_coverage}

    echo -e "${task.process}\tpureCLIP\t\$(pureclip --version | head -1 | cut -d' ' -f 3)" > ${task.process}.version.txt
    echo -e "${task.process}\tSeqAn\t\$(pureclip --version | tail -1 | cut -d' ' -f 3)" >> ${task.process}.version.txt
    """
}

/*
 * Transforms output from pureCLIP to WIG format
 */

process pureCLIP_to_wig{
    tag {query.simpleName}
    publishDir "${params.output}/cross_link_sites_peak_called/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

    input:
    tuple path(query), path(chrom_sizes)

    output:
    path("${query.simpleName}.wig2"), emit: wig2_peak_called_cl_sites, optional: true
    tuple val("cross_link_sites_peak_called"), path("${query.simpleName}_forward.wig"), emit: wig_peak_called_cl_sites_forward, optional: true
    tuple val("cross_link_sites_peak_called"), path("${query.simpleName}_reverse.wig"), emit: wig_peak_called_cl_sites_reverse, optional: true
    path('empty_sample,txt'), emit: txt_empty_sample, optional: true

    """
    if [[ "${query.size()}" > 0 ]]; then
        # Filter out of bounds coordinates, then split into strands
        awk -v sizes_file="${chrom_sizes}" '
        BEGIN {
            FS= "\\t";
            while((getline < sizes_file) > 0){
                chrom_sizes[\$1] = \$2;
            }
            close(sizes_file);
        }
        {
            chrom = \$1;
            pos = \$3;
            if(pos > 0 && <= chrom_sizes[chrom]){
                if (\$6 == "+") {
                    print \$0 > "${query.simpleName}_forward.bed";
                } else {
                    print \$0 > "${query.simpleName}_reverse.bed";
                }
            }
        }' ${query}
        # Add 1 into the 4th column (since it\'s from peak calling the actual size of the peak does not matter anymore)
        #  Then extracts important columns and converts the data into a valid WIG file
        if [[ -e "${query.simpleName}_forward.bed" ]]; then
            awk ' BEGIN {
                FS = "\\t";
                OFS = "\\t";
            }
            {
                \$NF = \$NF "\\t1";
                print \$0;
            }' ${query.simpleName}_forward.bed > ${query.simpleName}_forward.1.bed
            awk ' BEGIN {
                FS = "\\t";
            }
            {
                if (prev_header != \$1) {
                    prev_header = \$1;
                    printf("variableStep chrom=%s span=1\\n", \$1);
                }
                printf("%s %s\\n", \$3, \$8);
            }' ${query.simpleName}_forward.1.bed > ${query.simpleName}_forward.wig
        fi
        if [[ -e "${query.simpleName}_reverse.bed" ]]; then
            awk ' BEGIN {
                FS = "\\t";
                OFS = "\\t";
            }
            {
                \$NF = \$NF "\\t-1";
                print \$0;
            }' ${query.simpleName}_reverse.bed > ${query.simpleName}_reverse.1.bed
            awk 'BEGIN {
                FS = "\\t";
            }
            {
                if (prev_header != \$1) {
                    prev_header = \$1; 
                    printf("variableStep chrom=%s span=1\\n", \$1);
                }
                printf("%s %s\\n", \$3, \$8);
            }' ${query.simpleName}_reverse.1.bed > ${query.simpleName}_reverse.wig
        fi
        wig-to-wig2.py --wig ${query.simpleName}_forward.wig ${query.simpleName}_reverse.wig --output ${query.simpleName}.wig2
    else
        echo "${query.simpleName}" > empty_sample.txt
    fi
    """
}

/*
 * Merges several WIG2 files into a single representative form 
 */
process merge_wigs{
    tag {name}
    publishDir "${params.output}/cross_link_sites_merged/wig", mode: 'copy', pattern: "${name}_{forward,reverse}.wig"
    input:
    tuple val(name), path(query)

    output:
    path("${name}.wig2"),                   emit: wig2_merged
    tuple val("cross_link_sites_merged"), path("${name}_forward.wig"), path("${name}_reverse.wig"), emit: wig_merged_cross_link_sites
    tuple val("cross_link_sites_merged"), path("${name}_forward.wig"), emit: wig_merged_cross_link_sites_forward
    tuple val("cross_link_sites_merged"), path("${name}_reverse.wig"), emit: wig_merged_cross_link_sites_reverse
    path("${name}_{forward,reverse}.wig")
    path("${task.process}.version.txt"), 	emit: version

    script:
    def local_input     = name !=~ "unmatched*" ? name + '.wig2' : 'unmatched.wig2'
    def local_output    = name !=~ "unmatched*" ? name : 'unmatched'

    """
    merge-wig.py \
        --wig ${query} \
        --output ${local_input}

    wig2-to-wig.py \
        --input ${local_input} \
        --output ${local_output}

    echo -e "${task.process}\tmerge-wig.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
    """
}

/*
 * Splits WIG2 files into 2 WIG files - one for each strand
 */
process split_wig2_for_correlation{
    tag {query.simpleName}

    input:
    path(query)

    output:
    path("${query.simpleName}_forward.wig"),            emit: wig_split_forward, optional: true
    path("${query.simpleName}_reverse.wig"),            emit: wig_split_reverse, optional: true
    path("${query.simpleName}_{forward,reverse}.wig"),  emit: wig_split_both_strands, optional: true
    path("${task.process}.version.txt"), 	            emit: version

    """
    wig2-to-wig.py \
        --input ${query} \
        --output ${query.simpleName}

    echo -e "${task.process}\tmerge-wig.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
    """
}

/*
 * Calculates correlation for merging of WIG files. Both strands can get their correlation analysed separately or combined
 */
process calc_wig_correlation{

    tag {name}
    publishDir "${params.output}/correlation_of_replicates", mode: 'copy', pattern: "${name}_${strand}_correlation.{png,csv}"

    input:
    tuple val(name),path(query),val(strand),path(chrom_sizes)

    output:
    path("${name}_${strand}_correlation.png"),  emit: correlation_heatmap, optional: true
    path("${name}_${strand}_correlation.csv"),  emit: correlation_matrix, optional: true
    path("${task.process}.version.txt"), 	    emit: version

    script:
    String[] test_size = query

    def local_min_size          = params.combine_strands_correlation ? 2 : 1
    def local_combine_strands   = params.combine_strands_correlation ? '--both_strands' : ''

    """
    if [[ ${test_size.size()} > ${local_min_size} ]]; then
        wig_file_statistics.R \
            --input_path . \
            --chrom_length ${chrom_sizes} \
            --output ${name}_${strand} \
            --type png \
            ${local_combine_strands}
    fi

    echo -e "${task.process}\twig_file_statistics.R\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tR\t\$(R --version | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
    """
}

/*
 * Transforms WIG files to bigWig
 */
process wig_to_bigWig{
	tag {query.simpleName}
	publishDir "${params.output}/${out_dir}/bigWig", mode: 'copy', pattern: "*.bw"

	input:
	tuple val(out_dir), path(query), path(chrom_sizes)

	output:
	path("*.bw"), optional: true
	tuple val(out_dir), path("${query.baseName}.bw"),   emit: bigWig, optional: true
    path("${task.process}.version.txt"), 	            emit: version

	"""
	if [[ \$(cat ${query} | wc -l) > 1 ]]; then
		wigToBigWig \
            -clip \
            ${query} \
            ${chrom_sizes} \
            ${query.baseName}.bw
	fi

    echo -e "${task.process}\twigToBigWig\t\$(wigToBigWig 2>&1 | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
	"""
}

/*
 * Transforms bigWig files into the bedGraph format 
 */
process bigWig_to_bedgraph{
	tag {bigWig.simpleName}
	publishDir "${params.output}/${out_dir}/bedgraph", mode: 'copy', pattern: "*.bedgraph"

	input:
	tuple val(out_dir), path(bigWig)

	output:
	path("*.bedgraph"),                     emit: bedgraph
    path("${task.process}.version.txt"), 	emit: version

	"""
	bigWigToBedGraph \
        ${bigWig} \
        ${bigWig.baseName}.bedgraph

    echo -e "${task.process}\twigToBedGraph\tno_version_available" >> ${task.process}.version.txt
	"""
}

/*
 * Sorts and indexes BAM files before peak calling
 */
process index_for_peak_calling {
    tag{query.simpleName}

    input:
    path(query)

    output:
    tuple path("${query.simpleName}.sorted.bam"), path("${query.simpleName}.sorted.bam.bai"), emit: index
    path("${task.process}.version.txt"), 	emit: version

    """
    samtools sort ${query} > ${query.simpleName}.sorted.bam
    samtools index ${query.simpleName}.sorted.bam

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
    """
}

/*
 * Splits WIG2 into 2 WIG files - one for each strand
 */
process split_wig_2_for_peak_height_hist {
	tag {query.simpleName}

	input:
	path(query)

	output:
	tuple val("${query.simpleName}"), path("${query.simpleName}_forward.wig"), emit: wig_split_forward, optional: true
    tuple val("${query.simpleName}"), path("${query.simpleName}_reverse.wig"), emit: wig_split_reverse, optional: true
    path("${task.process}.version.txt"), 	emit: version

	"""
	wig2-to-wig.py \
        --input ${query} \
        --output ${query.simpleName}

    if [ ! -f ${query.simpleName}_forward.wig ]; then
        touch ${query.simpleName}_forward.wig
    fi
    if [[ ! -f ${query.simpleName}_reverse.wig ]]; then
        touch ${query.simpleName}_reverse.wig
    fi

    echo -e "${task.process}\twig2-to-wig.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
	"""
}

/*
 * Generates a peak height histogram with a bar showing the percentile cutoff being used by several other analyses 
 */
process generate_peak_height_histogram {
	tag {query}
	publishDir "${params.output}/peak_height_distribution", mode: 'copy'

	input:
	tuple val(query), path(forward), val(percentile)

	output:
	path("${query}.png"),                   emit: png_to_output_dir
    path("${task.process}.version.txt"), 	emit: version

	"""
	generate-peak-height-histogram.R \
        --input . \
        --output ${query} \
        --type png \
        --color "${params.color_barplot}" \
        --percentile ${percentile}

    echo -e "${task.process}\tgenerate-peak-height-histogram.R\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tR\t\$(R --version | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
	"""
}
