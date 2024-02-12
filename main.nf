#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import processes for preprocessing
include{
    quality_control
    quality_control_2
    adapter_removal
    quality_filter
} from './modules/preprocessing.nf'

//import processes for barcode handling
include{
    extract_rnd_barcode
    check_barcode_file
    split_exp_barcode
    get_length_exp_barcode
    remove_exp_barcode
    merge_preprocessed_reads
    generate_barcode_barplot
} from './modules/barcode_handling.nf'

//import processes used for aligments
include{
    build_index_bowtie
    mapping_bowtie
    build_index_STAR
    mapping_STAR
    filter_empty_bams
    collect_experiments_without_alignments
} from './modules/alignment.nf'

//import processes used for deduplication
include{
    sort_bam
    deduplicate
    merge_deduplicated_bam
} from './modules/deduplication.nf'

//import processes used for transcript analysis
if(params.map_to_transcripts == true){
    include{
        count_hits
        get_top_hits
        index_alignments
        extract_top_alignments
        remove_newlines
        extract_top_transcript_sequences
    } from './modules/transcript_handling.nf'
}

//import processes used to generate peaks
include{
    get_chromosome_sizes
    pureCLIP_to_wig
    calculate_crosslink_sites
    split_wig2_for_correlation
    calc_wig_correlation
    merge_wigs
    wig_to_bigWig
    bigWig_to_bedgraph
    index_for_peak_calling
    pureCLIP
    split_wig_2_for_peak_height_hist
    generate_peak_height_histogram
} from './modules/peak_generation.nf'

//import processes used to determine strand preference
include{
    sort_bam_before_strand_pref
    determine_strand_preference
    visualize_strand_preference
} from './modules/strand_preference.nf'

//import processes used for RNA subtype analysis
if(params.annotation != 'NO_FILE'){
    include{
        wig_to_bam
        feature_counts
        get_RNA_subtypes_distribution
        generate_RNA_subtypes_barplot
        collect_subtype_analysis_errors
    } from './modules/RNA_subtype_distribution.nf'
}

//import processes used for motif analysis
include{
    sequence_extraction
    motif_search
} from './modules/motif_analysis.nf'

//import processes used for peak distance analysis
include{
    calculate_peak_distance
    plot_peak_distance
} from './modules/peak_distance_analysis.nf'

//import processes used to generate the IGV session
include{
    prepare_annotation_for_igv
    generate_igv_session
} from './modules/generate_IGV_session.nf'

//import general processes
include{
    split_bam_by_chromosome
    sort_and_index_alignment
    output_reference
    multiqc
    collect_workflow_metrics
    get_md5sum
    collect_versions
} from './modules/general_processes.nf'

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        input reads         : ${params.reads}
        input reference     : ${params.reference}
        input barcodes      : ${params.barcodes}
        input annotation    : ${params.annotation}
        output to           : ${params.output}
        --
        Barcode pattern     : ${params.barcode_pattern}
        Barcode Mismatches  : ${params.barcode_mismatches}
        --
        Aligner                 : ${params.domain}
        MapQ-Score              : ${params.mapq}
        Max alignments          : ${params.max_alignments}
        Report all alignments   : ${params.report_all_alignments}
        --
        Merge replicates        : ${params.merge_replicates}
        Correlation analysis    : ${params.correlation_analysis}
        --
        RNA subtypes            : ${params.rna_subtypes}
        Gene identifier         : ${params.gene_id}
        --
        Omit peak calling       : ${params.omit_peak_calling}
        Peak calling high cov   : ${params.peak_calling_for_high_coverage}
        Peak calling regions    : ${params.peak_calling_regions}
        Region width            : ${params.peak_calling_regions_width}
        --
        Transcript analysis     : ${params.map_to_transcripts}
        Max top transcripts     : ${params.number_top_transcripts}
        --
        Omit sequence analysis          : ${params.omit_sequence_extraction}
        Extracted sequence length       : ${params.seq_len}
        Omit cl nucleotide              : ${params.omit_cl_nucleotide}
        Width around cl to be omitted   : ${params.omit_cl_width}
        Remove overlaps                 : ${params.remove_overlaps}
        Maximum number of motifs        : ${params.max_motif_num}
        Minimum motif length            : ${params.min_motif_width}
        Maximum motif length            : ${params.min_motif_width}
        --
        Omit peak distance analysis : ${params.omit_peak_distance}
        Distance                    : ${params.distance}
        --
        Color of barplots       : ${params.color_barplot}
        Percentiles             : ${params.percentile}
        --
        run as              : ${workflow.commandLine}
        started at          : ${workflow.start}
        config files        : ${workflow.configFiles}
        """
        .stripIndent()


//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
reference       = Channel.fromPath( params.reference )		//FASTA file containing reference sequence(s)
barcode_file    = Channel.fromPath( params.barcodes )		//TSV file containing experiment names and the corresponding experiemental barcode sequence
collect_input_files = input_reads
                    .concat(reference)
                    .concat(barcode_file)

//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
    collect_input_files = collect_input_files.concat(annotation_file)
}
annotation = file(params.annotation)

if(params.omit_peak_calling == false){
    percentile = Channel.from(0)
} else {
    percentile = Channel.from(params.percentile)
}

//preparation for RNA subtype analysis
rna_subtypes = Channel.from(params.rna_subtypes
                        .split(','))

workflow preprocessing {
    take: 
        input_reads
    main:
        if(params.speed) {
            input_reads.splitFastq(by: params.split_fastq_by, file:true ).set{input_reads}
        }
        //preprocessing
        quality_control(input_reads)
        adapter_removal(input_reads)
        quality_filter(adapter_removal.out.fastq_trimmed)
        quality_control_2(quality_filter.out.fastq_quality_filtered)

        // Collect versions
        versions = quality_control.out.version.first()
                    .concat(quality_control_2.out.version.first())
                    .concat(adapter_removal.out.version.first())
                    .concat(quality_filter.out.version.first())

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out.summary
        multiqc_quality_control_post_preprocessing  = quality_control_2.out.summary
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_quality_filter                      = quality_filter.out.report_quality_filter

        // data for downstream processes
        fastq_reads_quality_filtered                = quality_filter.out.fastq_quality_filtered

        versions    = versions
}

workflow barcode_handling {
    take: 
        fastq_reads_quality_filtered
        barcode_file
    main:
        extract_rnd_barcode(fastq_reads_quality_filtered)
        check_barcode_file(barcode_file)
        split_exp_barcode(extract_rnd_barcode.out.fastq_rnd_barcode_extracted
                            .combine(check_barcode_file.out.barcodes))
        get_length_exp_barcode(params.barcode_pattern)
        remove_exp_barcode(split_exp_barcode.out.fastq_split_experimental_barcode
                            .flatten()
                            .filter{ it.size() > 0 }
                            .combine(get_length_exp_barcode.out.var_length))
        remove_exp_barcode.out.fastq_trimmed
            .map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/,file)}
            .groupTuple()
            .set{fastq_collect_preprocessed_to_merge}
        merge_preprocessed_reads(fastq_collect_preprocessed_to_merge)
        generate_barcode_barplot(split_exp_barcode.out.report_split_experimental_barcode.first()) //TODO: Instead of getting the first emitted file get the input from a specific dir and use all inputs

        // Collect versions
        versions = extract_rnd_barcode.out.version.first()
                    .concat(check_barcode_file.out.version.first())
                    .concat(split_exp_barcode.out.version.first())
                    .concat(remove_exp_barcode.out.version.first())
                    .concat(generate_barcode_barplot.out.version.first())

    emit:
        // reports
        multiqc_rnd_barcode_extraction  = extract_rnd_barcode.out.report_rnd_barcode_extraction
        multiqc_exp_barcode_splitting   = split_exp_barcode.out.report_split_experimental_barcode
        versions = versions

        // data for downstream processes
        fastq_preprocessed_reads        = merge_preprocessed_reads.out.fastq_merged


}

workflow alignment {
    take:
        reference
        fastq_preprocessed_reads
        annotation
    main:
        if( params.domain == 'pro' || params.map_to_transcripts == true ){
            build_index_bowtie(reference)
            mapping_bowtie(build_index_bowtie.out.index.first(),
                fastq_preprocessed_reads)

            alignments          = mapping_bowtie.out.alignments
            report_alignments   = mapping_bowtie.out.report
            version_index       = build_index_bowtie.out.version.first()
            version_alignment   = mapping_bowtie.out.version.first()
        } else if ( params.domain == 'eu' ) {
            build_index_STAR(reference,
                annotation)
            mapping_STAR(fastq_preprocessed_reads
                            .combine(build_index_STAR.out.index))

            alignments          = mapping_STAR.out.alignments
            report_alignments   = mapping_STAR.out.report
            version_index       = build_index_STAR.out.version.first()
            version_alignment   = mapping_STAR.out.version.first()
        }
        filter_empty_bams(alignments)
        collect_experiments_without_alignments(filter_empty_bams.out.report_empty.flatten().toList())

        // Collect versions
        versions = version_index.concat(version_alignment)

    emit:
        // reports
        report_empty_alignments         = collect_experiments_without_alignments.out.output
        report_alignments               = report_alignments
        versions                        = versions

        // data for downstream processes
        bam_filtered_empty_alignments   = filter_empty_bams.out.bam_filtered
}

workflow deduplication {
    take:
        bam_input
    main:
        sort_bam(bam_input.flatten())
        deduplicate(sort_bam.out.index_and_alignment)
        deduplicate.out.bam_deduplicated
            .map{file -> tuple(file.name - ~/\.[\w.]+.bam$/,file)}
            .groupTuple()
            .set{ bam_dedup_sort_to_merge }
        merge_deduplicated_bam(bam_dedup_sort_to_merge)

        // Collect versions
        versions = sort_bam.out.version.first()
                    .concat(deduplicate.out.version.first())
                    .concat(merge_deduplicated_bam.out.version.first())

    emit:
        // reports
        report_deduplication    = deduplicate.out.report_deduplicated
        versions                = versions

        // data for downstream processes
        deduplicated_alignments = merge_deduplicated_bam.out.bam_merged
}

workflow transcript_analysis {
    take:
        bam_deduplicated
        reference
    main:
        count_hits(bam_deduplicated)
        get_top_hits(count_hits.out.hit_counts.flatten().toList())
        index_alignments(bam_deduplicated)
        extract_top_alignments(index_alignments.out.alignment_with_index
            .combine(get_top_hits.out.top_hits))
        remove_newlines(reference)
        extract_top_transcript_sequences(get_top_hits.out.top_hits
            .combine(remove_newlines.out.fasta_no_newlines))

        // Collect versions
        versions = count_hits.out.version.first()
                    .concat(index_alignments.out.version.first())
                    .concat(extract_top_alignments.out.version.first())
    emit:
        bam_extracted_transcript_alignments = extract_top_alignments.out.top_alignments
        fasta_top_transcripts               = extract_top_transcript_sequences.out.fasta_top_sequences
        versions                            = versions
}

workflow strand_preference {
    take:
        bam_collected
        reference
    main:
        sort_bam_before_strand_pref(bam_collected)
        bam_sorted = sort_bam_before_strand_pref.out.bam_sorted 
        if(params.merge_replicates == true){
            bam_sorted
                .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*(.sorted)?.bam$/,file)} 
                .groupTuple()
                .set{grouped_bam_to_strand_preference}
        } else {
            bam_sorted
                .map{file -> tuple(file.name - ~/(_filtered_top\d*)?(.sorted)?.bam$/,file)}
                .set{grouped_bam_to_strand_preference}
        }
        determine_strand_preference(grouped_bam_to_strand_preference
            .combine(reference))
        visualize_strand_preference(determine_strand_preference.out.txt_strand_proportions)

        // Collect versions
        versions = sort_bam_before_strand_pref.out.version.first()
                    .concat(determine_strand_preference.out.version.first())
                    .concat(visualize_strand_preference.out.version.first())

    emit:
        versions = versions
}

workflow peak_generation {
    take:
        bam_collected
        reference

    main:
        get_chromosome_sizes(reference)
        //calculate raw cross-link sites
        calculate_crosslink_sites(bam_collected
            .combine(get_chromosome_sizes.out.chrom_sizes))
        collect_cl_sites_to_transform = calculate_crosslink_sites.out.wig_cross_link_sites_split_forward
                                            .concat(calculate_crosslink_sites.out.wig_cross_link_sites_split_reverse)
        
        //Handle peak calling when it is not omitted. If pureCLIP is executed all further analyses are performed on the peaks determined by pureCLIP
        if(params.omit_peak_calling == false){
            index_for_peak_calling(bam_collected)
            pureCLIP(index_for_peak_calling.out.index
                .combine(reference))
            pureCLIP_to_wig(pureCLIP.out.bed_crosslink_sites)

            report_pureCLIP                 = pureCLIP.out.report_pureCLIP
            bed_pureCLIP_peaks              = pureCLIP.out.bed_crosslink_sites
            wig2_cross_link_sites           = pureCLIP_to_wig.out.wig2_peak_called_cl_sites
            collect_cl_sites_to_transform   = collect_cl_sites_to_transform
                                                .concat(pureCLIP_to_wig.out.wig_peak_called_cl_sites_forward)
                                                .concat(pureCLIP_to_wig.out.wig_peak_called_cl_sites_reverse)
        } else { 
            report_pureCLIP         = Channel.empty() 
            bed_pureCLIP_peaks      = Channel.empty()
            wig2_cross_link_sites   = calculate_crosslink_sites.out.wig2_cross_link_sites
        }
        
        if(params.merge_replicates == true){
            //group replicates according to their experiment and merge them
            wig2_cross_link_sites
                .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*.wig2$/,file)} 
                .groupTuple()
                .set{wig2_grouped_samples}
            merge_wigs(wig2_grouped_samples)

            if(params.correlation_analysis == true){
                split_wig2_for_correlation(wig2_cross_link_sites)
                // combines replicates for correlation analysis. if both strands are to be checked together they are combined
                //  if the strands are to be checked separately they are only the according strands are combined. 
                //  Code for grouping goes until the next comment
                if(params.combine_strands_correlation == true){
                    value_both_strands = Channel.from("both_strands")
                    split_wig2_for_correlation.out.wig_split_both_strands
                        .flatten()
                        .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*_(forward|reverse)?.wig$/,file)}
                        .groupTuple()
                        .combine(value_both_strands)
                        .combine(get_chromosome_sizes.out.chrom_sizes)
                        .set{prepared_input_correlation}
                } else {
                    value_forward = Channel.from("forward")
                    split_wig2_for_correlation.out.wig_split_forward
                        .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?_forward.wig$/,file)} 
                        .groupTuple()
                        .combine(value_forward)
                        .set{group_forward}

                    value_reverse = Channel.from("reverse")
                    split_wig2_for_correlation.out.wig_split_reverse
                        .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?_reverse.wig$/,file)} 
                        .groupTuple()
                        .combine(value_reverse)
                        .set{group_reverse}

                    group_forward
                        .concat(group_reverse)
                        .combine(get_chromosome_sizes.out.chrom_sizes)
                        .set{prepared_input_correlation}
                }
                // actual calculation of correlation
                calc_wig_correlation(prepared_input_correlation)
            }
            // Generates output channels
            merge_wigs.out.wig2_merged
                .set{wig2_cross_link_sites_collected}
            collect_cl_sites_to_transform = collect_cl_sites_to_transform
                .concat(merge_wigs.out.wig_merged_cross_link_sites_forward)
                .concat(merge_wigs.out.wig_merged_cross_link_sites_reverse)
        } else {
            wig2_cross_link_sites
                .set{wig2_cross_link_sites_collected}
        }
        //Transformation of peaks into different formats (bigWig and bedgraph)
        wig_to_bigWig(collect_cl_sites_to_transform
            .combine(get_chromosome_sizes.out.chrom_sizes))
        bigWig_to_bedgraph(wig_to_bigWig.out.bigWig)

        // Get correct bigWigs to display via IGV
        if(params.merge_replicates == true){
            wig_to_bigWig.out.bigWig
                .filter{ it[0] == 'cross-link-sites-merged' }
                .set{collect_bigWig_to_IGV}
        } else if(params.omit_peak_calling == false){
            wig_to_bigWig.out.bigWig
                .filter{ it[0] == 'cross-link-sites-peak-called' }
                .set{collect_bigWig_to_IGV}
        } else {
            wig_to_bigWig.out.bigWig
                .filter{ it[0] == 'cross-link-sites-raw' }
                .set{collect_bigWig_to_IGV}
        }

        //Generate peak height histogram
        split_wig_2_for_peak_height_hist(wig2_cross_link_sites_collected)
        split_wig_2_for_peak_height_hist.out.wig_split_forward
            .mix(split_wig_2_for_peak_height_hist.out.wig_split_reverse)
            .groupTuple()
            .set{combined_wig_for_peak_height_histogram}
        generate_peak_height_histogram(combined_wig_for_peak_height_histogram
            .combine(percentile))

        // Collect versions
        versions = get_chromosome_sizes.out.version.first()
                    .concat(calculate_crosslink_sites.out.version.first())
                    .concat(wig_to_bigWig.out.version.first())
                    .concat(bigWig_to_bedgraph.out.version.first())
                    .concat(split_wig_2_for_peak_height_hist.out.version.first())
                    .concat(generate_peak_height_histogram.out.version.first())

        versions = !params.omit_peak_calling ? versions.concat(index_for_peak_calling.out.version.first()).concat(pureCLIP.out.version.first()) : versions
        versions = params.merge_replicates ? versions.concat(merge_wigs.out.version.first()) : versions
        versions = params.merge_replicates && params.correlation_analysis ? versions.concat(split_wig2_for_correlation.out.version.first()).concat(calc_wig_correlation.out.version.first()) : versions

    emit:
        // reports
        report_pureCLIP     = report_pureCLIP
        versions            = versions

        // data for downstream analysis
        wig2_collected      = wig2_cross_link_sites_collected
        bigWig_collected    = collect_bigWig_to_IGV
        bed_pureCLIP_peaks  = bed_pureCLIP_peaks
}

workflow rna_subtype_analysis {
    take:
        wig2_cross_link_sites
        rna_subtypes
        annotation
    main:
        wig_to_bam(wig2_cross_link_sites)
        feature_counts(wig_to_bam.out.bam
            .combine(rna_subtypes)
            .combine(annotation))
        feature_counts.out.features
            .map{file -> tuple(file.name - ~/\.[\w.]+.tsv$/,file)}
            .groupTuple()
            .set{tsv_sort_to_calculate_distribution}
        get_RNA_subtypes_distribution(rna_subtypes.collect(),
                                    tsv_sort_to_calculate_distribution)
        generate_RNA_subtypes_barplot(get_RNA_subtypes_distribution.out.tsv_subtype_distribution)
        collect_subtype_analysis_errors(get_RNA_subtypes_distribution.out.report_errors)

        // Collect versions
        versions = wig_to_bam.out.version.first()
                    .concat(feature_counts.out.version.first())
                    .concat(get_RNA_subtypes_distribution.out.version.first())
                    .concat(generate_RNA_subtypes_barplot.out.version.first())

    emit:
        versions = versions
}

workflow motif_analysis {
    take:
        peaks
        reference
    main:
        sequence_extraction(peaks
            .combine(reference)
            .combine(percentile))
        if((2*params.seq_len)+1 >= params.min_motif_width){
            motif_search(sequence_extraction.out.extracted_sequences)
        }

        // Collect versions
        versions = sequence_extraction.out.version.first()
        versions = (2*params.seq_len)+1 >= params.min_motif_width ? versions.concat(motif_search.out.version.first()) : versions

    emit:
        versions = versions
}

workflow peak_distance_analysis {
    take:
        wig_peaks
    main:
        calculate_peak_distance(wig_peaks
            .combine(percentile))
        plot_peak_distance(calculate_peak_distance.out.tsv_distances)

        // Collect versions
        versions = calculate_peak_distance.out.version.first()
                    .concat(plot_peak_distance.out.version.first())

    emit:
        versions = versions
}
workflow igv_session {
    take:
        bigWig_tracks
        bam_alignments
        path_track_directory
        reference
        annotation
    main:
        if(params.annotation != 'NO_FILE'){
            prepare_annotation_for_igv(annotation)
            collect_annotation = prepare_annotation_for_igv.out.name_annotation
        } else {
            collect_annotation = annotation
        }
        generate_igv_session(bigWig_tracks.flatten().toList(),
            bam_alignments.flatten().toList(),
            path_track_directory,
            reference,
            collect_annotation)

        // Collect versions
        versions = Channel.empty()
        versions = params.annotation != 'NO_FILE' ? versions.concat(prepare_annotation_for_igv.out.version.first()) : versions
        versions = versions.concat(generate_igv_session.out.version.first())

    emit:
        versions = versions
}

params.version = false
if(params.version){
	println(workflow.manifest.version)
} else {
    workflow {
        preprocessing(input_reads)
        barcode_handling(preprocessing.out.fastq_reads_quality_filtered,barcode_file)
        alignment(reference,barcode_handling.out.fastq_preprocessed_reads,annotation)
        if(params.map_to_transcripts == false && params.speed){
            split_bam_by_chromosome(alignment.out.bam_filtered_empty_alignments)
            bam_split_alignments = split_bam_by_chromosome.out.bam_split 
        } else { bam_split_alignments = alignment.out.bam_filtered_empty_alignments }
        deduplication(bam_split_alignments)
        if(params.map_to_transcripts == true){
            transcript_analysis(deduplication.out.deduplicated_alignments,reference)
            bam_collected               = transcript_analysis.out.bam_extracted_transcript_alignments
            reference_for_downstream    = transcript_analysis.out.fasta_top_transcripts
        } else { 
            bam_collected               = deduplication.out.deduplicated_alignments
            reference_for_downstream    = reference
        }
        strand_preference(bam_collected,reference_for_downstream)
        peak_generation(bam_collected,reference)
        if(params.annotation != 'NO_FILE'){
            rna_subtype_analysis(peak_generation.out.wig2_collected, rna_subtypes, annotation_file)
        }
        if(params.omit_sequence_extraction == false){
            motif_analysis(peak_generation.out.wig2_collected,reference)
        }
        if(params.omit_peak_distance == false){
            peak_distance_analysis(peak_generation.out.wig2_collected)
        }
        sort_and_index_alignment(deduplication.out.deduplicated_alignments)
        if (params.merge_replicates == true){
            track_path_dir = Channel.value('cross-link-sites-merged/bigWig')
            peak_generation.out.bigWig_collected
                .filter{it[0] == 'cross-link-sites-merged'}
                .map{it[1]}.set{collected_bigWig}
        } else {
            track_path_dir = Channel.value('cross-link-sites/bigWig')
            peak_generation.out.bigWig_collected
                .map{it[1]}.set{collected_bigWig}
        }
        igv_session(collected_bigWig.flatten().toList(),
                    sort_and_index_alignment.out.bam_sorted_alignments.flatten().toList(),
                    track_path_dir,
                    reference,
                    annotation)

        multiqc(preprocessing.out.multiqc_adapter_removal.flatten().toList(),
                preprocessing.out.multiqc_quality_filter.flatten().toList(),
                preprocessing.out.multiqc_quality_control,
                preprocessing.out.multiqc_quality_control_post_preprocessing,
                barcode_handling.out.multiqc_exp_barcode_splitting.flatten().toList(),
                alignment.out.report_alignments.flatten().toList(),
                deduplication.out.report_deduplication.flatten().toList())

        output_reference(reference)
        collect_workflow_metrics()
        get_md5sum(collect_input_files.flatten().toList())

        versions = preprocessing.out.versions
            .concat(barcode_handling.out.versions)
            .concat(alignment.out.versions)
            .concat(deduplication.out.versions)
            .concat(strand_preference.out.versions)
            .concat(peak_generation.out.versions)
            .concat(igv_session.out.versions)

        versions = params.map_to_transcripts ? versions.concat(transcript_analysis.out.versions) : versions
        versions = params.annotation != 'NO_FILE' ? versions.concat(rna_subtype_analysis.out.versions) : versions
        versions = !params.omit_sequence_extraction ? versions.concat(motif_analysis.out.versions) : versions
        versions = !params.omit_peak_distance ? versions.concat(peak_distance_analysis.out.versions) : versions

        collect_versions(versions
                            .flatten()
                            .toList())
    }
}

workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}