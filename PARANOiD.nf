#!/usr/bin/env nextflow

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
} from './modules/general_processes.nf'


//essential input files
input_reads     = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
reference       = Channel.fromPath( params.reference )		//FASTA file containing reference sequence(s)
barcode_file    = Channel.fromPath( params.barcodes )		//TSV file containing experiment names and the corresponding experiemental barcode sequence
//non essential input files
if(params.annotation != 'NO_FILE'){
    annotation_file = Channel.fromPath( params.annotation )
}
annotation = file(params.annotation)

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

    emit:
        //data for multiqc
        multiqc_quality_control                     = quality_control.out
        multiqc_quality_control_post_preprocessing  = quality_control_2.out
        multiqc_adapter_removal                     = adapter_removal.out.report_trimming
        multiqc_quality_filter                      = quality_filter.out.report_quality_filter

        // data for downstream processes
        fastq_reads_quality_filtered                = quality_filter.out.fastq_quality_filtered
}

workflow barcode_handling {
    take: 
        fastq_reads_quality_filtered
        barcode_file
    main:
        extract_rnd_barcode(fastq_reads_quality_filtered)
        check_barcode_file(barcode_file)
        split_exp_barcode(extract_rnd_barcode.out.fastq_rnd_barcode_extracted
                            .combine(check_barcode_file.out))
        get_length_exp_barcode(params.barcode_pattern)
        remove_exp_barcode(split_exp_barcode.out.fastq_split_experimental_barcode
                            .flatten()
                            .filter{ it.size() > 0 }
                            .combine(get_length_exp_barcode.out))
        remove_exp_barcode.out
            .map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/,file)}
            .groupTuple()
            .set{fastq_collect_preprocessed_to_merge}
        merge_preprocessed_reads(fastq_collect_preprocessed_to_merge)
        generate_barcode_barplot(split_exp_barcode.out.report_split_experimental_barcode.first()) //TODO: Instead of getting the first emitted file get the input from a specific dir and use all inputs

    emit:
        // reports
        multiqc_rnd_barcode_extraction  = extract_rnd_barcode.out.report_rnd_barcode_extraction
        multiqc_exp_barcode_splitting   = split_exp_barcode.out.report_split_experimental_barcode

        // data for downstream processes
        fastq_preprocessed_reads        = merge_preprocessed_reads.out
}

workflow alignment {
    take:
        reference
        fastq_preprocessed_reads
        annotation
    main:
        if( params.domain == 'pro' || params.map_to_transcripts == true ){
            build_index_bowtie(reference)
            mapping_bowtie(build_index_bowtie.out.first(),
                fastq_preprocessed_reads)
            alignments = mapping_bowtie.out.bam_alignments
            report_alignments = mapping_bowtie.out.report_alignments
        } else if ( params.domain == 'eu' ) {
            build_index_STAR(reference,
                annotation)
            mapping_STAR(fastq_preprocessed_reads
                            .combine(build_index_STAR.out))
            alignments = mapping_STAR.out.bam_alignments
            report_alignments = mapping_STAR.out.report_alignments
        }
        filter_empty_bams(alignments)
        collect_experiments_without_alignments(filter_empty_bams.out.report_empty_alignments.flatten().toList())

    emit:
        // reports
        report_empty_alignments         = collect_experiments_without_alignments.out
        report_alignments               = report_alignments

        // data for downstream processes
        bam_filtered_empty_alignments   = filter_empty_bams.out.bam_filtered_empty
}

workflow deduplication {
    take:
        bam_input
    main:
        sort_bam(bam_input.flatten())
        deduplicate(sort_bam.out)
        deduplicate.out.bam_deduplicated
            .map{file -> tuple(file.name - ~/\.[\w.]+.bam$/,file)}
            .groupTuple()
            .set{ bam_dedup_sort_to_merge }
        merge_deduplicated_bam(bam_dedup_sort_to_merge)

    emit:
        // reports
        report_deduplication    = deduplicate.out.report_deduplicated

        // data for downstream processes
        deduplicated_alignments = merge_deduplicated_bam.out
}

workflow transcript_analysis {
    take:
        bam_deduplicated
        reference
    main:
        count_hits(bam_deduplicated)
        get_top_hits(count_hits.out.flatten().toList())
        index_alignments(bam_deduplicated)
        extract_top_alignments(index_alignments.out
            .combine(get_top_hits.out))
        remove_newlines(reference)
        extract_top_transcript_sequences(get_top_hits.out
            .combine(remove_newlines.out))
    emit:
        bam_extracted_transcript_alignments = extract_top_alignments.out
        fasta_top_transcripts               = extract_top_transcript_sequences.out
}

workflow strand_preference {
    take:
        bam_collected
        reference
    main:
        bam_sorted = sort_bam_before_strand_pref(bam_collected)
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
        visualize_strand_preference(determine_strand_preference.out)
}

workflow peak_generation {
    take:
        bam_collected
        reference

    main:
        get_chromosome_sizes(reference)
        calculate_crosslink_sites(bam_collected
            .combine(get_chromosome_sizes.out))
        wig2_cross_link_sites = calculate_crosslink_sites.out.wig2_cross_link_sites
        if(params.merge_replicates == true){
            wig2_cross_link_sites
                .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*.wig2$/,file)} 
                .groupTuple()
                .set{wig2_grouped_samples}
            merge_wigs(wig2_grouped_samples)
            if(params.correlation_analysis == true){
                split_wig2_for_correlation(wig2_cross_link_sites)
                if(params.combine_strands_correlation == true){
                    value_both_strands = Channel.from("both_strands")
                    split_wig2_for_correlation.out.wig_split_both_strands
                        .flatten()
                        .map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*_(forward|reverse)?.wig$/,file)}
                        .groupTuple()
                        .combine(value_both_strands)
                        .combine(get_chromosome_sizes.out)
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
                        .combine(get_chromosome_sizes.out)
                        .set{prepared_input_correlation}
                }
                calc_wig_correlation(prepared_input_correlation)
            }
            merge_wigs.out.wig2_merged
                .set{wig2_cross_link_sites_collected}
            merge_wigs.out.wig_merged_cross_link_sites
                .concat(calculate_crosslink_sites.out.wig_cross_link_sites_split)
                .set{wig_cross_link_sites_to_bigWig}
        } else {
            wig2_cross_link_sites
                .set{wig2_cross_link_sites_collected}
            calculate_crosslink_sites.out.wig_cross_link_sites_split
                .set{wig_cross_link_sites_to_bigWig}
        }
        wig_to_bigWig(wig_cross_link_sites_to_bigWig
            .combine(get_chromosome_sizes.out))
        bigWig_to_bedgraph(wig_to_bigWig.out.bigWig_forward
            .mix(wig_to_bigWig.out.bigWig_reverse))
        if(params.omit_peak_calling == false){
            index_for_peak_calling(bam_collected)
            pureCLIP(index_for_peak_calling.out
                .combine(reference))
            report_pureCLIP     = pureCLIP.out.report_pureCLIP
            bed_pureCLIP_peaks  = pureCLIP.out.bed_crosslink_sites 
        } else { 
            report_pureCLIP     = Channel.empty() 
            bed_pureCLIP_peaks  = Channel.empty()
        }
        split_wig_2_for_peak_height_hist(wig2_cross_link_sites_collected)
        generate_peak_height_histogram(split_wig_2_for_peak_height_hist.out)
    emit:
        // reports
        report_pureCLIP     = report_pureCLIP

        // data for downstream analysis
        wig2_collected      = wig2_cross_link_sites_collected
        bigWig_collected    = wig_to_bigWig.out.bigWig_both_strands
        bed_pureCLIP_peaks  = bed_pureCLIP_peaks
}

workflow rna_subtype_analysis {
    take:
        wig2_cross_link_sites
        rna_subtypes
        annotation
    main:
        wig_to_bam(wig2_cross_link_sites)
        feature_counts(wig_to_bam.out
            .combine(rna_subtypes)
            .combine(annotation))
        feature_counts.out
            .map{file -> tuple(file.name - ~/\.[\w.]+.tsv$/,file)}
            .groupTuple()
            .set{tsv_sort_to_calculate_distribution}
        get_RNA_subtypes_distribution(rna_subtypes.collect(),
                                    tsv_sort_to_calculate_distribution)
        generate_RNA_subtypes_barplot(get_RNA_subtypes_distribution.out.tsv_subtype_distribution)
        collect_subtype_analysis_errors(get_RNA_subtypes_distribution.out.report_errors)
}

workflow motif_analysis {
    take:
        peaks
        reference
    main:
        sequence_extraction(peaks
            .combine(reference))
        if((2*params.seq_len)+1 >= params.min_motif_width){
            motif_search(sequence_extraction.out.extracted_sequences)
        }
}

workflow peak_distance_analysis {
    take:
        wig_peaks
    main:
        calculate_peak_distance(wig_peaks)
        plot_peak_distance(calculate_peak_distance.out)
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
            bam_split_alignments = split_bam_by_chromosome(alignment.out.bam_filtered_empty_alignments)
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

        output_reference(reference)
        command_line = Channel.from(workflow.commandLine)
        collect_workflow_metrics(command_line)
        multiqc(preprocessing.out.multiqc_adapter_removal.flatten().toList(),
                preprocessing.out.multiqc_quality_filter.flatten().toList(),
                preprocessing.out.multiqc_quality_control,
                preprocessing.out.multiqc_quality_control_post_preprocessing,
                barcode_handling.out.multiqc_exp_barcode_splitting.flatten().toList(),
                alignment.out.report_alignments.flatten().toList(),
                deduplication.out.report_deduplication.flatten().toList())
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