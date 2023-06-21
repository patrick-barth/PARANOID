/*
 * Counts and orders the amount of alignment each sequence in the reference got assigned
 * Input: [BAM] Allignment file 
 * Output: [TSV] File showing the amount of alignment each sequence in the reference git assigned to  
 */
process count_hits {
    tag {bam.simpleName}

    publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

    input:
    path(bam)

    output:
    path("${bam.baseName}.hits.tsv")

    """
    samtools view ${bam} | cut -f3 | sort | uniq -c | sort -nr > ${bam.baseName}.hits.tsv
    """
}

/*
 * Extracts the top transcripts
 * Input: [TSV] Transcript names and the amount of reads aligned to each
 * Params: [INT] Amount of top transcripts extracted
 * Output: [TXT] Names of transcripts with most alignments
 */
process get_top_hits {

    publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

    input:
    path(tsv)

    output:
    path("transcript-targets-top${params.number_top_transcripts}.txt")

    """
    head -${params.number_top_transcripts} -q *.tsv  | rev | cut -f1 -d' ' | rev | sort | uniq > transcript-targets-top${params.number_top_transcripts}.txt
    """
}

/*
 * Sorts and indexes BAM files
 * Input: [BAM] Alignment file
 * Output: Tuple of [BAM] Sorted alignment file and [BAI] index file
 */
process index_alignments {
    tag {bam.simpleName}

    input:
    path(bam)

    output:
    tuple path("${bam.baseName}.sorted.bam"), path("${bam.baseName}.sorted.bam.bai")

    """
    samtools sort ${bam} -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    """
}

/*
 * Extractes alignment entries of transcripts with most alignments from BAM file
 * Input: Tuple of [BAM] sorted alignment file, [BAI] index file, [TXT] names of alignments to be extracted
 * Params: [INT] Amount of top transcripts extracted
 * Output: [BAM] Filtered alignments
 */
process extract_top_alignments {
    tag {bam.simpleName}

    input:
    tuple path(bam), path(bai), path(txt_alignments)

    output:
    path("${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam")

    """
    samtools view -hb ${bam} `cat ${txt_alignments}` > ${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam
    """
}

/*
 * Removes newline from within the sequences of the reference file
 * Input: [FASTA] Reference file 
 * Output: [FASTA] Reference file without newline in the sequences 
 */
process remove_newlines {

    input:
    path(ref)

    output:
    path("${ref.baseName}.removed_newlines.fna")

    """
    awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
    """
}

/*
 * Extarctes the transcripts with most alignments from the reference
 * Input: Tuple of [TXT] names of top transcripts, [FASTA] reference file without newlines in the sequences 
 * Params: [INT] Amount of top transcripts extracted
 * Output: [FASTA] Transcripts with most alignments 
 */
process extract_top_transcript_sequences {
    tag {txt_sequences.simpleName}
    publishDir "${params.output}/transcripts", mode: 'copy'

    input:
    tuple path(txt_sequences), path(ref)

    output:
    path("${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna")

    """
    egrep -A1 --no-group-separator -f ${txt_sequences} ${ref} > ${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna
    """
}