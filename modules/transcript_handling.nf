process count_hits {
    tag {bam.simpleName}

    publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

    input:
    path(bam)

    output:
    path("${bam.baseName}.hits.tsv")// into tsv_count_to_get_hits

    """
    samtools view ${bam} | cut -f3 | sort | uniq -c | sort -nr > ${bam.baseName}.hits.tsv
    """
}

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

process remove_newlines {

    input:
    path(ref)

    output:
    path("${ref.baseName}.removed_newlines.fna")

    """
    awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
    """
}

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