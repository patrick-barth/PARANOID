process count_hits {
    tag {bam.simpleName}

    publishDir "${params.output}/transcripts/overview_hits", mode: 'copy', pattern: "${bam.baseName}.hits.tsv"

    input:
    path(bam)

    output:
    path("${bam.baseName}.hits.tsv"),       emit: hit_counts
    path("${task.process}.version.txt"), 	emit: version

    """
    samtools view ${bam} \
        | cut -f3 \
        | sort \
        | uniq -c \
        | sort -nr \
        > ${bam.baseName}.hits.tsv

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
    """
}

process get_top_hits {

    publishDir "${params.output}/transcripts/overview_hits", mode: 'copy', pattern: "transcript-targets-top${params.number_top_transcripts}.txt"

    input:
    path(tsv)

    output:
    path("transcript-targets-top${params.number_top_transcripts}.txt"), emit: top_hits

    """
    head -${params.number_top_transcripts} -q *.tsv  \
        | rev \
        | cut -f1 -d' ' \
        | rev \
        | sort \
        | uniq \
        > transcript-targets-top${params.number_top_transcripts}.txt
    """
}

process index_alignments {
    tag {bam.simpleName}

    input:
    path(bam)

    output:
    tuple path("${bam.baseName}.sorted.bam"), path("${bam.baseName}.sorted.bam.bai"), emit: alignment_with_index
    path("${task.process}.version.txt"), 	emit: version

    """
    samtools sort ${bam} -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
    """
}

process extract_top_alignments {
    tag {bam.simpleName}

    input:
    tuple path(bam), path(bai), path(txt_alignments)

    output:
    path("${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam"), emit: top_alignments
    path("${task.process}.version.txt"), 	emit: version

    """
    samtools view -hb ${bam} `cat ${txt_alignments}` > ${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam

    echo -e "${task.process}\tsamtools\t\$(samtools --version | head -1 | rev | cut -f1 -d' ' | rev)" >> ${task.process}.version.txt
    """
}

process remove_newlines {

    input:
    path(ref)

    output:
    path("${ref.baseName}.removed_newlines.fna"),   emit: fasta_no_newlines

    """
    awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
    """
}

process extract_top_transcript_sequences {
    tag {txt_sequences.simpleName}
    publishDir "${params.output}/transcripts", mode: 'copy', pattern: "${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna"

    input:
    tuple path(txt_sequences), path(ref)

    output:
    path("${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna"), emit: fasta_top_sequences

    """
    egrep \
        -A1 \
        --no-group-separator \
        -f ${txt_sequences} \
        ${ref} \
        > ${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna
    """
}