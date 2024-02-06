/*
 * Extracts sequences with a specific length around detected peaks
 */
process sequence_extraction {
	tag {query.simpleName}

	publishDir "${params.output}/extracted_sequences", mode: 'copy', pattern: "*.extracted-sequences.*"

	input:
	tuple path(query),path(reference), val(percentile)

	output:
	path("*.extracted-sequences.fasta"), emit: extracted_sequences, optional: true
	path("*.extracted-sequences.*"), emit: complete_output, optional: true
	
	script:
	def local_remove_overlaps 	= params.remove_overlaps ? '--remove_overlaps' : ''
	def local_omit_cl			= params.omit_cl_nucleotide ? '--omit_cl' : ''
	def local_omit_width		= params.omit_cl_nucleotide ? '--omit_width ' + params.omit_cl_width : ''

	"""
	wig2-to-wig.py \
		--input ${query} \
		--output ${query.baseName}
	
	sequence-extraction.py \
		--input ${query.baseName}*.wig \
		--reference ${reference} \
		--output ${query.baseName}.extracted-sequences.fasta \
		--length ${params.seq_len} \
		--percentile ${percentile} \
		${local_omit_cl} \
		${local_omit_width} \
		${local_remove_overlaps}
	"""
}

/*
 * Extracts motifs from provided sequences
 */
process motif_search {
	tag {fasta.simpleName}
	publishDir "${params.output}/motif_search/", mode: 'copy'

	input:
	path(fasta)

	output:
	path("${fasta.baseName}_motif"), optional: true, emit: motifs

	script:
	"""
	if [[ \$(wc -l ${fasta} | cut -f1 -d' ') -ge 4 ]]; then
		streme \
			--oc ${fasta.baseName}_motif \
			--p ${fasta} \
			--dna \
			--seed 0 \
			--nmotifs ${params.max_motif_num} \
			--minw ${params.min_motif_width} \
			--maxw ${params.max_motif_width}
	fi
	"""
}