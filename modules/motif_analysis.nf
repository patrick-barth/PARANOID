/*
 * Extracts sequences with a specific length around detected peaks
 * Params: params.omit_cl_nucleotide -> if TRUE the nucleotide at the cl site is replaced with an N
 *		params.omit_peak_calling -> if TRUE no percentile cutoff is performed as PureCLIP should already filter out background noise
 *		params.seq_len -> Length of sequence being extracted from each site around the cross-link site
 * Input: Tuple of [WIG2] peak files, [FASTA] a reference genome and [FLOAT] percentile being used
 * Output: [FASTA] Sequences extracted around peaks
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
	if(params.omit_cl_nucleotide == true && params.omit_peak_calling == true)
		"""
		wig2-to-wig.py --input ${query} --output ${query.baseName}
		sequence-extraction.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${percentile} --omit_cl --omit_width ${params.omit_cl_width}
		"""
	else if(params.omit_cl_nucleotide == true && params.omit_peak_calling == false)
		"""
		sequence-extraction.py --input ${query} --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${percentile} --omit_cl --omit_width ${params.omit_cl_width}
		"""
	else if(params.omit_cl_nucleotide == false && params.omit_peak_calling == true)
		"""
		wig2-to-wig.py --input ${query} --output ${query.baseName}
		sequence-extraction.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${percentile}
		"""
	else if(params.omit_cl_nucleotide == false && params.omit_peak_calling == false)
		"""
		sequence-extraction.py --input ${query} --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${percentile}
		"""
}

/*
 * Extracts motifs from provided sequences
 * Params: params.max_motif_num -> Maximum number of motifs to be returned
 *		params.min_motif_width -> Minimum length extracted motifs are supposed to have
 *		params.max_motif_width -> Maximum length extracted motifs are supposed to have
 * Input: [FASTA] (Short) Sequences from which motifs are supposed to be extracted
 * Output: 
 */
process motif_search {
	tag {fasta.simpleName}
	publishDir "${params.output}/motif_search/", mode: 'copy'

	input:
	path(fasta)

	output:
	path("${fasta.baseName}_motif"), optional: true

	script:
	"""
	if [[ \$(wc -l ${fasta} | cut -f1 -d' ') -ge 4 ]]; then
		streme --oc ${fasta.baseName}_motif --p ${fasta} --dna --seed 0 --nmotifs ${params.max_motif_num} --minw ${params.min_motif_width} --maxw ${params.max_motif_width}
	fi
	"""
}