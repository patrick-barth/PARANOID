/*
 * Calculates distances between peaks
 */
process calculate_peak_distance {
	tag {query.simpleName}

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

	input:
	tuple path(query), val(percentile)

	output:
	path("${query.baseName}.peak-distance.tsv"), emit: tsv_distances

	"""
	wig2-to-wig.py \
		--input ${query} \
		--output ${query.baseName}

	peak-distance.py \
		--input ${query.baseName}_*.wig \
		--output ${query.baseName}.peak-distance.tsv \
		--percentile ${percentile} \
		--distance ${params.distance}
	"""
}

/*
 * Plots determined peak distances
 */
process plot_peak_distance {

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.simpleName}.peak-distance{_full,}.png"

	input:
	path(query)

	output:
	path("${query.simpleName}.peak-distance{_full,}.png"), emit: png_to_output_dir

	"""
	plot-distances.R \
		--input ${query} \
		--output ${query.simpleName}.peak-distance \
		--type png
	"""
}