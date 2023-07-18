/*
 * Calculates distances between peaks
 * Params: params.distance -> Maximum distance between peaks allowed
 * Input: Tuple of [WIG2] Peak files and [FLOAT] percentile being used
 * Output: [TSV] Tab separated files showing the amounts of found distances
 */
process calculate_peak_distance {
	tag {query.simpleName}

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

	input:
	tuple path(query), val(percentile)

	output:
	path("${query.baseName}.peak-distance.tsv")

	"""
	wig2-to-wig.py --input ${query} --output ${query.baseName}
	peak-distance.py --input ${query.baseName}_*.wig --output ${query.baseName}.peak-distance.tsv --percentile ${percentile} --distance ${params.distance}
	"""
}

/*
 * Plots determined peak distances
 * Input: [TSV] Tab separated files showing the amounts of found distances
 * Output: [PNG] 2 figures showing the amount of found distances as a plot
 */
process plot_peak_distance {

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.simpleName}.peak-distance{_full,}.png"

	input:
	path(query)

	output:
	path("${query.simpleName}.peak-distance{_full,}.png")

	"""
	plot-distances.R --input ${query} --output ${query.simpleName}.peak-distance --type png
	"""
}