/*
 * Calculates distances between peaks
 * Params: params.percentile -> Percentile of peaks being omited in order to filter out background noise. Omits peaks according to their height
 *		params.distance -> Maximum distance between peaks allowed
 * Input: [WIG2] Peak files
 * Output: [TSV] Tab separated files showing the amounts of found distances
 */
process calculate_peak_distance {
	tag {query.simpleName}

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

	input:
	path(query)

	output:
	path("${query.baseName}.peak-distance.tsv")

	"""
	wig2-to-wig.py --input ${query} --output ${query.baseName}
	peak-distance.py --input ${query.baseName}_*.wig --output ${query.baseName}.peak-distance.tsv --percentile ${params.percentile} --distance ${params.distance}
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