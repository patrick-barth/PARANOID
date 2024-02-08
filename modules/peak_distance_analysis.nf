/*
 * Calculates distances between peaks
 */
process calculate_peak_distance {
	tag {query.simpleName}

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

	input:
	tuple path(query), val(percentile)

	output:
	path("${query.baseName}.peak-distance.tsv"), 	emit: tsv_distances
	path("${task.process}.version.txt"), 			emit: version

	"""
	wig2-to-wig.py \
		--input ${query} \
		--output ${query.baseName}

	peak-distance.py \
		--input ${query.baseName}_*.wig \
		--output ${query.baseName}.peak-distance.tsv \
		--percentile ${percentile} \
		--distance ${params.distance}

	echo -e "${task.process}\tpeak-distance.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" >> ${task.process}.version.txt
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
	path("${query.simpleName}.peak-distance{_full,}.png"), 	emit: png_to_output_dir
	path("${task.process}.version.txt"), 					emit: version

	"""
	plot-distances.R \
		--input ${query} \
		--output ${query.simpleName}.peak-distance \
		--type png

	echo -e "${task.process}\tplot-distances.R\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tR\t\$(R --version | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
	"""
}