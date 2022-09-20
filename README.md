# PARANOID
Pipeline for Automated Read ANalysis Of iCLIP Data

//TODO: Blabla ueber Background

## Overview
[Basic usage](#Basic-usage)  
[Inputs](#Inputs)  
[Parameters](#Parameters)  
[Outputs](#Outputs)  

## Basic-usage
nextflow PARANOiD.nf --reads \<reads.fastq\> --reference \<reference_sequence.fasta\> --barcodes \<barcodes.tsv\>

## Inputs

### Reads (essential)

#### Usage
```
--reads reads_file.fastq
```

### Reference (essential)

#### Usage
```
--reference reference_file.fastq
```

### Barcodes (essential)

Barcode sequences are used to assign reads to their experiment. The file is provided as TSV-file (tab separated value).
The first consists of the experiment name and the second of the nucleotide sequence representing the barcode of the experiment. 
One experiment is described per lane and the columns are divided by a tab
The experiment name should be named as follows:
	\<experiment_name\>\_rep_\<replicate-number\>

Example:
```
experiment1_rep_1	GCATTG  
experiment1_rep_2	CAGTAA  
experiment1_rep_3	GGCCTA  
experiment2_rep_1	AATCCG  
experiment2_rep_2	CCGTTA  
experiment2_rep_3	GTCATT  
```

#### Usage
```
--barcodes barcode_file.tsv
```

### Annotation

#### Usage
```
--annotation annotation_file.gff
```

## Parameters

### --barcode_pattern

A string that allows to adapt to other barcode patterns (default is iCLIP2). N represent the random barcodes and X represent the experimental barcode.  
#### Usage
Example for iCLIP:
```
--barcode_pattern NNXXXXNNN
```
Default:
```
--barcode_pattern NNNNNXXXXXXNNNN
```

### --domain

#### Usage
```
--domain eu
```
Default:
```
--domain pro
```

### --output

#### Usage
```
--output /path/to/output
```
Default:
```
--output ./output
```

### --min_length

#### Usage
```
--min_length 30
```

### --min_qual

#### Usage
```
--min_qual 20
```

### --min_percent_qual_filter

#### Usage
```
--min_percent_qual_filter 90
```

### --barcode_mismatches

#### Usage
```
--barcode_mismatches 1
```

### --mapq

#### Usage
```
--mapq 2
```

### --split_fastq_by

#### Usage
```
--split_fastq_by 1000000
```

## Optional analyses

### Transcript analysis

#### --map_to_transcripts

##### Usage
```
--map_to_transcripts
```

#### --number_top_transcripts

##### Usage
```
--number_top_transcripts 10
```

#### Usage
```
--mapq 2
```

### Peak calling

#### --peak_calling

##### Usage
```
--peak_calling
```

### Merging of replicates

#### --merge_replicates

##### Usage
```
--merge_replicates
```

### RNA subtype distribution

#### --rna_subtypes

##### Usage
```
--rna_subtypes 3_prime_UTR,transcript,5_prime_UTR
```

#### --gene_id

##### Usage
```
--gene_id ID
```

#### --color_barplot

##### Usage
```
--color_barplot #69b3a2
```

### Peak distance analysis

#### --peak_distance

##### Usage
```
--peak_distance
```

#### --percentile
Shared with sequence exatrction

##### Usage
```
--percentile 90
```

#### --distance

##### Usage
```
--distance 50
```

### Sequence extraction

#### --sequence_extarction

##### Usage
```
--sequence_extraction
```

#### --percentile
Shared with sequence exatrction

##### Usage
```
--percentile 90
```

#### --seq_len

##### Usage
```
--seq_len 10
```

#### --sequence_format_txt

##### Usage
```
--sequence_format_txt
```

## Outputs