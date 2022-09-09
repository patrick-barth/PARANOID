# PARANOID
Pipeline for Automated Read ANalysis Of iCLIP Data

//TODO: Blabla ueber Background

## Overview
[Basic usage](#Basic-usage)  
[Parameters](#Parameters)  
[Inputs](#Inputs)  
[Outputs](#Outputs)  

## Basic-usage
nextflow PARANOiD.nf --reads \<reads.fastq\> --reference \<reference_sequence.fasta\> --barcodes \<barcodes.tsv\>

## Parameters

## Inputs

### Barcodes

Barcode sequences are used to assign reads to their experiment. The file is provided as TSV-file (tab separated value).
The first consists of the experiment name and the second of the nucleotide sequence representing the barcode of the experiment. 
One experiment is described per lane and the columns are divided by a tab
The experiment name should be named as follows:
	\<experiment_name\>\_rep_\<replicate-number\>

example:
```
experiment1_rep_1	GCATTG  
experiment1_rep_2	CAGTAA  
experiment1_rep_3	GGCCTA  
experiment2_rep_1	AATCCG  
experiment2_rep_2	CCGTTA  
experiment2_rep_3	GTCATT  
```



## Outputs