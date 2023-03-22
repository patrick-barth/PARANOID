# PARANOID
Pipeline for Automated Read ANalysis Of iCLIP Data

PARANOiD is a versatile software for the fully automated analysis of iCLIP and iCLIP2 data. It contains all steps necessary for preprocessing, the determination of cross-link locations and several additional steps which can be used to detect specific characteristics, e.g. definite distances between cross-link events or binding motifs. The cross-link sites are presented as WIG files that can be easily visualized e.g. using IGV, for which a config file is offered. Additionally, results are offered as statistical plots for a quick overview and as standardized bioinformatics file formats or TSV files which can be used for further analysis steps. 


## Basic-usage
```
nextflow PARANOiD.nf --reads \<reads.fastq\> --reference \<reference_sequence.fasta\> --barcodes \<barcodes.tsv\>
```

## Inputs

### Reads (essential)

Reads generated by iCLIP experiments. Can be provided as one or more files. If providing more than one file, regular expressions can be used within quotation marks.  
Format: FASTQ

#### Usage
```
--reads reads_file.fastq
--reads "reads_{1,2}.fastq"
--reads "*.fastq"
```

### Reference (essential)

File containing the reference to which the reads will be mapped.  
Format: FASTA

#### Usage
```
--reference reference_file.fasta
```

### Barcodes (essential)

Barcode sequences are used to assign reads to their experiment. The file is provided as TSV-file (tab separated value).
The first consists of the experiment name and the second of the nucleotide sequence representing the barcode of the experiment. 
One experiment is described per lane and the columns are divided by a tab.  
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

File containing annotations of the reference provided. Advised when working with splicing capable organisms. Necessary for RNA subtype analysis.  
Formats: GFF GTF

#### Usage
```
--annotation annotation_file.gff
```

For more detailed descriptions of input files, parameters and implemented anaylses please visit the [Read the Docs manual](https://paranoid.readthedocs.io/en/readthedocs/)