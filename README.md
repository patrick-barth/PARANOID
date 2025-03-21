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


## Tools used in this workflow 

This section shows all tools and their versions required by PARANOiD.
However, Docker containers are provided for every step and it is highly recommended to use them instead. This way only a Nextflow installation and the [according container software](#container-usage) 
are required. Currently supported are Docker, Podman and Singularity. 

### Essential

These tools are essential to run PARANOiD

| Tool     | Version      |
|----------|--------------|
| Nextflow | 23.04.1.5866 |

### Non-essential

These tools are only required when running PARANOiD without the provided containers. However, running PARANOiD this way is not recommended.

| Tool				| Version	| Note																|
|-------------------|-----------|-------------------------------------------------------------------|
| FastQC			| 0.11.9	|																	|
| Cutadapt			| 4.2		|																	|
| Trim_Galore		| 0.6.7		|																	|
| fastx_toolkit		| 0.0.14	|																	|
| umi_tools			| 1.1.4		|																	|
| Python			| 3.11		|																	|
| Samtools			| 1.16.1	|																	|
| Bamtools			| 2.5.2		|																	|
| wigToBigWig		| 2.9		|																	|
| bigWigToBedGraph	|			|																	|
| Bowtie2			| 2.5.1		| Only when using --domain pro or running the transcript analysis	|
| STAR				| 2.7.10b	| Only when using --domain eu										|
| subread			| 2.0.3		|																	|
| pureCLIP			| 1.3.1		| Only when using peak calling										|
| multiqc			| 1.13		|																	|
| pysam				| 0.19.1	|																	|
| R					| 4.0.3		|																	|
| optparse			|			| R package															|
| wig				|			| R package															|
| reshape2			|			| R package															|
| ggplot2			|			| R package															|
| numpy				|			| Python package													|
| biopython			|			| Python package													|
| ggf3sort.pl		| 1.0.0		| Only when providing an [annotation file](#annotation)				|
| bgzip				| 1.16		| Only when providing an [annotation file](#annotation)				|
| tabix				| 1.16		| Only when providing an [annotation file](#annotation)				|
| meme				| 5.4.1		| Only when using motif detection									|

## Container usage

Containers offer an environment that is separated from your current working environment while providing all tools and dependencies in their correct version necessary to execute a specific task.
This highly reduces the software a user needs to install and lessens the error potential that comes from installing wrong versions, contributing to the reproducibility of the results. Therefore,
it is highly recommended to use containers when running PARANOiD.
Currently supported are Docker, Podman and Singularity. They can be used by adding the following at the end of the command line command:

### Use Docker
```
-profile docker
```
### Use Podman
```
-profile podman
```
### Use Singularity
```
-profile singularity
```


For more detailed descriptions of input files, parameters and implemented analyses please visit the [Read the Docs manual](https://paranoid.readthedocs.io/en/latest/)
