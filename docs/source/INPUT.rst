PARANOiD Inputs
===============

Detailed description of all input files

.. _read-file:

Reads
-----

``FASTQ`` file containing all reads. Each read is represented by 4 lines:

1. **Sequence identifier** and optional description. Starts with a ``@``
2. Actual **nucleotide sequence** of the read
3. **Delimiter line**. Starts with a ``+``
4. **Quality values** of nucleotide sequence (line2). Must contain same number of symbols as line 2

Example:

1. @SEQ_ID
2. GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
3. \+
4. !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65


.. _barcodes:

Barcodes
--------

| ``TSV`` file containing experiment names and the corresponding barcode sequences. Reads from the input ``FASTQ`` file are split according to the detected barcode sequence and assigned to the appropriate experiment. This results in one ``FASTQ`` file per experiment. 
| When choosing the option to merge replicates :ref:`-\-merge_replicates <merge-replicates>` the experiment names have to be chosen appropriately indicating which experiments belong together. In order to do that the appendix ``_rep_<number>`` has to be added to the experiment names, exchanging ``<number>`` with the replicate number. 
| Barcode files consist of 2 columns separated by a single ``tab``: 

1. experiment name 
2. barcode sequence present in reads 

| Experiment names should only consist of Letters ``{a-zA-Z}``, numbers ``{1-9}`` and underscores ``_``. Any whitespaces (e.g. space, tab) will result in errors and thus the termination of the pipeline execution. The length of the barcode sequence is dependant on the protocol used an can be adapted via :ref:`-\-barcode_pattern <barcode-pattern>`.

Example:

| knockdown_N_rep_1    TGATAG 
| knockdown_N_rep_2    AGTGGA 
| knockdown_N_rep_3    GCTCGA 
| mock_N_rep_1    TAAGTA 
| mock_N_rep_2    GCAGTC 
| mock_N_rep_3    CCTAGG

.. _reference:

Reference
---------

``FASTA`` file containing nucleotide data of interest. Is used to align reads to and thus find the location of cross-link sites. Can contain genomic or transcriptomic sequences of an organism or completely artificial sequences.
Every sequence consists of at least 2 lines:
1. Header
2-n. Nucleotide sequence
The header starts with a ``>`` and is followed by a description of the sequence
The sequence consists of nucleotides ``{ACGTN}`` and can span an arbitrary amount of lines

Example:

>NW_024429180.1 Mesocricetus auratus isolate SY011 unplaced genomic scaffold
AACTCTGTTGtaaaaaggctttcccacattcattcCATTCATAAGGTTTCTGTACATTATGGATTCTTTCATGCCTTTTA
AGATGATTATGATATACATAGACTTTAACACCTCAAGAATAttcaggtttctctccagtatgacaATTTGGTCTAATTAT
AAAGAAGAATCAGATATTAAGGTTTTATCACTGTTTACACTCATGCTGTTCCCCTTCATTAAGGTTGGTTTGGATCTTTG
AATATACCTGGGTTCCTATAGTCTCCACCATCACATCTTTATGGAGATTCTTCTGGGAGGGATCCAGCAAATCCCACTCT
\.\.\.

.. _annotation:

Annotation
----------


``GFF`` or ``GTF`` file. Contains annotation information belonging to the reference used in the input. Describes features and their positions. PARANOiD does not rely on the annotation for it's analysis, however it is highly recommended to provid it when working with splicing capable organisms (:ref:`-\-domain eu <domain>`) as annotation files typically contain information about intron-exon structures which highly improve the mapping capability.
Furthermore, providing an annotation file enables the :ref:`RNA subtype analysis<RNA-subtype-analysis>`.
Consists of several header lines followed by one line feature.
Header lines start with a ``#`` and contain general information about the annotation.

| Feature lines consist of 9 columns which are separated by tabs:
| 1. **seqname**: name of the chromosome or scaffold on which the feature is located
| 2. **source**: name of the program that generated this feature or the source
| 3. **feature type**: type of the current feature - e.g. exon - intron - CDS - mRNA - 3_prime_UTR - transcript - 5_prime_UTR
| 4. **start**: Start position of the feature (1-based)
| 5. **end**: End position of the feature (1-based)
| 6. **score**: Float point value (can also simply be a ``.``)
| 7. **strand**: Stradn on which the feature is present. ``+`` for forward; ``-`` for reverse
| 8. **frame**: Indicates which base of the feature is actually the first base of a codon. 0 -> the first base of the feature is the first abse of a codon; 1 -> the second base of the feature is the first base of a codon .... (can slo simply be a ``.``)
| 9. **attributes**: Semicolon separated list of tag-value pairs providing additional information

Example:

| ##gff-version 3 
| #!gff-spec-version 1.21 
| #!processor NCBI annotwriter 
| #!genome-build BCM_Maur_2.0 
| #!genome-build-accession NCBI_Assembly:GCF_017639785.1 
| #!annotation-source NCBI Mesocricetus auratus Annotation Release 103 
| ##sequence-region NW_024429180.1 1 52462669 
| ##species \https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10036
| NW_024429180.1	RefSeq	region	1	52462669	.	+	.	ID=NW_024429180.1:1..52462669;Dbxref=taxon:10036;Name=Unknown;chromosome=Unknown;dev-stage=adult;gbkey=Src;genome=genomic;isolate=SY011;mol_type=genomic DNA;sex=female;tissue-type=liver 
| NW_024429180.1	Gnomon	pseudogene	37366	38359	.	+	.	ID=gene-LOC101842720;Dbxref=GeneID:101842720;Name=LOC101842720;gbkey=Gene;gene=LOC101842720;gene_biotype=pseudogene;pseudo=true