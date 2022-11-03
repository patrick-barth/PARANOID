PARANOiD Inputs
===============

Detailed description of all input files

.. _read file:

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

``TSV`` file containing experiment names and the corresponding barcode sequences. Reads from the input ``FASTQ`` file are split according to the detected barcode sequence and assigned to the appropriate experiment. This results in one ``FASTQ`` file per experiment.
When choosing the option to merge datasets :[--merge_replicates]:`merge-replicates`

.. _reference:

Reference
---------

.. _annotation:

Annotation
----------