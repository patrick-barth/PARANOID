Basic usage of PARANOiD
=======================

Description of the minimal setup required for running PARANOiD.

The minimal setup includes: :ref:`FASTQ file <read-file>` containing reads, a :ref:`reference genome <reference>` and
a :ref:`barcode file <barcodes>`. Generated outputs include :ref:`alignments <output-alignments>`, :ref:`cross-link sites <output-cross-link-sites-raw>` in three different file formats, an overview of the :ref:`peak height distribution <output-peak-height-distribution>`, :ref:` processing statistics <output-statistics>`, :ref:`strand distributions <output-strand-distribution>` and an :ref:`IGV session file <output-igv-session>` which can be loaded directly into the IGV to visualize first results. A directory named output (:ref:`unless stated otherwise <output-dir>`) containing all results  will be generated. 
All parts marked with ``<>`` are files that need to be specified by the user.


.. code-block:: shell

   nextflow /path/to/directory/PARANOiD.nf --reads <read-file> --reference <reference-file> --barcodes <barcode-file> --omit_peak_calling
