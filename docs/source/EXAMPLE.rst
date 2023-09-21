.. _section-example-run:

Example run
===========

This page shows the minimal execution of PARANOiD on example data.

.. _subsection-example-download-data:

Download test data
------------------

The example files consist of 2 different experiments and can be downloaded from Zenodo via following link:
https://zenodo.org/record/7733740

Alternatively they can be downloaded using the CLI with following commands:

.. code-block:: shell

    RVFV sample:
    curl "https://zenodo.org/record/7733740/files/barcodes-RVFV.tsv" -o barcodes-RVFV.tsv
    curl "https://zenodo.org/record/7733740/files/virion-reads-M-fragment-only.fastq.gz" -o virion-reads-M-fragment-only.fastq.gz
    curl "https://zenodo.org/record/7733740/files/reference_RVFV.fasta.gz" -o reference_RVFV.fasta.gz
    gzip -d reference_RVFV.fasta.gz
    gzip -d virion-reads-M-fragment-only.fastq.gz

    BHK sample:
    curl "https://zenodo.org/record/7733740/files/barcodes-BHK.tsv" -o barcodes-BHK.tsv
    curl "https://zenodo.org/record/7733740/files/BHK-reads-M-fragment-only.fastq.gz" -o BHK-reads-M-fragment-only.fastq.gz
    curl "https://zenodo.org/record/7733740/files/reference_RVFV.fasta.gz" -o reference_RVFV.fasta.gz
    gzip -d reference_RVFV.fasta.gz
    gzip -d BHK-reads-M-fragment-only.fastq.gz

.. _subsection-example-execute:

Run PARANOiD on test data
-------------------------

To automatically download and then execute PARANOiD the following commands can be used:

.. code-block:: shell
    
    RVFV sample:
    nextflow run patrick-barth/PARANOiD -r main --reads virion-reads-M-fragment-only.fastq --reference reference_RVFV.fasta --barcodes barcodes-RVFV.tsv --output output-RVFV --omit_peak_calling --omit_peak_distance --omit_sequence_extraction  -profile podman

    BHK sample:
    nextflow run patrick-barth/PARANOiD -r main --reads BHK-reads-M-fragment-only.fastq --reference reference_RVFV.fasta --barcodes barcodes-BHK.tsv --output output-BHK --omit_peak_calling --omit_peak_distance --omit_sequence_extraction  -profile podman

In case the resource 
To manually download and execute PARANOiD following commands can be used:

.. code-block:: shell
    
    git clone git@github.com:patrick-barth/PARANOID.git

    RVFV sample:
    nextflow PARANOID/main.nf --reads virion-reads-M-fragment-only.fastq --reference reference_RVFV.fasta --barcodes barcodes-RVFV.tsv --output output-RVFV --omit_peak_calling --omit_peak_distance --omit_sequence_extraction  -profile podman

    BHK sample:
    nextflow PARANOID/main.nf --reads BHK-reads-M-fragment-only.fastq --reference reference_RVFV.fasta --barcodes barcodes-BHK.tsv --output output-BHK --omit_peak_calling --omit_peak_distance --omit_sequence_extraction  -profile podman


If another container execution system is to be used then `podman` can be displaced with `singularity` or `docker` as described :ref:`here <section-container>`.
If the jobs are supposed to be distributed to a cluster the distribution system can be added to the profile argument as described :ref:`here <section-cluster>`. 
Please note that without distributing jobs to a cluster all processes will be calculated locally. 
This currently uses a minimum of 8 cores and 100 GB memory which can exceed the available resources of typical computers. In this case resource usage can be adapted in the config file.

.. _subsection-example-output:

Output
------

The minimal execution of PARANOiD only includes the :ref:`basic analysis <basic-analysis>` and should provide the following outputs if executed correctly:

1. Directory containing alignments
2. Raw cross-link sites
3. Execution metrics
4. An IGV session
5. Distribution of peak heights
6. The reference sequence used for the run
7. Statistics and reports of the run and several processes
8. Strand distributions

More detailed information about the output files can be found :ref:`here <section-output>`.