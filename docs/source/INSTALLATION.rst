Installation
============

To run PARANOiD users need to have `Nextflow <https://www.nextflow.io/docs/latest/install.html>` and one of the container runtimes mentioned in the :ref:`container section <section-container>` installed. PARANOiD can then manually be downloaded from the `GitHub page <https://github.com/patrick-barth/PARANOID>` or via the command line:

.. code-block:: shell

    git clone https://github.com/patrick-barth/PARANOID

Alternatively PARANOiD can be directly started via the following command, which automatically downloads the workflow to ``<HOME>/.nextflow/assets/patrick-barth/PARANOiD``:

.. code-block:: shell

    nextflow run patrick-barth/PARANOiD -r main --reads reads.fastq --reference reference.fasta --barcodes barcodes.tsv --output output --omit_peak_calling --omit_peak_distance --omit_sequence_extraction  -profile apptainer




