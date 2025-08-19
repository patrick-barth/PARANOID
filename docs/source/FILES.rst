.. _section-files:

Files included in PARANOiD
=========================

This section provides an overview of the files and directories included in PARANOiD.

1. :ref:`bin <subsection-files-bin>`
2. :ref:`dockerfiles <subsection-files-dockerfiles>`
3. :ref:`docs <subsection-files-docs>`
4. :ref:`modules <subsection-files-modules>`
5. :ref:`build_docker.sh <subsection-files-build-docker>`
6. :ref:`featuretypes-from-gtfgff.awk <subsection-files-featuretypes-from-gffgtf>`
7. :ref:`LICENSE <subsection-files-license>`
8. :ref:`LICENSE.pybam <subsection-files-license-pybam>`
9. :ref:`main.nf <subsection-files-main>`
10. :ref:`PARANOiD-deprecated-DSL1.nf <subsection-files-paranoid-dsl1>`
11. :ref:`pull_images.sh <subsection-files-pull-images>`
12. :ref:`README.md <subsection-files-readme>`

.. _subsection-files-bin:

bin
---
This directory mainly contains custom scripts required for several PARANOiD steps.
This directory is only required if PARANOiD is run without containerization.
Users typically do not need to interact with the files in this directory.

.. _subsection-files-dockerfiles:

dockerfiles
-----------

This directory contains Dockerfiles used to build container images if necessary.
These images can be used to create containers for every PARANOiD step except PureCLIP.
Users typically do not need to interact with the files in this directory

.. _subsection-files-docs:

docs
----

This directory contains files required to build and display this documentation
Users typically do not need to interact with the files in this directory.

.. _subsection-files-modules:

modules
-------

This directory contains all Nextflow modules used by PARANOiD. 
These modules are collections of processes that can be included in the Nextflow workflow.
Each process defines a specific step, including its required and optional inputs, as well as its generated outputs.
Users typically do not need to interact with the files in this directory.

.. _subsection-files-build-docker:

build_docker.sh
---------------

Shell script to automatically build images from all Dockerfiles in the :ref:corresponding directory <subsection-files-dockerfiles> and upload them to Docker Hub.
Users typically do not need to interact with the files in this directory.

.. _subsection-files-featuretypes-from-gffgtf:

featuretypes-from-gtfgff.awk
----------------------------

AWK script to extract all feature types described in a GTF or GFF file.
This can  be useful for the :ref:`RNA subtype analysis <RNA-subtype-analysis>`, which requires exact subtype names.
Usage instructions can be found :ref:`here <determine-feature-types>`.

.. _subsection-files-license:

LICENSE
-------

MIT copyright declaration. Basically says that PARANOiD can be used however you please. You can copy, change and publish this software or parts of it as long as it is under MIT copyright.

.. _subsection-files-license-pybam:

LICENSE.pybam
-------------

Apache copyright declaration which is only valid for pybam, which is used in the process of generating cross-link pile ups from bam files after the alignment.
The Apache copyright allows you to use or change the software as much as you want, as long as you do it under the Apache copyright and make notices on all altered files.

.. _subsection-files-main:

main.nf
-------

Nextflow script to :ref:`run <section-example-run>` when starting a PARANOiD anaylsis.
Uses processes described within the :ref:`modules directory <subsection-files-modules>` and connects them in the right order and with the correct logic to form the pipeline.

.. _subsection-files-config:

nextflow.config
---------------

Config file that is automatically used by PARANOiD (given that it is present in the same directory as the :ref:`main.nf script <subsection-files-main>`).
Consists of 3 parts:

Parameters
^^^^^^^^^^

A list of all :ref:`parameters <section-parameters>` usable when running PARANOiD together and their default values. 
Default parameters can be adapted by users to better suit their needs.

Profiles
^^^^^^^^

Describes usage of :ref:`container executors <section-container>` and :ref:`cluster distribution <section-cluster>`.
The specifications should work on most systems but there is a possibility that they need to be adapted if errors related to the profiles arise.

Resource allocations
^^^^^^^^^^^^^^^^^^^^

Describes the computational resources that will be required to run each process. The current resource requirements are chosen in order to work for most datasets and might not 
be necessary for all datasets. In some cases they might even be set too low; it depends on the size of the :ref:`read file <read-file>` and the :ref:`reference <reference>`.
However, they can (and in some cases should) be adapted if the used system does not meet the required resources which are currently set to 8 cores and 100 GB RAM.
If PARANOiD will be executed on a local computer with less resources available than necessary, the resource requirements can be adapted in this file.
Lowering the required resources can also increase the computing speed as more processes are allowed to be run in parallel.
In this case the file *nextflow.config* can be opened via a text editor and the relevant resource requirements changed.
The most relevant processes will be 'build_index_STAR|mapping_STAR' as they require the highest amount of resources. When opening the config file the relevant entry looks like this::

    withName: 'build_index_STAR|mapping_STAR' {
		cpus = 8
		memory = '100 GB'
		container = 'docker://pbarth/star:1.0'
	}

To change the required cores the number after **cpus = ** needs to be changed - to lower it to 4 cores it should be **cpus = 4**. 
To change the required memory the number after **memory = ** needs to be changed - to lower it to 50 GB it should be **memory = '50 GB'**.

.. _subsection-files-paranoid-dsl1:

PARANOiD-deprecated-DSL1.nf
---------------------------

An older version of PARANOiD that uses DSL1 instead of the later DSL2. Should not be used as it is already deprecated and will not receive any updates in future.

.. _subsection-files-paranoid-galaxy:

PARANOiD_galaxy.xml
-------------------

XML file used to integrate PARANOiD into `Galaxy <https://galaxyproject.org/>`_. Describes all usable parameters, input and output files in order for Galaxy to correctly display 
and translate them into a valid CLI command.

.. _subsection-files-pull-images:

pull_images.sh
--------------

Shell script that can be used to download all images used to build containers by PARANOiD into a specific directory.
Can be used as preparation if PARANOiD is supposed to be run without internet connection.
Additional information on how t run the script can be found :ref:`here <pull-images>`.

.. _subsection-files-readme:

README.md
---------

Readme displayed on `github <https://github.com/patrick-barth/PARANOID>`_. 
Typically there is no need for users to interact with this file.
