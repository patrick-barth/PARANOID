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

MIT copyright declaration. Basically says that PARANOiD can be used however you please. You can copy, change and publish this software or parts of it as long as the MIT license terms are retained.

.. _subsection-files-license-pybam:

LICENSE.pybam
-------------

Apache copyright declaration which applies only to pybam, which is used to generate cross-link pileups from BAM files after alignment
The Apache License allows you to use or modify the software freely, as long as you comply with the terms of the Apache License and include notices on all modified files

.. _subsection-files-main:

main.nf
-------

Nextflow script that is used to :ref:`run <section-example-run>` a PARANOiD analysis.
It uses processes described in the :ref:`modules directory <subsection-files-modules>` and connects them in the correct order and logical structure to form the pipeline.

.. _subsection-files-config:

nextflow.config
---------------

Configuration file that is automatically used by PARANOiD if it is located in the same directory as the :ref:`main.nf script <subsection-files-main>`.
It consists of three parts:

Parameters
^^^^^^^^^^

A list of all :ref:`parameters <section-parameters>` that can be used when running PARANOiD, along with their default values
Users can modify the default parameters to better suit their needs.
Profiles
^^^^^^^^

Describes the usage of :ref:`container executors <section-container>` and :ref:`cluster distribution <section-cluster>`.
The specifications should work on most systems, but they may need to be adapted if profile-related errors occur.

Resource allocations
^^^^^^^^^^^^^^^^^^^^

Describes the computational resources required to run each process. The current resource requirements are chosen in order to work for most datasets and may not 
be necessary for all datasets. In some cases they might even be set too low; it depends on the size of the :ref:`read file <read-file>` and the :ref:`reference <reference>`.
However, they can (and in some cases should) be adapted if the system in use does not meet the required resources which are currently set to 8 cores and 100 GB of RAM.
When running PARANOiD on a system with limited resources, you may need to adjust the resource settings defined in this file.
Lowering the required resources can also increase overall processing speed as more processes are allowed to be run in parallel.
In this case the file *nextflow.config* can be opened in a text editor and the relevant resource requirements can be adjusted.
The most resource-intensive processes are. 'build_index_STAR|mapping_STAR' as they require the highest amount of resources. When opening the config file the relevant entry looks like this:

    withName: 'build_index_STAR|mapping_STAR' {
		cpus = 8
		memory = '100 GB'
		container = 'docker://pbarth/star:1.0'
	}

To change the number of required CPU cores, adjust the value after **cpus =**  - For example, to reduce it to 4 cores: **cpus = 4**. 
To change the required memory the value after **memory =**  needs to be changed - For example, to reduce it to 50 GB: **memory = '50 GB'**.

.. _subsection-files-paranoid-dsl1:

PARANOiD-deprecated-DSL1.nf
---------------------------

An older version of PARANOiD that uses DSL1 instead of the newer DSL2. It should not be used as it is already deprecated and will not receive any updates in the future.
.. _subsection-files-paranoid-galaxy:

PARANOiD_galaxy.xml
-------------------

XML file used to integrate PARANOiD into `Galaxy <https://galaxyproject.org/>`_.
It describes all usable parameters, input and output files to enable Galaxy to correctly display and translate them into a valid CLI command.

.. _subsection-files-pull-images:

pull_images.sh
--------------

Shell script that can be used to download all images used to build containers for PARANOiD into a specific directory.
Can be used as preparation if PARANOiD is supposed to run without internet connection.
Additional information on how to run the script can be found :ref:`here <pull-images>`.

.. _subsection-files-readme:

README.md
---------

README displayed on `GitHub <https://github.com/patrick-barth/PARANOID>`_. 
Typically, users do not need to interact with this file.
