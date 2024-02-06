.. _section-files:

Files present in PARANOiD
=========================

This is an overview of the files and directories that come with PARANOiD.

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

Directory that mainly consists of custom scripts that are needed for several PARANOiD steps. 
This directory is only necessary if no containers are used to execute PARANOiD.
Typically there is no need for users to interact with files in this directory.

.. _subsection-files-dockerfiles:

dockerfiles
-----------

Directory that contains dockerfiles from which container images can be built if necessary. Images built from these dockerfiles can be used to generate containers for every step executed by 
PARANOiD (except PureCLIP).
Typically there is no need for users to interact with files in this directory.

.. _subsection-files-docs:

docs
----

Directory that contains files necessary to build and display this documentation.
Typically there is no need for users to interact with files in this directory.

.. _subsection-files-modules:

modules
-------

Directory that contains all nextflow modules included by PARANOiD. These modules are a collection of processes that can be included in nextflow. 
Each process describes the implementation of a specific step together with the necessary and optional inputs and the generated outputs.
Typically there is no need for users to interact with files in this directory.

.. _subsection-files-build-docker:

build_docker.sh
---------------

Shell script that can be used to automatically build images from all docker files included in the :ref:`correspondent directory <subsection-files-dockerfiles>` 
and upload them to docker hub.
Typically there is no need for users to interact with files in this directory.

.. _subsection-files-featuretypes-from-gffgtf:

featuretypes-from-gtfgff.awk
----------------------------

Short awk script that can be used to get all feature-types described within a gtf or gff file. Can be useful for the :ref:`RNA subtype analysis <RNA-subtype-analysis>` as it needs the exact
subtype names. Usage can be found :ref:`here <determine-feature-types>`.

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