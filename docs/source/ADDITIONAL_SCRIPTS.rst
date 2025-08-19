Supplementary scripts for PARANOiD
==================================

.. _determine-feature-types:

Determine valid RNA subtypes
----------------------------

| Script to determine valid RNA subtypes for the :ref:`RNA subtype analysis <RNA-subtype-analysis>`.
| Script name: featuretypes-from-gtfgff.awk

Usage:

.. code-block:: shell

   featuretypes-from-gtfgff.awk /path/to/annotation_file.gff


.. _pull-images:

Pull images
-----------

| Pulls images via Singularity. These images are used to build the containers executed by individual processes.
| This script should be used to pull all images before starting PARANOiD.

Usage:

.. code-block:: shell

   pull_images.sh /path/to/PARANOiD/dockerfiles /path/to/image_directory
