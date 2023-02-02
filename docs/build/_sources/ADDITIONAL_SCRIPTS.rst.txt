Supplementary scripts for PARANOiD
==================================

Supplementary scripts for PARANOiD

.. _determine-feature-types:

Determine valid RNA subtypes
----------------------------

Script added to determine valid RNA subtypes for the :ref:`[RNA subtype analysis]<RNA-subtype>`.
Script name: featuretypes-from-gtfgff.awk

Usage:
```
featuretypes-from-gtfgff.awk /path/to/anntotation_file.gff
```

.. _pull-images:

Pull images
-----------

Pulls images via singularity. Images are used to build the containers used by processes.
Should be used to pull all images before starting PARANOiD.

```
pull_images.sh /path/to/PARANOiD/dockerfiles /path/to/image_directory
```