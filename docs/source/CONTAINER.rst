.. _section-container:

Container usage
===============

| PARANOiD offers the usage of docker containers via several different executors.
| Currently supported executors are Docker, Podman, and Singularity.
| Since containers are handled via Nextflow, the download should start automatically when running the pipeline.
| The default directory for images is inside the work directory created by Nextflow.
| If problems occur while downloading images via Nextflow, the script :ref:`pull_images.sh <pull-images>` can be used to preload them.
| It is recommended to run PARANOiD with containers to ensure correct versioning of tools.
| If several profiles are used (e.g. using singularity together with SLURM), they are separated by a single comma (e.g. ``-profile apptainer,slurm``).

.. _container-docker:

Docker
------

Uses Docker to run processes within containers

.. code-block:: shell

    -profile docker


.. _container-podman:

Podman
------

Uses Podman to run processes within containers

.. code-block:: shell

    -profile podman


.. _container-singularity:

Singularity
-----------

Uses Singularity to run processes within containers

.. code-block:: shell

    -profile singularity

Apptainer
---------

Uses Apptainer to run processes within containers

.. code-block:: shell

    -profile apptainer
