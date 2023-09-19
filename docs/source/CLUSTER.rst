.. _section-cluster:

Cluster usage
=============

| PARANOiD supports the distribution of jobs to clusters via job schedulers.
| Currently supported are SGE and SLURM, but the list can potentially expanded upon request. 
| Using job scheduling systems to distribute jobs can immensly shorten execution time.  
| If several profiles are used (e.g. using singularity together with SLURM) they are separated by a single ``,`` (``-profile singularity,slurm``)

.. _cluster-sge:

SGE
------

Uses SGE to distribute jobs.

.. code-block:: shell

    -profile sge


.. _cluster-slurm:

SLURM
-----------

Uses SLURM to distribute jobs.

.. code-block:: shell

    -profile slurm