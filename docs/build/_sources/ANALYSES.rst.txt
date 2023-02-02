Included Analyses
=================

Short overview of all analyses implemented in PARANOiD

.. _merge_replicates:

Merge replicates
----------------

Merges several replicates into a single representative version which can be used for publications, posters or presentations. This version shows the mean hit count for every position. Additionally, a correllation analysis is performed to give the user an evaluation of the sample similarity and therefore a rationale for this analysis.
Is deactivated by default. Can be activated via :ref:`[--merge_replicates]<merge_replicates>`

.. _RNA_subtype:

RNA subtypes
------------

Analysis to determine if the protein of interest is prone to bind to specific RNA subtypes or regions. As this is determined via the :ref:`[annotation file]<annotation>` only subtypes included there can be determined (shown in column 3). To see which RNA subtypes are included in the annotation file a :ref:`[script]<determine_feature_types>` was added. When choosing RNA subtypes one has to be careful not to use subtypes that are hierarchically higher or lower to each other as these will at least partially cover the same reference regoins making hits in these regions ambiguous. The `SO ontologies <https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/Ontology_Files/subsets/SOFA.obo>` can be used to get an overview of the official hierarchical structures of annotation files. 
Is activated when an :ref:`[annotation file]<annotation>` is provided.

.. _transcript-analysis:

Transcript analysis
-------------------

Analysis to show if specific RNAs are more prone to interact with the the protein of interest. 
If choosing this analysis a file containing all RNAs of interest should be used as input reference instead of the genome. Here all RNAs of interest (or artificial RNAs present in the sample) can be combined to a single fasta file. If the general transcriptome of an organism shall be examined, they can often be accessed next to the genome and annotation of the organism. If not a FASTA file containing the transcripts can be generated as follows (needs the genome and an annotation file):

```
--map_to_transcripts
```

.. _peak_calling:

Peak calling
------------

.. _peak_distance:

Peak distance analysis
----------------------

.. _motif_detection:

Motif detection
---------------