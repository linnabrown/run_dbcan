Run with CGCFinder
==================

A CAZyme gene cluster (CGC) refers to a group of genes co-located on the genome that are collectively involved in the metabolism of carbohydrates. These gene clusters encode enzymes and other proteins that work together to perform specific functions related to carbohydrate processing. The concept of a CAZyme gene cluster is particularly relevant in the context of microbial genomes, where such clusters often play crucial roles in the utilization of diverse carbohydrate sources.

Here is an example of how to use run_dbcan to look for CGCs from `Escherichia coli Strain MG1655`_:

Use `-c cluster` to turn on CGCFinder function for complete genome file:

.. code-block:: shell

    run_dbcan EscheriaColiK12MG1655.fna prok -c cluster --out_dir output_EscheriaColiK12MG1655

Or use `-c EscheriaColiK12MG1655.gff` to turn on CGCFinder function for protein sequence since A GFF or BED format file with gene position information is required to run CGCFinder when using a protein input.

.. code-block:: shell

    run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --out_dir output_EscheriaColiK12MG1655

.. _Escherichia coli Strain MG1655: https://www.ncbi.nlm.nih.gov/nuccore/U00096.2
