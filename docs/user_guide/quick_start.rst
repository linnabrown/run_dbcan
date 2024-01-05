Quick Start
===========

This section provides a quick guide to running the run_dbcan tool suite with example data and explains the output files generated.

1. Running Example Data
-----------------------

To run the dbCAN tool suite on the `Escherichia coli Strain MG1655`_ example data, use the following command. The input file `EscheriaColiK12MG1655.fna` represents the FASTA format complete genome DNA sequence, and `prok` specifies that the organism is a prokaryote.

.. code-block:: shell

    run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655

.. _Escherichia coli Strain MG1655: https://www.ncbi.nlm.nih.gov/nuccore/U00096.2

2. Understanding the Output
---------------------------

After running the tool, several output files are generated in `output_EscheriaColiK12MG1655`, each with specific information:

**uniInput**
  The unified input file for subsequent tools, created by Prodigal if a nucleotide sequence is used.

**dbsub.out**
  Output from the dbCAN_sub run.

**diamond.out**
  Results from the Diamond BLAST.

**hmmer.out**
  Output from the HMMER run.

**tf.out**
  Diamond BLAST output predicting Transcription Factors (TFs) for CGCFinder.

**tc.out**
  Diamond BLAST output predicting Transporter Classifications (TCs) for CGCFinder.

**cgc.gff**
  GFF input file for CGCFinder.

**cgc.out**
  Output from the CGCFinder run.

**cgc_standard.out**
  Simplified version of cgc.out, containing columns like CGC_id, Type, Contig_id, Gene_id, Start, End, Strand, and Annotation.

**overview.txt**
  Summarizes CAZyme predictions across tools, including SignalP results.
