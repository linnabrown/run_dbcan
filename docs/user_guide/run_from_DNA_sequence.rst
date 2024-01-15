Run from DNA Sequence
=========================

This section provides a quick guide to running the run_dbcan tool suite with example data and explains the output files generated.


To run the dbCAN tool suite on the `Escherichia coli Strain MG1655`_ example data, use the following command. The input file `EscheriaColiK12MG1655.fna` represents the FASTA format complete genome DNA sequence, and `prok` specifies that the organism is a prokaryote.

.. code-block:: shell

    run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655

.. _Escherichia coli Strain MG1655: https://www.ncbi.nlm.nih.gov/nuccore/U00096.2
