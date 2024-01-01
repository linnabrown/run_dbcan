Run from Protein Sequence
=========================

This section provides an example to run the run_dbcan tool suite with protein sequence data.

To run the dbCAN tool suite on the `Escherichia coli Strain MG1655`_ example data, use the following command. The input file `EscheriaColiK12MG1655.faa` represents the FASTA format complete genome protein sequence, and `prok` specifies that the organism is a prokaryote.

.. code-block:: shell

    run_dbcan EscheriaColiK12MG1655.faa protein --out_dir output_EscheriaColiK12MG1655

.. _Escherichia coli Strain MG1655: https://www.ncbi.nlm.nih.gov/nuccore/U00096.2
