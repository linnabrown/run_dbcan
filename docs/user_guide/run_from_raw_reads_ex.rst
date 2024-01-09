More examples: Run from Raw Reads
=================================

.. _wastyk_2021:

Example 2: Wastyk2021 Dataset :cite:`2021:Wastyk`
-------------------------------------------------

The Wastyk2021 dataset :cite:`2021:Wastyk` was published in 20211 from a human dietary intervention study. In the published paper, researchers studied how high-fermented and high-fiber diets influence the human microbiome metabolism and modulate the human immune status. Among various data analyses conducted in the paper1, CAZymes were mined from shotgun metagenomic reads of 18 healthy human participants, and each participant had four time points of stool samples for metagenome sequencing. CAZyme abundance profiles were compared before and after the high-fiber intervention (baseline vs high-fiber). One of the main findings from their CAZyme analysis was that high-fiber consumption increased the CAZyme abundance. For this protocol, we will select two samples (paired-end 2x146bp reads) of two time points (day 2 before high-fiber diet as baseline, and 10 weeks after high-fiber diet as intervention) from one participant (Table S2). The protocol is for the individual sample route (Fig. 3).

The raw read data, intermediate data from each analysis step, and final result data and visualization files are organized in nested folders available on our website https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/, Fig. 5) and https://dbcan.readthedocs.io.

Procedure
---------

Module 1: Reads processing (Fig. 3) to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the Wastyk2021 dataset:

.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1__shotgun_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1__shotgun_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_7__shotgun_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_7__shotgun_2.fastq.gz

Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report fefifo_8022_1.kreport --output fefifo_8022_1.kraken.output fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report fefifo_8022_7.kreport --output fefifo_8022_7.kraken.output fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz

Kraken2 found very little contamination in the data. Consequently, there was no need for the contamination removal step.

If contamination is identified, users can align the reads to the reference genomes of potential contamination source organisms to remove
the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts.

Box 1: Example to remove contamination reads from human
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Kraken2 will produce the following output files:

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 1.1G Sep 21 23:17 fefifo_8022_1.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 991K Sep 21 23:19 fefifo_8022_1.kreport
        -rw-rw-r-- 1 jinfang jinfang 574M Sep 21 23:21 fefifo_8022_7.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 949K Sep 21 23:22 fefifo_8022_7.kreport


    Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

    .. code-block:: shell

        wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa mem hg38 fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz -t 32 -o fefifo_8022_1.hg38.sam
        bwa mem hg38 fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz -t 32 -o fefifo_8022_7.hg38.sam
        samtools view -f 12 fefifo_8022_1.hg38.sam > fefifo_8022_1.hg38.unmap.bam
        samtools view -f 12 fefifo_8022_7.hg38.sam > fefifo_8022_7.hg38.unmap.bam
        samtools fastq -1 fefifo_8022_1_1.clean.fq.gz -2 fefifo_8022_1_2.clean.fq.gz fefifo_8022_1.hg38.unmap.bam
        samtools fastq -1 fefifo_8022_7_1.clean.fq.gz -2 fefifo_8022_7_2.clean.fq.gz fefifo_8022_7.hg38.unmap.bam

P2. Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired fefifo_8022_1__shotgun_1.fastq.gz fefifo_8022_1__shotgun_2.fastq.gz --illumina -j 36
    trim_galore --paired fefifo_8022_7__shotgun_1.fastq.gz fefifo_8022_7__shotgun_2.fastq.gz --illumina -j 36

We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform.
Nonetheless, trim_galore can automatically detect adapters, providing flexibility for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang  429M Oct 30 22:44 fefifo_8022_1__shotgun_1.fastq.gz
        -rw-rw-r-- 1 jinfang jinfang  4.1K Oct 31 05:15 fefifo_8022_1__shotgun_1.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang  390M Oct 31 05:16 fefifo_8022_1__shotgun_1_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang  540M Oct 30 22:44 fefifo_8022_1__shotgun_2.fastq.gz
        -rw-rw-r-- 1 jinfang jinfang  4.2K Oct 31 05:16 fefifo_8022_1__shotgun_2.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang  499M Oct 31 05:16 fefifo_8022_1__shotgun_2_val_2.fq.gz
        -rw-rw-r-- 1 jinfang jinfang  931M Oct 30 22:34 fefifo_8022_7__shotgun_1.fastq.gz
        -rw-rw-r-- 1 jinfang jinfang  4.2K Oct 31 05:17 fefifo_8022_7__shotgun_1.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang  861M Oct 31 05:20 fefifo_8022_7__shotgun_1_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang  1.1G Oct 30 22:34 fefifo_8022_7__shotgun_2.fastq.gz
        -rw-rw-r-- 1 jinfang jinfang  4.4K Oct 31 05:20 fefifo_8022_7__shotgun_2.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 1003M Oct 31 05:20 fefifo_8022_7__shotgun_2_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the ``--retain_unpaired`` parameter in ``trim_galore`` allows for the preservation of single-end reads.
    In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.

P3. Assemble reads into contigs (TIMING ~84min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Megahit for assembling reads into contigs:

.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_fefifo_8022_1 -1 fefifo_8022_1__shotgun_1_val_1.fq.gz -2 fefifo_8022_1__shotgun_2_val_2.fq.gz --out-prefix fefifo_8022_1 --min-contig-len 1000
    megahit -m 0.5 -t 32 -o megahit_fefifo_8022_7 -1 fefifo_8022_7__shotgun_1_val_1.fq.gz -2 fefifo_8022_7__shotgun_2_val_2.fq.gz --out-prefix fefifo_8022_7 --min-contig-len 1000


`MEGAHIT` generates two output folders `megahit_fefifo_8022_1` and `megahit_fefifo_8022_7`.
Each contains five files and one sub-folder (Box 3). `fefifo_8022_1.contigs.fa` is the final contig sequence file.
We set `--min-contig-len 1000`, a common practice to retain all contigs longer than 1,000 base pairs.

Box 3: Example output of `MEGAHIT`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  262 Oct 31 05:49 checkpoints.txt
        -rw-rw-r--  1 jinfang jinfang    0 Oct 31 05:49 done
        -rw-rw-r--  1 jinfang jinfang  97M Oct 31 05:49 fefifo_8022_1.contigs.fa
        -rw-rw-r--  1 jinfang jinfang 149K Oct 31 05:49 fefifo_8022_1.log
        drwxrwxr-x  2 jinfang jinfang 4.0K Oct 31 05:49 intermediate_contigs
        -rw-rw-r--  1 jinfang jinfang 1.1K Oct 31 05:27 options.json


P4. Predict genes by `Prokka` (TIMING ~40h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    prokka --kingdom Bacteria --cpus 36 --outdir prokka_fefifo_8022_1 --prefix fefifo_8022_1 --addgenes --addmrna --locustag fefifo_8022_1 megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa
    prokka --kingdom Bacteria --cpus 36 --outdir prokka_fefifo_8022_7 --prefix fefifo_8022_7 --addgenes --addmrna --locustag fefifo_8022_7 megahit_fefifo_8022_7/fefifo_8022_7.contigs.fa

The parameter `--kingdom Bacteria` is required for bacterial gene prediction.
To optimize performance, `--CPU 36` instructs the utilization of 36 computer processors. The output files comprise of both protein and CDS sequences in Fasta format (e.g., `fefifo_8022_1.faa` and `fefifo_8022_1.ffn` in Box 4).

Box 4: Example output of `Prokka`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  181 Oct 31 22:28 errorsummary.val
        -rw-rw-r--  1 jinfang jinfang 2.6M Oct 31 22:28 fefifo_8022_1.err
        -rw-rw-r--  1 jinfang jinfang  31M Oct 31 22:09 fefifo_8022_1.faa
        -rw-rw-r--  1 jinfang jinfang  83M Oct 31 22:09 fefifo_8022_1.ffn
        -rw-rw-r--  1 jinfang jinfang  22K Oct 31 22:21 fefifo_8022_1.fixedproducts
        -rw-rw-r--  1 jinfang jinfang  98M Oct 31 21:50 fefifo_8022_1.fna
        -rw-rw-r--  1 jinfang jinfang  99M Oct 31 22:09 fefifo_8022_1.fsa
        -rw-rw-r--  1 jinfang jinfang 222M Oct 31 22:24 fefifo_8022_1.gbf
        -rw-rw-r--  1 jinfang jinfang 142M Oct 31 22:09 fefifo_8022_1.gff
        -rw-rw-r--  1 jinfang jinfang 692K Oct 31 22:29 fefifo_8022_1.log
        -rw-rw-r--  1 jinfang jinfang 406M Oct 31 22:22 fefifo_8022_1.sqn
        -rw-rw-r--  1 jinfang jinfang  26M Oct 31 22:09 fefifo_8022_1.tbl
        -rw-rw-r--  1 jinfang jinfang  12M Oct 31 22:09 fefifo_8022_1.tsv
        -rw-rw-r--  1 jinfang jinfang  131 Oct 31 22:09 fefifo_8022_1.txt
        -rw-rw-r--  1 jinfang jinfang 145K Oct 31 22:24 fefifo_8022_1.val

Module 2. run_dbcan annotation (Fig. 3) to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

Users can skip P5 and P6, and directly run P7 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates.

P5. CAZyme annotation at the CAZyme family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.faa protein --hmm_cpu 32 --out_dir fefifo_8022_1.CAZyme --tools hmmer --db_dir db
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.faa protein --hmm_cpu 32 --out_dir fefifo_8022_7.CAZyme --tools hmmer --db_dir db

Two arguments are required for ``run_dbcan``: the input sequence file (faa files) and the sequence type (protein).
By default, ``run_dbcan`` will use three methods (``HMMER`` vs ``dbCAN HMMdb``, ``DIAMOND`` vs ``CAZy``, ``HMMER`` vs ``dbCAN-sub HMMdb``) for
CAZyme annotation (Table 1, Fig. 2). This default setting is equivalent to the use ``--tools all`` parameter (Box 5).
Here we only invoke the ``HMMER`` vs ``dbCAN HMMdb`` for CAZyme annotation at the family level.

Box 5: CAZyme annotation with default setting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the ``--tools`` parameter is not set, it defaults to the equivalent of ``--tools all``.
This setting will take a much longer time to finish (~5 hours) due to the large size of ``dbCAN-sub HMMdb``
(used for substrate prediction for CAZymes, see Table 1).

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.faa protein --out_dir fefifo_8022_1.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 [--tools all]
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.faa protein --out_dir fefifo_8022_7.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 [--tools all]


The sequence type can be `protein`, `prok`, `meta`. If the input sequence file contains metagenomic contig sequences (`fna` file),
the sequence type has to be `meta`, and `prodigal` will be called to predict genes.

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.fna meta --out_dir fefifo_8022_1.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.fna meta --out_dir fefifo_8022_7.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32


P6. CGC prediction (TIMING ~15 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein `faa` and gene location `gff` files.

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_fefifo_8022_1/fefifo_8022_1.gff --out_dir fefifo_8022_1.PUL --dia_cpu 32 --hmm_cpu 32
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_fefifo_8022_7/fefifo_8022_7.gff --out_dir fefifo_8022_7.PUL --dia_cpu 32 --hmm_cpu 32

As mentioned above (see Table 1, Fig. 2), CGC prediction is a featured function added into dbCAN2 in 2018.
To identify CGCs with the protein sequence type, a gene location file (``gff``) must be provided together. If the input sequence type
is ``prok`` or ``meta``, meaning users only have contig ``fna`` files, the CGC prediction can be activated by setting the ``-c cluster`` parameter.

.. warning::

    **Creating own gff file**
    If the users would like to create their own ``gff`` file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the ``gff`` file matches the protein ID in the protein ``faa`` file.

    **[Troubleshooting]CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in ``gff`` file and ``faa`` file do not match.
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.


P7. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs, and their substrates with the `--cgc_substrate` parameter.

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.faa protein --dbcan_thread 32 --tf_cpu 32 --stp_cpu 32 -c prokka_fefifo_8022_1/fefifo_8022_1.gff --cgc_substrate --hmm_cpu 32 --out_dir fefifo_8022_1.dbCAN --dia_cpu 32
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.faa protein --dbcan_thread 32 --tf_cpu 32 --stp_cpu 32 -c prokka_fefifo_8022_7/fefifo_8022_7.gff --cgc_substrate --hmm_cpu 32 --out_dir fefifo_8022_7.dbCAN --dia_cpu 32

The above commands do not set the `--tools` parameter, which means all three methods for CAZyme annotation will be activated (Box 5). Because dbCAN-sub HMMdb (for CAZyme substrate prediction) is 200 times larger than dbCAN HMMdb, the runtime will be much longer. Users can specify `--tools hmmer`, so that the HMMER search against dbCAN-sub will be disabled. However, this will turn off the substrate prediction for CAZymes and CGCs based on CAZyme substrate majority voting. Consequently, the substrate prediction will be solely based on homology search against PULs in dbCAN-PUL (Fig. 1, Table 1).

.. code-block:: shell

    run_dbcan prokka_fefifo_8022_1/fefifo_8022_1.faa protein --tools hmmer --stp_cpu 32 -c prokka_fefifo_8022_1/fefifo_8022_1.gff --cgc_substrate --out_dir fefifo_8022_1.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32
    run_dbcan prokka_fefifo_8022_7/fefifo_8022_7.faa protein --tools hmmer --stp_cpu 32 -c prokka_fefifo_8022_7/fefifo_8022_7.gff --cgc_substrate --out_dir fefifo_8022_7.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32


Box 6. Example output folder content of run_dbcan substrate prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In the `fefifo_8022_1.dbCAN <https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1.dbCAN/>_` directory, a total of 17 files and 1 folder are generated:

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  39M Nov  1 22:18 PUL_blast.out
        -rw-rw-r--  1 jinfang jinfang 3.1M Nov  1 22:15 CGC.faa
        -rw-rw-r--  1 jinfang jinfang 6.9M Nov  1 22:15 cgc.gff
        -rw-rw-r--  1 jinfang jinfang 702K Nov  1 22:15 cgc.out
        -rw-rw-r--  1 jinfang jinfang 321K Nov  1 22:15 cgc_standard.out
        -rw-rw-r--  1 jinfang jinfang 1.5M Nov  1 22:15 cgc_standard.out.json
        -rw-rw-r--  1 jinfang jinfang 556K Nov  1 22:14 dbcan-sub.hmm.out
        -rw-rw-r--  1 jinfang jinfang 345K Nov  1 22:14 diamond.out
        -rw-rw-r--  1 jinfang jinfang 455K Nov  1 22:14 dtemp.out
        -rw-rw-r--  1 jinfang jinfang 298K Nov  1 22:14 hmmer.out
        -rw-rw-r--  1 jinfang jinfang 270K Nov  1 22:15 overview.txt
        -rw-rw-r--  1 jinfang jinfang 1.1M Nov  1 22:15 stp.out
        -rw-rw-r--  1 jinfang jinfang  54K Nov  1 22:18 substrate.out
        drwxrwxr-x  2 jinfang jinfang  32K Nov  2 09:48 synteny.pdf
        -rw-rw-r--  1 jinfang jinfang 288K Nov  1 22:14 tf-1.out
        -rw-rw-r--  1 jinfang jinfang 237K Nov  1 22:14 tf-2.out
        -rw-rw-r--  1 jinfang jinfang 804K Nov  1 22:15 tp.out
        -rw-rw-r--  1 jinfang jinfang  31M Nov  1 21:07 uniInput


    Descriptions of Output Files:

    - ``PUL_blast.out``: BLAST results between CGCs and PULs.
    - ``CGC.faa``: CGC Fasta sequences.
    - ``cgc.gff``: reformatted from the user input gff file by marking CAZymes, TFs, TCs, and STPs.
    - ``cgc.out``: raw output of CGC predictions.

        1.	CGC_id: CGC1
        2.	type: CAZyme
        3.	contig_id: k141_32617
        4.	gene_id: fefifo_8022_1_00137
        5.	start: 1755
        6.	end: 3332
        7.	strand: -
        8.	annotation: GH13

    **Explanation**: Explanation: the gene fefifo_8022_1_00137 encodes a GH13 CAZyme in the CGC1 of the contig k141_32617. CGC1 also has other genes, which are provided in other rows. fefifo_8022_1_00137 is on the negative strand of k141_32617 from 1755 to 3332. The type can be one of the four signature gene types (CAZymes, TCs, TFs, STPs) or the null type (not annotated as one of the four signature genes).

    - ``cgc_standard.out.json``: JSON format of cgc_standard.out.
    - ``dbcan-sub.hmm.out``: HMMER search result against dbCAN-sub HMMdb, including a column with CAZyme substrates extracted from `fam-substrate-mapping-08012023.tsv`.
    - ``diamond.out``: DIAMOND search result against the CAZy annotated protein sequences (`CAZyDB.07262023.fa`).
    - ``dtemp.out``: temporary file.
    - ``hmmer.out``: HMMER search result against dbCAN HMMdb.
    - ``overview.txt``: summary of CAZyme annotation from three methods in TSV format. An example row has the following columns:

        1. ``Gene_ID``: fefifo_8022_1_00719
        2. ``EC#``: PL8_e13:2
        3. ``dbCAN``: PL8_2(368-612)
        4. ``dbCAN_sub``: PL8_e13
        5. ``DIAMOND``: PL8_2
        6. ``#ofTools``: 3

    **Explanation**: Explanation: the protein fefifo_8022_1_00719 is annotated by 3 tools to be a CAZyme: (1) PL8_2 (CAZy defined subfamily 2 of PL8) by HMMER vs dbCAN HMMdb with a domain range from aa position 368 to 612, (2) PL8_e13 (eCAMI defined subfamily e13; e indicates it is from eCAMI not CAZy) by HMMER vs dbCAN-sub HMMdb (derived from eCAMI subfamilies), and (3) PL8_2 by DIAMOND vs CAZy annotated protein sequences. The second column 4.2.2.20:2 is extracted from eCAMI, meaning that the eCAMI subfamily PL8_e13 contains two member proteins which have an EC 4.2.2.20 according to CAZy. In most cases, the 3 tools will have the same CAZyme family assignment. When they give different assignment. We recommend a preference order: dbCAN > eCAMI/dbCAN-sub > DIAMOND. See our dbCAN2 paper2, dbCAN3 paper3, and eCAMI4 for more details.

    **Note**: If users invoked the ``--use_signalP`` parameter when running run_dbcan, there will be an additional column called ``signalP`` in the overview.txt.

    - ``stp.out``: HMMER search result against the MiST5 compiled signal transduction protein HMMs from Pfam.
    - ``tf-1.out``: HMMER search result against the DBD6 compiled transcription factor HMMs from Pfam 7.
    - ``tf-2.out``: HMMER search result against the DBD compiled transcription factor HMMs from Superfamily 8.
    - ``tp.out``: DIAMOND search result against the TCDB 9 annotated protein sequences.
    - ``substrate.out``: summary of substrate prediction results for CGCs in TSV format from two approaches3 (dbCAN-PUL blast search and dbCAN-sub majority voting). An example row has the following columns:

        1. ``CGC_ID``: k141_31366|CGC2
        2. ``Best hit PUL_ID in dbCAN-PUL``: PUL0008
        3. ``Substrate of the hit PUL``: fructan
        4. ``Sum of bitscores for homologous gene pairs between CGC and PUL``: 6132.0
        5. ``Types of homologous gene pairs``: CAZyme-CAZyme;CAZyme-CAZyme;TC-TC;CAZyme-CAZyme;CAZyme-CAZyme;TC-TC
        6. ``Substrate predicted by majority voting of CAZymes in CGC``: fructan
        7. ``Voting score``: 2.0

    **Explanation**: The CGC1 of contig ``k141_31366`` has its best hit ``PUL0008`` (from ``PUL_blast.out``) with fructan as substrate (from ``dbCAN-PUL_12-12-2023.xlsx``). Six signature genes are matched between ``k141_31366|CGC2 and PUL0008 (from PUL_blast.out)``: four are CAZymes and the other two are TCs. The sum of blast bitscores of the six homologous pairs (``CAZyme-CAZyme, CAZyme-CAZyme, TC-TC, CAZyme-CAZyme, CAZyme-CAZyme and TC-TC``) is 6132.0. Hence, the substrate of ``k141_31366|CGC2`` is predicted to be fructan according to dbCAN-PUL blast search. The last two columns are based on the dbCAN-sub result (``dbcan-sub.hmm.out``), according to which two CAZymes in ``k141_31366|CGC2`` are predicted to have fructan substrate. The voting score is thus 2.0, so that according to the majority voting rule, ``k141_31366|CGC2`` is predicted to have a fructan substrate.

    *Note*: : for many CGCs, only one of the two approaches produces substrate prediction. In some cases, the two approaches produce different substrate assignments. We recommend a preference order: ``dbCAN-PUL blast search > dbCAN-sub`` majority voting. See our `dbCAN3 <https://academic.oup.com/nar/article/51/W1/W115/7147496>_` :cite:`2023:dbCAN3` paper3 for more details.

    - ``synteny.pdf``: a folder with syntenic block alignment plots between all CGCs and PULs.
    - ``uniInput``: renamed Fasta file from input protein sequence file.


Module 3. Read mapping (Fig. 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
``````````````````````````````````````````````````````````````````````````````````````````````````````````````

P8. Read mapping to all CDS of each sample (TIMING ~10 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    bwa index prokka_fefifo_8022_1/fefifo_8022_1.ffn
    bwa index prokka_fefifo_8022_7/fefifo_8022_7.ffn
    mkdir samfiles
    bwa mem -t 32 -o samfiles/fefifo_8022_1.CDS.sam prokka_fefifo_8022_1/fefifo_8022_1.ffn fefifo_8022_1__shotgun_1_val_1.fq.gz fefifo_8022_1__shotgun_2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/fefifo_8022_7.CDS.sam prokka_fefifo_8022_7/fefifo_8022_7.ffn fefifo_8022_7__shotgun_1_val_1.fq.gz fefifo_8022_7__shotgun_2_val_2.fq.gz


Reads are mapped to the ``ffn`` files from Prokka.


P9. Read mapping to all contigs of each sample (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    bwa index megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa
    bwa index megahit_fefifo_8022_7/fefifo_8022_7.contigs.fa
    bwa mem -t 32 -o samfiles/fefifo_8022_1.sam megahit_fefifo_8022_1/fefifo_8022_1.contigs.fa fefifo_8022_1__shotgun_1_val_1.fq.gz fefifo_8022_1__shotgun_2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/fefifo_8022_7.sam megahit_fefifo_8022_7/fefifo_8022_7.contigs.fa fefifo_8022_7__shotgun_1_val_1.fq.gz fefifo_8022_7__shotgun_2_val_2.fq.gz


Reads are mapped to the `contig` files from MEGAHIT.

P10. Sort SAM files by coordinates (TIMING ~6min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd samfiles
    samtools sort -@ 32 -o fefifo_8022_1.CDS.bam fefifo_8022_1.CDS.sam
    samtools sort -@ 32 -o fefifo_8022_7.CDS.bam fefifo_8022_7.CDS.sam
    samtools sort -@ 32 -o fefifo_8022_1.bam fefifo_8022_1.sam
    samtools sort -@ 32 -o fefifo_8022_7.bam fefifo_8022_7.sam
    rm -rf *sam
    cd ..


P11. Read count calculation for all proteins of each sample using Bedtools (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    mkdir fefifo_8022_1_abund && cd fefifo_8022_1_abund
    seqkit fx2tab -l -n -i ../prokka_fefifo_8022_1/fefifo_8022_1.ffn | awk '{print $1"\t"$2}' > fefifo_8022_1.length
    seqkit fx2tab -l -n -i ../prokka_fefifo_8022_1/fefifo_8022_1.ffn | awk '{print $1"\t"0"\t"$2}' > fefifo_8022_1.bed
    bedtools coverage -g fefifo_8022_1.length -sorted -a fefifo_8022_1.bed -counts -b ../samfiles/fefifo_8022_1.CDS.bam > fefifo_8022_1.depth.txt

    cd .. && mkdir fefifo_8022_7_abund && cd fefifo_8022_7_abund
    seqkit fx2tab -l -n -i ../prokka_fefifo_8022_7/fefifo_8022_7.ffn | awk '{print $1"\t"$2}' > fefifo_8022_7.length
    seqkit fx2tab -l -n -i ../prokka_fefifo_8022_7/fefifo_8022_7.ffn | awk '{print $1"\t"0"\t"$2}' > fefifo_8022_7.bed
    bedtools coverage -g fefifo_8022_7.length -sorted -a fefifo_8022_7.bed -counts -b ../samfiles/fefifo_8022_7.CDS.bam > fefifo_8022_7.depth.txt
    cd ..


Read counts are saved in ``depth.txt`` files of each sample.


P12. Read count calculation for a given region of contigs using Samtools (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    cd fefifo_8022_1_abund
    samtools index ../samfiles/fefifo_8022_1.bam
    samtools depth -r k141_2168:4235-19858 ../samfiles/fefifo_8022_1.bam > fefifo_8022_1.cgc.depth.txt
    cd ..


The parameter ``-r k141_2168:4235-19858`` specifies a region in a contig. For any CGC, its positional range can be found in the file ``cgc_standard.out`` produced by run_dbcan (Box 6). The ``depth.txt`` files contain the raw read counts for the specified region.


.. warning::

    The contig IDs are automatically generated by MEGAHIT. There is a small chance that a same contig ID appears in both samples. However, the two contigs in the two samples do not match each other even the ID is the same. For example, the contig ID ``k141_2168`` is most likely only found in the ``fefifo_8022_1`` sample. Even if there is a ``k141_2168`` in ``fefifo_8022_7``, the actual contigs in two samples are different.

P13. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_utils fam_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils CGC_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt fefifo_8022_1.depth.txt -i ../fefifo_8022_1.dbCAN -a TPM

    cd .. && cd fefifo_8022_7_abund
    dbcan_utils fam_abundfam_substrate_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils fam_substrate_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils CGC_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    dbcan_utils CGC_substrate_abund -bt fefifo_8022_7.depth.txt -i ../fefifo_8022_7.dbCAN -a TPM
    cd ..


We developed a set of Python scripts as ``dbcan_utils`` (included in the ``run_dbcan`` package) to take the raw read counts for all genes as input and output the normalized abundances (refer to Box 7) of CAZyme families, subfamilies, CGCs, and substrates (see Fig. 4). The parameter ``-a TPM`` can also be set to two other metrics: RPM, or RPKM61.

- **RPKM** is calculated as the number of mapped reads to a gene G divided by [(total number of mapped reads to all genes / 10^6) x (gene G length / 1000)].
- **RPM** is the number of mapped reads to a gene G divided by (total number of mapped reads to all genes / 10^6).
- **TPM** is calculated as [number of mapped reads to a gene G / (gene G length / 1000)] divided by the sum of [number of mapped reads to each gene / (the gene length / 1000)].


Box 7. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example, `fefifo_8022_1_abund <https://bcb.unl.edu/dbCAN_tutorial/dataset2-Wastyk2021/fefifo_8022_1_abund/>_` folder has 7 TSV files:

.. code-block:: shell

    -rw-rw-r--  1 jinfang jinfang 178K Jan  2 04:08 CGC_abund.out
    -rw-rw-r--  1 jinfang jinfang 3.3K Jan  2 04:08 CGC_substrate_majority_voting.out
    -rw-rw-r--  1 jinfang jinfang  12K Jan  2 04:08 CGC_substrate_PUL_homology.out
    -rw-rw-r--  1 jinfang jinfang 2.5K Jan  2 04:08 EC_abund.out
    -rw-rw-r--  1 jinfang jinfang 4.1K Jan  2 04:08 fam_abund.out
    -rw-rw-r--  1 jinfang jinfang  42K Jan  2 04:08 fam_substrate_abund.out
    -rw-rw-r--  1 jinfang jinfang  26K Jan  2 04:08 subfam_abund.out

Explanation of columns in these TSV files is as follows:

    - ``fam_abund.out``: CAZy family (from HMMER vs dbCAN HMMdb), sum of TPM, number of CAZymes in the family.
    - ``subfam_abund.out``: eCAMI subfamily (from HMMER vs dbCAN-sub HMMdb), sum of TPM, number of CAZymes in the subfamily.
    - ``EC_abund.out``: EC number (extracted from dbCAN-sub subfamily), sum of TPM, number of CAZymes with the EC.
    - ``fam_substrate_abund.out``: Substrate (from HMMER vs dbCAN-sub HMMdb), sum of TPM (all CAZymes in this substrate group), GeneID (all CAZyme IDs in this substrate group).
    - ``CGC_abund.out``: CGC_ID (e.g., k141_338400|CGC1), mean of TPM (all genes in the CGC), Seq_IDs (IDs of all genes in the CGC), TPM (of all genes in the CGC), Families (CAZyme family or other signature gene type of all genes in the CGC).
    - ``CGC_substrate_PUL_homology.out``: Substrate (from dbCAN-PUL blast search), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-PUL blast search), TPM (of CGCs in this substrate group).
    - ``CGC_substrate_majority_voting.out``: Substrate (from dbCAN-sub majority voting), sum of TPM, CGC_IDs (all CGCs predicted to have the substrate from dbCAN-sub majority voting), TPM (of CGCs in this substrate group).


Module 4: dbcan_plot for data visualization (Fig. 3) of abundances of CAZymes, CGCs, and substrates (TIMING variable)
`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````

**CRITICAL STEP**

To visualize the CAZyme annotation result, we provide a set of Python scripts as dbcan_plot to make publication quality plots with the dbcan_utils results as the input. The dbcan_plot scripts are included in the run_dbcan package. Once the plots are made in PDF format, they can be transferred to users' Windows or Mac computers for visualization.

Five data folders will be needed as the input for ``dbcan_plot``:

1. two abundance folders ``fefifo_8022_1_abund`` and ``fefifo_8022_7_abund``,
2. two CAZyme annotation ``folders fefifo_8022_1.dbCAN`` and ``fefifo_8022_7.dbCAN``, and
3. the ``dbCAN-PUL folder`` (under the db folder, released from ``dbCAN-PUL.tar.gz``).

P14. Heatmap for CAZyme substrate abundance across samples (Fig. S4B) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot heatmap_plot --samples fefifo_8022_1,fefifo_8022_7 -i fefifo_8022_1_abund/ fam_substrate_abund.out, fefifo_8022_7_abund/fam_substrate_abund.out --show_abund --top 20


Here we plot the top 20 substrates in the two samples. The input files are the two CAZyme substrate abundance files calculated based on dbCAN-sub result. The default heatmap is ranked by substrate abundances. To rank the heatmap according to abundance profile using the function clustermap of seaborn package, users can invoke the ``--cluster_map`` parameter.

P15. Barplot for CAZyme family/subfamily/EC abundance across samples (Fig. S4C) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot bar_plot --samples fefifo_8022_1,fefifo_8022_7 --vertical_bar --top 20 -i fefifo_8022_1_abund/fam_abund.out,fefifo_8022_7_abund/fam_abund.out

Users can choose to generate a barplot instead of heatmap using the ``bar_plot`` method.

P16. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (Fig. S4A) (TIMING 1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    dbcan_plot CGC_synteny_coverage_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1' --readscount fefifo_8022_1_abund/fefifo_8022_1.cgc.depth.txt


The ``fefifo_8022_1.dbCAN`` folder contains the ``PUL_blast.out`` file. Using this file, the ``cgc_standard.out`` file,
and the best PUL's ``gff`` file in ``dbCAN-PUL.tar.gz``, the CGC_synteny_plot method will create the ``CGC-PUL synteny plot``.
The ``-cgcid`` parameter is required to specify which CGC to be plotted (``'k141_2168|CGC1'`` in this example).
The ``fefifo_8022_1.cgc.depth.txt`` file is used to plot the read mapping coverage.

If users only want to plot the CGC structure:

.. code-block:: shell

    dbcan_plot CGC_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1'

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    dbcan_plot CGC_coverage_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1' --readscount fefifo_8022_1_abund/fefifo_8022_1.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    dbcan_plot CGC_synteny_plot -i fefifo_8022_1.dbCAN --cgcid 'k141_2168|CGC1'


.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying ``-i fefifo_8022_1.dbCAN`` is to plot
    the ``'k141_2168|CGC1'`` in the fefifo_8022_1 sample. The ``'k141_2168|CGC1'`` in the fefifo_8022_7 sample most likely does not exist,
    and even it does, the CGC has a different sequence even if the ID is the same.


.. _priest_2023:

Example 3: Priest2023 Dataset :cite:`2023:Priest`
-------------------------------------------------
