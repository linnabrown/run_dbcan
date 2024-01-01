Run from Raw Reads
==================

Introduction
------------

dbCAN and run_dbcan require assembled contigs for CAZyme annotation.
Typically, microbiome researchers begin with raw sequencing reads (metagenomic or metatranscriptomic) from various samples.
These reads must be pre-processed and assembled prior to annotation.
Additionally, there's often a need for CAZyme abundance comparison
and visualization across multiple samples. To address these requirements,
this protocol paper provides a comprehensive guide on CAZyme annotation.
It includes steps from initial sequencing reads to the visualization of CAZyme occurrence and abundance across samples.
Key topics covered are software setup, read pre-processing, metagenome assembly, gene prediction,
CAZyme and CGC prediction, glycan substrate prediction, and data visualization.

.. image:: ../_static/img/Picture1.png
   :alt: workflow figure
   :width: 800px
   :align: center


For this tutorial, we provide a comprehensive pipeline to teach users how to run CAZyme annotations from raw reads to generate abundance information.
We use Carter2023 and the individual sample assembly route of the figure above. The procedure has 4 modules and 16 steps (P1-P16).
First, we need to create the environment.

Installation and Data Preparation
---------------------------------


1. Downloading Carter2023 (Table 2) Raw Reads


To download the required raw reads, use the following wget commands:

.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_toturial/raw_reads/Dry2014_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_toturial/raw_reads/Dry2014_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_toturial/raw_reads/Wet2014_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_toturial/raw_reads/Wet2014_2.fastq.gz

2. Create Anaconda Environment


Create and activate a new Anaconda environment with the following steps:

.. code-block:: shell

    conda create -n CAZyme_annotation python=3.9
    conda activate CAZyme_annotation

3. Installing Bioinformatics Dependencies and dbCAN

Install all necessary bioinformatics tools either with a single command or individually:

.. code-block:: shell

    conda install -f dbcan.configure

Alternatively, install the tools one by one:

.. code-block:: shell

    conda install -c conda-forge -c bioconda -c defaults prokka -y
    conda install -c bioconda megahit trim-galore -y
    conda install -c bioconda blast bwa diamond -y
    conda install -c bioconda hmmer -y
    conda install -c bioconda samtools bedtools seqkit -y
    conda install -c bioconda kraken2 -y
    conda install -c agbiome bbtools
    conda install -c bioconda seqtk flye minimap2
    conda install -c conda-forge -c bioconda mmseqs2
    conda install dbcan -c conda-forge -c bioconda


4. Database Installation

To install the databases, execute the following commands:

.. code-block:: shell

    test -d db || mkdir db
        cd db \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08012023.tsv \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL_12112023.faa && mv PUL_12112023.faa PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.xlsx \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_12-12-2023.txt \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
            && wget https://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
            && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/CAZyDB.07262023.fa && diamond makedb --in CAZyDB.07262023.fa -d CAZy \
            && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/dbCAN-HMMdb-V12.txt && mv dbCAN-HMMdb-V12.txt dbCAN.txt && hmmpress dbCAN.txt \
            && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-1.hmm && hmmpress tf-1.hmm \
            && wget http://bcb.unl.edu/dbCAN2/download/Databases/V12/tf-2.hmm && hmmpress tf-2.hmm \
            && wget https://bcb.unl.edu/dbCAN2/download/Databases/V12/stp.hmm && hmmpress stp.hmm \
            && kraken2-build --standard --db K2

The downloaded files must be all in the right location (the db folder).
The CAZyDB.08062022.fa file is needed for DIAMOND search (Table 1).
The dbCAN-HMMdb-V11.txt and dbCAN_sub.hmm files are for HMMER search.
The tcdb.fa, tf-1.hmm, tf-2.hmm, and stp.hmm files are for CGC prediction.
The PUL.faa file consists of protein sequences from experimentally validated PULs for BLAST search to predict substrates for CGCs.
The dbCAN-PUL_07-01-2022.txt and dbCAN-PUL_07-01-2022.xlsx files contain PUL-substrate mapping curated from literature.
Lastly, the fam-substrate-mapping-08252022.tsv file is the family-EC-substrate mapping table for the prediction of CAZyme substrates.

.. warning::
    The conda installation and configuration step may experience prolonged time while resolving environment dependencies. Users should be patient during this process. Alternatively, users consider "mamba",
    another Python package manager that offers similar functionality to Anaconda.
    Information and access to mamba software can be found at https://github.com/mamba-org/mamba.



Module 1: Reads processing to obtain contigs
--------------------------------------------


P1. Contamination Check
^^^^^^^^^^^^^^^^^^^^^^^

Use `kraken2` to check for contaminated reads:

.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report Wet2014.kreport --output Wet2014. kraken.output Wet2014_1.fastq.gz Wet2014_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report Dry2014.kreport --output Dry2014. kraken.output Dry2014_1.fastq.gz Dry2014_2.fastq.gz

Kraken2 found very little contamination in the Carter2023 data. Consequently, there was no need for the contamination removal step.

If contamination is identified, users can align the reads to the reference genomes of potential
contamination source organisms to remove the aligned reads (Box 1). The most common source in human microbiome studies is from human hosts.

Box 1: Removing Contamination Reads from Humans
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Kraken2 will produce the following output files.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 2.0G Dec 12 10:24 Dry2014.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 1.2M Dec 12 10:25 Dry2014.kreport
        -rw-rw-r-- 1 jinfang jinfang 5.1G Dec 12 09:47 Wet2014.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 1.1M Dec 12 09:48 Wet2014.kreport

    Suppose from these files, we have identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

    .. code-block:: shell

        wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa index -p hg38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        bwa mem hg38 Wet2014_1.fastq.gz Wet2014_2.fastq.gz -t 32 -o Wet2014.hg38.sam
        bwa mem hg38 Dry2014_1.fastq.gz Dry2014_2.fastq.gz -t 32 -o Dry2014.hg38.sam
        samtools view -f 12 Wet2014.hg38.sam > Wet2014.hg38.unmap.bam
        samtools view -f 12 Dry2014.hg38.sam > Dry2014.hg38.unmap.bam
        samtools fastq -1 Wet2014_1.clean.fq.gz -2 Wet2014_2.clean.fq.gz Wet2014.hg38.unmap.bam
        samtools fastq -1 Dry2014_1.clean.fq.gz -2 Dry2014_2.clean.fq.gz Dry2014.hg38.unmap.bam

P2. Trimming Adapters and Low-Quality Reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    trim_galore --paired Wet2014_1.fastq.gz Wet2014_2.fastq.gz --illumina -j 36
    trim_galore --paired Dry2014_1.fastq.gz Dry2014_2.fastq.gz --illumina -j 36


Trim_galore is used to trim adapters and low-quality reads.
We specified `--illumina` to indicate that the reads were generated using the Illumina sequencing platform. Nonetheless, trim_galore possesses the ability to automatically detect the adapter,
providing flexibility in adapter handling for users who may know the specific sequencing platform.
Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of `trim_galore`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    In addition to the trimmed read files, `Trim_galore`` also generates a trimming report file.
    The trimming report contains details on read trimming, such as the number of trimmed reads.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 4.2K Dec 13 01:48 Dry2014_1.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.0G Dec 13 01:55 Dry2014_1_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.4K Dec 13 01:55 Dry2014_2.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.4G Dec 13 01:55 Dry2014_2_val_2.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.4K Dec 13 01:30 Wet2014_1.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.4G Dec 13 01:46 Wet2014_1_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.6K Dec 13 01:46 Wet2014_2.fastq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.7G Dec 13 01:46 Wet2014_2_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety.
    Using the --retain_unpaired parameter in trim_galore allows for the preservation of single-end reads.
    In this protocol, this option was not select, so that both reads of a forward-revise pair were removed.

P3. Assemble reads into contigs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use Megahit for assembling reads into contigs:


.. code-block:: shell

    megahit -m 0.5 -t 32 -o megahit_ Wet2014 -1 Wet2014_1_val_1.fq.gz -2 Wet2014_2_val_2.fq.gz --out-prefix Wet2014 --min-contig-len 1000
    megahit -m 0.5 -t 32 -o megahit_ Dry2014 -1 Dry2014_1_val_1.fq.gz -2 Dry2014_2_val_2.fq.gz --out-prefix Dry2014 --min-contig-len 1000


MEGAHIT generates two output folders. Each contains five files and one sub-folder (Box 3).
Wet2014.contigs.fa is the final contig sequence file. We set --min-contig-len 1000,
a common practice to retain all contigs longer than 1,000 base pairs.

Box 3: Example output of `megahit`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  262 Dec 13 04:19 checkpoints.txt
        -rw-rw-r--  1 jinfang jinfang    0 Dec 13 04:19 done
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 13 04:19 intermediate_contigs
        -rw-rw-r--  1 jinfang jinfang 1.1K Dec 13 02:22 options.json
        -rw-rw-r--  1 jinfang jinfang 258M Dec 13 04:19 Wet2014.contigs.fa

P4. Predict Genes with Prokka
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    prokka --kingdom Bacteria --cpus 36 --outdir prokka_Wet2014 --prefix Wet2014 --addgenes --addmrna --locustag Wet2014 megahit_Wet2014/Wet2014.contigs.fa
    prokka --kingdom Bacteria --cpus 36 --outdir prokka_Dry2014 --prefix Dry2014 --addgenes --addmrna --locustag Dry2014 megahit_Dry2014/Dry2014.contigs.fa

The parameter --kingdom Bacteria is required for bacterial gene prediction.
To optimize performance, --CPU 36 instructs the utilization of 36 computer processors.
The output files comprise of both protein and CDS sequences in Fasta format (e.g., Wet2014.faa and Wet2014.ffn in Box 4).


Box 4: Example output of `Prokka`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang 8.4M Dec 14 00:51 Wet2014.err
        -rw-rw-r--  1 jinfang jinfang  75M Dec 13 21:38 Wet2014.faa
        -rw-rw-r--  1 jinfang jinfang 204M Dec 13 21:38 Wet2014.ffn
        -rw-rw-r--  1 jinfang jinfang 259M Dec 13 20:47 Wet2014.fna
        -rw-rw-r--  1 jinfang jinfang 264M Dec 13 21:38 Wet2014.fsa
        -rw-rw-r--  1 jinfang jinfang 599M Dec 14 00:52 Wet2014.gbk
        -rw-rw-r--  1 jinfang jinfang 372M Dec 13 21:38 Wet2014.gff
        -rw-rw-r--  1 jinfang jinfang 2.2M Dec 14 00:52 Wet2014.log
        -rw-rw-r--  1 jinfang jinfang 1.2G Dec 14 00:52 Wet2014.sqn
        -rw-rw-r--  1 jinfang jinfang  68M Dec 13 21:38 Wet2014.tbl
        -rw-rw-r--  1 jinfang jinfang  30M Dec 13 21:38 Wet2014.tsv
        -rw-rw-r--  1 jinfang jinfang  152 Dec 13 21:38 Wet2014.txt

Module 2. run_dbcan annotation to obtain CAZymes, CGCs, and substrates
----------------------------------------------------------------------

P5. CAZyme annotation at family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.faa protein --hmm_cpu 32 --out_dir Wet2014.CAZyme --tools hmmer --db_dir db
    run_dbcan prokka_Dry2014/Dry2014.faa protein --hmm_cpu 32 --out_dir Dry2014.CAZyme --tools hmmer --db_dir db

Two arguments are required for run_dbcan: the input sequence file (faa files) and the sequence type (protein).
By default, run_dbcan will use three methods (HMMER vs dbCAN HMMdb, DIAMOND vs CAZy, HMMER vs dbCAN-sub HMMdb) for CAZyme annotation (Table 1, Figure 2).
This default setting is equivalent to the use --tools all parameter (Box 5).
Here we only invoke the HMMER vs dbCAN HMMdb for CAZyme annotation at the family level.

Box 5: CAZyme annotation with default setting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the --tools parameter is not set, it is the default setting, which is the same as --tools all.
This will take much longer time to finish (~5h) due to the large size of dbCAN-sub HMMdb (used for substrate prediction for CAZymes, see Table 1).

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.faa protein --out_dir Wet2014.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 --tools all
    run_dbcan prokka_Dry2014/Dry2014.faa protein --out_dir Dry2014.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 --tools all

The sequence type can be `protein`, `prok`, `meta`. If the input sequence file contains metagenomic contig sequences (`fna` file),
the sequence type has to be meta, and prodigal will be called to predict genes.

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.fna meta --out_dir Wet2014.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32
    run_dbcan prokka_Dry2014/Dry2014.fna meta --out_dir Dry2014.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32

5.1. Combine proteins from multiple samples

.. warning::
    As shown in Figure 3 (step3), proteins from multiple samples can be combined to generate a non-redundant set of proteins.
    This will reduce the runtime for the run_dbcan step (step4), as only one faa file will be processed.
    However, this does not work for the CGC prediction, as contigs (fna files) from each sample will be needed.
    Therefore, this step (5.1) is recommended if users only want the CAZyme annotation, and not recommended if CGCs are also to be predicted.


This protein sequence clustering step will create a mapping table with sequence cluster ID and protein IDs from each sample.

.. code-block:: shell

    mkdir mmseqs_cluster && cd mmseqs_cluster
    ln -s ../db .
    cat ../prokka_Wet2014/Wet2014.faa ../prokka_Dry2014/Dry2014.faa > Dry_Wet.faa
    mmseqs easy-cluster --threads 32 -c 0.95 --min-seq-id 0.95 --cov-mode 2 Dry_Wet.faa Dry_Wet_cluster tmp
    mv Dry_Wet_cluster_cluster_rep.fasta Dry_Wet.cluster.faa

This `Dry_Wet.cluster.faa` file now contains the non-redundant set of proteins from the two samples.

.. code-block:: shell

    grep "^>" Dry_Wet.cluster.faa | tr ">" " " |awk '{print $1}' > Dry_Wet.geneids
    seqkit grep -f Dry_Wet.geneids ../prokka_Dry2014/Wet2014.ffn > Dry_Wet.ffn
    seqkit grep -f Dry_Wet.geneids ../prokka_Dry2014/Dry2014.ffn >> Dry_Wet.ffn

This `Dry_Wet.ffn file` now contains the CDS sequences of the non-redundant set of proteins from the two samples.

.. code-block:: shell

    bwa index Dry_Wet.ffn
    ln -s ../Dry2014_1_val_1.fq.gz . && ln -s ../Dry2014_2_val_2.fq.gz . && ln -s ../Wet2014_2_val_2.fq.gz . && ln -s ../Wet2014_1_val_1.fq.gz .
    bwa mem -t 32 -o samfiles/Wet2014.CDS.sam Dry_Wet.ffn Wet2014_1_val_1.fq.gz Wet2014 _2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/Dry2014.CDS.sam Dry_Wet.ffn Dry2014_1_val_1.fq.gz Dry2014_2_val_2.fq.gz

The two sam files now contain the read mapping result from each sample to the `Dry_Wet.ffn` file.

P6. CGC prediction.
^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein `faa` and gene location `gff` files.

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_Wet2014/Wet2014.gff --out_dir Wet2014.PUL --dia_cpu 32 --hmm_cpu 32
    run_dbcan prokka_Dry2014/Dry2014.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_ Dry2014/Dry2014.gff --out_dir Dry2014.PUL --dia_cpu 32 --hmm_cpu 32

As mentioned above (Table 1, Figure 2),
CGC prediction is a featured function added into dbCAN2 in 2018.
To identify CGCs with the protein sequence type,
a gene location file (gff) must be provided together.
If the input sequence type is prok or meta, meaning users only have contig fna files, the CGC prediction can be activated by setting -c cluster.

.. warning::

    **Creating own gff file**
    If the users would like to create their own gff file (instead of using Prokka or Prodigal),
    it is important to make sure the value of ID attribute in the gff file matches the protein ID in the protein faa file.

    **CGC not found**
    If no result is found in CGC output file, it is most likely because the sequence IDs in gff file and faa file do not match. Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.

P7. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs, and their substrates with the `--cgc_substrate` parameter.

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.faa protein --dbcan_thread 32 --tf_cpu 32 --stp_cpu 32 -c prokka_Wet2014/Wet2014.gff --cgc_substrate --hmm_cpu 32 --out_dir Wet2014.dbCAN --dia_cpu 32
    run_dbcan prokka_Dry2014/Dry2014.faa protein --dbcan_thread 32 --stp_cpu 32 -c prokka_Dry2014/Dry2014.gff --cgc_substrate --out_dir Dry2014.dbCAN --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32

.. warning::
    The above commands do not set the --tools parameter,
    which means all three methods for CAZyme annotation will be activated (Box 5).
    Because dbCAN-sub HMMdb (for CAZyme substrate prediction) is 200 times larger than dbCAN HMMdb,
    the runtime will be much longer. Users can specify --tools hmmer, so that the HMMER search against dbCAN-sub will be disabled.
    However, this will turn off the substrate prediction for CAZymes and CGCs based on CAZyme substrate majority voting.
    Consequently, the substrate prediction will be solely based on homology search against PULs in dbCAN-PUL

.. code-block:: shell

    run_dbcan prokka_Wet2014/Wet2014.faa protein --tools hmmer --stp_cpu 32 -c prokka_Wet2014/Wet2014.gff --cgc_substrate --out_dir Wet2014.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32
    run_dbcan prokka_Dry2014/Dry2014.faa protein --tools hmmer --stp_cpu 32 -c prokka_Dry2014/Dry2014.gff --cgc_substrate --out_dir Dry2014.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32

.. warning::
    The above commands do not set the --tools parameter, which means all three methods for CAZyme annotation will be activated (Box 5).
    Because dbCAN-sub HMMdb (for CAZyme substrate prediction) is 200 times larger than dbCAN HMMdb, the runtime will be much longer.
    Users can specify --tools hmmer, so that the HMMER search against dbCAN-sub will be disabled.
    However, this will turn off the substrate prediction for CAZymes and CGCs based on CAZyme substrate majority voting.
    Consequently, the substrate prediction will be solely based on homology search against PULs in dbCAN-PUL (Figure 1, Table 1).

    .. code-block:: shell

        run_dbcan prokka_Wet2014/Wet2014.faa protein --tools hmmer --stp_cpu 32 -c prokka_Wet2014/Wet2014.gff --cgc_substrate --out_dir Wet2014.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32
        run_dbcan prokka_Dry2014/Dry2014.faa protein --tools hmmer --stp_cpu 32 -c prokka_Dry2014/Dry2014.gff --cgc_substrate --out_dir Dry2014.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32


Box 6. Example Output Folder Content of run_dbcan Substrate Prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The output directory of run_dbcan substrate prediction typically contains 17 files and 1 folder:

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  33M Dec 17 09:36 blastp.out
        -rw-rw-r--  1 jinfang jinfang 3.3M Dec 17 09:35 CAZyme.pep
        -rw-rw-r--  1 jinfang jinfang  18M Dec 17 09:35 cgc.gff
        -rw-rw-r--  1 jinfang jinfang 836K Dec 17 09:35 cgc.out
        -rw-rw-r--  1 jinfang jinfang 374K Dec 17 09:35 cgc_standard.out
        -rw-rw-r--  1 jinfang jinfang 1.8M Dec 17 09:35 cgc_standard.out.json
        -rw-rw-r--  1 jinfang jinfang 785K Dec 17 09:31 dbsub.out
        -rw-rw-r--  1 jinfang jinfang 511K Dec 17 09:31 diamond.out
        -rw-rw-r--  1 jinfang jinfang 638K Dec 17 09:31 dtemp.out
        -rw-rw-r--  1 jinfang jinfang 414K Dec 17 09:31 hmmer.out
        -rw-rw-r--  1 jinfang jinfang 386K Dec 17 09:35 overview.txt
        -rw-rw-r--  1 jinfang jinfang 2.8M Dec 17 09:35 stp.out
        -rw-rw-r--  1 jinfang jinfang  63K Dec 17 09:36 sub.prediction.out
        drwxrwxr-x  2 jinfang jinfang  36K Dec 17 09:39 syntenic.svg
        -rw-rw-r--  1 jinfang jinfang 799K Dec 17 09:32 tf-1.out
        -rw-rw-r--  1 jinfang jinfang 645K Dec 17 09:34 tf-2.out
        -rw-rw-r--  1 jinfang jinfang 2.3M Dec 17 09:35 tp.out
        -rw-rw-r--  1 jinfang jinfang  75M Dec 17 02:07 uniInput

    Descriptions of Key Output Files:

    - `blastp.out`: BLAST results between CGCs and PULs.
    - `CAZyme.pep`: Fasta sequences of CAZymes.
    - `cgc.gff`: Reformatted user input GFF file, marking CAZymes, TFs, TCs, and STPs.
    - `cgc.out`: Raw output of CGC predictions.
    - `cgc_standard.out`: Simplified version of `cgc.out` in TSV format for easy parsing (refer to Box 7 for columns).
    - `cgc_standard.out.json`: JSON format of `cgc_standard.out`.
    - `dbsub.out`: HMMER search result against dbCAN-sub HMMdb, with CAZyme substrates extracted from fam-substrate-mapping-08252022.tsv.
    - `diamond.out`: DIAMOND search result against the CAZy annotated protein sequences (CAZyDB.08062022.fa).
    - `dtemp.out`: Temporary file.
    - `hmmer.out`: HMMER search result against dbCAN HMMdb.
    - `overview.txt`: Summary of CAZyme annotation from three methods in TSV format (refer to Box 7 for columns).
    - `stp.out`: HMMER search result against the MiST65 compiled signal transduction protein HMMs from Pfam.
    - `tf-1.out` and `tf-2.out`: HMMER search results against transcription factor HMMs from Pfam and Superfamily databases.
    - `tp.out`: DIAMOND search result against the TCDB annotated protein sequences.
    - `sub.prediction.out`: Summary of substrate prediction results (refer to Box 7) for CGCs.
    - `syntenic.svg`: Syntenic block alignment plots between all CGCs and PULs.
    - `uniInput`: Renamed Fasta file from input protein sequence file.



Module 3. Read mapping (Figure 3) to calculate abundance for CAZyme families, subfamilies, CGCs, and substrates
---------------------------------------------------------------------------------------------------------------

P8. Read mapping to all CDS of each sample (TIMING ~20 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    bwa index prokka_Wet2014/Wet2014.ffn
    bwa index prokka_Dry2014/Dry2014.ffn
    mkdir samfiles
    bwa mem -t 32 -o samfiles/Wet2014.CDS.sam prokka_Wet2014/Wet2014.ffn Wet2014_1_val_1.fq.gz Wet2014 _2_val_2.fq.gz
    bwa mem -t 32 -o samfiles/Dry2014.CDS.sam prokka_Dry2014/Dry2014.ffn Dry2014_1_val_1.fq.gz Dry2014_2_val_2.fq.gz

Reads are mapped to the ffn files from Prokka.


P9. Read mapping to all contigs of each sample (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ bwa index megahit_Wet2014/Wet2014.contigs.fa
    $ bwa index megahit_Dry2014/Dry2014.contigs.fa
    $ bwa mem -t 32 -o samfiles/Wet2014.sam megahit_Wet2014/Wet2014.contigs.fa Wet2014_1_val_1.fq.gz Wet2014_2_val_2.fq.gz
    $ bwa mem -t 32 -o samfiles/Dry2014.sam megahit_Dry2014/Dry2014.contigs.fa Dry2014_1_val_1.fq.gz Dry2014_2_val_2.fq.gz

Reads are mapped to the contig files from MEGAHIT.


P10. Sort SAM files by coordinates (TIMING ~8min)

.. code-block:: shell

    $ cd samfiles
    $ samtools sort -@ 32 -o Wet2014.CDS.bam Wet2014.CDS.sam
    $ samtools sort -@ 32 -o Dry2014.CDS.bam Dry2014.CDS.sam
    $ samtools sort -@ 32 -o Wet2014.bam Wet2014.sam
    $ samtools sort -@ 32 -o Dry2014.bam Dry2014.sam
    $ rm -rf *sam
    $ cd ..

P11. Read count calculation for all proteins of each sample using Bedtools (TIMING ~2min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ mkdir Wet2014_abund && cd Wet2014_abund
    $ seqkit fx2tab -l -n -i ../prokka_Wet2014/Wet2014.ffn | awk '{print $1"\t"$2}' > Wet2014.length
    $ seqkit fx2tab -l -n -i ../prokka_Wet2014/Wet2014.ffn | awk '{print $1"\t"0"\t"$2}' > Wet2014.bed
    $ bedtools coverage -g Wet2014.length -sorted -a Wet2014.bed -counts -b ../samfiles/Wet2014.CDS.bam > Wet2014.depth.txt

    $ cd .. && mkdir Dry2014_abund && cd Dry2014_abund
    $ seqkit fx2tab -l -n -i ../prokka_Dry2014/Dry2014.ffn | awk '{print $1"\t"$2}' > Dry2014.length
    $ seqkit fx2tab -l -n -i ../prokka_Dry2014/Dry2014.ffn | awk '{print $1"\t"0"\t"$2}' > Dry2014.bed
    $ bedtools coverage -g Dry2014.length -sorted -a Dry2014.bed  -counts -b ../samfiles/Dry2014.CDS.bam > Dry2014.depth.txt
    $ cd ..

Read counts are saved in depth.txt files of each sample.

P12. Read count calculation for a given region of contigs using Samtools (TIMING ~2min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ cd Wet2014_abund
    $ samtools index ../samfiles/Wet2014.bam
    $ samtools depth -r k141_41392:152403-165349 ../samfiles/Wet2014.bam > Wet2014.cgc.depth.txt
    $ cd ..
    $ cd Dry2014_abund
    $ samtools index ../samfiles/Dry2014.bam
    $ samtools depth -r k141_41392:152403-165349 ../samfiles/Dry2014.bam > Dry2014.cgc.depth.txt

The parameter -r k141_41392:152403-165349 specifies a region in a contig. For any CGC, its positional range can be found in the file cgc_standard.out produced by run_dbcan (Box 6). The depth.txt files contain the raw read counts for the specified region.

P13. dbcan_utils to calculate the abundance of CAZyme families, subfamilies, CGCs, and substrates (TIMING ~1min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ dbcan_utils CAZyme_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    $ dbcan_utils CAZymeSub_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    $ dbcan_utils PUL_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM
    $ dbcan_utils PULSub_abund -bt Wet2014.depth.txt -i ../Wet2014.dbCAN -a TPM

    $ cd .. && cd Dry2014_abund
    $ dbcan_utils CAZyme_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    $ dbcan_utils CAZymeSub_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    $ dbcan_utils PUL_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    $ dbcan_utils PULSub_abund -bt Dry2014.depth.txt -i ../Dry2014.dbCAN -a TPM
    cd ..

We developed a set of Python scripts as dbcan_utils to take the raw read counts for all CDS as input and output the normalized abundances (Box 8) of CAZyme families, subfamilies, CGCs, and substrates (Figure 4). The parameter -a TPM can also be two other metrics: RPM, or FPKM.

Box 8. Example output of dbcan_utils
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Executing these commands will yield five distinct files for each sample: CAZyme_abund_output,
    PUL_abund_output, CAZymeSub_abund_output, PULSub_abund_output.major_voting, and PULSub_abund_output.homo.
    These files encompass the abundance of CAZyme, CGC, and substrates within the respective samples,
    providing detailed information. Notably, users can conveniently trace back the abundance of these entities.
    The abundance calculations adhere to the TPM definition.

Module 4: dbcan_plot for data visualization (Figure 3) of abundances of CAZymes, CGCs, and substrates (TIMING variable)
-----------------------------------------------------------------------------------------------------------------------

To visualize the CAZyme annotation result, we provide a set of Python scripts as dbcan_plot to make publication quality plots with the dbcan_utils results as the input. The dbcan_plot scripts can be installed with commands as follows:

.. code-block:: shell

    $ python3 setup.py install

In addition to the two abundance folders Wet2014_abund and Dry2014_abund, the two CAZyme annotation folders Wet2014.dbCAN and Dry2014.dbCAN, are also needed well as two abundance folders Wet2014_abund.

P14. Heatmap for CAZyme substrate abundance across samples (Figure 6A) (TIMING ~xx)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ dbcan_plot heatmap_plot --samples Wet2014,Dry2014 -i Wet2014_abund/CAZymeSub_abund_output,Dry2014_abund/CAZymeSub_abund_output --show_abund --top 20

Here we plot the top 20 substrates in the two samples.
The input files are the two CAZyme substrate abundance files calculated based on dbCAN-sub result.
The default heatmap is ranked by substrate abundances.
To rank the heatmap according to abundance profile using the xxx clustering algorithm,
users can invoke the `--cluster_map` parameter.

P15. Barplot for CAZyme substrate abundance across samples (Figure 6B) (TIMING ~xx)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: shell

    $ dbcan_plot bar_plot --samples Wet2014,Dry2014 --vertical_bar --top 20 -i Wet2014_abund/CAZyme_abund_output,Dry2014_abund/CAZyme_abund_output

Users can choose to generate a barplot instead of heatmap using the bar_plot method.


P16. Synteny plot between a CGC and its best PUL hit with read mapping coverage to CGC (Figure 6C) (TIMING ~xx)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    $ dbcan_plot CGC_syntenic_with_PUL_abund -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3' --readscount Wet2014_abund/Wet2014.cgc.depth.txt

The Wet2014.dbCAN folder contains the PUL.out file. Using this file, the cgc_standard.out file, and the best PUL’s gff file in dbCAN-PUL.tar.gz, the CGC_synteny_plot method will create the CGC-PUL synteny plot. The –cgcid parameter is required to specify which CGC to be plotted (‘k141_41392|CGC3' in this example). The Wet2014.cgc.depth.txt file is used to plot the read mapping coverage.

If users only want to plot the CGC structure:

.. code-block:: shell

    $ dbcan_plot CGC -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3'

If users only want to plot the CGC structure plus the read mapping coverage:

.. code-block:: shell

    $ dbcan_plot CGC_abund -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3' --readscount Wet2014_abund/Wet2014.cgc.depth.txt

If users only want to plot the synteny between the CGC and PUL:

.. code-block:: shell

    $ dbcan_plot CGC_syntenic_with_PUL -i Wet2014.dbCAN --cgcid 'k141_41392|CGC3'

.. warning::

    The CGC IDs in different samples do not match each other. For example, specifying -i Wet2014.dbCAN is to plot the `k141_41392|CGC3`` in the Wet2014 sample. The `k141_41392|CGC3`` in the Dry2014 sample will be different.
