More examples: Run from Raw Reads
=================================

.. _priest_2023:

Example 3: Priest2023 Dataset :cite:`2023:Priest`
-------------------------------------------------

The Wastyk2021 dataset :cite:`2023:Priest` was published in 2021 from a human dietary intervention study. In the published paper, researchers studied how high-fermented and high-fiber diets influence the human microbiome metabolism and modulate the human immune status. 
Among various data analyses conducted in the paper, 
CAZymes were mined from shotgun metagenomic reads of 18 healthy human participants, and each participant had four time points of stool samples for metagenome sequencing. 
CAZyme abundance profiles were compared before and after the high-fiber intervention (baseline vs high-fiber). One of the main findings from their CAZyme analysis was that high-fiber consumption increased the CAZyme abundance. 
For this protocol, we will select two samples (paired-end 2x146bp reads) of two time points (day 2 before high-fiber diet as baseline, and 10 weeks after high-fiber diet as intervention) from one participant
The protocol is for the individual sample route.

============  ===========  ===========
Header 1      Header 2     Header 3
============  ===========  ===========
row 1, col 1  row 1, col 2  row 1, col 3
row 2, col 1  row 2, col 2  row 2, col 3
============  ===========  ===========


The Priest2023 dataset :cite:`2021:Wastyk` was published in 20211 from a human dietary intervention study. In the published paper, researchers studied how high-fermented and high-fiber diets influence the human microbiome metabolism and modulate the human immune status. Among various data analyses conducted in the paper1, CAZymes were mined from shotgun metagenomic reads of 18 healthy human participants, and each participant had four time points of stool samples for metagenome sequencing. CAZyme abundance profiles were compared before and after the high-fiber intervention (baseline vs high-fiber). One of the main findings from their CAZyme analysis was that high-fiber consumption increased the CAZyme abundance. For this protocol, we will select two samples (paired-end 2x146bp reads) of two time points (day 2 before high-fiber diet as baseline, and 10 weeks after high-fiber diet as intervention) from one participant (Table S2). The protocol is for the individual sample route (Fig. 3).

Procedure
---------

Module 1: Reads processing (Fig. 3) to obtain contigs
`````````````````````````````````````````````````````

P1. Contamination Check (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the Priest2023 dataset:
.. code-block:: shell

    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MG.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MG.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MT_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MT_2.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MT_1.fastq.gz
    wget https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/SRF_MT_2.fastq.gz

Use `kraken2` to check for contaminated reads:


.. code-block:: shell

    kraken2 --threads 32 --quick --paired --db K2 --report SRF_MT.kreport --output SRF_MT.kraken.output SRF_MT_1.fastq.gz SRF_MT_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report BML_MT.kreport --output BML_MT.kraken.output BML_MT_1.fastq.gz BML_MT_2.fastq.gz

    kraken2 --threads 32 --quick --paired --db K2 --report SRF_MG.kreport --output SRF_MG.kraken.output SRF_MG_1.fastq.gz SRF_MG_2.fastq.gz
    kraken2 --threads 32 --quick --paired --db K2 --report BML_MG.kreport --output BML_MG.kraken.output BML_MG_1.fastq.gz BML_MG_2.fastq.gz

Kraken2 found much contamination (Box 1) from human in the Priest2023 data. Consequently, human reads need to be removed before assembly. 

Reads can be aligned to the reference genomes of potential contamination source organisms to remove the aligned reads. The most common source in microbiome studies is from human. 

Box 1: Example of Kraken2 output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    The `kreport` files can be examined to identify potential contamination source organisms.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang   54M Dec 27 07:42 BML_MG.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.2M Dec 27 07:42 BML_MG.kreport
        -rw-rw-r-- 1 jinfang jinfang  3.4G Dec 27 08:01 BML_MT.kraken.output
        -rw-rw-r-- 1 jinfang jinfang 1023K Dec 27 08:02 BML_MT.kreport
        -rw-rw-r-- 1 jinfang jinfang   61M Dec 27 07:39 SRF_MG.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.2M Dec 27 07:39 SRF_MG.kreport
        -rw-rw-r-- 1 jinfang jinfang  2.6G Dec 27 07:50 SRF_MT.kraken.output
        -rw-rw-r-- 1 jinfang jinfang  1.1M Dec 27 07:51 SRF_MT.kreport 


P2. Remove contamination reads from human (TIMING ~40min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the Kraken2 output files, we identified humans as the contamination source, we can use the following commands to remove the contamination reads by aligning reads to the human reference genome.

.. code-block:: shell

    mkdir hg38 && cd hg38 && wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    cd .. && mkdir contamination && cd contamination
    minimap2 -a -x map-hifi -MD -t 32 -o SRF_MG.hg38.sam ../hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ../SRF_MG.fastq.gz
    minimap2 -a -x map-hifi -MD -t 32 -o BML_MG.hg38.sam ../hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ../BML_MG.fastq.gz
    samtools fastq -f 4 -@ 32 -0 ../SRF_MG.clean.fq.gz SRF_MG.hg38.sam
    samtools fastq -f 4 -@ 32 -0 ../BML_MG.clean.fq.gz BML_MG.hg38.sam
    bwa mem ../hg38/hg38 ../SRF_MT_1.fastq.gz ../SRF_MT_2.fastq.gz -t 32 -o SRF_MT.hg38.sam
    bwa mem ../hg38/hg38 ../BML_MT_1.fastq.gz ../BML_MT_2.fastq.gz -t 32 -o BML_MT.hg38.sam
    samtools fastq -f 12 -@ 32 -1 ../SRF_MT_1.clean.fq.gz -2 ../SRF_MT_2.clean.fq.gz SRF_MT.hg38.sam
    samtools fastq -f 12 -@ 32 -1 ../BML_MT_1.clean.fq.gz -2 ../BML_MT_2.clean.fq.gz BML_MT.hg38.sam
    cd ..

P3| Trim adapter and low-quality reads (TIMING ~20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. code-block:: shell

    trim_galore --illumina -j 8 --paired BML_MT_1.clean.fastq.gz BML_MT_2.clean.fastq.gz
    trim_galore --illumina -j 8 --paired SRF_MT_1.clean.fastq.gz SRF_MT_2.clean.fastq.gz

The HiFi long reads do not need to be trimmed. Hence, this step only applies to MT illumina short read data. We specified --illumina to indicate that the reads were generated using the Illumina sequencing platform. Nonetheless, trim_galore possesses the ability to automatically detect the adapter, providing flexibility in adapter handling for users who may know the specific sequencing platform. Details of trimming are available in the trimming report file (Box 2).

Box 2: Example output of trim_galore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    In addition to the trimmed read files, Trim_galore also generates a trimming report file. The trimming report contains details on read trimming, such as the number of trimmed reads.

    .. code-block:: shell

        -rw-rw-r-- 1 jinfang jinfang 4.2K Dec 28 21:56 BML_MT_1.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.3G Dec 28 22:05 BML_MT_1.clean_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.7K Dec 28 22:05 BML_MT_2.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.0G Dec 28 22:05 BML_MT_2.clean_val_2.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 4.9K Dec 28 10:07 SRF_MT_1.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 2.7G Dec 28 10:19 SRF_MT_1.clean_val_1.fq.gz
        -rw-rw-r-- 1 jinfang jinfang 5.1K Dec 28 10:19 SRF_MT_2.clean.fq.gz_trimming_report.txt
        -rw-rw-r-- 1 jinfang jinfang 3.3G Dec 28 10:19 SRF_MT_2.clean_val_2.fq.gz

.. warning::

    During the trimming process, certain reads may be entirely removed due to low quality in its entirety. Using the `--retain_unpaired` parameter in trim_galore allows for the preservation of single-end reads. In this protocol, this option was not selected, so that both reads of a forward-revise pair were removed.



P4. Assemble HiFi reads into metagenome (TIMING ~4h20min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flye was used to assemble the HiFi long reads into contigs. 

.. code-block:: shell

    flye --threads 32 --meta --pacbio-hifi BML_MG.clean.fq.gz --hifi-error 0.01 --keep-haplotypes --out-dir flye_BML_MG
    flye --threads 32 --meta --pacbio-hifi SRF_MG.clean.fq.gz --hifi-error 0.01 --keep-haplotypes --out-dir flye_SRF_MG

Flye generates two folders `flye_BML_MG` and `flye_SRF_MG`. Each folder 
contains 6 files and 5 sub-folders (Box 3), among them `assembly.fasta` is the final contig sequence file. 
We set `--hifi-error` 0.01, a generally accepted error rate of HiFi sequencing. 
Parameter `--meta` is set to assemble reads into metagenomes.

Box 3: Example output of Flye
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. code-block:: shell

        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 20:15 00-assembly
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 20:43 10-consensus
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 21:14 20-repeat
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 21:16 30-contigger
        drwxrwxr-x  2 jinfang jinfang 4.0K Dec 27 22:06 40-polishing
        -rw-rw-r--  1 jinfang jinfang 314M Dec 27 22:06 assembly.fasta
        -rw-rw-r--  1 jinfang jinfang 311M Dec 27 22:06 assembly_graph.gfa
        -rw-rw-r--  1 jinfang jinfang 6.6M Dec 27 22:06 assembly_graph.gv
        -rw-rw-r--  1 jinfang jinfang 867K Dec 27 22:06 assembly_info.txt
        -rw-rw-r--  1 jinfang jinfang  61M Dec 27 22:06 flye.log
        -rw-rw-r--  1 jinfang jinfang   92 Dec 27 22:06 params.json


P5. Predict genes by Prokka (~21h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   prokka --outdir prokka_BML_MG --prefix BML_MG --addgenes --addmrna --locustag BML_MG --kingdom Bacteria --cpus 36 flye_BML_MG/assembly.fasta
   prokka --outdir prokka_SRF_MG --prefix SRF_MG --addgenes --addmrna --locustag SRF_MG --kingdom Bacteria --cpus 36 flye_SRF_MG/assembly.fasta

The parameter `--kingdom` Bacteria is required for bacterial gene prediction. 
To optimize performance, `--CPU` 36 instructs the utilization of 36 computer processors. 
The output files comprise of both protein and CDS sequences in Fasta format (e.g., `BML_MG.faa` and `SRF_MG.ffn` in Box 4).

Box 3: Example output of Prokka 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. code-block:: shell
        
        -rw-rw-r--  1 jinfang jinfang 2.2M Dec 28 05:38 BML_MG.err
        -rw-rw-r--  1 jinfang jinfang 105M Dec 27 23:26 BML_MG.faa
        -rw-rw-r--  1 jinfang jinfang 288M Dec 27 23:26 BML_MG.ffn
        -rw-rw-r--  1 jinfang jinfang 314M Dec 27 22:06 BML_MG.fna
        -rw-rw-r--  1 jinfang jinfang 315M Dec 27 23:26 BML_MG.fsa
        -rw-rw-r--  1 jinfang jinfang 724M Dec 28 05:39 BML_MG.gbk
        -rw-rw-r--  1 jinfang jinfang 467M Dec 27 23:26 BML_MG.gff
        -rw-rw-r--  1 jinfang jinfang 1.9M Dec 28 05:39 BML_MG.log
        -rw-rw-r--  1 jinfang jinfang 1.5G Dec 28 05:39 BML_MG.sqn
        -rw-rw-r--  1 jinfang jinfang  89M Dec 27 23:26 BML_MG.tbl
        -rw-rw-r--  1 jinfang jinfang  40M Dec 27 23:26 BML_MG.tsv
        -rw-rw-r--  1 jinfang jinfang  152 Dec 27 23:26 BML_MG.txt


Module 1: run_dbcan annotation (Fig. 3) to obtain CAZymes, CGCs, and substrates
```````````````````````````````````````````````````````````````````````````````````````````````

Users can skip P6 and P7, and directly run P8 (much slower though), if they want to predict not only CAZymes and CGCs, but also substrates. 

P6. CAZyme annotation at family level (TIMING ~10min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

   run_dbcan prokka_BML_MG/BML_MG.faa protein --hmm_cpu 32 --out_dir BML_MG.CAZyme --tools hmmer --db_dir db
   run_dbcan prokka_SRF_MG/SRF_MG.faa protein --hmm_cpu 32 --out_dir SRF_MG.CAZyme --tools hmmer --db_dir db

Two arguments are required for run_dbcan: the input sequence file (faa) and the sequence type (protein). By default, run_dbcan will use three methods (HMMER vs dbCAN HMMdb, DIAMOND vs CAZy, HMMER vs dbCAN-sub HMMdb) for CAZyme annotation (Table 1, Fig. 2). This default setting is equivalent to the use --tools all parameter (Box 5). Here we only invoke the HMMER vs dbCAN HMMdb for CAZyme annotation at the family level. 


Box 3: CAZyme annotation with default setting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    If the `--tools` parameter is not set, it is the default setting, which is the same as `--tools` all. 
    This will take much longer time to finish (~5h) due to the large size of dbCAN-sub HMMdb (used for substrate prediction for CAZymes, see Table 1).

    .. code-block:: shell

       run_dbcan prokka_BML_MG/BML_MG.faa protein --out_dir BML_MG.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 --tools all
       run_dbcan prokka_SRF_MG/SRF_MG.faa protein --out_dir SRF_MG.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32 --tools all

    The sequence type can be protein, prok, meta. If the input sequence file contains metagenomic contig sequences (fna file), 
    the sequence type has to be meta, and prodigal will be called to predict genes. 

    .. code-block:: shell

        run_dbcan prokka_BML_MG/BML_MG.fna meta --out_dir BML_MG.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32
        run_dbcan prokka_SRF_MG/SRF_MG.fna meta --out_dir SRF_MG.CAZyme --dia_cpu 32 --hmm_cpu 32 --dbcan_thread 32
    
P7. CGC prediction (TIMING ~15 min)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to not only predict CAZymes but also CGCs with protein faa and gene location gff files.

 .. code-block:: shell

    run_dbcan prokka_BML_MG/BML_MG.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_BML_MG/BML_MG.gff --out_dir BML_MG.PUL --dia_cpu 32 --hmm_cpu 32 
    run_dbcan prokka_SRF_MG/SRF_MG.faa protein --tools hmmer --tf_cpu 32 --stp_cpu 32 -c prokka_SRF_MG/SRF_MG.gff --out_dir SRF_MG.PUL --dia_cpu 32 --hmm_cpu 32 

As mentioned above (Table 1, Fig. 2), CGC prediction is a featured function added into dbCAN2 in 2018. 
To identify CGCs with the protein sequence type, a gene location file (gff) must be provided together. 
If the input sequence type is prok or meta, meaning users only have contig fna files, 
the CGC prediction can be activated by setting `-c cluster`.

.. warning::

    **CAUTION **

    If the users would like to create their own gff file (instead of using Prokka or Prodigal), 
    it is important to make sure the value of ID attribute in the gff file matches the protein ID in the protein faa file. 

    **Troubleshooting**

    If no result is found in CGC output file, it is most likely because the sequence IDs in gff file and faa file do not match. 
    Another less likely reason is that the contigs are too short and fragmented and not suitable for CGC prediction.

P8. Substrate prediction for CAZymes and CGCs (TIMING ~5h)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following commands will re-run run_dbcan to predict CAZymes, CGCs, 
and their substrates with the `--cgc_substrate` parameter.

.. code-block:: shell
    run_dbcan prokka_BML_MG/BML_MG.faa protein --dbcan_thread 32 --tf_cpu 32 --stp_cpu 32 -c prokka_BML_MG/BML_MG.gff --cgc_substrate --hmm_cpu 32 --out_dir BML_MG.dbCAN --dia_cpu 32 
    run_dbcan prokka_SRF_MG/SRF_MG.faa protein --dbcan_thread 32 --stp_cpu 32 -c prokka_SRF_MG/SRF_MG.gff --cgc_substrate --out_dir SRF_MG.dbCAN --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32 

.. warning::

    The above commands do not set the `--tools` parameter, 
    which means all three methods for CAZyme annotation will be activated (Box 5). 
    Because dbCAN-sub HMMdb (for CAZyme substrate prediction) is 200 times larger than dbCAN HMMdb, 
    the runtime will be much longer. Users can specify `--tools` hmmer, 
    so that the HMMER search against dbCAN-sub will be disabled. 
    However, this will turn off the substrate prediction for CAZymes and CGCs based on CAZyme substrate majority voting. 
    Consequently, 
    the substrate prediction will be solely based on homology search against PULs in dbCAN-PUL (Fig. 1, Table 1). 

.. code-block:: shell

    run_dbcan prokka_BML_MG/BML_MG.faa protein --tools hmmer --stp_cpu 32 -c prokka_BML_MG/BML_MG.gff --cgc_substrate --out_dir BML_MG.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32
    run_dbcan prokka_SRF_MG/SRF_MG.faa protein --tools hmmer --stp_cpu 32 -c prokka_SRF_MG/SRF_MG.gff --cgc_substrate --out_dir SRF_MG.PUL.Sub --dia_cpu 32 --hmm_cpu 32 --tf_cpu 32 

Box 6: Example output folder content of run_dbcan substrate prediction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    In the output directory (https://bcb.unl.edu/dbCAN_tutorial/dataset3-Priest2023/BML_MG.dbCAN/), a total of 17 files and 1 folder are generated:

    .. code-block:: shell

        -rw-rw-r--  1 jinfang jinfang  9.6M Dec 28 10:18 PUL_blast.out
        -rw-rw-r--  1 jinfang jinfang  1.8M Dec 28 10:18 CGC.faa
        -rw-rw-r--  1 jinfang jinfang   26M Dec 28 10:18 cgc.gff
        -rw-rw-r--  1 jinfang jinfang  450K Dec 28 10:18 cgc.out
        -rw-rw-r--  1 jinfang jinfang  212K Dec 28 10:18 cgc_standard.out
        -rw-rw-r--  1 jinfang jinfang 1005K Dec 28 10:18 cgc_standard.out.json
        -rw-rw-r--  1 jinfang jinfang  406K Dec 28 10:11 dbcan-sub.hmm.out
        -rw-rw-r--  1 jinfang jinfang  325K Dec 28 10:11 diamond.out
        -rw-rw-r--  1 jinfang jinfang  332K Dec 28 10:11 dtemp.out
        -rw-rw-r--  1 jinfang jinfang  220K Dec 28 10:11 hmmer.out
        -rw-rw-r--  1 jinfang jinfang  240K Dec 28 10:18 overview.txt
        -rw-rw-r--  1 jinfang jinfang  1.7M Dec 28 10:17 stp.out
        -rw-rw-r--  1 jinfang jinfang   17K Dec 28 10:18 substrate.out
        drwxrwxr-x  2 jinfang jinfang   12K Dec 28 10:19 synteny.pdf
        -rw-rw-r--  1 jinfang jinfang  293K Dec 28 10:13 tf-1.out
        -rw-rw-r--  1 jinfang jinfang  222K Dec 28 10:15 tf-2.out
        -rw-rw-r--  1 jinfang jinfang  1.7M Dec 28 10:17 tp.out
        -rw-rw-r--  1 jinfang jinfang  105M Dec 28 05:57 uniInput


Descriptions of Output Files:

    - ``PUL_blast.out``: BLAST results between CGCs and PULs.
    - ``CGC.faa``: CGC Fasta sequences.
    - ``cgc.gff``: reformatted from the user input gff file by marking CAZymes, TFs, TCs, and STPs.
    - ``cgc.out``: raw output of CGC predictions.
Each entry in cgc.out includes:


    1.	CGC_id: CGC1
    2.	type: CAZyme
    3.	contig_id: contig_10157
    4.	gene_id: BML_MG_01992
    5.	start: 33003
    6.	end: 36077
    7.	strand: +
    8.	annotation: GH2

Explanation: the gene BML_MG_01992 encodes a GH2 CAZyme in the CGC1 of the contig contig_10157. CGC1 also has other genes, which are provided in other rows. BML_MG_01992 is on the positive strand of contig_10157 from 33003 to 36077. The type can be one of the four signature gene types (CAZymes, TCs, TFs, STPs) or the null type (not annotated as one of the four signature genes).

`cgc_standard.out.json`: JSON format of cgc_standard.out.
`dbcan-sub.hmm.out`: HMMER search result against dbCAN-sub HMMdb, including a column with CAZyme substrates extracted from fam-substrate-mapping-08012023.tsv. 
`diamond.out`: DIAMOND search result against the CAZy annotated protein sequences (CAZyDB.07262023.fa).
`dtemp.out`: temporary file.
`hmmer.out`: HMMER search result against dbCAN HMMdb.
`overview.txt`: summary of CAZyme annotation from three methods in TSV format. An example row has the following columns:
    1.	Gene_ID: BML_MG_01761
    2.	EC#: 2.4.99.-:5
    3.	dbCAN: GT112(19-370)
    4.	dbCAN_sub: GT112_e0
    5.	DIAMOND: GT112
    6.	#ofTools: 3
Explanation: the protein BML_MG_01761 is annotated by 3 tools to be a CAZyme: (1) GT112 (CAZy defined family GT112) by HMMER vs dbCAN HMMdb with a domain range from aa position 19 to 370, (2) GT112_e0 (eCAMI defined subfamily e0; e indicates it is from eCAMI not CAZy) by HMMER vs dbCAN-sub HMMdb (derived from eCAMI subfamilies), and (3) GT112 by DIAMOND vs CAZy annotated protein sequences. The second column 2.4.99.-:5 is extracted from eCAMI, meaning that the eCAMI subfamily GT112_e0 contains 5 member proteins which have an EC 2.4.99.- according to CAZy. In most cases, the 3 tools will have the same CAZyme family assignment. When they give different assignment. We recommend a preference order: dbCAN > eCAMI/dbCAN-sub > DIAMOND. See our dbCAN2 paper2, dbCAN3 paper3, and eCAMI4 for more details.
Note: If users invoked the --use_signalP parameter when running run_dbcan, there will be an additional column called signal in the overview.txt. 
stp.out: HMMER search result against the MiST5 compiled signal transduction protein HMMs from Pfam.
tf-1.out: HMMER search result against the DBD6 compiled transcription factor HMMs from Pfam 7.
tf-2.out: HMMER search result against the DBD compiled transcription factor HMMs from Superfamily 8.
tp.out: DIAMOND search result against the TCDB 9 annotated protein sequences.
substrate.out: summary of substrate prediction results for CGCs in TSV format from two approaches3 (dbCAN-PUL blast search and dbCAN-sub majority voting). An example row has the following columns:
    1.	CGC_ID: contig_10778|CGC2
    2.	Best hit PUL_ID in dbCAN-PUL: PUL0400 
    3.	Substrate of the hit PUL: alginate
    4.	Sum of bitscores for homologous gene pairs between CGC and PUL: 851.0
    5.	Types of homologous gene pairs: CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme 
    6.	Substrate predicted by majority voting of CAZymes in CGC: alginate
    7.	Voting score: 2.0
Explanation: The CGC2 of contig_10778 has its best hit PUL0400 (from PUL_blast.out) with alginate as substrate (from dbCAN-PUL_12-12-2023.xlsx). Four signature genes are matched between contig_10778|CGC2 and PUL0400 (from PUL_blast.out): all the four are CAZymes. The sum of blast bitscores of the 4 homologous pairs (CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme) is 851.0. Hence, the substrate of contig_10778|CGC2 is predicted to be alginate according to dbCAN-PUL blast search. The last two columns are based on the dbCAN-sub result (dbcan-sub.hmm.out), according to which two CAZymes in contig_10778|CGC2 are predicted to have alginate substrate. The voting score is thus 2.0, so that according to the majority voting rule, contig_10778|CGC2 is predicted to have an alginate substrate.
Note: for many CGCs, only one of the two approaches produces substrate prediction. In some cases, the two approaches produce different substrate assignments. We recommend a preference order: dbCAN-PUL blast search > dbCAN-sub majority voting. See our dbCAN3 paper3 for more details.
synteny.pdf: a folder with syntenic block alignment plots between all CGCs and PULs.
uniInput: renamed Fasta file from input protein sequence file.













