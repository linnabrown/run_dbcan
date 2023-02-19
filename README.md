
run_dbcan4
========================

Status
----
[![dbcan](https://github.com/linnabrown/run_dbcan/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/linnabrown/run_dbcan/actions/workflows/ci.yml)
[![Package status](https://anaconda.org/bioconda/dbcan/badges/version.svg)](https://anaconda.org/bioconda/dbcan)
[![GitHub license](https://img.shields.io/badge/license-GUN3.0-blue.svg)](https://github.com/linnabrown/run_dbcan/blob/master/LICENSE)
[![Platform](https://anaconda.org/bioconda/dbcan/badges/platforms.svg)](https://anaconda.org/bioconda/dbcan/badges/platforms.svg)
[![GitHub downloads](https://anaconda.org/bioconda/dbcan/badges/downloads.svg)](https://anaconda.org/bioconda/dbcan)
[![GitHub versions](https://img.shields.io/pypi/pyversions/dbcan.svg)](https://anaconda.org/bioconda/dbcan)

<!-- [![Package version](https://img.shields.io/pypi/v/dbcan.svg)](https://pypi.org/project/dbcan/#files) -->

A standalone tool of [dbCAN3 web server](http://bcb.unl.edu/dbCAN2/).

Update Info
---
run_dbcan 4.0.0 is released.
1. CAZyme substrate prediction based on dbCAN-sub ; 

2. CGC substrate prediction based on dbCAN-PUL searching and [dbCAN-sub](https://bcb.unl.edu/dbCAN_sub/) majority voting. For CGC substrate prediction, please see our [dbCAN-seq update paper](https://academic.oup.com/nar/article/51/D1/D557/6833251?login=false) for details. With these new functions (esp. the dbCAN-sub search), run_dbcan4.0 is now slower to get the result back to you. Please be patience!


**Please update all of the databases**.

[Previous update information](https://github.com/linnabrown/run_dbcan/wiki/Update-information-Archive)

Function
----
- Accepts user input
- Predicts genes if needed
- Runs input against HMMER, DIAMOND, and dbCAN_sub
- Optionally predicts CGCs with CGCFinder

Support Platform
-----
Linux(Ubuntu, CentOS), MacOS

Installation via Bioconda
-----
1. Please install [Anoconda](https://www.anaconda.com) first.

2. Install [NCBI Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download).

3. Create virtual environment with dbcan and activate the virtual environment.

```
conda create -n run_dbcan python=3.8 dbcan -c conda-forge -c bioconda
conda activate run_dbcan
```

4. Database Installation.
```
test -d db || mkdir db
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/fam-substrate-mapping-08252022.tsv \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa && makeblastdb -in PUL.faa -dbtype prot \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.xlsx \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL_07-01-2022.txt \
	&& wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz && tar xvf dbCAN-PUL.tar.gz \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm && hmmpress dbCAN_sub.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/CAZyDB.08062022.fa && diamond makedb --in CAZyDB.08062022.fa -d CAZy \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt && mv dbCAN-HMMdb-V11.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/V11/tf-2.hmm && hmmpress tf-2.hmm \
    && wget https://bcb.unl.edu/dbCAN2/download/Databases/V11/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff

```
5. (Optional) SignalP Installation.
Our program include Signalp Petitide prediction with SignalP. Make sure to set `use_signalP=True` and *have to* obtain your own academic license of SignalP and download it from [here](https://services.healthtech.dtu.dk/service.php?SignalP-4.1), and then move the perl file from the tarball file (signalp-4.1g.Linux.tar.gz) into `/usr/bin/signalp` by yourself. Following statement is singalP-4.1 installation instruction.

Decompress signalp-4.1g.Linux.tar.gz than open the directory

```
tar -xvf signalp-4.1g.Linux.tar.gz && cd signalp-4.1
```

Then you can find those files/directories located in `signalp-4.1` directory
```
(base) lehuang@lehuang:~/Downloads/signalp-4.1$ ls
bin  lib  signalp  signalp.1  signalp-4.1.readme  syn  test
```

*signalp* is the perl file that you will use in your program
Edit the paragraph labeled  "GENERAL SETTINGS, CUSTOMIZE ..." in the top of
   the file 'signalp'. The following twovmandatory variables need to be set:

   	**SIGNALP**		full path to the signalp-4.1 directory on your system
	**outputDir**	where to store temporary files (writable to all users)
	**MAX_ALLOWED_ENTRIES** the number of input sequences allowed per run.
```
Here is the example for me to change line 13, line 17 and line 20 in `singalp` file. I suggest you to set MAX_ALLOWED_ENTRIES as 100000
###############################################################################
#               GENERAL SETTINGS: CUSTOMIZE TO YOUR SITE
###############################################################################

# full path to the signalp-4.1 directory on your system (mandatory)
BEGIN {
    $ENV{SIGNALP} = '/home/lehuang/Downloads/signalp-4.1';
}

# determine where to store temporary files (must be writable to all users)
my $outputDir = "/home/lehuang/Downloads/signalp-4.1/output";

# max number of sequences per run (any number can be handled)
my $MAX_ALLOWED_ENTRIES=100000;
```

And then, use this command:

```
sudo cp signalp /usr/bin/signalp
sudo chmod 755 /usr/bin/signalp
```
If you don't have the permission to access `/usr/bin`, you can use the parameter `-sp` or `--signalP_path` to indicate your `signalp` file path in the run_dbcan program. Please see the step 6.
6. Check Program.
```
run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655
```

If you want to run the code with SignalP
```
run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655 --use_signalP=TRUE

```
If you don't have the permission to access `/usr/bin` when running with signalP, you can use the parameter `-sp` or `--signalP_path` to indicate your `signalp` file path in the run_dbcan program. 

```
run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655 --use_signalP=TRUE -sp /home/lehuang/Downloads/signalp-4.1/signalp 
```
Installation via Docker
----
1. Make sure docker is installed on your computer successfully.
2. Docker pull image
```
docker pull haidyi/run_dbcan:latest
```
3. Run. Mount `input sequence file` and `output directory` to the container.
```
docker run --name <preferred_name> -v <host-path>:<container-path> -it haidyi/run_dbcan:latest run_dbcan <input_file> [params] --out_dir <output_dir>
```



REQUIREMENTS
----

**TOOLS**

----
P.S.: You do not need to download `CGCFinder` and `hmmscan-parser` because they are included in run_dbcan V4. If you use python package or docker, you don't need to download Prodigal because they includes these denpendencies. Otherwise we recommend you to install and copy them into `/usr/bin` as system application or add their path into system envrionmental profile.


[Python3]--Be sure to use python3, not python2

[DIAMOND](https://github.com/bbuchfink/diamond)-- Included in run_dbcan4.

[HMMER](hmmer.org)--Included in run_dbcan4.

[hmmscan-parser](https://github.com/linnabrown/run_dbcan/blob/master/hmmscan-parser.py)--This is included in run_dbcan4.

dbCAN_sub--Included in run_dbcan4.

[signalp](http://www.cbs.dtu.dk/services/SignalP/)--please download and install if you need.

[Prodigal](https://github.com/hyattpd/Prodigal)--Included in run_dbcan4.

[CGCFinder](https://github.com/linnabrown/run_dbcan/blob/master/dbcan/utils/CGCFinder.py)--Included in run_dbcan4.

**DATABASES Installation (those are included in step4 Database Installation)**

----
[Databse](http://bcb.unl.edu/dbCAN2/download/Databases) -- Database Folder

[CAZy.fa](http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.09242021.fa)--use `diamond makedb --in CAZyDB.09242021.fa -d CAZy`

[dbCAN_sub](http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN_sub.hmm) --use `hmmpress dbCAN_sub.hmm` .

[dbCAN-PUL](http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-PUL.tar.gz) The substrates files from dbCAN-PUL.

[PUL](http://bcb.unl.edu/dbCAN2/download/Databases/PUL.faa)--The PUL sequences, use `makeblastdb -in PUL.faa -dbtype prot`.
[dbCAN-HMMdb-V11.txt](http://bcb.unl.edu/dbCAN2/download/Databases/V11/dbCAN-HMMdb-V11.txt)--First use `mv dbCAN-HMMdb-V11.txt dbCAN.txt`, then use `hmmpress dbCAN.txt`

[tcdb.fa](http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa)--use `diamond makedb --in tcdb.fa -d tcdb`

[tf-1.hmm](http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm)--use `hmmpress tf-1.hmm`

[tf-2.hmm](http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm)--use `hmmpress tf-2.hmm`

[stp.hmm](http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm)--use `hmmpress stp.hmm`


Params
----
```
Required arguments:
  inputFile             User input file. Must be in FASTA format.
  {protein,prok,meta}   Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome

optional arguments:
  -h, --help            show this help message and exit
  --dbCANFile DBCANFILE
                        Indicate the file name of HMM database such as dbCAN.txt, please use the newest one from dbCAN2 website.
  --dia_eval DIA_EVAL   DIAMOND E Value
  --dia_cpu DIA_CPU     Number of CPU cores that DIAMOND is allowed to use
  --hmm_eval HMM_EVAL   HMMER E Value
  --hmm_cov HMM_COV     HMMER Coverage val
  --hmm_cpu HMM_CPU     Number of CPU cores that HMMER is allowed to use
  --out_pre OUT_PRE     Output files prefix
  --out_dir OUT_DIR     Output directory
  --db_dir DB_DIR       Database directory
  --tools {hmmer,diamond,dbcansub,all} [{hmmer,diamond,dbcansub,all} ...], -t {hmmer,diamond,dbcansub,all} [{hmmer,diamond,dbcansub,all} ...]
                        Choose a combination of tools to run
  --use_signalP USE_SIGNALP
                        Use signalP or not, remember, you need to setup signalP tool first. Because of signalP license, Docker version does not have signalP.
  --signalP_path SIGNALP_PATH, -sp SIGNALP_PATH
                        The path for signalp. Default location is signalp
  --gram {p,n,all}, -g {p,n,all}
                        Choose gram+(p) or gram-(n) for proteome/prokaryote nucleotide, which are params of SingalP, only if user use singalP
  -v VERSION, --version VERSION

dbCAN-sub parameters:
  --dbcan_thread DBCAN_THREAD, -dt DBCAN_THREAD
  --tf_eval TF_EVAL     tf.hmm HMMER E Value
  --tf_cov TF_COV       tf.hmm HMMER Coverage val
  --tf_cpu TF_CPU       tf.hmm Number of CPU cores that HMMER is allowed to use
  --stp_eval STP_EVAL   stp.hmm HMMER E Value
  --stp_cov STP_COV     stp.hmm HMMER Coverage val
  --stp_cpu STP_CPU     stp.hmm Number of CPU cores that HMMER is allowed to use

CGC_Finder parameters:
  --cluster CLUSTER, -c CLUSTER
                        Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used
  --cgc_dis CGC_DIS     CGCFinder Distance value
  --cgc_sig_genes {tf,tp,stp,tp+tf,tp+stp,tf+stp,all}
                        CGCFinder Signature Genes value

CGC_Substrate parameters:
  --cgc_substrate       run cgc substrate prediction?
  --pul PUL             dbCAN-PUL PUL.faa
  -o OUT, --out OUT
  -w WORKDIR, --workdir WORKDIR
  -env ENV, --env ENV
  -oecami, --oecami     out eCAMI prediction intermediate result?
  -odbcanpul, --odbcanpul
                        output dbCAN-PUL prediction intermediate result?

dbCAN-PUL homologous searching parameters:
  how to define homologous gene hits and PUL hits

  -upghn UNIQ_PUL_GENE_HIT_NUM, --uniq_pul_gene_hit_num UNIQ_PUL_GENE_HIT_NUM
  -uqcgn UNIQ_QUERY_CGC_GENE_NUM, --uniq_query_cgc_gene_num UNIQ_QUERY_CGC_GENE_NUM
  -cpn CAZYME_PAIR_NUM, --CAZyme_pair_num CAZYME_PAIR_NUM
  -tpn TOTAL_PAIR_NUM, --total_pair_num TOTAL_PAIR_NUM
  -ept EXTRA_PAIR_TYPE, --extra_pair_type EXTRA_PAIR_TYPE
                        None[TC-TC,STP-STP]. Some like sigunature hits
  -eptn EXTRA_PAIR_TYPE_NUM, --extra_pair_type_num EXTRA_PAIR_TYPE_NUM
                        specify signature pair cutoff.1,2
  -iden IDENTITY_CUTOFF, --identity_cutoff IDENTITY_CUTOFF
                        identity to identify a homologous hit
  -cov COVERAGE_CUTOFF, --coverage_cutoff COVERAGE_CUTOFF
                        query coverage cutoff to identify a homologous hit
  -bsc BITSCORE_CUTOFF, --bitscore_cutoff BITSCORE_CUTOFF
                        bitscore cutoff to identify a homologous hit
  -evalue EVALUE_CUTOFF, --evalue_cutoff EVALUE_CUTOFF
                        evalue cutoff to identify a homologous hit

dbCAN-sub major voting parameters:
  how to define dbsub hits and dbCAN-sub subfamily substrate

  -hmmcov HMMCOV, --hmmcov HMMCOV
  -hmmevalue HMMEVALUE, --hmmevalue HMMEVALUE
  -ndsc NUM_OF_DOMAINS_SUBSTRATE_CUTOFF, --num_of_domains_substrate_cutoff NUM_OF_DOMAINS_SUBSTRATE_CUTOFF
                        define how many domains share substrates in a CGC, one protein may include several subfamily domains.
  -npsc NUM_OF_PROTEIN_SUBSTRATE_CUTOFF, --num_of_protein_substrate_cutoff NUM_OF_PROTEIN_SUBSTRATE_CUTOFF
                        define how many sequences share substrates in a CGC, one protein may include several subfamily domains.
  -subs SUBSTRATE_SCORS, --substrate_scors SUBSTRATE_SCORS
                        each cgc contains with substrate must more than this value
```

RUN & OUTPUT
----
Use following command to run the program.
```
run_dbcan [inputFile] [inputType] [-c AuxillaryFile] [-t Tools] etc.
```

Several files will be produced via `run_dbcan`. They are as follows:

	uniInput - The unified input file for the rest of the tools
			(created by prodigal if a nucleotide sequence was used)

	dbsub.out - the output from the dbCAN_sub run

	diamond.out - the output from the diamond blast

	hmmer.out - the output from the hmmer run

	tf.out - the output from the diamond blast predicting TF's for CGCFinder

	tc.out - the output from the diamond blast predicting TC's for CGCFinder

	cgc.gff - GFF input file for CGCFinder

	cgc.out - ouput from the CGCFinder run

	overview.txt - Details the CAZyme predictions across the three tools with signalp results

EXAMPLE
----

An example setup is available in the example directory. Included in this directory are two FASTA sequences (one protein, one nucleotide).

To run this example type, run:

```
run_dbcan EscheriaColiK12MG1655.fna prok --out_dir output_EscheriaColiK12MG1655
```

or

```
run_dbcan EscheriaColiK12MG1655.faa protein --out_dir output_EscheriaColiK12MG1655
```

To run the examples with CGCFinder turned on, run:
```
run_dbcan EscheriaColiK12MG1655.fna prok -c cluster --out_dir output_EscheriaColiK12MG1655
```

or

```
run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --out_dir output_EscheriaColiK12MG1655
```

Notice that the protein command has a GFF file following the -c option. A GFF or BED format file with gene position information is required to run CGCFinder when using a protein input.

If you have any questions, please feel free to contact with Dr. Yin (yanbin.yin@gmail.com or yyin@unl.edu) or me (Le Huang) on [Issue Dashboard](https://github.com/linnabrown/run_dbcan/issues).


Reference
----

This is the standalone version of dbCAN annotation tool for automated CAZyme annotation (known as run_dbCAN), written by Le Huang and Tanner Yohe.

If you want to use our dbCAN3 webserver, please go to http://bcb.unl.edu/dbCAN2/.

If you use dbCAN standalone tool (run_dbcan) or/and our web server for publication, please cite us:

*Han Zhang, Tanner Yohe, **Le Huang**, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin;
dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research,
Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418*
```
@article{doi:10.1093/nar/gky418,
author = {Zhang, Han and Yohe, Tanner and Huang, Le and Entwistle, Sarah and Wu, Peizhi and Yang, Zhenglu and Busk, Peter K and Xu, Ying and Yin, Yanbin},
title = {dbCAN2: a meta server for automated carbohydrate-active enzyme annotation},
journal = {Nucleic Acids Research},
volume = {46},
number = {W1},
pages = {W95-W101},
year = {2018},
doi = {10.1093/nar/gky418},
URL = {http://dx.doi.org/10.1093/nar/gky418},
eprint = {/oup/backfile/content_public/journal/nar/46/w1/10.1093_nar_gky418/1/gky418.pdf}
}
```

If you want to use pre-computed bacterial CAZyme sequences/annotations directly, please go to http://bcb.unl.edu/dbCAN_seq/ and cite us:

**Le Huang**, Han Zhang, Peizhi Wu, Sarah Entwistle, Xueqiong Li, Tanner Yohe, Haidong Yi, Zhenglu Yang, Yanbin Yin;
dbCAN-seq: a database of carbohydrate-active enzyme (CAZyme) sequence and annotation, Nucleic Acids Research,
Volume 46, Issue D1, 4 January 2018, Pages D516–D521, https://doi.org/10.1093/nar/gkx894*
```
@article{doi:10.1093/nar/gkx894,
author = {Huang, Le and Zhang, Han and Wu, Peizhi and Entwistle, Sarah and Li, Xueqiong and Yohe, Tanner and Yi, Haidong and Yang, Zhenglu and Yin, Yanbin},
title = {dbCAN-seq: a database of carbohydrate-active enzyme (CAZyme) sequence and annotation},
journal = {Nucleic Acids Research},
volume = {46},
number = {D1},
pages = {D516-D521},
year = {2018},
doi = {10.1093/nar/gkx894},
URL = {http://dx.doi.org/10.1093/nar/gkx894},
eprint = {/oup/backfile/content_public/journal/nar/46/d1/10.1093_nar_gkx894/2/gkx894.pdf}
}
```


## Commit History

[![Commit History Chart](https://commit-history-api.herokuapp.com/svg?repos=linnabrown/run_dbcan&type=Date)](https://the-commit-history.vercel.app/#linnabrown/run_dbcan&Date)


