# run_dbcan 2.0

If you use webserver, please cite us:

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

If you want to use microbial sequence annotations directly, please cite us:

***Le Huang**, Han Zhang, Peizhi Wu, Sarah Entwistle, Xueqiong Li, Tanner Yohe, Haidong Yi, Zhenglu Yang, Yanbin Yin; 
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



## run_dbcan.py Stand Alone Version2.0 User Mannual

Rewritten by Huang Le in the Zhang Lab at NKU

Last updated 12/24/18

**update info**
- more user friendly
- added stp signature gene in CGC_Finder.py

### Function

- Accepts user input
- Predicts genes if needed
- Runs input against HMMER, DIAMOND, and Hotpep
- Optionally predicts CGCs with CGCFinder


### REQUIREMENTS

#### TOOLS
[DIAMOND](https://github.com/bbuchfink/diamond)

[HMMER](hmmer.org)

[hmmscan-parser](http://cys.bios.niu.edu/dbCAN2/download/Tools/)

[Hotpep-Python](http://cys.bios.niu.edu/dbCAN2/download/Tools/)

[signalp](http://www.cbs.dtu.dk/services/SignalP/)

[Prodigal](https://github.com/hyattpd/Prodigal)

[FragGeneScan](https://github.com/COL-IU/FragGeneScan)

[CGCFinder](http://cys.bios.niu.edu/dbCAN2/download/Tools/)

#### DATABASES [Link](http://cys.bios.niu.edu/dbCAN2/download/Databases)

	CAZY.dbCAN2_07202017.fa 		-- diamond makedb

	PPR					-- included in Hotpep

	dbCAN-HMMdb-V6.txt			-- use `hmmpress`

	tcdb.fa				        -- use `diamond makedb --in <inputFile> -d <dbName>`
	
	tf.fa					-- use `diamond makedb --in <inputFile> -d <dbName>`
	
#### PYTHON MODULE

	natsort                  -- use `pip install natsort`


### DIRECTORY STRUCTURE

- We recommend that all the databases are put in a seperate directory. This directory can be defined with the `--db_dir` option or leave it as directory with default name `db`. Otherwise, all databases must be in the same directory as run_dbcan.py.

- The Hotpep directory must also be in the same directory as run_dbcan.py. 

- When you use **FragGeneScan**, you should download this file, unpress, and add its path to user enviroenment variable. The concrete method is written on [this website(http://www.baidu/com).

- The "example" directory contains an example of a working directory setup.

### INPUT

```
python run_dbcan.py [inputFile] [inputType] [-c AuxillaryFile] [-t Tools] etc.
```

	[inputFile] - FASTA format file of either nucleotide or protein sequences
		
	[inputType] - protein=proteome, prok=prokaryote, meta=metagenome/mRNA/CDSs/short DNA seqs
		
	[--out_dir] - REQUIRED, user specifies an output directory.
		
	[-c AuxillaryFile]- optional, include to enable CGCFinder. If using a proteome input,
	the AuxillaryFile must be a GFF or BED format file containing gene positioning
	information. Otherwise, the AuxillaryFile may be left blank.
					
	[-t Tools] - optional, allows user to select a combination of tools to run. The options are any
					combination of 'diamond', 'hmmer', and 'hotpep'. The default value is 'all' which runs all three tools.
	[--dia_eval] - optional, allows user to set the DIAMOND E Value. Default = 1e-121.
		
	[--dia_cpu] - optional, allows user to set how many CPU cores DIAMOND can use. Default = 5.
		
	[--hmm_eval] - optional, allows user to set the HMMER E Value. Default = 1e-35.
		
	[--hmm_cov] - optional, allows user to set the HMMER Coverage value. Default = 0.35.
		
	[--hmm_cpu] - optional, allows user to set how many CPU cores HMMER can use. Default = 1.
		
	[--hot_hits] - optional, allows user to set the Hotpep Hits value. Default = 4.
		
	[--hot_freq] - optional, allows user to set the Hotpep Frequency value. Default = 2.0.
		
	[--hot_cpu] - optional, allows user to set how many CPU cores Hotpep can use. Default = 4.
		
	[--out_pre] - optional, allows user to set a prefix for all output files.
		
	[--db_dir] - optional, allows user to specify a database directory. Default = db/
		
	[--cgc_dis] - optional, allows user to specify CGCFinder Distance value. Allowed values are integers between 0-10. Default = 2.
		
	[--cgc_sig_genes] - optional, allows user to specify CGCFinder Signature Genes. The options are, 'tp': TP and CAZymes, 'tf': TF and CAZymes, and 'all': TF, TP, and CAZymes. Default = 'tp'.



### OUTPUT

Several files will be outputted. they are as follows:
	
	uniInput - The unified input file for the rest of the tools 
			(created by prodigal or FragGeneScan if a nucleotide sequence was used)
			
	Hotpep.out - the output from the Hotpep run
	
	diamond.out - the output from the diamond blast
		
	hmmer.out - the output from the hmmer run
		
	tf.out - the output from the diamond blast predicting TF's for CGCFinder
		
	tc.out - the output from the diamond blast predicting TC's for CGCFinder
		
	cgc.gff - GFF input file for CGCFinder
		
	cgc.out - ouput from the CGCFinder run

### EXAMPLE

An example setup is available in the example directory. Included in this directory are two FASTA sequences (one protein, one nucleotide).
	
To run this example type, run:

```
python run_dbcan.py EscheriaColiK12MG1655.fna prok
```
or 

```	
python run_dbcan.py EscheriaColiK12MG1655.faa protein 
```

While this example directory contains all the databases you will need (already formatted) and the Hotpep and FragGeneScan programs, you will still need to have the remaining programs installed on your machine (DIAMOND, HMMER, etc.).
	
To run the examples with CGCFinder turned on, run:
```	
python run_dbcan.py EscheriaColiK12MG1655.fna prok -c
```

or 

```
python run_dbcan.py EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff
```

Notice that the protein command has a GFF file following the -c option. A GFF or BED format file with gene position information is required to run CGCFinder when using a protein input.

If you have any questions, please feel free to contact with me on [Issue Dashboard](https://github.com/linnabrown/run_dbcan/issues) or Dr. Yin (yyin@niu.edu).
