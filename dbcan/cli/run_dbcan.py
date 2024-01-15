#!/usr/bin/env python3
#########################################################
# dbCAN3 (Stand Alone Version)
#
# Written by Tanner Yohe in the Yin Lab at NIU
# Revised by Qiwei Ge in Yin Lab at UNL && Le Huang at NKU
# Updated by Le Huang at UNC, Mohamad Majd Raslan in the Yin Lab at NIU, Wei Li, Qiwei Ge in Dr.Yin's Lab at UNL
# Updated by Haidong Yi for reconstructing codes,  Alex Fraser for adding functions.
# Updated by Jinfang Zheng in Yinlab at UNL, new function, substrate prediciton based on dbCAN-PUL and dbCAN-sub database.

import argparse
import os
import time

# Recent updated information:
#   Jan/01/23: Add doc code [Haidong Yi, Le Huang]
#   Oct/10/23: Recontructed the run_dbcan [Haidong Yi]
#   Sep/07/23: Replace hmmscan with hmmsearch. Update perl code [Le Huang, Yanbin Yin]
#   Dec/15/22: 1.adding function to convert cgc_standard.out to json format. 2. adding function cgc_[Jinfang Zheng]
#   Dec/06/22: fix gene ID in CGCfinder output file cgc.out[Jinfang Zheng]
#   Nov/06/22: Using dbCAN_sub, eCAMI has been removed [Qiwei Ge]
#   Jun/13/22: Allowing direct calls to main function from other scripts [Alex Fraser]
#   Sep/29/22: Hotpep has been removed, added eCAMI tool. 2. cgc out reformatting. 3. Fixed multiple GT2s [Qiwei Ge]
#
# Accepts user input
# Predicts genes if needed
# Runs input against HMMER, DIAMOND, and dbCAN_sub
# Optionally predicts CGCs with CGCFinder
# Creats an overview table using output files from core
# tools from dbcan-sub.hmm.out,hmmer.out and diamond.out if they exist.
##########################################################
from subprocess import Popen, call

from dbcan.cli import hmmer_parser
from dbcan.utils.cgc_substrate_prediction import cgc_substrate_prediction
from dbcan.utils.CGCFinder import cgc_finder, simplify_cgc_output



def runHmmer(outPath, hmm_cpu, dbDir, hmm_eval, hmm_cov, db_name):
    """
    Executes HMMER's hmmsearch tool to search sequence databases for sequence homologs.

    This function runs the hmmsearch command from the HMMER suite, directing the output to specified files. It then parses the hmmsearch output, writes the parsed output to a file, and cleans up intermediate files.

    Parameters
    ----------
        outPath (str): The path where output files will be saved.
        hmm_cpu (int): The number of CPUs to use for hmmsearch.
        dbDir (str): The directory containing the HMM database files.
        hmm_eval (float): The E-value threshold for considering hits as significant in hmmsearch.
        hmm_cov (float): The coverage threshold for considering hits as significant.
        db_name (str): The name of the database to use for hmmsearch.

    Side Effects:
        - Creates and writes to files in the specified output path.
        - Deletes intermediate output files after processing is complete.
    """
    domtblout_file = f"{outPath}h{db_name}.out"
    hmm_file = f"{dbDir}{db_name}.hmm"
    uniInput_file = f"{outPath}uniInput"

    hmmer = Popen(
        [
            "hmmsearch",
            "--domtblout",
            domtblout_file,
            "--cpu",
            str(hmm_cpu),
            "-o",
            "/dev/null",
            hmm_file,
            uniInput_file,
        ]
    )
    hmmer.wait()
    parsed_hmm_output = hmmer_parser.run(input_file=f"{outPath}h{db_name}.out", eval_num=hmm_eval, coverage=hmm_cov)
    with open(f"{outPath}{db_name}.out", "w") as f:
        f.write(parsed_hmm_output)
    if os.path.exists(domtblout_file):
        call(["rm", domtblout_file])


def split_uniInput(uniInput, dbcan_thread, outPath, dbDir, hmm_eval, hmm_cov, hmm_cpu, dbcan_offset):
    """
    Splits a large input file into smaller files and runs HMMER hmmsearch on each part.

    This function is designed to handle large genomic data files. It splits the input file into smaller parts based on the specified offset and then runs hmmsearch on each part in parallel. The function manages file creation and deletion, ensuring that temporary files are removed after processing.

    Parameters
    ----------
        uniInput (str): Path to the input file containing genomic data.
        dbcan_thread (int): Number of threads to use for parallel processing.
        outPath (str): Path where output files and temporary files will be stored.
        dbDir (str): Directory containing the HMM database files.
        hmm_eval (float): E-value threshold for considering hits as significant in hmmsearch.
        hmm_cov (float): Coverage threshold for considering hits as significant.
        dbcan_offset (float): Factor to determine the size of each split file based on the input file size.

    Side Effects:
        - Creates multiple temporary files during processing.
        - Prints the count of identifiers and processing time to the standard output.
        - Deletes all temporary files created during the process.

    Notes
    -----
        - The function assumes a specific format for the input data.
        - Time taken for execution is printed at the end of the function's run.
    """
    ticks = time.time()
    
    dbsub = Popen(
        [
            "hmmsearch",
            "--domtblout",
            f"{outPath}d.txt",
            "--cpu",
            str(hmm_cpu),
            "-o",
            "/dev/null",
            f"{dbDir}dbCAN_sub.hmm",
            f"{outPath}uniInput",
        ]
    )
    dbsub.wait()

    hmm_parser_output = hmmer_parser.run(f"{outPath}d.txt", eval_num=hmm_eval, coverage=hmm_cov)
    with open(f"{outPath}dtemp.out", "w") as temp_hmmer_file:
        temp_hmmer_file.write(hmm_parser_output)

    print("total time:", time.time() - ticks)


def run_dbCAN(
    inputFile,
    inputType,
    cluster=None,
    dbCANFile="dbCAN.txt",
    dia_eval=1e-102,
    dia_cpu=8,
    hmm_eval=1e-15,
    hmm_cov=0.35,
    hmm_cpu=8,
    dbcan_thread=5,
    dbcan_offset=2,
    tf_eval=1e-4,
    tf_cov=0.35,
    tf_cpu=1,
    stp_eval=1e-4,
    stp_cov=0.3,
    stp_cpu=1,
    prefix="",
    outDir="output",
    dbDir="db",
    cgc_dis=2,
    cgc_sig_genes="tp",
    tool_arg="all",
    use_signalP=False,
    signalP_path="signalp",
    gram="all",
):
    """
    Executes the dbCAN tool for analyzing genomic sequences.

    This function runs the dbCAN tool, which involves several steps: setting up paths and directories, gene prediction, running SignalP, core analysis tools (DIAMOND, HMMER, dbCAN_sub), and parsing the results. The function also combines results from different analysis tools and prepares a summary overview.

    Parameters
    ----------
        inputFile (str): Path to the input file in FASTA format.
        inputType (str): Type of sequence input (protein, prok, meta).
        cluster (str, optional): Path to the cluster file if available.
        dbCANFile (str): Filename of the HMM database.
        dia_eval (float): E-value threshold for DIAMOND.
        dia_cpu, hmm_cpu, dbcan_thread, tf_cpu, stp_cpu (int): CPU cores for respective tools.
        hmm_eval (float): E-value threshold for HMMER.
        hmm_cov, tf_cov, stp_cov (float): Coverage values for HMMER and related tools.
        dbcan_offset (int): Offset parameter for dbCAN.
        prefix (str): Prefix for output files.
        outDir (str): Directory for output files.
        dbDir (str): Directory containing the databases.
        cgc_dis (int): Distance parameter for CGC Finder.
        cgc_sig_genes (str): Signature genes for CGC Finder.
        tool_arg (str): Tools to run ('all', 'diamond', 'hmmer', etc.).
        use_signalP (bool): Whether to use SignalP tool.
        signalP_path (str): Path to the SignalP tool.
        gram (str): Gram type for proteome/prokaryote nucleotide (p, n, all).

    """
    # Begin Setup and Input Checks
    if not dbDir.endswith("/") and len(dbDir) > 0:
        dbDir += "/"

    if not outDir.endswith("/") and len(outDir) > 0:
        outDir += "/"

    outPath = outDir + prefix
    auxFile = ""

    find_clusters = False
    if cluster is not None:
        find_clusters = True
        if inputType == "protein":
            auxFile = cluster
        else:
            auxFile = "%sprodigal.gff" % outPath

    if not os.path.isdir(dbDir):
        print(dbDir, "ERROR: The database directory does not exist")
        exit()

    if not os.path.isfile(os.path.join(dbDir, "CAZy.dmnd")):
        print(
            "ERROR: No CAZy DIAMOND database found. \
        Please make sure that your CAZy DIAMOND databased is named 'CAZy.dmnd' and is located in your database directory"
        )
        exit()

    if not os.path.isfile(os.path.join(dbDir, dbCANFile)):
        print(
            "ERROR: No dbCAN HMM database found. \
        Please make sure that your dbCAN HMM database is named 'dbCAN-HMMdb-V11.txt' or the newest one, has been through hmmpress, and is located in your database directory"
        )
        exit()

    if not os.path.isfile(os.path.join(dbDir, "dbCAN_sub.hmm")):
        print(
            "ERROR: No dbCAN_sub HMM database found. \
        Please make sure that your dbCAN_sub HMM databased is named 'dbCAN_sub.hmm' or has been through hmmpress, and is located in your database directory"
        )
        exit()

    if not os.path.isdir(outDir):
        call(["mkdir", outDir])

    if find_clusters and inputType == "protein":
        if len(auxFile) > 0:
            print(auxFile)
            if not os.path.isfile(auxFile):
                print("ERROR: It seems that the auxillary filename that you provided does not exist, or is not a file")
                exit()
        else:
            print(
                "ERROR: Please provide an auxillary input file with the position of each gene. This file can either be in BED or GFF format"
            )
            exit()

    tools = [True, True, True]  # DIAMOND, HMMER, dbCAN_sub
    if "all" not in tool_arg:
        if "diamond" not in tool_arg:
            tools[0] = False
        if "hmmer" not in tool_arg:
            tools[1] = False
        if "dbcansub" not in tool_arg:
            tools[2] = False

    # End Setup and Input Checks
    #########################
    #########################
    # Begin Gene Prediction Tools
    if inputType == "prok":
        call(
            [
                "prodigal",
                "-i",
                inputFile,
                "-a",
                "%suniInput" % outPath,
                "-o",
                "%sprodigal.gff" % outPath,
                "-f",
                "gff",
                "-q",
            ]
        )
    if inputType == "meta":
        call(
            [
                "prodigal",
                "-i",
                inputFile,
                "-a",
                "%suniInput" % outPath,
                "-o",
                "%sprodigal.gff" % outPath,
                "-f",
                "gff",
                "-p",
                "meta",
                "-q",
            ]
        )
    # Proteome
    if inputType == "protein":
        call(["cp", inputFile, "%suniInput" % outPath])

    # End Gene Prediction Tools
    #######################
    # Begin SignalP
    if use_signalP:
        print("\n\n***************************0. SIGNALP start*************************************************\n\n")
        if gram == "p" or gram == "all":
            signalpos = Popen(f"{signalP_path} -t gram+ {outPath}uniInput > {outPath}signalp.pos", shell=True)
        if gram == "n" or gram == "all":
            signalpneg = Popen(f"{signalP_path} -t gram- {outPath}uniInput > {outPath}signalp.neg", shell=True)
        if gram == "euk" or gram == "all":
            signalpeuk = Popen(f"{signalP_path} -t euk {outPath}uniInput > {outPath}signalp.euk", shell=True)

    # End SignalP
    #######################
    # Begin Core Tools

    if tools[0]:  ### run diamond
        print("\n\n***************************1. DIAMOND start*************************************************\n\n")
        os.system(
            "diamond blastp -d %s -e %s -q %suniInput -k 1 -p %d -o %sdiamond.out -f 6"
            % (os.path.join(dbDir, "CAZy"), str(dia_eval), outPath, dia_cpu, outPath)
        )
        print("\n\n***************************1. DIAMOND end***************************************************\n\n")

    if tools[1]:  ### run hmmsearch (hmmer)
        print("\n\n***************************2. HMMER start*************************************************\n\n")
        os.system(
            f"hmmsearch --domtblout {outPath}h.out --cpu {hmm_cpu} -o /dev/null {os.path.join(dbDir, dbCANFile)} {outPath}uniInput "
        )
        print("\n\n***************************2. HMMER end***************************************************\n\n")

        hmm_parser_output = hmmer_parser.run(f"{outPath}h.out", eval_num=hmm_eval, coverage=hmm_cov)
        with open(f"{outPath}hmmer.out", "w") as hmmer_file:
            hmmer_file.write(hmm_parser_output)
        # could clean this up and manipulate hmm_parser_output data directly instead of passing it into a temp file
        with open(f"{outPath}hmmer.out", "r+") as f:
            text = f.read()
            f.close()
            call(["rm", f"{outPath}hmmer.out"])
            text = text.split("\n")
            if "" in text:
                text.remove("")
            for i in range(len(text)):
                if "GT2_" in text[i]:
                    profile = text[i].split("\t")[0].split(".")[0]
                    text[i] = text[i].replace(profile, "GT2")
                with open(f"{outPath}hmmer.out", "a") as f:
                    f.write(text[i] + "\n")
                    f.close()
        if os.path.exists(f"{outPath}h.out"):
            call(["rm", f"{outPath}h.out"])

    if tools[2]:
        print(
            "\n\n***************************3. dbCAN_sub start***************************************************\n\n"
        )
        split_uniInput("%suniInput" % outPath, dbcan_thread, outPath, dbDir, hmm_eval, hmm_cov, hmm_cpu, dbcan_offset)
        print("\n\n***************************3. dbCAN_sub end***************************************************\n\n")

        # Process dtemp.out and create dbcan-sub.hmm.out
        with open(f"{outPath}dtemp.out") as f:
            with open("%sdbcan-sub.hmm.out" % outPath, "w") as out:
                processed_lines = []
                for line in f:
                    row = line.rstrip().split("\t")
                    row.append(str(float(int(row[6]) - int(row[5])) / int(row[1])))
                    if float(row[4]) <= 1e-15 and float(row[-1]) >= 0.35:
                        processed_lines.append("\t".join(row))
        
        # Process dbcan-sub.hmm.out content
        updated_lines = [line for line in processed_lines if line.strip()]
                
        for i in range(len(updated_lines)):
            if "GT2_" in updated_lines[i]:
                profile = updated_lines[i].split("\t")[0].split(".")[0]
                updated_lines[i] = updated_lines[i].replace(profile, "GT2")
        if updated_lines:
            with open(f"{outPath}dbcan-sub.hmm.out", "w") as f2:
                f2.write("\n".join(updated_lines))

        
    # End Core Tools
    ########################
    # Begin Parse Results

    # parse dbCAN_sub result
    if tools[2]:
        subs_dict = {}
        with open(f"{dbDir}fam-substrate-mapping.tsv") as f:
            next(f)
            for line in f:
                r = line.split("\t")
                key = (r[2], "-") if len(r[4]) == 1 else (r[2], r[4].strip())
                subs_dict[key] = r[0]
        # Process dbcan-sub.hmm.out if it exists
        if os.path.exists(f"{outPath}dbcan-sub.hmm.out"):
            if os.path.exists(f"{outPath}temp"):
                os.system(f"rm {outPath}temp")
            with open(f"{outPath}dbcan-sub.hmm.out") as f, open(f"{outPath}temp", "w") as out:
                    out.write(
                        "dbCAN subfam\tSubfam Composition\tSubfam EC\tSubstrate\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n"
                    )

                    for line in f:
                        profile = line.split("\t")
                        subfam = []
                        sub_composition = []
                        sub_ec = []
                        newline = []
                        substrate = []
                        key1 = "-"
                        key2 = ["-"]

                        for p in profile[0].split("|"):
                            if ".hmm" in p:
                                subfam.append(p.split(".")[0])
                                key1 = p.split(".")[0].split("_")[0]
                            elif len(p.split(".")) == 4:
                                sub_ec.append(p)
                                key2.append(p.split(":")[0])
                            else:
                                sub_composition.append(p)

                        for i in range(len(key2)):
                            try:
                                substrate.append(subs_dict[key1, key2[i]])
                            except KeyError as e:
                                print("No substrate for it", e)

                        subfam = "|".join(subfam)

                        if sub_composition:
                            sub_composition = "|".join(sub_composition)
                        else:
                            sub_composition = "-"

                        if sub_ec:
                            sub_ec = "|".join(sub_ec)
                        else:
                            sub_ec = "-"

                        if substrate:
                            substrate = ", ".join(substrate)
                        else:
                            substrate = "-"

                        rest = "\t".join(profile[1:])

                        newline = subfam + "\t" + sub_composition + "\t" + sub_ec + "\t" + substrate + "\t" + rest
                        out.write(newline)
            call(["mv", f"{outPath}temp", f"{outPath}dbcan-sub.hmm.out"])
        else:
            print(f"File not found: {outPath}dbcan-sub.hmm.out")


    # parse hmmer result
    if tools[1]:
        try:
            with open(outDir + prefix + "hmmer.out") as f:
                with open(outDir + prefix + "temp", "w") as out:
                    out.write(
                        "HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n"
                    )
                    for line in f:
                        out.write(line)
            call(["mv", outDir + prefix + "temp", outDir + prefix + "hmmer.out"])
        except FileNotFoundError as e:
            print(f"File not found: {e}")
            # Optionally, create the file if it doesn't exist
            with open(outDir + prefix + "temp", "w") as out:
                out.write(
                    "HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n"
                )
            call(["mv", outDir + prefix + "temp", outDir + prefix + "hmmer.out"])
        except OSError as e:
            print(f"IO error occurred: {e}")

    # parse diamond result
    if tools[0]:
        with open(outDir + prefix + "diamond.out") as f:
            with open(outDir + prefix + "temp", "w") as out:
                out.write(
                    "Gene ID\tCAZy ID\t% Identical\tLength\tMismatches\tGap Open\tGene Start\tGene End\tCAZy Start\tCAZy End\tE Value\tBit Score\n"
                )
                for line in f:
                    out.write(line)
        call(["mv", outDir + prefix + "temp", outDir + prefix + "diamond.out"])

    # End Parse Results
    ########################
    # Begin CGCFinder

    if find_clusters:  ### run cgc_finder or not
        print("*****************************CGC-Finder start************************************")

        ########################
        # Begin TF,TP, STP prediction
        """
        tf hmmer
        """
        runHmmer(outPath, str(tf_cpu), dbDir, str(tf_eval), str(tf_cov), "tf-1")
        runHmmer(outPath, str(tf_cpu), dbDir, str(tf_eval), str(tf_cov), "tf-2")
        """
        stp hmmer
        """
        runHmmer(outPath, str(stp_cpu), dbDir, str(stp_eval), str(stp_cov), "stp")

        """
        tp diamond
        """
        call(
            [
                "diamond",
                "blastp",
                "-d",
                dbDir + "tcdb.dmnd",
                "-e",
                "1e-10",
                "-q",
                "%suniInput" % outPath,
                "-k",
                "1",
                "-p",
                "1",
                "-o",
                outPath + "tp.out",
                "-f",
                "6",
            ]
        )

        # Initialize set and dictionaries for TF,  STP, TP,
        tf, stp, tp = set(), set(), set()
        tf_genes, stp_genes, tp_genes = {}, {}, {}
        
        with open("%stf-1.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split("\t")
                tf.add(row[2])
                row[0] = "DBD-Pfam|" + row[0]
                if row[2] not in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += "," + row[0]

        with open("%stf-2.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split("\t")
                tf.add(row[2])
                row[0] = "DBD-SUPERFAMILY|" + row[0]
                if row[2] not in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += "," + row[0]

        with open("%sstp.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split("\t")
                stp.add(row[2])
                row[0] = "STP|" + row[0]
                if row[2] not in stp_genes:
                    stp_genes[row[2]] = row[0]
                else:
                    stp_genes[row[2]] += "," + row[0]

        with open(outDir + prefix + "tp.out") as f:
            for line in f:
                row = line.rstrip().split("\t")
                tp.add(row[0])
                if row[0] not in tp_genes:
                    tp_genes[row[0]] = row[1]
                else:
                    tp_genes[row[0]] += "," + row[1]

        
        # End TF and TP prediction
        ##########################
        # Begine CAZyme Extraction
        
        #define dictionary for cazyme_genes: key -> gene; val -> fams
        cazyme_genes = {}

        dia = set()
        hmm = set()
        dbs = set()

        if tools[0]:  # Deal with diamond result
            with open(outDir + prefix + "diamond.out") as f:
                next(f)  # Skip the first line
                for line in f:
                    row = line.rstrip().split("\t")
                    gene_key = row[0]
                    fam_vals = row[1].strip("|").split("|")[1:]
                    dia.add(gene_key)
                    cazyme_genes.setdefault(gene_key, set()).update(fam_vals)

        if tools[1]:  # Deal with hmmer result
            with open(outDir + prefix + "hmmer.out") as f:
                next(f)  # Skip the first line
                for line in f:
                    row = line.rstrip().split("\t")
                    gene_key = row[2]
                    fam_val = row[0].split(".hmm")[0]
                    hmm.add(gene_key)
                    cazyme_genes.setdefault(gene_key, set()).add(fam_val)

        if tools[2]:  ### deal with dbcan_sub result
            with open(outDir + prefix + "dbcan-sub.hmm.out") as f:
                next(f)
                for line in f:
                    row = line.rstrip().split("\t")
                    dbs.add(row[5])
                    cazyme_genes.setdefault(row[5], set()).add(row[0])
        cazyme_genes = {key: '|'.join(values) for key, values in cazyme_genes.items()}
        if tools.count(True) > 1:
            temp1 = hmm.intersection(dbs)
            temp2 = hmm.intersection(dia)
            temp3 = dia.intersection(dbs)
            cazyme = temp1.union(temp2, temp3)
        else:
            cazyme = hmm.union(dia, dbs)
        # End CAZyme Extraction
        ######################
        # Begin GFF preperation

        if inputType in ["prok", "meta"]:  # use Prodigal GFF output
            with open(outDir + prefix + "prodigal.gff") as f:
                with open(outDir + prefix + "cgc.gff", "w") as out:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.rstrip().rstrip(";").split("\t")
                            num = row[-1].split(";")[0].split("_")[-1]
                            gene = row[0] + "_" + num
                            row[8] = ""
                            if gene in cazyme:
                                row[2] = "CAZyme"
                                row[8] = "DB=" + cazyme_genes[gene]
                            elif gene in tf:
                                row[2] = "TF"
                                row[8] = "DB=" + tf_genes[gene]
                            elif gene in tp:
                                row[2] = "TC"
                                row[8] = "DB=" + tp_genes[gene]
                            elif gene in stp:
                                row[2] = "STP"
                                row[8] = "DB=" + stp_genes[gene]
                            row[8] += ";ID=" + gene
                            out.write("\t".join(row) + "\n")
        else:  # user provided GFF/BED file
            gff = False
            with open(auxFile) as f:
                for line in f:
                    if not line.startswith("#"):
                        if len(line.split("\t")) == 9:
                            gff = True
                            break
            if gff:  # user file was in GFF format
                with open(auxFile) as f:
                    with open(outDir + prefix + "cgc.gff", "w") as out:
                        for line in f:
                            if not line.startswith("#"):
                                row = line.rstrip().split("\t")
                                if row[2] == "CDS":
                                    note = row[8].strip().rstrip(";").split(";")
                                    gene = ""
                                    notes = {}
                                    for x in note:
                                        temp = x.split("=")
                                        notes[temp[0]] = temp[1]
                                    if "ID" in notes:
                                        gene = notes["ID"]
                                    else:
                                        continue
                                    if gene in cazyme:
                                        row[2] = "CAZyme"
                                        row[8] = "DB=" + cazyme_genes[gene]
                                    elif gene in tf:
                                        row[2] = "TF"
                                        row[8] = "DB=" + tf_genes[gene]
                                    elif gene in tp:
                                        row[2] = "TC"
                                        row[8] = "DB=" + tp_genes[gene]
                                    elif gene in stp:
                                        row[2] = "STP"
                                        row[8] = "DB=" + stp_genes[gene]
                                    else:
                                        row[8] = ""
                                    row[8] += ";ID=" + gene
                                    out.write("\t".join(row) + "\n")
            else:  # user file was in BED format
                with open(auxFile) as f:
                    with open(outDir + prefix + "cgc.gff", "w") as out:
                        for line in f:
                            if line.startswith("track"):
                                continue
                            row = line.rstrip().rstrip(";").split("\t")
                            outrow = ["."] * 8 + [""]
                            gene = row[1]
                            if gene in cazyme:
                                outrow[2] = "CAZyme"
                                outrow[8] = "DB=" + cazyme_genes[gene]
                            elif gene in tf:
                                outrow[2] = "TF"
                                outrow[8] = "DB=" + tf_genes[gene]
                            elif gene in tp:
                                outrow[2] = "TC"
                                outrow[8] = "DB=" + tp_genes[gene]
                            elif gene in stp:
                                outrow[2] = "STP"
                                outrow[8] = "DB=" + stp_genes[gene]
                            else:
                                outrow[2] = "CDS"
                            outrow[0] = row[0]
                            outrow[3] = row[2]
                            outrow[4] = row[3]
                            outrow[6] = row[4]
                            outrow[8] += ";ID=" + gene
                            out.write("\t".join(outrow) + "\n")
        # End GFF
        ####################
        # Begin CGCFinder call
        print("**************************************CGC-Finder start***********************************************")
        cgc_finder(outDir + prefix + "cgc.gff", cgc_dis, cgc_sig_genes, outDir + prefix + "cgc.out")
        simplify_cgc_output(outDir + prefix + "cgc.out")
        print("**************************************CGC-Finder end***********************************************")
        # End CGCFinder call
        # End CGCFinder
        ####################
    # Begin SignalP combination
    if use_signalP:  ### signalP
        print("Waiting on signalP")
        with open(outDir + prefix + "temp", "w") as out:
            if gram == "all" or gram == "p":
                signalpos.wait()
                print("SignalP pos complete")

                with open(outDir + prefix + "signalp.pos") as f:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.split(" ")
                            row = [x for x in row if x != ""]
                            if row[9] == "Y":
                                out.write(line)
                call(["rm", outDir + prefix + "signalp.pos"])
            if gram == "all" or gram == "n":
                signalpneg.wait()
                print("SignalP neg complete")
                with open(outDir + prefix + "signalp.neg") as f:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.split(" ")
                            row = [x for x in row if x != ""]
                            if row[9] == "Y":
                                out.write(line)
                call(["rm", outDir + prefix + "signalp.neg"])
            if gram == "all" or gram == "euk":
                signalpeuk.wait()
                print("SignalP euk complete")
                with open(outDir + prefix + "signalp.euk") as f:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.split(" ")
                            row = [x for x in row if x != ""]
                            if row[9] == "Y":
                                out.write(line)
                call(["rm", outDir + prefix + "signalp.euk"])
        call("sort -u " + outDir + prefix + "temp > " + outDir + prefix + "signalp.out", shell=True)
        call(["rm", outDir + prefix + "temp"])

    # End SignalP combination
    #######################
    #######################
    # start overview Generation
    def unique(sequence):
        """ Remove duplicates from a list while keeping the original order. """
        seen = set()
        return [x for x in sequence if not (x in seen or seen.add(x))]

    def read_lines_if_file_exists(filepath):
        """ Read non-empty lines from a file if it exists, return an empty list otherwise. """
        if os.path.exists(filepath):
            with open(filepath, 'r') as file:
                return [line.strip() for line in file if line.strip()]
        return []
    
    print("Preparing overview table from hmmer, dbCAN_sub and diamond output...")

    arr_diamond = read_lines_if_file_exists(outPath + "diamond.out")
    arr_hmmer = read_lines_if_file_exists(outPath + "hmmer.out")
    arr_dbsub = read_lines_if_file_exists(outPath + "dbcan-sub.hmm.out")

    diamond_genes = unique([line.split()[0] for line in arr_diamond[1:]])
    hmmer_genes = unique([line.split()[2] for line in arr_hmmer[1:]])
    dbsub_genes = unique([line.split("\t")[5] for line in arr_dbsub[1:]])

    sigp_genes = {}
    if use_signalP and os.path.exists(outPath + "signalp.out"):
        arr_sigp = read_lines_if_file_exists(outPath + "signalp.out")
        for line in arr_sigp:
            parts = line.split()
            sigp_genes[parts[0]] = parts[4]
            
    diamond_fams, hmmer_fams, dbsub_fams = {}, {}, {}

    if tools[0]:  # Diamond
        diamond_fams = {row.split("\t")[0]: row.split("\t")[1].strip("|").split("|")[1:] for row in arr_diamond[1:]}
    
    if tools[1]:  # Hmmer
        for row in arr_hmmer[1:]:
            parts = row.split("\t")
            fam = parts[0].split(".")[0] + "(" + parts[7] + "-" + parts[8] + ")"
            hmmer_fams.setdefault(parts[2], []).append(fam)
    
    if tools[2]:  # dbCAN_sub
        for row in arr_dbsub[1:]:
            parts = row.split("\t")
            gene_id = parts[5]
            dbsub_fams.setdefault(gene_id, {"fam_name": [], "ec_num": []})
            dbsub_fams[gene_id]["fam_name"].append(parts[0])
            dbsub_fams[gene_id]["ec_num"].append(parts[2])
    # Create the overview table
    all_genes = unique(hmmer_genes + dbsub_genes + diamond_genes)

    with open(outPath + "overview.txt", "w+") as fp:
        headers = ["Gene ID", "EC#", "HMMER", "dbCAN_sub", "DIAMOND", "#ofTools"]
        if use_signalP:
            headers.insert(5, "Signalp")
        fp.write("\t".join(headers) + "\n")
        
        for gene in all_genes:
            row = [gene]
            num_tools = 0

            # dbCAN_sub
            row.append("|".join(dbsub_fams[gene]["ec_num"]) if gene in dbsub_fams else "-")

            # HMMER
            if gene in hmmer_fams:
                num_tools += 1
                row.append("+".join(hmmer_fams[gene]))
            else:
                row.append("-")
            
            # dbCAN_sub (fam_name)
            if tools[2] and gene in dbsub_fams:
                row.append("+".join(dbsub_fams[gene]["fam_name"]))
                num_tools += 1
            else:
                row.append("-")

            # DIAMOND
            if gene in diamond_fams:
                num_tools += 1
                row.append("+".join(diamond_fams[gene]))
            else:
                row.append("-")

            # SignalP
            if use_signalP:
                row.append("Y(1-" + sigp_genes[gene] + ")" if gene in sigp_genes else "N")
            
            row.append(str(num_tools))
            fp.write("\t".join(row) + "\n")
    print(f"Overview table complete. Saved as {outPath}overview.txt")
    # End overview


def rundbCAN_parser():
    """
    Sets up and returns an argument parser for the dbCAN3 Command Line Tool.

    This function configures an `argparse.ArgumentParser` for parsing command-line arguments necessary to run the dbCAN3 tool. It defines the expected input file, input type, and various optional parameters for running dbCAN3, including database files, DIAMOND and HMMER settings, output preferences, and additional tools and parameters for CGC finder and substrate prediction.

    Returns
    -------
        argparse.ArgumentParser: Configured argument parser for dbCAN3 command-line tool.

    The parser includes arguments for:
        - inputFile (str): Path to the input file in FASTA format.
        - inputType (str): Type of sequence input (protein, prok, meta).
        - dbCANFile (str): Filename of the HMM database.
        - Various parameters for DIAMOND, HMMER, and other tools.
        - Output settings, database directory, and tool selections.
        - Parameters specific to CGC finder and substrate prediction modules.
    """

    '''

    Example usages of run_dbcan:
    1. CAZyme annotation with isolated genome sequence as input:
    run_dbcan EscheriaColiK12MG1655.fna prok
    
    2. CAZyme annotation with protein sequence as input:
    run_dbcan EscheriaColiK12MG1655.faa protein
    
    3. CAZyme annotation with meta genome as input:
    run_dbcan EscheriaColiK12MG1655.fna meta
    
    4. CAZyme and CGC annotation with meta genome as input:
    run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff
    
    5. CGC substrate prediction with the results from #4:
    run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff --cgc_substrate --only_sub
    
    6. CAZyme, CGC annotation and substrate prediction with mete genome as input:
    run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff --cgc_substrate
    
    7. CAZyme, CGC annotation and substrate prediction with protein and gff as input:
    run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --cgc_substrate
    
    8. CAZyme, CGC annotation and substrate prediction using hmmer against dbCAN-HMM datbase with protein and gff as input:
    run_dbcan EscheriaColiK12MG1655.faa protein -c EscheriaColiK12MG1655.gff --cgc_substrate --tools hmmer
    

    More usage can be found in our https://dbcan.readthedocs.io/en/latest.
    '''

    parser = argparse.ArgumentParser(description="dbCAN: Automatic CAZyme Annotation", prog='run_dbcan')

    parser.add_argument("inputFile", help="User input file. Must be in FASTA format.")
    parser.add_argument(
        "inputType",
        choices=["protein", "prok", "meta"],  # protein=proteome, prok=prokaryote nucleotide, meta=metagenome nucleotide
        help="Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome",
    )
    parser.add_argument(
        "--dbCANFile",
        default="dbCAN.txt",
        help="Indicate the file name of HMM database such as dbCAN.txt, please use the newest one from dbCAN2 website.",
    )
    parser.add_argument("--dia_eval", default=1e-102, type=float, help="DIAMOND E Value")
    parser.add_argument("--dia_cpu", default=8, type=int, help="Number of CPU cores that DIAMOND is allowed to use")
    parser.add_argument("--hmm_eval", default=1e-15, type=float, help="HMMER E Value")
    parser.add_argument("--hmm_cov", default=0.35, type=float, help="HMMER Coverage val")
    parser.add_argument("--hmm_cpu", default=8, type=int, help="Number of CPU cores that HMMER is allowed to use")
    parser.add_argument("--out_pre", default="", help="Output files prefix")
    parser.add_argument("--out_dir", default="output", help="Output directory")
    parser.add_argument("--db_dir", default="db", help="Database directory")
    parser.add_argument(
        "--tools",
        "-t",
        nargs="+",
        choices=["hmmer", "diamond", "dbcansub", "all"],
        default="all",
        help="Choose a combination of tools to run",
    )
    parser.add_argument(
        "--use_signalP",
        default=False,
        type=bool,
        help="Use signalP or not, remember, you need to setup signalP tool first. Because of signalP license, Docker version does not have signalP.",
    )
    parser.add_argument(
        "--signalP_path", "-sp", default="signalp", type=str, help="The path for signalp. Default location is signalp"
    )
    parser.add_argument(
        "--gram",
        "-g",
        choices=["p", "n", "all"],
        default="all",
        help="Choose gram+(p) or gram-(n) for proteome/prokaryote nucleotide, which are params of SingalP, only if user use singalP",
    )
    parser.add_argument("-v", "--version", default="4.1.1", type=str)
    # dbCAN-sub
    dbCAN_sub_group = parser.add_argument_group("dbCAN-sub parameters")
    dbCAN_sub_group.add_argument("--dbcan_thread", "-dt", default=12, type=int)
    dbCAN_sub_group.add_argument("--tf_eval", default=1e-4, type=float, help="tf.hmm HMMER E Value")
    dbCAN_sub_group.add_argument("--tf_cov", default=0.35, type=float, help="tf.hmm HMMER Coverage val")
    dbCAN_sub_group.add_argument(
        "--tf_cpu", default=8, type=int, help="tf.hmm Number of CPU cores that HMMER is allowed to use"
    )
    dbCAN_sub_group.add_argument("--stp_eval", default=1e-4, type=float, help="stp.hmm HMMER E Value")
    dbCAN_sub_group.add_argument("--stp_cov", default=0.3, type=float, help="stp.hmm HMMER Coverage val")
    dbCAN_sub_group.add_argument(
        "--stp_cpu", default=8, type=int, help="stp.hmm Number of CPU cores that HMMER is allowed to use"
    )

    ### cgc finder
    cgcfinder_group = parser.add_argument_group("CGC_Finder parameters")
    cgcfinder_group.add_argument(
        "--cluster",
        "-c",
        help="Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used",
    )
    cgcfinder_group.add_argument("--cgc_dis", default=2, type=int, help="CGCFinder Distance value")
    cgcfinder_group.add_argument(
        "--cgc_sig_genes",
        default="tp",
        choices=["tf", "tp", "stp", "tp+tf", "tp+stp", "tf+stp", "all"],
        help="CGCFinder Signature Genes value",
    )

    ### cgc substrate prediction
    cgcsubstrate_group = parser.add_argument_group("CGC_Substrate parameters")
    cgcsubstrate_group.add_argument('--only_sub',action='store_false',help="Only run substrate prediction for PUL. If this parameter is presented, dbcan will skip the CAZyme annotation and CGC prediction.")
    cgcsubstrate_group.add_argument("--cgc_substrate", action="store_true", help="run cgc substrate prediction?")
    cgcsubstrate_group.add_argument("--pul", help="dbCAN-PUL PUL.faa")
    cgcsubstrate_group.add_argument("-o", "--out", default="substrate.out")
    cgcsubstrate_group.add_argument("-w", "--workdir", type=str, default=".")
    cgcsubstrate_group.add_argument("-env", "--env", type=str, default="local")
    cgcsubstrate_group.add_argument(
        "-odbcan_sub", "--odbcan_sub", action="store_true", help="Output dbCAN-sub prediction intermediate result? for debug"
    )
    cgcsubstrate_group.add_argument(
        "-odbcanpul", "--odbcanpul", action="store_true", help="Output dbCAN-PUL prediction intermediate result? for debug."
    )

    ### cgc substrate prediction:dbCAN-PUL
    group1 = parser.add_argument_group(
        "dbCAN-PUL homologous searching parameters", "how to define homologous gene hits and PUL hits"
    )
    group1.add_argument("-upghn", "--uniq_pul_gene_hit_num", default=2, type=int)
    group1.add_argument("-uqcgn", "--uniq_query_cgc_gene_num", default=2, type=int)
    group1.add_argument("-cpn", "--CAZyme_pair_num", default=1, type=int)
    group1.add_argument("-tpn", "--total_pair_num", default=2, type=int)
    group1.add_argument(
        "-ept", "--extra_pair_type", default=None, type=str, help="None[TC-TC,STP-STP]. Some like sigunature hits"
    )
    group1.add_argument(
        "-eptn", "--extra_pair_type_num", default="0", type=str, help="specify signature pair cutoff.1,2"
    )
    group1.add_argument(
        "-iden", "--identity_cutoff", default=0.3, type=float, help="identity to identify a homologous hit"
    )
    group1.add_argument(
        "-cov", "--coverage_cutoff", default=0.3, type=float, help="query coverage cutoff to identify a homologous hit"
    )
    group1.add_argument(
        "-bsc", "--bitscore_cutoff", default=50, type=float, help="bitscore cutoff to identify a homologous hit"
    )
    group1.add_argument(
        "-evalue", "--evalue_cutoff", default=0.01, type=float, help="evalue cutoff to identify a homologous hit"
    )

    ### cgc substrate prediction:dbCAN-sub
    group2 = parser.add_argument_group(
        "dbCAN-sub major voting parameters", "how to define dbsub hits and dbCAN-sub subfamily substrate"
    )
    group2.add_argument("-hmmcov", "--hmmcov", default=0.3, type=float)
    group2.add_argument("-hmmevalue", "--hmmevalue", default=0.01, type=float)
    group2.add_argument(
        "-ndsc",
        "--num_of_domains_substrate_cutoff",
        default=2,
        type=int,
        help="define how many domains share substrates in a CGC, one protein may include several subfamily domains.",
    )
    group2.add_argument(
        "-npsc",
        "--num_of_protein_substrate_cutoff",
        default=2,
        type=int,
        help="define how many sequences share substrates in a CGC, one protein may include several subfamily domains.",
    )
    group2.add_argument(
        "-subs",
        "--substrate_scors",
        default=2,
        type=int,
        help="each cgc contains with substrate must more than this value",
    )

    return parser


# Putting the ArgumentParser in this block allows the script to be called from command line as before, while
# allowing the main function to be called directly from other scripts without invoking a subprocess. This prevents extra
# subprocesses or extra python interpreters being spawned, as well as simplifying python scripts which call run_dbcan.
def cli_main():
    
    parser = rundbCAN_parser()
    args = parser.parse_args()

    ### rundbCAN3
    if args.only_sub:
        run_dbCAN(
            inputFile=args.inputFile,
            inputType=args.inputType,
            cluster=args.cluster,
            dbCANFile=args.dbCANFile,
            dia_eval=args.dia_eval,
            dia_cpu=args.dia_cpu,
            hmm_eval=args.hmm_eval,
            hmm_cov=args.hmm_cov,
            hmm_cpu=args.hmm_cpu,
            dbcan_thread=args.dbcan_thread,
            tf_eval=args.tf_eval,
            tf_cov=args.tf_cov,
            tf_cpu=args.tf_cpu,
            stp_eval=args.stp_eval,
            stp_cov=args.stp_cov,
            stp_cpu=args.stp_cpu,
            prefix=args.out_pre,
            outDir=args.out_dir,
            dbDir=args.db_dir,
            cgc_dis=args.cgc_dis,
            cgc_sig_genes=args.cgc_sig_genes,
            tool_arg=args.tools,
            use_signalP=args.use_signalP,
            signalP_path=args.signalP_path,
            gram=args.gram,
        )

    ### convert cgc_standard.out to json format

    if args.cluster:  ### run cgc_finder
        os.system(f"cgc_standard2json -i {args.out_dir}/cgc_standard.out -o {args.out_dir}/cgc_standard.out.json")
    ### substarate prediction
    if args.cgc_substrate:
        cgc_substrate_prediction(args)


if __name__ == "__main__":
    cli_main()
