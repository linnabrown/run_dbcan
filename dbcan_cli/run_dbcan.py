#!/usr/bin/env python3
#########################################################
# dbCAN3 (Stand Alone Version)
#
# Written by Tanner Yohe in the Yin Lab at NIU
# Revised by Qiwei Ge in Yin Lab at UNL && Le Huang at NKU
# Updated by Le Huang at NKU, Mohamad Majd Raslan in the Yin Lab at NIU, Wei Li, Qiwei Ge in Dr.Yin's Lab at UNL, Alex Fraser.
# Updated by Jinfang Zheng in Yinlab at UNL, new function, substrate prediciton based on dbCAN-PUL and dbCAN-sub database.

# Recent updated information: 
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
# tools from dbsub.out,hmmer.out and diamond.out
##########################################################
from subprocess import Popen, call, check_output
import os
import argparse
import dbcan
from dbcan.utils.simplify_cgc import simplify_output
from dbcan.utils.CGCFinder import cgc_finder
from dbcan_cli import hmmscan_parser
import time
from dbcan.utils.cgc_substrate_prediction import cgc_substrate_prediction


def runHmmScan(outPath, hmm_cpu, dbDir, hmm_eval, hmm_cov, db_name):
    '''
    Run Hmmer
    '''
    hmmer = Popen(['hmmscan', '--domtblout', '%sh%s.out' % (outPath, db_name), '--cpu', hmm_cpu, '-o', '/dev/null', '%s%s.hmm' % (dbDir,db_name), '%suniInput' % outPath])
    hmmer.wait()
    parsed_hmm_output = hmmscan_parser.run(input_file=f"{outPath}h{db_name}.out", eval_num=hmm_eval, coverage=hmm_cov)
    with open(f"{outPath}{db_name}.out", 'w') as f:
        f.write(parsed_hmm_output)
    if os.path.exists('%sh%s.out' % (outPath, db_name)):
        call(['rm', '%sh%s.out' % (outPath, db_name)])

def split_uniInput(uniInput,dbcan_thread,outPath,dbDir,hmm_eval,hmm_cov):
    '''
    Run dbcan_sub
    '''
    ticks = time.time()
    file = open(uniInput, "r")
    uniInput_file = file.readlines()
    file.close()
    signal_count = 0
    split_size = 0
    min_files = dbcan_thread
    check_id = False
    file_number = None
    split_files = []
    off_set = 3
    fsize = int(os.path.getsize(uniInput)/float(1024*1024)*off_set)

    if fsize < 1:
        fsize = 1

    for line in uniInput_file:
        if ">" in line:
            signal_count+=1
    print("ID count: %s" % signal_count)

    if signal_count >= min_files:
        for i in range(fsize):
            f = open("%s%s.txt"%(outPath,i),"w")
            f.close()
            split_files.append("%s.txt"%i)
        for i in range(len(uniInput_file)):
            if ">" in uniInput_file[i]:
                file_number = i%fsize
                f = open('%s%s.txt'%(outPath,file_number), 'a')
                f.write(uniInput_file[i])
                f.close()
            else:
                f = open('%s%s.txt'%(outPath,file_number), 'a')
                f.write(uniInput_file[i])
                f.close()

        ths = []
        for j in split_files:
            ths.append(Popen(['hmmscan', '--domtblout', '%sd%s'%(outPath,j), '--cpu', '5', '-o', '/dev/null', '%sdbCAN_sub.hmm'%dbDir, "%s%s"%(outPath,j)]))
        for th in ths:
            th.wait()

        for m in split_files:
            hmm_parser_output = hmmscan_parser.run("%sd%s"%(outPath,m), eval_num=hmm_eval, coverage=hmm_cov)
            with open("%stemp_%s"%(outPath,m), 'w') as temp_hmmer_file:
                temp_hmmer_file.write(hmm_parser_output)
            call(['rm', '%sd%s'%(outPath,m)])
            call(['rm', '%s%s'%(outPath,m)]) #remove temporary files
            
        f = open("%sdtemp.out"%outPath,"w")
        f.close()

        for n in split_files:
            file_read = open("%stemp_%s"%(outPath,n),"r")
            files_lines = file_read.readlines()
            file_read.close()
            call(['rm', "%stemp_%s"%(outPath,n)]) #remove temporary files
            for j in range(len(files_lines)):
                f = open("%sdtemp.out"%outPath,"a")
                f.write(files_lines[j])
                f.close()
    else:
        dbsub = Popen(['hmmscan', '--domtblout', '%sd.txt'%outPath, '--cpu', '5', '-o', '/dev/null', '%sdbCAN_sub.hmm'%dbDir, '%suniInput'%outPath])
        dbsub.wait()

        hmm_parser_output = hmmscan_parser.run("%sd.txt"%outPath, eval_num=hmm_eval, coverage=hmm_cov)
        with open("%sdtemp.out"%outPath, 'w') as temp_hmmer_file:
            temp_hmmer_file.write(hmm_parser_output)

    print("total time:",time.time() - ticks)

def run(inputFile, inputType, cluster=None, dbCANFile="dbCAN.txt", dia_eval=1e-102, dia_cpu=4, hmm_eval=1e-15,
        hmm_cov=0.35, hmm_cpu=4, dbcan_thread=5, tf_eval=1e-4, tf_cov=0.35, tf_cpu=1, stp_eval=1e-4, stp_cov=0.3, stp_cpu=1, prefix="",
        outDir="output", dbDir="db", cgc_dis=2, cgc_sig_genes="tp", tool_arg="all", use_signalP=False,
        signalP_path="signalp", gram="all"):
    '''
    Run dbCAN
    '''

    # Begin Setup and Input Checks
    if not dbDir.endswith("/") and len(dbDir) > 0:
        dbDir += "/"

    if not outDir.endswith("/") and len(outDir) > 0:
        outDir += "/"

    outPath = outDir + prefix
    auxFile = ""

    find_clusters = False
    if cluster != None:
        find_clusters = True
        if inputType == "protein":
            auxFile = cluster
        else:
            auxFile = '%sprodigal.gff'%outPath

    if not os.path.isdir(dbDir):
        print(dbDir , "ERROR: The database directory does not exist")
        exit()

    if not os.path.isfile(os.path.join(dbDir,'CAZy.dmnd')):
        print("ERROR: No CAZy DIAMOND database found. \
        Please make sure that your CAZy DIAMOND databased is named 'CAZy.dmnd' and is located in your database directory")
        exit()

    if not os.path.isfile(os.path.join(dbDir, dbCANFile)):
        print("ERROR: No dbCAN HMM database found. \
        Please make sure that your dbCAN HMM database is named 'dbCAN-HMMdb-V11.txt' or the newest one, has been through hmmpress, and is located in your database directory")
        exit()

    if not os.path.isfile(os.path.join(dbDir,'dbCAN_sub.hmm')):
        print("ERROR: No dbCAN_sub HMM database found. \
        Please make sure that your dbCAN_sub HMM databased is named 'dbCAN_sub.hmm' or has been through hmmpress, and is located in your database directory")
        exit()

    if not os.path.isdir(outDir):
        call(['mkdir', outDir])

    if find_clusters and inputType == "protein":
        if len(auxFile) > 0:
            print(auxFile)
            if not os.path.isfile(auxFile):
                print("ERROR: It seems that the auxillary filename that you provided does not exist, or is not a file")
                exit()
        else:
            print("ERROR: Please provide an auxillary input file with the position of each gene. This file can either be in BED or GFF format")
            exit()

    tools = [True, True, True] #DIAMOND, HMMER, dbCAN_sub
    if 'all' not in tool_arg:
        if 'diamond' not in tool_arg:
            tools[0] = False
        if 'hmmer' not in tool_arg:
            tools[1] = False
        if 'dbcansub' not in tool_arg:
            tools[2] = False

    # End Setup and Input Checks
    #########################
    #########################
    # Begin Gene Prediction Tools
    if inputType == 'prok':
        call(['prodigal', '-i', inputFile, '-a', '%suniInput'%outPath, '-o', '%sprodigal.gff'%outPath, '-f', 'gff', '-q'])
    if inputType == 'meta':
        call(['prodigal', '-i', inputFile, '-a', '%suniInput'%outPath, '-o', '%sprodigal.gff'%outPath, '-f', 'gff', '-p', 'meta','-q'])
    #Proteome
    if inputType == 'protein':
        call(['cp', inputFile, '%suniInput'%outPath])

    # End Gene Prediction Tools
    #######################
    # Begin SignalP
    if use_signalP:
        print("\n\n***************************0. SIGNALP start*************************************************\n\n")
        if gram == "p" or gram=="all":
            signalpos = Popen('%s -t gram+ %suniInput > %ssignalp.pos' % (signalP_path, outPath, outPath), shell=True)
        if gram == "n" or gram == "all":
            signalpneg = Popen('%s -t gram- %suniInput > %ssignalp.neg' % (signalP_path, outPath, outPath), shell=True)
        if gram == "euk" or gram=="all":
            signalpeuk = Popen('%s -t euk %suniInput > %ssignalp.euk' % (signalP_path, outPath, outPath), shell=True)

    # End SignalP
    #######################
    # Begin Core Tools

    if tools[0]: ### run diamond
        # diamond blastp -d db/CAZy -e 1e-102 -q output_EscheriaColiK12MG1655/uniInput -k 1 -p 2 -o output_EscheriaColiK12MG1655/diamond1.out -f 6
        print("\n\n***************************1. DIAMOND start*************************************************\n\n")
        os.system('diamond blastp -d %s -e %s -q %suniInput -k 1 -p %d -o %sdiamond.out -f 6'%(os.path.join(dbDir, "CAZy"), str(dia_eval), outPath, dia_cpu, outPath))
        print("\n\n***************************1. DIAMOND end***************************************************\n\n")

    if tools[1]: ### run hmmscan (hmmer)
        print("\n\n***************************2. HMMER start*************************************************\n\n")
        os.system(f"hmmscan --domtblout {outPath}h.out --cpu {hmm_cpu} -o /dev/null {os.path.join(dbDir, dbCANFile)} {outPath}uniInput ")
        print("\n\n***************************2. HMMER end***************************************************\n\n")

        hmm_parser_output = hmmscan_parser.run(f"{outPath}h.out", eval_num=hmm_eval, coverage=hmm_cov)
        with open(f"{outPath}hmmer.out", 'w') as hmmer_file:
            hmmer_file.write(hmm_parser_output)
        # could clean this up and manipulate hmm_parser_output data directly instead of passing it into a temp file
        with open(f"{outPath}hmmer.out", "r+") as f:
            text = f.read()
            f.close()
            call(['rm', f"{outPath}hmmer.out"])
            text = text.split('\n')
            if '' in text:
                text.remove('')
            for i in range(len(text)):
                if 'GT2_' in text[i]:
                    profile = text[i].split('\t')[0].split('.')[0]
                    text[i] = text[i].replace(profile,'GT2')
                with open(f"{outPath}hmmer.out", 'a') as f:
                    f.write(text[i]+'\n')
                    f.close()
        if os.path.exists(f"{outPath}h.out"):
            call(['rm', f"{outPath}h.out"])

    if tools[2]:
        print("\n\n***************************3. dbCAN_sub start***************************************************\n\n")
        split_uniInput('%suniInput'%outPath,dbcan_thread,outPath,dbDir,hmm_eval,hmm_cov)
        print("\n\n***************************3. dbCAN_sub end***************************************************\n\n")
        with open(f"{outPath}dtemp.out", 'r') as f:
            with open('%sdbsub.out'%outPath, 'w') as out:
                for line in f:
                    row = line.rstrip().split('\t')
                    row.append(float(int(row[6])-int(row[5]))/int(row[1]))
                    if float(row[4]) <= 1e-15 and float(row[-1]) >= 0.35:
                        out.write('\t'.join([str(x) for x in row]) + '\n')
        with open(f"{outPath}dbsub.out", 'r+') as f: #formated GT2_ in hmmer.out
            text = f.read()
            f.close()
            call(['rm', f"{outPath}dbsub.out"])
            text = text.split('\n')
            if '' in text:
                text.remove('')
            for i in range(len(text)):
                if 'GT2_' in text[i]:
                    profile = text[i].split('\t')[0].split('.')[0]
                    text[i] = text[i].replace(profile,'GT2')
                with open(f"{outPath}dbsub.out", 'a') as f:
                    f.write(text[i]+'\n')
                    f.close()
    # End Core Tools
    ########################
    # Begin Parse Results

    # parse dbCAN_sub result
    if tools[2]: 
        subs_dict = {}
        with open(f"{dbDir}fam-substrate-mapping-08252022.tsv", 'r') as f:
            next(f)
            for line in f:
                r = line.split("\t")
                if len(r[4]) == 1:
                    subs_dict[r[2],"-"] = r[0]
                else:
                    subs_dict[r[2],r[4].strip()] = r[0]
        with open(f"{outPath}dbsub.out") as f:
            with open(f"{outPath}temp", 'w') as out:
                out.write('dbCAN subfam\tSubfam Composition\tSubfam EC\tSubstrate\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n')

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
                            # print(key1,key2[i])
                            substrate.append(subs_dict[key1,key2[i]])
                        except:
                            print("No substrate for it")

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
        call(['mv', outDir+prefix+'temp', outDir+prefix+'dbsub.out'])

    # parse hmmer result
    if tools[1]:
        try:
            with open(outDir+prefix+'hmmer.out') as f:
                with open(outDir+prefix+'temp', 'w') as out:
                    out.write('HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n')
                    for line in f:
                        out.write(line)
            call(['mv', outDir+prefix+'temp', outDir+prefix+'hmmer.out'])
        except:
            with open(outDir+prefix+'temp', 'w') as out:
                out.write('HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n')
            call(['mv', outDir+prefix+'temp', outDir+prefix+'hmmer.out'])

    # parse diamond result
    if tools[0]: 
        with open(outDir+prefix+'diamond.out') as f:
            with open(outDir+prefix+'temp', 'w') as out:
                out.write('Gene ID\tCAZy ID\t% Identical\tLength\tMismatches\tGap Open\tGene Start\tGene End\tCAZy Start\tCAZy End\tE Value\tBit Score\n')
                for line in f:
                    out.write(line)
        call(['mv', outDir+prefix+'temp', outDir+prefix+'diamond.out'])

    # End Parse Results
    ########################
    # Begin CGCFinder

    if find_clusters: ### run cgc_finder or not
        print("*****************************CGC-Finder start************************************")

        ########################
        # Begin TF,TP, STP prediction
        '''
        tf hmmer
        '''
        #call(['diamond', 'blastp', '-d', dbDir+'tf_v1/tf.dmnd', '-e', '1e-10', '-q', '%suniInput' % outPath, '-k', '1', '-p', '1', '-o', outDir+prefix+'tf.out', '-f', '6'])
        runHmmScan(outPath, str(tf_cpu), dbDir, str(tf_eval), str(tf_cov), "tf-1")
        runHmmScan(outPath, str(tf_cpu), dbDir, str(tf_eval), str(tf_cov), "tf-2")
        '''
        stp hmmer
        '''
        runHmmScan(outPath, str(stp_cpu), dbDir, str(stp_eval), str(stp_cov), "stp")

        '''
        tp diamond
        '''
        call(['diamond', 'blastp', '-d', dbDir+'tcdb.dmnd', '-e', '1e-10', '-q', '%suniInput' % outPath, '-k', '1', '-p', '1', '-o', outPath+'tp.out', '-f', '6'])


        tp = set()
        tf = set()
        stp = set()

        tp_genes = {}
        tf_genes = {}
        stp_genes = {}

        with open("%stf-1.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split('\t')
                tf.add(row[2])
                row[0] = "DBD-Pfam|" + row[0]
                if not row[2] in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += ',' + row[0]

        with open("%stf-2.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split('\t')
                tf.add(row[2])
                row[0] = "DBD-SUPERFAMILY|" + row[0]
                if not row[2] in tf_genes:
                    tf_genes[row[2]] = row[0]
                else:
                    tf_genes[row[2]] += ',' + row[0]

        with open(outDir+prefix+'tp.out') as f:
            for line in f:
                row = line.rstrip().split('\t')
                tp.add(row[0])
                if not row[0] in tp_genes:
                    tp_genes[row[0]] = row[1]
                else:
                    tp_genes[row[0]] += ','+row[1]

        with open("%sstp.out" % outPath) as f:
            for line in f:
                row = line.rstrip().split('\t')
                stp.add(row[2])
                row[0] = "STP|" + row[0]
                if not row[2] in stp_genes:
                    stp_genes[row[2]] = row[0]
                else:
                    stp_genes[row[2]] += ',' + row[0]
        # End TF and TP prediction
        ##########################
        # Begine CAZyme Extraction
        cazyme_genes = {}
        
        dia = set()
        hmm = set()
        dbs = set()

        if tools[0]: ### deal with diamond result
            with open(outDir+prefix+'diamond.out') as f:
                next(f)
                for line in f:
                    row = line.rstrip().split('\t')
                    dia.add(row[0])
                    if row[0] not in cazyme_genes:
                        cazyme_genes[row[0]] = set()
                        cazyme_genes[row[0]].update(set(row[1].strip("|").split('|')[1:]))

        if tools[1]: ### deal with hmmscan result
            with open(outDir+prefix+'hmmer.out') as f:
                next(f)
                for line in f:
                    row = line.rstrip().split('\t')
                    hmm.add(row[2])
                    if row[2] not in cazyme_genes:
                        cazyme_genes[row[2]] = set()
                    if row[0].split('.hmm')[0] in cazyme_genes[row[2]]:
                        cazyme_genes[row[2]].add(" "+row[0].split('.hmm')[0])
                    else:
                        cazyme_genes[row[2]].add(row[0].split('.hmm')[0])

        if tools[2]: ### deal with eCAMI result 
            with open(outDir+prefix+'dbsub.out') as f:
                next(f)
                for line in f:
                    row = line.rstrip().split('\t')
                    dbs.add(row[5])
                    if row[5] not in cazyme_genes:
                        cazyme_genes[row[5]] = set()
                        cazyme_genes[row[5]].add(row[0])

        
        if tools.count(True) > 1:
            temp1 = hmm.intersection(dbs)
            # print(hmm, 'This intersection  hmm')
            temp2 = hmm.intersection(dia)
            # print(dia, 'This intersection  dia')
            temp3 = dia.intersection(dbs)
            # print(dbs, 'This intersection  dbs')
            cazyme = temp1.union(temp2, temp3)
        else:
            cazyme = hmm.union(dia, dbs)
        # End CAZyme Extraction
        ######################
        # Begin GFF preperation

        if inputType == "prok" or inputType == "meta":   #use Prodigal GFF output
            with open(outDir+prefix+'prodigal.gff') as f:
                with open(outDir+prefix+'cgc.gff', 'w') as out:
                    for line in f:
                        if not line.startswith("#"):
                            row = line.rstrip().rstrip(";").split('\t')
                            num = row[-1].split(";")[0].split('_')[-1]
                            gene = row[0] + '_' + num
                            row[8] = ""
                            if gene in cazyme:
                                row[2] = "CAZyme"
                                # Uncomment this, if all CAZyme results need to be write into cgc.out
                                row[8] = "DB="+'|'.join(cazyme_genes[gene])
                                #
                                # cazyme_genes_list = list(cazyme_genes[gene])
                                # row[8] = "DB="+cazyme_genes_list[0]
                                #
                            elif gene in tf:
                                row[2] = "TF"
                                row[8] = "DB="+tf_genes[gene]
                            elif gene in tp:
                                row[2] = "TC"
                                row[8] = "DB="+tp_genes[gene]
                            elif gene in stp:
                                row[2] = "STP"
                                row[8] = "DB="+stp_genes[gene]
                            row[8] += ";ID="+gene
                            out.write('\t'.join(row)+'\n')
        else:  #user provided GFF/BED file
            gff = False
            with open(auxFile) as f:
                for line in f:
                    if not line.startswith('#'):
                        if len(line.split('\t')) == 9:
                            gff = True
                            break
            if gff:  #user file was in GFF format
                with open(auxFile) as f:
                    with open(outDir+prefix+'cgc.gff', 'w') as out:
                        for line in f:
                            if not line.startswith("#"):
                                row = line.rstrip().split('\t')
                                if row[2] == "CDS":
                                    note = row[8].strip().rstrip(";").split(";")
                                    gene = ""
                                    notes = {}
                                    for x in note:
                                        temp = x.split('=')
                                        notes[temp[0]] = temp[1]
                                    # if "Name" in notes:
                                    #     gene = notes["Name"]
                                    # elif "ID" in notes:
                                    #     gene = notes["ID"]
                                    if "ID" in notes:
                                        gene = notes["ID"]
                                    else:
                                        continue

                                    if gene in cazyme:
                                        row[2] = "CAZyme"
                                        # Uncomment this, if all CAZyme results need to be write into cgc.out
                                        row[8] = "DB="+'|'.join(cazyme_genes[gene])
                                        #
                                        # cazyme_genes_list = list(cazyme_genes[gene])
                                        # row[8] = "DB="+cazyme_genes_list[0]
                                        #
                                    elif gene in tf:
                                        row[2] = "TF"
                                        row[8] = "DB="+tf_genes[gene]
                                    elif gene in tp:
                                        row[2] = "TC"
                                        row[8] = "DB="+tp_genes[gene]
                                    elif gene in stp:
                                        row[2] = "STP"
                                        row[8] = "DB=" + stp_genes[gene]
                                    else:
                                        row[8] = ""
                                    row[8] += ";ID="+gene
                                    out.write('\t'.join(row)+'\n')
            else:  #user file was in BED format
                with open(auxFile) as f:
                    with open(outDir+prefix+'cgc.gff', 'w') as out:
                        for line in f:
                            if line.startswith("track"):
                                continue
                            row = line.rstrip().rstrip(";").split('\t')
                            outrow = ['.'] * 8 + ['']
                            gene = row[1]
                            if gene in cazyme:
                                outrow[2] = 'CAZyme'
                                # Uncomment this, if all CAZyme results need to be write into cgc.out
                                outrow[8] = "DB="+'|'.join(cazyme_genes[gene])
                                #
                                # cazyme_genes_list = list(cazyme_genes[gene])
                                # outrow[8] = "DB="+cazyme_genes_list[0]
                                #
                            elif gene in tf:
                                outrow[2] = 'TF'
                                outrow[8] =  "DB="+tf_genes[gene]
                            elif gene in tp:
                                outrow[2] = 'TC'
                                outrow[8] = "DB="+tp_genes[gene]
                            elif gene in stp:
                                outrow[2] = 'STP'
                                outrow[8] = "DB=" + stp_genes[gene]
                            else:
                                outrow[2] = 'CDS'
                            outrow[0] = row[0]
                            outrow[3] = row[2]
                            outrow[4] = row[3]
                            outrow[6] = row[4]
                            outrow[8] += ";ID="+gene
                            out.write('\t'.join(outrow)+'\n')
        # End GFF
        ####################
        # Begin CGCFinder call
        print("**************************************CGC-Finder start***********************************************")
        # call(['CGCFinder.py', outDir+prefix+'cgc.gff', '-o', outDir+prefix+'cgc.out', '-s', args.cgc_sig_genes, '-d', str(args.cgc_dis)])
        cgc_finder(outDir+prefix+'cgc.gff', cgc_dis, cgc_sig_genes, outDir+prefix+'cgc.out')
        simplify_output(outDir+prefix+'cgc.out')
        print("**************************************CGC-Finder end***********************************************")
        # End CGCFinder call
        # End CGCFinder
        ####################
    # Begin SignalP combination
    if use_signalP: ### signalP
        print("Waiting on signalP")
        with open(outDir+prefix+'temp', 'w') as out:
            if gram == "all" or gram =="p":
                signalpos.wait()
                print("SignalP pos complete")

                with open(outDir+prefix+'signalp.pos') as f:
                    for line in f:
                        if not line.startswith('#'):
                            row = line.split(' ')
                            row = [x for x in row if x != '']
                            if row[9] == 'Y':
                                out.write(line)
                call(['rm', outDir+prefix+'signalp.pos'])
            if gram == "all" or gram == "n":
                signalpneg.wait()
                print("SignalP neg complete")
                with open(outDir+prefix+'signalp.neg') as f:
                    for line in f:
                        if not line.startswith('#'):
                            row = line.split(' ')
                            row = [x for x in row if x != '']
                            if row[9] == 'Y':
                                out.write(line)
                call(['rm', outDir+prefix+'signalp.neg'])
            if gram == "all" or gram == "euk":
                signalpeuk.wait()
                print("SignalP euk complete")
                with open(outDir+prefix+'signalp.euk') as f:
                    for line in f:
                        if not line.startswith('#'):
                            row = line.split(' ')
                            row = [x for x in row if x != '']
                            if row[9] == 'Y':
                                out.write(line)
                call(['rm', outDir+prefix+'signalp.euk'])
        call('sort -u '+outDir+prefix+'temp > '+outDir+prefix+'signalp.out', shell=True)
        call(['rm', outDir+prefix+'temp'])

    # End SignalP combination
    #######################
    #######################
    # start Overview
    print("Preparing overview table from hmmer, dbCAN_sub and diamond output...")
    workdir = outDir+prefix
    # a function to remove duplicates from lists while keeping original order
    def unique(seq):
        exists = set()
        return [x for x in seq if not (x in exists or exists.add(x))]

    arr_dbsub = None
    arr_hmmer = None

    # check if files exist. if so, read files and get the gene numbers
    if tools[0]:
        arr_diamond = open(workdir+"diamond.out").readlines()
        diamond_genes = [arr_diamond[i].split()[0] for i in range(1, len(arr_diamond))] # or diamond_genes = []

    if tools[1]:
        arr_hmmer = open(workdir+"hmmer.out").readlines()
        hmmer_genes = [arr_hmmer[i].split()[2] for i in range(1, len(arr_hmmer))] # or hmmer_genes = []

    if tools[2]:
        arr_dbsub = open(workdir+"dbsub.out").readlines()
        dbsub_genes = [arr_dbsub[i].split("\t")[5] for i in range(1, len(arr_dbsub))]# or dbsub_genes = []

    if use_signalP and (os.path.exists(workdir + "signalp.out")):
        arr_sigp = open(workdir+"signalp.out").readlines()
        sigp_genes = {}
        for i in range (0,len(arr_sigp)):
            row = arr_sigp[i].split()
            sigp_genes[row[0]] = row[4] #previous one is row[2], use Y-score instead from suggestion of Dongyao Li

    ##Catie Ausland edits BEGIN, Le add variable exists or not, remove duplicates from input lists
    if not tools[0]:
        diamond_genes = []
    if not tools[1]:
        hmmer_genes   = []
    if not tools[2]:
        dbsub_genes  = []

    if len(dbsub_genes) > 0:
        if (dbsub_genes[-1] == None):
            #print('I am in &&&&&&&&&&&&&&&&&&&&&&')
            dbsub_genes.pop()
            dbsub_genes = unique(dbsub_genes)
            if 'hmmer_genes' in locals():
                hmmer_genes.pop()
                hmmer_genes = unique(hmmer_genes)
            if 'diamond_genes' in locals():
                diamond_genes.pop()
                diamond_genes = unique(diamond_genes)
    ## Catie edits END, Le add variable exists or not, remove duplicates from input lists

    # parse input, stroe needed variables
    if tools[0] and (len(arr_diamond) > 1):
        diamond_fams = {}
        for i in range (1,len(arr_diamond)):
            row = arr_diamond[i].split("\t")
            fam = row[1].strip("|").split("|")
            diamond_fams[row[0]] = fam[1:]

    if tools[1] and (len(arr_hmmer) > 1):
        hmmer_fams = {}
        for i in range (1, len(arr_hmmer)):
            row = arr_hmmer[i].split("\t")
            fam = row[0].split(".")
            fam = fam[0]+"("+row[7]+"-"+row[8]+")"
            if(row[2] not in hmmer_fams):
                hmmer_fams[row[2]] = []
            hmmer_fams[row[2]].append(fam)

    if tools[2] and (len(arr_dbsub) > 1) :
        dbsub_fams = {}
        for i in range (1,len(arr_dbsub)):
            row_ori = arr_dbsub[i].split("\t")
            fams_ID = row_ori[5]
            if fams_ID not in dbsub_fams:
                dbsub_fams[fams_ID] = {}
                dbsub_fams[fams_ID]["fam_name"] = []
                dbsub_fams[fams_ID]["ec_num"] = []

            dbsub_fams[fams_ID]["fam_name"].append(row_ori[0])
            dbsub_fams[fams_ID]["ec_num"].append(row_ori[2])

    #overall table

    all_genes = unique(hmmer_genes+dbsub_genes+diamond_genes)

    with open(workdir+"overview.txt", 'w+') as fp:
        if use_signalP:
            fp.write("Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\tSignalp\t#ofTools\n")
        else:
            fp.write("Gene ID\tEC#\tHMMER\tdbCAN_sub\tDIAMOND\t#ofTools\n")
        for gene in all_genes:
            csv=[gene]
            num_tools = 0

            if tools[2] and arr_dbsub != None and (gene in dbsub_genes):
                if dbsub_fams[gene]["ec_num"] == []:
                    csv.append("-")
                else:
                    csv.append("|".join(dbsub_fams[gene]["ec_num"]))
            else:
                csv.append("-")

            if tools[1] and arr_hmmer != None and (gene in hmmer_genes):
                num_tools += 1
                csv.append("+".join(hmmer_fams[gene]))
            else:
                csv.append("-")

            if tools[2] and arr_dbsub!= None and (gene in dbsub_genes):
                num_tools += 1
                csv.append("+".join(dbsub_fams[gene]["fam_name"]))
            else:
                csv.append("-")

            if tools[0] and arr_diamond != None and (gene in diamond_genes):
                num_tools += 1
                csv.append("+".join(diamond_fams[gene]))
            else:
                csv.append("-")
            if use_signalP:
                if (gene in sigp_genes):
                    csv.append("Y(1-"+sigp_genes[gene]+")")
                else:
                    csv.append("N")
            csv.append(str(num_tools))
            temp = "\t".join(csv) + "\n"
            fp.write(temp)
    print("overview table complete. Saved as "+workdir+"overview.txt")
    # End overview


# Putting the ArgumentParser in this block allows the script to be called from command line as before, while
# allowing the main function to be called directly from other scripts without invoking a subprocess. This prevents extra
# subprocesses or extra python interpreters being spawned, as well as simplifying python scripts which call run_dbcan.
def cli_main():
    example_command='''
    example command:
    1. CAZyme annotation with isolated genome sequence as input
    run_dbcan EscheriaColiK12MG1655.fna prok
    2. CAZyme annotation with isolated protein sequence as input
    run_dbcan EscheriaColiK12MG1655.faa protein
    3. CAZyme annotation with meta genome as input
    run_dbcan EscheriaColiK12MG1655.fna meta
    4. CAZyme and CGC annotation with mete genome as input
    run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff
    5. CAZyme, CGC annotation and substrate prediction with mete genome as input
    run_dbcan EscheriaColiK12MG1655.fna meta -c EscheriaColiK12MG1655.gff --cgc_substrate
    '''
    parser = argparse.ArgumentParser(description='dbCAN4 Driver Script')
    parser.add_argument('inputFile', help='User input file. Must be in FASTA format.')
    parser.add_argument('inputType', choices=['protein', 'prok', 'meta'], #protein=proteome, prok=prokaryote nucleotide, meta=metagenome nucleotide
                        help='Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome')
    parser.add_argument('--dbCANFile',default="dbCAN.txt", help='Indicate the file name of HMM database such as dbCAN.txt, please use the newest one from dbCAN2 website.')
    parser.add_argument('--dia_eval', default=1e-102,type=float, help='DIAMOND E Value')
    parser.add_argument('--dia_cpu', default=4, type=int, help='Number of CPU cores that DIAMOND is allowed to use')
    parser.add_argument('--hmm_eval', default=1e-15, type=float, help='HMMER E Value')
    parser.add_argument('--hmm_cov', default=0.35, type=float, help='HMMER Coverage val')
    parser.add_argument('--hmm_cpu', default=4, type=int, help='Number of CPU cores that HMMER is allowed to use')
    parser.add_argument('--out_pre', default="", help='Output files prefix')
    parser.add_argument('--out_dir', default="output", help='Output directory')
    parser.add_argument('--db_dir', default="db", help='Database directory')
    parser.add_argument('--tools', '-t', nargs='+', choices=['hmmer', 'diamond', 'dbcansub', 'all'], default='all', help='Choose a combination of tools to run')
    parser.add_argument('--use_signalP', default=False, type=bool, help='Use signalP or not, remember, you need to setup signalP tool first. Because of signalP license, Docker version does not have signalP.')
    parser.add_argument('--signalP_path', '-sp',default="signalp", type=str, help='The path for signalp. Default location is signalp')
    parser.add_argument('--gram', '-g', choices=["p","n","all"], default="all", help="Choose gram+(p) or gram-(n) for proteome/prokaryote nucleotide, which are params of SingalP, only if user use singalP")
    parser.add_argument('-v', '--version',default="3.0.0", type=str)
    # dbCAN-sub
    dbCAN_sub_group = parser.add_argument_group('dbCAN-sub parameters')
    dbCAN_sub_group.add_argument('--dbcan_thread', '-dt', default=5,type=int, help='number of cpu for dbcan-sub')
    ### cgc finder 
    cgcfinder_group = parser.add_argument_group('CGC_Finder parameters')
    cgcfinder_group.add_argument('--cluster', '-c', help='Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used')
    cgcfinder_group.add_argument('--cgc_dis', default=2, type=int, help='CGCFinder Distance value')
    cgcfinder_group.add_argument('--cgc_sig_genes', default='tp', choices=['tf', 'tp', 'stp', 'tp+tf', 'tp+stp', 'tf+stp', 'all'], help='CGCFinder Signature Genes value')
    cgcfinder_group.add_argument('--tf_eval', default=1e-4, type=float, help='tf.hmm HMMER E Value')
    cgcfinder_group.add_argument('--tf_cov', default=0.35, type=float, help='tf.hmm HMMER Coverage val')
    cgcfinder_group.add_argument('--tf_cpu', default=1, type=int, help='tf.hmm Number of CPU cores that HMMER is allowed to use')
    cgcfinder_group.add_argument('--stp_eval', default=1e-4, type=float, help='stp.hmm HMMER E Value')
    cgcfinder_group.add_argument('--stp_cov', default=0.3, type=float, help='stp.hmm HMMER Coverage val')
    cgcfinder_group.add_argument('--stp_cpu', default=1, type=int, help='stp.hmm Number of CPU cores that HMMER is allowed to use')
    ### cgc substrate prediction 
    cgcsubstrate_group = parser.add_argument_group('CGC_Substrate parameters')
    cgcsubstrate_group.add_argument('--cgc_substrate',action='store_true',help="run cgc substrate prediction?")
    cgcsubstrate_group.add_argument('--pul',help="dbCAN-PUL PUL.faa")
    cgcsubstrate_group.add_argument('-o','--out',default="sub.prediction.out")
    cgcsubstrate_group.add_argument('-w','--workdir',type=str,default=".")
    cgcsubstrate_group.add_argument('-env','--env',type=str,default="local")
    cgcsubstrate_group.add_argument('-oecami','--oecami',action='store_true',help="out eCAMI prediction intermediate result?")
    cgcsubstrate_group.add_argument('-odbcanpul','--odbcanpul',action='store_true',help="output dbCAN-PUL prediction intermediate result?")

    ### cgc substrate prediction:dbCAN-PUL
    group1 = parser.add_argument_group('dbCAN-PUL homologous searching parameters', 'how to define homologous gene hits and PUL hits')
    group1.add_argument('-upghn','--uniq_pul_gene_hit_num',default = 2,type=int)
    group1.add_argument('-uqcgn','--uniq_query_cgc_gene_num',default = 2,type=int)
    group1.add_argument('-cpn','--CAZyme_pair_num',default = 1,type=int)
    group1.add_argument('-tpn','--total_pair_num',default = 2,type=int)
    group1.add_argument('-ept','--extra_pair_type',default = None,type=str,help="None[TC-TC,STP-STP]. Some like sigunature hits")
    group1.add_argument('-eptn','--extra_pair_type_num',default ="0",type=str,help="specify signature pair cutoff.1,2")
    group1.add_argument('-iden','--identity_cutoff',default = 0.3,type=float,help="identity to identify a homologous hit")
    group1.add_argument('-cov','--coverage_cutoff',default = 0.3,type=float,help="query coverage cutoff to identify a homologous hit")
    group1.add_argument('-bsc','--bitscore_cutoff',default = 50,type=float,help="bitscore cutoff to identify a homologous hit")
    group1.add_argument('-evalue','--evalue_cutoff',default = 0.01,type=float,help="evalue cutoff to identify a homologous hit")

    ### cgc substrate prediction:dbCAN-sub
    group2 = parser.add_argument_group('dbCAN-sub major voting parameters', 'how to define dbsub hits and dbCAN-sub subfamily substrate')
    group2.add_argument('-hmmcov','--hmmcov',default = 0.3,type=float)
    group2.add_argument('-hmmevalue','--hmmevalue',default = 0.01,type=float)
    group2.add_argument('-ndsc','--num_of_domains_substrate_cutoff',default = 2,type=int,help="define how many domains share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-npsc','--num_of_protein_substrate_cutoff',default = 2,type=int,help="define how many sequences share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-subs','--substrate_scors',default = 2,type=int,help="each cgc contains with substrate must more than this value")

    args = parser.parse_args()
    
    ### rundbCAN3
    run(inputFile=args.inputFile, inputType=args.inputType, cluster=args.cluster, dbCANFile=args.dbCANFile,
        dia_eval=args.dia_eval, dia_cpu=args.dia_cpu, hmm_eval=args.hmm_eval, hmm_cov=args.hmm_cov,
        hmm_cpu=args.hmm_cpu, dbcan_thread=args.dbcan_thread, tf_eval=args.tf_eval, tf_cov=args.tf_cov, tf_cpu=args.tf_cpu,
        stp_eval=args.stp_eval, stp_cov=args.stp_cov, stp_cpu=args.stp_cpu, prefix=args.out_pre, outDir=args.out_dir,
        dbDir=args.db_dir, cgc_dis=args.cgc_dis, cgc_sig_genes=args.cgc_sig_genes, tool_arg=args.tools,
        use_signalP=args.use_signalP, signalP_path=args.signalP_path, gram=args.gram)
    
    ### convert cgc_standard.out to json format
    
    if args.cluster: ### run cgc_finder
        os.system(f"cgc_standard2json -i {args.out_dir}/cgc_standard.out -o {args.out_dir}/cgc_standard.out.json")
    ### substarate prediction 
    if args.cgc_substrate:
        cgc_substrate_prediction(args)

if __name__ == '__main__':
    cli_main()
