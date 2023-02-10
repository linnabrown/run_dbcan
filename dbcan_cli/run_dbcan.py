#!/usr/bin/env python3
#########################################################
# dbCAN2 Driver Script (Stand Alone Version)
#
# Written by Tanner Yohe in the Yin Lab at NIU
# Revised by Le Huang in the Zhang Lab at NKU
# Updated by Mohamad Majd Raslan in the Yin Lab at NIU
# Updated by Wei Li created table
# Updated by Le Huang at NKU
# Updated by Qiwei Ge in Dr.Yin's Lab at UNL
# Updated by Alex Fraser to allow direct calls to main function from other scripts on 13/06/22
# updated information[Qiwei Ge]: 1. Hotpep has been removed, added eCAMI tool. 2. cgc out reformatting. 3. Fixed issues multiple GT2s.
# Accepts user input
# Predicts genes if needed
# Runs input against HMMER, DIAMOND, and Hotpep
# Optionally predicts CGCs with CGCFinder
# Creats an overview table using output files from core
# tools from Hotpep.out,hmmer.out and diamond.out
##########################################################
from subprocess import Popen, call, check_output
import os
import argparse
import dbcan
from dbcan.utils.simplify_cgc import simplify_output
from dbcan.utils.CGCFinder import cgc_finder
from dbcan.eCAMI import eCAMI_config, eCAMI_main
from dbcan_cli import hmmscan_parser

'''
def some functions
'''


def runHmmScan(outPath, hmm_cpu, dbDir, hmm_eval, hmm_cov, db_name):
    hmmer = Popen(['hmmscan', '--domtblout', '%sh%s.out' % (outPath, db_name), '--cpu', hmm_cpu, '-o', '/dev/null', '%s%s.hmm' % (dbDir,db_name), '%suniInput' % outPath])
    hmmer.wait()
    # call('hmmscan_parser.py %sh%s.out %s %s > %s%s.out'%(outPath, db_name, hmm_eval, hmm_cov, outPath, db_name), shell=True)
    parsed_hmm_output = hmmscan_parser.run(input_file=f"{outPath}h{db_name}.out", eval_num=hmm_eval, coverage=hmm_cov)
    with open(f"{outPath}{db_name}.out", 'w') as f:
        f.write(parsed_hmm_output)

    if os.path.exists('%sh%s.out' % (outPath, db_name)):
        call(['rm', '%sh%s.out' % (outPath, db_name)])


def run(inputFile, inputType, cluster=None, dbCANFile="dbCAN.txt", dia_eval=1e-102, dia_cpu=4, hmm_eval=1e-15,
        hmm_cov=0.35, hmm_cpu=4, eCAMI_kmer_db="CAZyme", eCAMI_k_mer=8, eCAMI_jobs=8, eCAMI_important_k_mer_number=5,
        eCAMI_beta=2, tf_eval=1e-4, tf_cov=0.35, tf_cpu=1, stp_eval=1e-4, stp_cov=0.3, stp_cpu=1, prefix="",
        outDir="output", dbDir="db", cgc_dis=2, cgc_sig_genes="tp", tool_arg="all", use_signalP=False,
        signalP_path="signalp", gram="all"):

    ####
    #  run_dbcan.py [inputFile] [inputType]
    ####

    ##########################
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
        Please make sure that your dbCAN HMM database is named 'dbCAN-HMMdb-V8.txt' or the newest one, has been through hmmpress, and is located in your database directory")
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
    tools = [True, True, True] #DIAMOND, HMMER, eCAMI
    if 'all' not in tool_arg:
        if 'diamond' not in tool_arg:
            tools[0] = False
        if 'hmmer' not in tool_arg:
            tools[1] = False
        if 'eCAMI' not in tool_arg:
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

    if tools[0]:
        # diamond blastp -d db/CAZy -e 1e-102 -q output_EscheriaColiK12MG1655/uniInput -k 1 -p 2 -o output_EscheriaColiK12MG1655/diamond1.out -f 6
        print("\n\n***************************1. DIAMOND start*************************************************\n\n")
        os.system('diamond blastp -d %s -e %s -q %suniInput -k 1 -p %d -o %sdiamond.out -f 6'%(os.path.join(dbDir, "CAZy"), str(dia_eval), outPath, dia_cpu, outPath))
        # diamond = Popen(['diamond', 'blastp', '-d', '%sCAZy.dmnd' % dbDir, '-e', str(args.dia_eval), '-q', '%suniInput' % outPath, '-k', '1', '-p', str(args.dia_cpu), '-o', '%sdiamond.out'%outPath, '-f', '6'])
        print("\n\n***************************1. DIAMOND end***************************************************\n\n")

    if tools[1]:
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
        print("\n\n***************************3. eCAMI start***************************************************\n\n")
        print("Using "+eCAMI_kmer_db+" db in eCAMI")
        ecami_config = eCAMI_config(
            db_type = eCAMI_kmer_db,
            input = f"{outPath}uniInput",
            output= f"{outPath}eCAMI.out",
            k_mer = eCAMI_k_mer,
            jobs = eCAMI_jobs,
            important_k_mer_number = eCAMI_important_k_mer_number,
            beta = eCAMI_beta
        )
        eCAMI_main(ecami_config)
        # os.system('python eCAMI/prediction.py -input %suniInput -kmer_db eCAMI/%s -output %seCAMI.out -k_mer %s -jobs %s -important_k_mer_number %s -beta %s' % (outPath, str(args.eCAMI_kmer_db), outPath, str(args.eCAMI_k_mer), str(args.eCAMI_jobs), str(args.eCAMI_important_k_mer_number),str(args.eCAMI_beta)))

        print("\n\n***************************3. eCAMI end***************************************************\n\n")
    # End Core Tools
    ########################
    # Begin Adding Column Headers

    if tools[2]:
        with open(outPath+'eCAMI.out') as f:
            with open(outPath+'temp', 'w') as out:
                out.write('protein_name\tfam_name:group_number\tsubfam_name_of_the_group:subfam_name_count\n')
                # for line in f:
                for count, line in enumerate(f):
                    if count % 2 == 0:
                        more_information = line.split(">")
                        out.write(more_information[1])
        call(['mv', outPath+'temp', outPath+'eCAMI.out'])

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

    if tools[0]:
        with open(outDir+prefix+'diamond.out') as f:
            with open(outDir+prefix+'temp', 'w') as out:
                out.write('Gene ID\tCAZy ID\t% Identical\tLength\tMismatches\tGap Open\tGene Start\tGene End\tCAZy Start\tCAZy End\tE Value\tBit Score\n')
                for line in f:
                    out.write(line)
        call(['mv', outDir+prefix+'temp', outDir+prefix+'diamond.out'])

    # End Adding Column Headers
    ########################
    # Begin CGCFinder

    if find_clusters:
        print("*****************************CGC-Finder start************************************")

        ########################
        # Begin TF,TP, STP prediction
        '''
        previous tf uses diamond, Now tf-1 and tf-2 uses hmmer
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
        # hot = set()
        hmm = set()
        eca = set()
        if tools[0]:
            with open(outDir+prefix+'diamond.out') as f:
                next(f)
                for line in f:
                    row = line.rstrip().split('\t')
                    dia.add(row[0])
                    if row[0] not in cazyme_genes:
                        cazyme_genes[row[0]] = set()
                    cazyme_genes[row[0]].update(set(row[1].strip("|").split('|')[1:]))
        if tools[1]:
            with open(outDir+prefix+'hmmer.out') as f:
                next(f)
                for line in f:
                    row = line.rstrip().split('\t')
                    hmm.add(row[2])
                    if row[2] not in cazyme_genes:
                        cazyme_genes[row[2]] = set()
                    cazyme_genes[row[2]].add(row[0].split('.hmm')[0])
        if tools[2]:
            with open(outDir+prefix+'eCAMI.out') as f:
                next(f)
                for line in f:
                    row_ori = line.rstrip().split('\t')
                    if ' ' in row_ori[0]:
                        fams_ID = row_ori[0].split(' ')[0]
                    else:
                        fams_ID = row_ori[0]
                    eca.add(fams_ID)
                    if fams_ID not in cazyme_genes:
                        cazyme_genes[fams_ID] = set()
                    cazyme_genes[fams_ID].add(row_ori[1].split(':')[0])

        if tools.count(True) > 1:
            temp1 = hmm.intersection(eca)
            # print(hmm, 'This intersection  hmm')
            temp2 = hmm.intersection(dia)
            # print(dia, 'This intersection  dia')
            temp3 = dia.intersection(eca)
            # print(eca, 'This intersection  eca')
            cazyme = temp1.union(temp2, temp3)
        else:
            cazyme = hmm.union(dia, eca)
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
                                # row[8] = "DB="+'|'.join(cazyme_genes[gene])
                                #
                                cazyme_genes_list = list(cazyme_genes[gene])
                                row[8] = "DB="+cazyme_genes_list[0]
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
                                        # row[8] = "DB="+'|'.join(cazyme_genes[gene])
                                        #
                                        cazyme_genes_list = list(cazyme_genes[gene])
                                        row[8] = "DB="+cazyme_genes_list[0]
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
                            outrow = ['.','.','.','.','.','.','.','.','']
                            gene = row[1]
                            if gene in cazyme:
                                outrow[2] = 'CAZyme'
                                # Uncomment this, if all CAZyme results need to be write into cgc.out
                                # outrow[8] = "DB="+'|'.join(cazyme_genes[gene])
                                #
                                cazyme_genes_list = list(cazyme_genes[gene])
                                outrow[8] = "DB="+cazyme_genes_list[0]
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

        # call(['CGCFinder.py', outDir+prefix+'cgc.gff', '-o', outDir+prefix+'cgc.out', '-s', args.cgc_sig_genes, '-d', str(args.cgc_dis)])
        cgc_finder(outDir+prefix+'cgc.gff', cgc_dis, cgc_sig_genes, outDir+prefix+'cgc.out')
        simplify_output(outDir+prefix+'cgc.out')
        print("**************************************CGC-Finder end***********************************************")
        # End CGCFinder call
        # End CGCFinder
        ####################
    # Begin SignalP combination
    if use_signalP:
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
    print("Preparing overview table from hmmer, eCAMI and diamond output...")
    workdir = outDir+prefix
    # a function to remove duplicates from lists while keeping original order
    def unique(seq):
        exists = set()
        return [x for x in seq if not (x in exists or exists.add(x))]

    arr_eCAMI = None
    arr_hmmer = None

    # check if files exist. if so, read files and get the gene numbers
    if tools[0]:
        arr_diamond = open(workdir+"diamond.out").readlines()
        diamond_genes = [arr_diamond[i].split()[0] for i in range(1, len(arr_diamond))] # or diamond_genes = []

    if tools[1]:
        arr_hmmer = open(workdir+"hmmer.out").readlines()
        hmmer_genes = [arr_hmmer[i].split()[2] for i in range(1, len(arr_hmmer))] # or hmmer_genes = []

    if tools[2]:
        arr_eCAMI = open(workdir+"eCAMI.out").readlines()
        eCAMI_genes = [arr_eCAMI[i].split()[0] for i in range(1, len(arr_eCAMI))]# or eCAMI_genes = []

    if use_signalP and (os.path.exists(workdir + "signalp.out")):
        arr_sigp = open(workdir+"signalp.out").readlines()
        sigp_genes = {}
        for i in range (0,len(arr_sigp)):
            row = arr_sigp[i].split()
            sigp_genes[row[0]] = row[4] #previous one is row[2], use Y-score instead from suggestion of Dongyao Li

    ##Catie Ausland edits BEGIN, Le add variable exists or not, remove duplicates from input lists
    if not tools[0]:
        diamond_genes =[]
    if not tools[1]:
        hmmer_genes   = []
    if not tools[2]:
        eCAMI_genes  =[]

    if len(eCAMI_genes) > 0:
        if (eCAMI_genes[-1] == None):
            #print('I am in &&&&&&&&&&&&&&&&&&&&&&')
            eCAMI_genes.pop()
            eCAMI_genes = unique(eCAMI_genes)
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


    if tools[2] and (len(arr_eCAMI) > 1) :
        eCAMI_fams = {}
        for i in range (1,len(arr_eCAMI)):
            row_ori = arr_eCAMI[i].split("\t")
            subfam_names = row_ori[2].split('|')
            fam = row_ori[1].split(':')
            if ' ' in row_ori[0]:
                fams_ID = row_ori[0].split(' ')[0]
            else:
                fams_ID = row_ori[0]

            diam_name = []
            for name in subfam_names:
                if '.' in name:
                    diam_name.append(name.split(":")[0])

            if(fams_ID not in eCAMI_fams):
                eCAMI_fams[fams_ID] = {}
                eCAMI_fams[fams_ID]["fam_name"] = []
                eCAMI_fams[fams_ID]["ec_num"] = []

            eCAMI_fams[fams_ID]["fam_name"].append(fam[0])
            eCAMI_fams[fams_ID]["ec_num"] = diam_name

    #overall table

    all_genes = unique(hmmer_genes+eCAMI_genes+diamond_genes)

    with open(workdir+"overview.txt", 'w+') as fp:
        if use_signalP:
            fp.write("Gene ID\tEC#\tHMMER\teCAMI\tDIAMOND\tSignalp\t#ofTools\n")
        else:
            fp.write("Gene ID\tEC#\tHMMER\teCAMI\tDIAMOND\t#ofTools\n")
        for gene in all_genes:
            csv=[gene]
            num_tools = 0

            if tools[2] and arr_eCAMI != None and (gene in eCAMI_genes):
                if eCAMI_fams[gene]["ec_num"] == []:
                    csv.append("-")
                else:
                    csv.append("|".join(eCAMI_fams[gene]["ec_num"]))
            else:
                csv.append("-")

            if tools[1] and arr_hmmer != None and (gene in hmmer_genes):
                num_tools += 1
                csv.append("+".join(hmmer_fams[gene]))
            else:
                csv.append("-")

            if tools[2] and arr_eCAMI != None and (gene in eCAMI_genes):
                num_tools += 1
                csv.append("+".join(eCAMI_fams[gene]["fam_name"]))
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
    parser = argparse.ArgumentParser(description='dbCAN2 Driver Script')
    parser.add_argument('inputFile', help='User input file. Must be in FASTA format.')
    parser.add_argument('inputType', choices=['protein', 'prok', 'meta'], #protein=proteome, prok=prokaryote nucleotide, meta=metagenome nucleotide
                        help='Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome')
    parser.add_argument('--cluster', '-c', help='Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used')
    parser.add_argument('--dbCANFile',default="dbCAN.txt", help='Indicate the file name of HMM database such as dbCAN.txt, please use the newest one from dbCAN2 website.')
    parser.add_argument('--dia_eval', default=1e-102,type=float, help='DIAMOND E Value')
    parser.add_argument('--dia_cpu', default=4, type=int, help='Number of CPU cores that DIAMOND is allowed to use')
    parser.add_argument('--hmm_eval', default=1e-15, type=float, help='HMMER E Value')
    parser.add_argument('--hmm_cov', default=0.35, type=float, help='HMMER Coverage val')
    parser.add_argument('--hmm_cpu', default=4, type=int, help='Number of CPU cores that HMMER is allowed to use')
    # eCAMI
    parser.add_argument('--eCAMI_kmer_db', default="CAZyme",type=str, help="Change n_mer directories path for prediction")
    parser.add_argument('--eCAMI_k_mer', default=8, type=int, help="Peptide length for prediction")
    parser.add_argument('--eCAMI_jobs', default=8, type=int, help='Number of processor for use for prediction')
    parser.add_argument('--eCAMI_important_k_mer_number', default=5, type=int, help="Minimum number of n_mer for prediction")
    parser.add_argument('--eCAMI_beta', default=2, type=float, help="Minimum sum of percentage of frequency of n_mer for prediction")
    # eCAMI
    parser.add_argument('--tf_eval', default=1e-4, type=float, help='tf.hmm HMMER E Value')
    parser.add_argument('--tf_cov', default=0.35, type=float, help='tf.hmm HMMER Coverage val')
    parser.add_argument('--tf_cpu', default=1, type=int, help='tf.hmm Number of CPU cores that HMMER is allowed to use')
    parser.add_argument('--stp_eval', default=1e-4, type=float, help='stp.hmm HMMER E Value')
    parser.add_argument('--stp_cov', default=0.3, type=float, help='stp.hmm HMMER Coverage val')
    parser.add_argument('--stp_cpu', default=1, type=int, help='stp.hmm Number of CPU cores that HMMER is allowed to use')
    parser.add_argument('--out_pre', default="", help='Output files prefix')
    parser.add_argument('--out_dir', default="output", help='Output directory')
    parser.add_argument('--db_dir', default="db", help='Database directory')
    parser.add_argument('--cgc_dis', default=2, type=int, help='CGCFinder Distance value')
    parser.add_argument('--cgc_sig_genes', default='tp', choices=['tf', 'tp', 'stp', 'tp+tf', 'tp+stp', 'tf+stp', 'all'], help='CGCFinder Signature Genes value')
    parser.add_argument('--tools', '-t', nargs='+', choices=['hmmer', 'diamond', 'eCAMI', 'all'], default='all', help='Choose a combination of tools to run')
    parser.add_argument('--use_signalP', default=False, type=bool, help='Use signalP or not, remember, you need to setup signalP tool first. Because of signalP license, Docker version does not have signalP.')
    parser.add_argument('--signalP_path', '-sp',default="signalp", type=str, help='The path for signalp. Default location is signalp')
    parser.add_argument('--gram', '-g', choices=["p","n","euk","all"], default="all", help="Choose gram+(p) or gram-(n) for proteome/prokaryote nucleotide, or euk(euk) for which are params of signalP, only if user use signalP")
    args = parser.parse_args()

    run(inputFile=args.inputFile, inputType=args.inputType, cluster=args.cluster, dbCANFile=args.dbCANFile,
        dia_eval=args.dia_eval, dia_cpu=args.dia_cpu, hmm_eval=args.hmm_eval, hmm_cov=args.hmm_cov,
        hmm_cpu=args.hmm_cpu, eCAMI_kmer_db=args.eCAMI_kmer_db, eCAMI_k_mer=args.eCAMI_k_mer,
        eCAMI_jobs=args.eCAMI_jobs, eCAMI_important_k_mer_number=args.eCAMI_important_k_mer_number,
        eCAMI_beta=args.eCAMI_beta, tf_eval=args.tf_eval, tf_cov=args.tf_cov, tf_cpu=args.tf_cpu,
        stp_eval=args.stp_eval, stp_cov=args.stp_cov, stp_cpu=args.stp_cpu, prefix=args.out_pre, outDir=args.out_dir,
        dbDir=args.db_dir, cgc_dis=args.cgc_dis, cgc_sig_genes=args.cgc_sig_genes, tool_arg=args.tools,
        use_signalP=args.use_signalP, signalP_path=args.signalP_path, gram=args.gram)


if __name__ == '__main__':
    cli_main()
