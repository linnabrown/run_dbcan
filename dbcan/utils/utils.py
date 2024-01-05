from Bio import SeqIO
from itertools import (takewhile, repeat)
import time
from subprocess import Popen, call, check_output
from dbcan.cli import hmmer_parser
import argparse,os
import numpy as np
'''
design some functions and classes for dbCAN tutorial

'''

def fq_file_line_count(file_name):
    if not file_name.endswith(".gz"): ### not gz files
        buffer = 1024 * 1024
        with open(file_name) as f:
            buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
            return sum(buf.count('\n') for buf in buf_gen)/4
    else:
        r = os.popen("zcat " + file_name + " | echo $((`wc -l`/4))")
        text = r.read()
        r.close()
        return int(text)

def total_mapped_reads_count(file_name):
    ### sum the total mapped reads count
    total_mapped_reads = 0 
    for line in open(file_name):
        lines = line.split()
        total_mapped_reads += int(lines[-1])
    return total_mapped_reads

class abund_parameters():
    def __init__(self,args):
        self.input = args.input if args.input.endswith("/") else args.input +"/"
        self.R1 = args.R1
        self.R2 = args.R2
        self.bedtools = args.bedtools
        #self.output = args.function + "_" + args.output
        self.CAZyme_annotation  = self.input + "overview.txt"
        self.dbCANsub_substrate_annotation  = self.input + "dbcan-sub.hmm.out"
        self.PUL_substrate_annotation  = self.input + "substrate.out"
        self.PUL_annotation  = self.input + "cgc_standard.out"
        self.function = args.function
        self.parameters_check()
    
    def parameters_check(self):
        if self.function == "fam_abund":
            print("You are estimating the abundance of CAZyme!")
            self.output = "fam_abund.out"
            if not os.path.exists(self.CAZyme_annotation):
                print(f"CAZyme annotation file {self.CAZyme_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
        if self.function == "fam_substrate_abund":
            self.output = "fam_substrate_abund.out"
            print("You are estimating the abundance of Substrate according to dbCAN-sub!")
            if not os.path.exists(self.dbCANsub_substrate_annotation):
                print(f"CAZyme annotation file {self.dbCANsub_substrate_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
        if self.function == "CGC_abund":
            self.output = "CGC_abund.out"
            print("You are estimating the abundance of CGC!")
            if not os.path.exists(self.PUL_annotation):
                print(f"PUL annotation file {self.PUL_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
        if self.function == "CGC_substrate_abund":
            print("You are estimating the abundance of CGC substrate!")
            if not os.path.exists(self.PUL_substrate_annotation):
                print(f"CGC substrate prediction file {self.PUL_substrate_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
            if not os.path.exists(self.PUL_annotation):
                print(f"CGC annotation file {self.PUL_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
        if self.R2:
            self.ngs_pe = True
            print("Reads are pair end!")
        else:
            self.ngs_pe = False
            print("Reads are single end!")
        if not os.path.exists(self.bedtools):
            print(f"Reads count file {self.bedtools} dose not exit!")
            exit(1)

class bedtools_read_count():
    def __init__(self,lines):
        self.seqid = lines[0]
        self.length = int(lines[2])
        self.read_count = int(lines[3])
    def __repr__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self)])

def ReadBedtoos(filename):
    lines = open(filename).readlines()
    seqid2info = {line.split()[0]:bedtools_read_count(line.split()) for line in lines[1:]}
    normalized_tpm = 0.
    for seqid in seqid2info:
        seqid_depth = seqid2info[seqid]
        normalized_tpm += seqid_depth.read_count/seqid_depth.length
    return seqid2info,normalized_tpm

#Gene ID	EC#	HMMER	dbCAN_sub	DIAMOND	#ofTools
#fefifo_8002_1_00016	-	GH25(447-633)	GH25_e123	GH25	3
#fefifo_8002_1_00040	3.2.1.-:3|3.2.1.63:2	GH95(1-376)	GH95_e1	GH95	3

def Is_EC(ec):
    return ec[0].isdigit()

def Clean_Hmmer_sub(pred):
    preds = pred.split("+")
    clean_preds = []
    for f in preds:
        f = f.split("(")[0]
        clean_preds.append(f)
        #if "_" in f: ### also include subfamily
        #    f_s = f.split("_")[0] ## get family
        #    clean_preds.append(f_s)
    return list(set(clean_preds))

def Clean_diamond(pred):
    preds = pred.split("+")
    clean_preds = []
    for f in preds:
        if Is_EC(f): ### exclude EC
            continue
        f = f.split("(")[0]
        clean_preds.append(f)
        #if "_" in f: ### also include subfamily
        #    f_s = f.split("_")[0] ## get family
        #    clean_preds.append(f_s)
    return list(set(clean_preds))

class OverView():
    def __init__(self,lines):
        self.seqid = lines[0]
        self.ECs = list(set([i.split(":")[0] for i in lines[1].split("|") if i.split(":")[0] != '-']))
        self.hmmer = Clean_Hmmer_sub(lines[2]) if lines[2] != "-" else []
        self.dbcan_sub = Clean_Hmmer_sub(lines[3]) if lines[3] != "-" else []
        self.diamond = Clean_diamond(lines[4]) if lines[4] != "-" else []
        self.justify_final_pred()
    def __repr__(self):
        outline = []
        if not self.preds:
            return ""
        for value in vars(self):
            attr = getattr(self, value)
            if isinstance(attr,list):
                outline.append("+".join(attr))
            else:
                outline.append(str(attr))
        return "\t".join(outline)
    ### dbcan_sub is for substrate prediction
    def justify_final_pred(self):
        union_fam = list(set(self.hmmer)|set(self.dbcan_sub)) + self.ECs ## union
        if union_fam:
            self.preds = union_fam
        else:
            self.preds = []

def ReadOverView(filename):
    return {line.split()[0]:OverView(line.split()) for line in open(filename).readlines()[1:]}

#dbCAN subfam	Subfam Composition	Subfam EC	Substrate	Profile Length	Gene ID	Gene Length	E Value	Profile Start	Profile End	Gene Start	Gene End	Coverage
#GH25_e123	GH25:8	-	-	181	fefifo_8002_1_00016	644	7.5e-68	1	181	448	634	0.994475138121547	0.994475138121547

class DbcanSub_line():
    def __init__(self,lines):
        self.dbcan_subfam = lines[0]
        self.subfam_comp = lines[1]
        self.subfam_EC = lines[2]
        self.substrate = lines[3]
        self.hmmlen  = lines[4]
        self.seqid = lines[5]
        self.protlen = lines[6] 
        self.evalue = lines[7]
    def __repr__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self)])

def Read_dbcansub_out(filename):
    return {line.split("\t")[5]:DbcanSub_line(line.rstrip().split("\t")) for line in open(filename).readlines()[1:]}

#CGC#	Gene Type	Contig ID	Protein ID	Gene Start	Gene Stop	Direction	Protein Family
#CGC1	TC	k141_70330	fefifo_8002_1_00592	100	1263	-	2.A.1.7.16

class cgc_standard_line():
    def __init__(self,lines):
        self.cgcid = lines[2] + "|" + lines[0]
        self.cgc_order = lines[0]
        self.gene_type = lines[1] ### CAZyme
        self.contig_id = lines[2]
        self.seqid = lines[3]
        self.gene_start = int(lines[4])
        self.gene_end = int(lines[5])
        self.strand = lines[6]
        self.protfam = lines[7]
    def __repr__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self)])

def Read_cgc_standard_out(filename):
    seqid2records = {};cgcid2records ={}
    for line in open(filename).readlines()[1:]:
        lines = line.rstrip("\n").split("\t")
        cgcid = lines[2] + "|" + lines[0]
        tmp_record = cgc_standard_line(line.rstrip().split("\t"))
        seqid2records[lines[3]] = tmp_record
        cgcid2records.setdefault(cgcid,[]).append(tmp_record)
    return seqid2records,cgcid2records

## the prediction of cgc substrate includes two results, homologous searach and major votting
#cgcid	PULID	dbCAN-PUL substrate	bitscore	signature pairs	dbCAN-sub substrate	dbCAN-sub substrate score
#k141_145965|CGC1	PUL0538	galactomannan	3075.0	TC-TC;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;CAZyme-CAZyme;TC-TC		

class cgc_substrate():
    def __init__(self,lines):
        self.cgcid = lines[0]
        self.homo_pul = lines[1]
        self.homo_sub = lines[2]
        self.bitscore = lines[3]
        self.signature_pairs = lines[4]
        self.major_voting_sub = lines[5]
        self.major_voting_score = lines[6]
    def __repr__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self)])

def Read_cgc_substrate(filename):
    return {line.split()[0]:cgc_substrate(line.rstrip("\n").split("\t")) for line in open(filename).readlines()[1:]}

def get_length_readcount(seqid2dbcan_annotation,seqid2readcount):
    for seqid in seqid2dbcan_annotation:
        read_count = seqid2readcount.get(seqid,"")
        if not read_count:
            print(f"Can not find read count information for CAZyme: {seqid}")
            exit(1)
        seqid_annotation = seqid2dbcan_annotation[seqid]
        seqid_annotation.length = read_count.length
        seqid_annotation.read_count = read_count.read_count

class CAZyme_Abundance_estimate():
    def __init__(self,parameters):
        self.pars = parameters
        #print("Counting reads!")
        #self.fq_reads_count = fq_file_line_count(self.pars.R1)
        self.fq_reads_count = total_mapped_reads_count(self.pars.bedtools)
        print(f"Total reads count: {self.fq_reads_count}!")
        seqid2readcount,normalized_tpm = ReadBedtoos(parameters.bedtools)
        self.normalized_tpm = normalized_tpm
        ### read overview to 
        if parameters.function == "fam_abund":
            seqid2dbcan_annotation = ReadOverView(parameters.CAZyme_annotation)
        ### read dbsub.out to 
        if parameters.function == "fam_substrate_abund":
            seqid2dbcan_annotation = Read_dbcansub_out(parameters.dbCANsub_substrate_annotation)
        ### read cgc_standard.out
        if parameters.function == "CGC_abund": ### seqid = cgcid
            seqid2dbcan_annotation,cgcid2cgc_standard = Read_cgc_standard_out(parameters.PUL_annotation)
            self.cgcid2cgc_standard = cgcid2cgc_standard
        ### read cgc_standard.out and sub.prediction.out
        if parameters.function == "CGC_substrate_abund": ### seqid = cgcid
            seqid2dbcan_annotation,cgcid2cgc_standard = Read_cgc_standard_out(parameters.PUL_annotation)
            cgcid2cgc_substrate = Read_cgc_substrate(parameters.PUL_substrate_annotation)
            self.cgcid2cgc_standard = cgcid2cgc_standard
            self.cgcid2cgc_substrate = cgcid2cgc_substrate
        
        get_length_readcount(seqid2dbcan_annotation,seqid2readcount)
        self.seqid2dbcan_annotation = seqid2dbcan_annotation
    
    ###             reads per transcript
    ### FPKM  =  -------------------------------
    ###          total reads   transcript length
    ###          ----------- X -----------------
    ###          10E6          1000

    def Cal_Seq_Abundance(self,method="FPKM"):
        if method == "FPKM":
            for seqid in self.seqid2dbcan_annotation: ### for each protein
                annotation = self.seqid2dbcan_annotation[seqid]
                normalized_total_reads_counts = self.fq_reads_count/pow(10,6)
                normalized_seq_length = annotation.length/1000.0 ## CDS length 
                annotation.abund = annotation.read_count/normalized_total_reads_counts/normalized_seq_length
        if method == "RPM":
            for seqid in self.seqid2dbcan_annotation:
                annotation = self.seqid2dbcan_annotation[seqid]
                normalized_total_reads_counts = self.fq_reads_count/pow(10,6)
                annotation.abund = annotation.read_count/normalized_total_reads_counts
        if method  == "TPM":
            for seqid in self.seqid2dbcan_annotation:
                annotation = self.seqid2dbcan_annotation[seqid]
                normalized_total_reads_counts = annotation.read_count/annotation.length*pow(10,6)
                annotation.abund = normalized_total_reads_counts/self.normalized_tpm

    def Cal_Famliy_Abundance(self):
        family2seqid = {}
        for seqid in self.seqid2dbcan_annotation:
            annotation = self.seqid2dbcan_annotation[seqid]
            if annotation.preds: ### if dbcan has predecitions
                for family in annotation.preds:
                    family2seqid.setdefault(family,[]).append(seqid) 
        self.family2seqid = family2seqid
        family2abund = {familyid:0.0 for familyid in family2seqid}
        for familyid in family2seqid:
            for seqid in family2seqid[familyid]:
                family2abund[familyid] += self.seqid2dbcan_annotation[seqid].abund  ### sum of all seqs
        self.family2abund = family2abund 
    
    def Cal_Substrate_Abundance(self):
        substrate2seqid = {}
        for seqid in self.seqid2dbcan_annotation:
            annotation = self.seqid2dbcan_annotation[seqid]
            substrates_tmp = annotation.substrate.replace("and",",").split(",")
            substrates = list(set([tmp.strip() for tmp in substrates_tmp if tmp != "-" and tmp]))
            for sub in substrates:
                substrate2seqid.setdefault(sub,[]).append(seqid)
        
        substrate2abund = {sub:0.0 for sub in substrate2seqid}
        for sub in substrate2seqid:
            for seqid in substrate2seqid[sub]:
                substrate2abund[sub] += self.seqid2dbcan_annotation[seqid].abund
        self.substrate2abund = substrate2abund
        self.substrate2seqid = substrate2seqid
    
    def Cal_PUL_Abundance(self):
        cgcid2seqid = {}; cgcid2_standard_cgc_lines = {};
        for seqid in self.seqid2dbcan_annotation:
            annotation = self.seqid2dbcan_annotation[seqid]
            cgcid = annotation.cgcid
            cgcid2seqid.setdefault(cgcid,[]).append(seqid)
        cgcid2abund = {cgcid:0.0 for cgcid in cgcid2seqid}
        cgcid2seqabund = {}
        for cgcid in cgcid2seqid:
            for seqid in cgcid2seqid[cgcid]:
                #cgcid2abund[cgcid] += self.seqid2dbcan_annotation[seqid].abund
                cgcid2seqabund.setdefault(cgcid,[]).append(self.seqid2dbcan_annotation[seqid].abund)

        ### calculate mean abund
        for cgcid in cgcid2seqid:
            cgcid2abund[cgcid] = np.mean(cgcid2seqabund[cgcid])

        self.cgcid2abund = cgcid2abund
        self.cgcid2seqid = cgcid2seqid
        self.cgcid2seqabund = cgcid2seqabund
    
    ### self.seqid2dbcan_annotation,self.cgcid2cgc_standard
    ### self.cgcid2seqid, cgcid2seqabund
    ### self.cgcid2cgc_substrate

    def Cal_PUL_Substrate_Abundance(self):
        ''' two substates prediction'''
        cgcsubstrate2cgcid_homo = {}; cgcsubstrate2cgcid_major_votting = {}
        for cgcid in self.cgcid2cgc_substrate:
            cgc_substrate_line = self.cgcid2cgc_substrate[cgcid]
            if cgc_substrate_line.homo_sub and cgc_substrate_line.homo_sub != "X":### homo has prediction 
                cgcsubstrate2cgcid_homo.setdefault(cgc_substrate_line.homo_sub,[]).append(cgcid)
            
            if cgc_substrate_line.major_voting_sub: ### homo has prediction
                substrates = cgc_substrate_line.major_voting_sub.split(",")
                for tmp_sub in substrates:
                    cgcsubstrate2cgcid_major_votting.setdefault(tmp_sub,[]).append(cgcid)
        
        ### for homologous search substrate 
        cgcsubstrate2abunds_homo = {}; self.cgcsubstrate2cgcid_homo = cgcsubstrate2cgcid_homo
        for substrate in cgcsubstrate2cgcid_homo:
            cgcids = cgcsubstrate2cgcid_homo[substrate]
            for cgcid in cgcids:
                cgc_abunds = self.cgcid2seqabund[cgcid]  ### list constis of sequence abundance
                #cgcsubstrate2abunds_homo.setdefault(substrate,[]).extend(cgc_abunds)
                cgcsubstrate2abunds_homo.setdefault(substrate,[]).append(np.mean(cgc_abunds))
        self.cgcsubstrate2abunds_homo = cgcsubstrate2abunds_homo
        ### for major votting substrate
        cgcsubstrate2abunds_major_votting = {}; 
        self.cgcsubstrate2cgcid_major_votting = cgcsubstrate2cgcid_major_votting
        for substrate in cgcsubstrate2cgcid_major_votting:
            cgcids = cgcsubstrate2cgcid_major_votting[substrate]
            for cgcid in cgcids:
                cgc_abunds = self.cgcid2seqabund[cgcid]  ### list constis of sequence abundance
                #cgcsubstrate2abunds_major_votting.setdefault(substrate,[]).extend(cgc_abunds)
                cgcsubstrate2abunds_major_votting.setdefault(substrate,[]).append(np.mean(cgc_abunds))
        self.cgcsubstrate2abunds_major_votting = cgcsubstrate2abunds_major_votting
    
    def output_cgcsubstrate_abund(self):
        ### for cgc substrate homologous
        cgc_substrates = []; abunds = [] ; cgcids = [] ; cgcids_abund = [] # list
        for cgc_substrate in self.cgcsubstrate2abunds_homo:
            abunds.append(np.sum(self.cgcsubstrate2abunds_homo[cgc_substrate]))
            cgc_substrates.append(cgc_substrate)
            cgcids.append(self.cgcsubstrate2cgcid_homo[cgc_substrate])
            cgcids_abund.append(self.cgcsubstrate2abunds_homo[cgc_substrate])
        abund_sortidx = np.argsort(abunds)[::-1]
        print(f"Writring abundance of substrate predicted on dbCAN-PUL to file CGC_substrate_PUL_homology.out!")
        with open("CGC_substrate_PUL_homology.out",'w') as f:
            f.write("Substrate\tAbundance(sum of CGC)\tcgcs\tcgcs_abunds\n")
            for idx in abund_sortidx:
                #if cgc_substrates[idx] == "pectin":
                #    print(len(cgcids[idx]),len(cgcids_abund[idx]))
                cgc = ";".join(cgcids[idx])
                abunds_tmp = ";".join([str(round(abund,3)) for abund in cgcids_abund[idx]])
                f.write(f"{cgc_substrates[idx]}\t{round(abunds[idx],3)}\t{cgc}\t{abunds_tmp}\n")
        
        ### for cgc substrate major votting
        cgc_substrates = []; abunds = [] ; cgcids = []; cgcids_abund = [] # list
        for cgc_substrate in self.cgcsubstrate2abunds_major_votting:
            abunds.append(np.sum(self.cgcsubstrate2abunds_major_votting[cgc_substrate]))
            cgc_substrates.append(cgc_substrate)
            cgcids.append(self.cgcsubstrate2cgcid_major_votting[cgc_substrate])
            cgcids_abund.append(self.cgcsubstrate2abunds_major_votting[cgc_substrate])
        abund_sortidx = np.argsort(abunds)[::-1]
        print(f"Writring abundance of substrate predicted by major votting to file CGC_substrate_majority_voting.out!")
        with open("CGC_substrate_majority_voting.out",'w') as f:
            f.write("Substrate\tAbundance(sum of CGC)\tcgcs\tcgcs_abunds\n")
            for idx in abund_sortidx:
                cgc = ";".join(cgcids[idx])
                abunds_tmp = ";".join([str(round(abund,3)) for abund in cgcids_abund[idx]])
                f.write(f"{cgc_substrates[idx]}\t{round(abunds[idx],3)}\t{cgc}\t{abunds_tmp}\n")
    

    ### need to consider HMM model, subfamily and EC
    def output_family_abund(self,method="family"):
        fams = [] ; abunds = [];seqs = []
        for familyid in self.family2abund:
            fams.append(familyid)
            abunds.append(self.family2abund[familyid])
            seqs.append(self.family2seqid[familyid])
        abund_sortidx = np.argsort(abunds)[::-1]
        print(f"Writring family abundance to file fam_abund.out!")
        print(f"Writring subfamily(dbsub) abundance to file subfam_abund.out!")
        print(f"Writring EC abundance to file EC_abund.out!")
        fam_file = open(self.pars.output,'w')
        subfam_file = open("subfam_abund.out",'w')
        EC_file = open("EC_abund.out",'w')

        fam_file.write("Family\tAbundance\tSeqNum\n")
        subfam_file.write("Subfamily\tAbundance\tSeqNum\n")
        EC_file.write("EC\tAbundance\tSeqNum\n")

        for idx in abund_sortidx:
            famid = fams[idx]
            if Is_EC(famid):
                EC_file.write(f"{fams[idx]}\t{round(abunds[idx],3)}\t{len(seqs[idx])}\n")
            elif "_e" in famid:
                subfam_file.write(f"{fams[idx]}\t{round(abunds[idx],3)}\t{len(seqs[idx])}\n")
            else:
                fam_file.write(f"{fams[idx]}\t{round(abunds[idx],3)}\t{len(seqs[idx])}\n")
    
    def output_substrate_abund(self):
        subs = [] ; abunds = [] ; genes = []
        for sub in self.substrate2abund:
            subs.append(sub)
            abunds.append(self.substrate2abund[sub])
            genes.append(self.substrate2seqid[sub])
        abund_sortidx = np.argsort(abunds)[::-1]
        print(f"Writring substrate abundance to file {self.pars.output}!")
        with open(self.pars.output,'w') as f:
            f.write("Substrate\tAbundance\tGeneID\n")
            for idx in abund_sortidx:
                if subs[idx]:
                    f.write(f"{subs[idx]}\t{round(abunds[idx],3)}\t{';'.join(genes[idx])}\n")
    
    def output_cgc_abund(self):
        cgcids = [] ; abunds = []; seqids = []; seq_abunds = []
        cgc_standard_records = []
        for cgcid in self.cgcid2abund:
            cgcids.append(cgcid)
            abunds.append(self.cgcid2abund[cgcid])
            seqids.append(self.cgcid2seqid[cgcid])
            seq_abunds.append(self.cgcid2seqabund[cgcid])
            cgc_standard_records.append(self.cgcid2cgc_standard[cgcid])

        abund_sortidx = np.argsort(abunds)[::-1]
        print(f"Writring CGC abundance to file {self.pars.output}!")
        with open(self.pars.output,'w') as f:
            f.write("#CGCID\tAbundance(mean)\tSeqid\tSeq_abund\tFams\n")
            for idx in abund_sortidx:
                seqs_tmp_abunds = ";".join([str(round(i,3)) for i in seq_abunds[idx]])
                seqs_tmp = ";".join(seqids[idx])
                fams = ";".join(record.protfam if record.gene_type== "CAZyme" else record.gene_type for record in cgc_standard_records[idx])
                f.write(f"{cgcids[idx]}\t{round(abunds[idx],3)}\t{seqs_tmp}\t{seqs_tmp_abunds}\t{fams}\n")

def CAZyme_abundance(args):
    paras = abund_parameters(args)
    CAZyme_abund = CAZyme_Abundance_estimate(paras)
    CAZyme_abund.Cal_Seq_Abundance(args.abundance) ### calculate abundance 
    CAZyme_abund.Cal_Famliy_Abundance() ### calculate abundance 
    CAZyme_abund.output_family_abund()

def CAZymeSub_abundance(args):
    paras = abund_parameters(args)
    CAZymeSub_abund = CAZyme_Abundance_estimate(paras)
    CAZymeSub_abund.Cal_Seq_Abundance(args.abundance)
    CAZymeSub_abund.Cal_Substrate_Abundance()
    CAZymeSub_abund.output_substrate_abund()

def PUL_abundance(args):
    paras = abund_parameters(args)
    PUL_abund = CAZyme_Abundance_estimate(paras)
    PUL_abund.Cal_Seq_Abundance(args.abundance)
    PUL_abund.Cal_PUL_Abundance()
    PUL_abund.output_cgc_abund()

def PUL_Substrate_abundance(args):
    paras = abund_parameters(args)
    PUL_abund = CAZyme_Abundance_estimate(paras)
    PUL_abund.Cal_Seq_Abundance(args.abundance)
    PUL_abund.Cal_PUL_Abundance()
    PUL_abund.Cal_PUL_Substrate_Abundance()
    PUL_abund.output_cgcsubstrate_abund()

def Get_GFF_attri(field,ID="ID"):
    fields = field.split(";")
    for field in fields:
        if field.startswith(ID+"="):
            return field.split("=")[-1]
    return ""

def read_PULgff(filename):
    geneID2feature = {}
    for line in open(filename):
        if not line.startswith("#"):
            lines = line.rstrip("\n").split("\t")
            if lines[2] == "CDS":
                ID = Get_GFF_attri(lines[-1]).split(".")[0]
                gene_name = Get_GFF_attri(lines[-1],"Name")
                locus_tag =  Get_GFF_attri(lines[-1],"locus_tag")
                protein_id =  Get_GFF_attri(lines[-1],"protein_id")
                geneID2feature[ID] = (gene_name,locus_tag,protein_id)
    return geneID2feature

def read_cgcgff(filename):
    geneID2feature = {} ; line_order = 0
    for line in open(filename):
        if not line.startswith("#"):
            line_order += 1
            lines = line.rstrip("\n").split("\t")
            ID = Get_GFF_attri(lines[-1])
            signature =  Get_GFF_attri(lines[-1],"DB")
            family = lines[2]
            geneID2feature[ID] = (line_order,family,signature)
    return geneID2feature

def generate_PULfaa(args):
    ''' read protein sequence from PUL.faa,cgc.gff and PUL.gff ''' 
    ''' cgc.gff provide coordinate information and gene, locus and protid'''
    
    if not args.input.endswith("/"):
        args.input = args.input +"/"
    
    ### get all dbCAN_PUL folders
    dbCAN_PUL_folders = os.listdir(args.input+"dbCAN-PUL")
    seqs = []
    for folder in dbCAN_PUL_folders:
        if folder.endswith(".out"): ### PULID.out
            PULID = folder.split(".")[0]
            PULID_num = int(PULID[3:len(PULID)])
            if PULID_num < 602 or PULID_num > 656:
                continue
            print(f"deal with PUL {PULID}")
            faa_file = args.input + "dbCAN-PUL/" + folder + "/" + PULID +".faa"
            cgc_gff  = args.input + "dbCAN-PUL/" + folder + "/" + "cgc.gff"
            gff  = args.input + "dbCAN-PUL/" + folder + "/" + PULID +".gff"
            if not os.path.exists(faa_file):
                print(f"File {faa_file} does not exit!")
                continue
            if not os.path.exists(cgc_gff):
                print(f"File {cgc_gff} does not exit!")
                continue
            if not os.path.exists(gff):
                print(f"File {gff} does not exit!")
                continue
            PUL_geneid2feature = read_PULgff(gff)
            cgcgff_geneid2feature = read_cgcgff(cgc_gff)
            seqid2seqs = {}
            for seq in SeqIO.parse(faa_file,'fasta'):
                oldseqid = seq.id
                #print(oldseqid,PUL_geneid2feature,cgcgff_geneid2feature)
                if oldseqid in cgcgff_geneid2feature and oldseqid in PUL_geneid2feature:
                    #pulid_order:pulid:gene:locus_tag:protein_id:type:family
                    pulid_order = PULID + "_" +str(cgcgff_geneid2feature[oldseqid][0])
                    gene = PUL_geneid2feature[oldseqid][0]
                    locus_tag  = PUL_geneid2feature[oldseqid][1]
                    protein_id = PUL_geneid2feature[oldseqid][2]
                    type1 = cgcgff_geneid2feature[oldseqid][1]
                    family  = cgcgff_geneid2feature[oldseqid][2]
                    new_seqid = f"{pulid_order}:{PULID}:{gene}:{locus_tag}:{protein_id}:{type1}:{family}"
                    seq.id = new_seqid
                    seq.description = ""
                    seqs.append(seq)
    SeqIO.write(seqs,"PUL_updated.faa",'fasta')

def parse_argv():

    usage = '''
    %(prog)s [positional arguments] [options]
    -----------------------------------------
    fam_abund [calculate the abundance of CAZyme]. Example usage: dbcan_utils fam_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
    -----------------------------------------
    fam_substrate_abund [calculate the abundance of CAZyme substrate]. Example usage: dbcan_utils fam_substrate_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
    -----------------------------------------
    CGC_abund [calculate the abundance of PUL]. Example usage: dbcan_utils CGC_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
    -----------------------------------------
    CGC_substrate_abund [calculate the abundance of PUL substrate]. Example usage: dbcan_utils CGC_substrate_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
    '''
    parser = argparse.ArgumentParser(description='Calculate the abundance CAZyme, PUL and substrate.',prog='dbcan_utils',usage=usage)
    parser.add_argument('function', help='What function will be used to analyze?',choices=["fam_abund","fam_substrate_abund","CGC_abund","CGC_substrate_abund","generate_PULfaa"])
    parser.add_argument('-i','--input',help='dbCAN CAZyme annotation output folder.',default="output",required=True)
    parser.add_argument('-bt','--bedtools',help='bedtools gene reads count results.')
    parser.add_argument('-1','--R1',help='R1 reads, support gz compress file')
    parser.add_argument('-2','--R2',help='R2 reads, support gz compress file, None for single end sequencing',default=None)
    parser.add_argument('-o','--output',help='output files',default="output")
    parser.add_argument('-a','--abundance',default="RPM",help='normalized method',choices=["FPKM","RPM","TPM"])
    parser.add_argument('--db_dir', default="db", help='dbCAN database directory')
    args = parser.parse_args()
    return args

def main():
    args = parse_argv()
    if args.function == "fam_abund":
        # dbcan_utils CAZyme_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
        CAZyme_abundance(args)
    if args.function == "fam_substrate_abund":
        # dbcan_utils CAZymeSub_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
        CAZymeSub_abundance(args)
    if args.function == "CGC_abund":
        # dbcan_utils PUL_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
        PUL_abundance(args)
    if args.function == "CGC_substrate_abund":
        # dbcan_utils PULSub_abund -bt samfiles/fefifo_8002_7.depth.txt -i fefifo_8002_7.dbCAN
        PUL_Substrate_abundance(args)
    if args.function == "generate_PULfaa": ### update dbCAN-PUL faa proteins
        # dbcan_utils generate_PULfaa -i db
        generate_PULfaa(args)
if __name__ =="__main__": ### for test
    pass
