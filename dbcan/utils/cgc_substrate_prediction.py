import argparse,os
from Bio import SeqIO
import uuid,sys
import shutil
import pandas as pd
import numpy as np
import math,json
import time

#ROOT_FOLDR = "/mnt/raid5-1/jinfang/dbCAN3/db/"

def Sum_bitscore(genes):
    return sum([gene.bitscore for gene in genes])

class blastp_hit(object):
    def __init__(self,lines):
        self.qseqid = lines[0]  ### NC_000913.3|CGC30|NC_000913.3_2025|CAZyme|GT2
        self.sseqid = lines[1]  ### PUL0090_11:PUL0090:wbuB::ACA24857.1:CAZyme:GT4;PULID_Protorder:PULID:Gene name:locus:protein ID:sigature:family/feature
        self.pident = lines[2]
        self.length = int(lines[3])
        self.mismatch = int(lines[4])
        self.gapopen  = int(lines[5])
        self.qstart   = int(lines[6])
        self.qend     = int(lines[7])
        self.sstart   = int(lines[8])
        self.send     = int(lines[9])
        self.evalue   = float(lines[10])
        self.bitscore = float(lines[11])
        if len(lines) >= 13:
            self.qlen = int(lines[12])
        if len(lines) >= 14:
            self.slen = int(lines[13])

    def __repr__(self):
        return "\t".join([str(self.__dict__[attr]) for attr in self.__dict__])
    
    def __eq__(self, other):
        if self.evalue == self.evalue:
            if self.bitscore == self.bitscore:
                if self.pident == other.pident:
                    return 1
    
    def __le__(self,other):
        if self.evalue > other.evalue:
            return 1
        elif self.evalue == other.evalue:
            if self.bitscore < other.bitscore:
                return 1
            elif self.bitscore == other.bitscore:
                if self.pident < other.pident:
                    return 1
    
    def format_str(self):
        qseqids = self.qseqid.split("|")
        sseqids = self.sseqid.split(":")
        qtype = qseqids[3]
        stype = sseqids[5]
        if qtype == "CAZyme":
            families = ";".join(qseqids[4:])
            qseqid = qseqids[2] + "|" + qseqids[3] + "|" + families
        else:
            qseqid = qseqids[2] + "|" + qseqids[3]
        
        if stype == "CAZyme":
            sseqid = sseqids[0] + "|" + sseqids[5] + "|" +  sseqids[6].replace("|",";")
        else:
            sseqid = sseqids[0] + "|" + sseqids[5]
        cgcid = qseqids[0] + "|" + qseqids[1]
        pulid = sseqids[1]
        return "\t".join([qseqid,sseqid,cgcid,pulid,self.pident,str(self.length),str(self.mismatch),str(self.gapopen),str(self.qstart),str(self.qend),str(self.sstart),str(self.send),str(self.evalue),str(self.bitscore),str(self.qlen),str(self.slen)])
class Gene(object):
    '''
    design for each line in cgc_standard.out 
    '''
    def __init__(self,lines,cluster_number=0):
        self.clusterid = lines[0]
        #self.order =    cluster_number + 1
        self.type  =    lines[1]
        self.contig =   lines[2]
        self.Protein_ID = lines[3]
        self.start =    int(lines[4])
        self.end =      int(lines[5])
        self.strand =   lines[6]
        ### for top3 if type == "TC"
        if self.type == "TC":
            self.Protein_Fam = ".".join(lines[7].split(".")[0:3])
        else:
            self.Protein_Fam =    lines[7]
        
        self.CGC_ID = self.contig + "|" +self.clusterid

    def gatherAttrs(self):
        return "".join("{} = {}\n".format(k, getattr(self, k))for k in self.__dict__.keys())
    def __repr__(self):
        return "{}{}".format("", self.gatherAttrs())
    
    def format_out(self):
        return "\t".join([self.clusterid,self.type,self.contig,self.Protein_ID,str(self.start),str(self.end),self.strand,self.Protein_Fam])
    
    def __getattr__(self,name):
        return ("-")
    def Get_CAzyID(self):
        return self.contig+"|"+self.clusterid+"|"+self.Protein_ID+"|"+self.type + "|" +self.Protein_Fam

class dbCAN_Out(object):
    '''
    design for the whole cgc_standard.out file
    '''
    def __init__(self,filename):
        hits = open(filename).readlines()[1:]
        self.genes= []
        for line in hits:
            if line.startswith("CGC#"):
                continue
            lines = line.split()
            self.genes.append(Gene(lines))
    
    def __iter__(self):
        return iter(self.genes)

    def Out_Null_gene(self,label="null",filename="out"):
        with open(filename,'w') as f:
            for gene in self:
                if gene.type == label:
                    f.write(gene.Protein_ID+"\n")
    
    def Out_Null_gene_array(self,label="null"):
        return [gene.Protein_ID for gene in self if gene.type == label]

    def CGCID2genes(self):
        cgcdict = {}
        for gene in self:
            cgcdict.setdefault(gene.CGC_ID,[]).append(gene)
        return cgcdict

    def ProteinID2genes(self):
        ProteinIDdict = {}
        for gene in self:
            ProteinIDdict[gene.Protein_ID] = gene
        return ProteinIDdict


class CGC(object):
    def __init__(self,genes):
        self.genes = genes
        self.ID = genes[0].CGC_ID ### get cgc id
        self.start = min([gene.start for gene in genes])
        self.end = max([gene.end for gene in genes])
        self.gene_num = len(genes)

    def __iter__(self):
        return iter(self.genes)
    
    def __repr__(self):
        #return " ".join([self.ID,str(self.start),str(self.end),str(self.gene_num)])
        return "\t".join([self.ID,str(self.start),str(self.end),str(self.gene_num)])

    def clean_signature(self,prot2domain):
        for gene in self:
            if gene.type == "CAZyme":
                gene.Protein_Fam = "|".join(prot2domain.get(gene.Protein_ID,[]))

    def Out2file(self,fasta):
        seqs = []
        for gene in self:
            seq = fasta[gene.Protein_ID]
            #description = seq.description
            #seq.description = ""
            seq.id = gene.contig + "|" + gene.clusterid + "|" + gene.Protein_ID + "|" + gene.type
            seqs.append(seq)
        SeqIO.write(seqs,self.ID+".fasta",'fasta')
    
    def get_cgc_CAZyme(self):
        return [gene.type for gene in self]
    
    def get_proteinfam(self):
        return [gene.Protein_Fam for gene in self]

    def get_proteinfam_assign_null(self,pfam):
        contents = []
        for gene in self:
            if gene.type == "null":
                domains = pfam.get(gene.Protein_ID,"null")
                contents.append(domains)
            else:
                contents.append(gene.Protein_Fam)
        return ",".join(contents)

    def get_cgc_CAZyme_sig(self,type):
        return [gene.Protein_Fam for gene in self if gene.type ==type]


    def get_protein_null(self):
        return [gene.Protein_ID for gene in self if gene.type == "null"]

    def get_proteinID(self):
        return [gene.Protein_ID for gene in self ]


    def __len__(self):
        return len(self.genes)
    
    def get_CAZyme_num(self,orders=["CAZyme","TC","TF","STP","null"]):
        types = self.get_cgc_CAZyme()  
        return [str(types.count(tmp)) for tmp in orders]
    
    def get_positions(self):
        starts = [] ; ends = [] ; strands = []
        for gene in self:
            starts.append(gene.start)
            ends.append(gene.end)
            strands.append(gene.strand)
        return starts,ends,strands

class CGC_hub(object):
    '''
    design for the instance of dbCAN_Out
    '''
    def __init__(self,dbcan):
        self.CGCs = []
        cgcdict = dbcan.CGCID2genes()
        for cgc in cgcdict:
            self.CGCs.append(CGC(cgcdict[cgc]))
    
    def __iter__(self):
        return iter(self.CGCs)
    
    def Out2file(self,fasta):
        for cgc in self:
            cgc.Out2file(fasta)
    
    def CGCID2CGC(self):
        return {cgc.ID:cgc for cgc in self}

class dbSub(object):
    '''
    design for dbCAN3, dbCAN_sub output
    '''
    def __init__(self,filename,dbsub_parameters):
        self.Genes = []
        for line in open(filename).readlines()[1:]: ### ignored the first line
            lines = line.rstrip("\n").split("\t")
            hmmevalue = float(lines[7])
            hmmcov    = float(lines[12])
            if hmmevalue <= dbsub_parameters.hmmevalue and hmmcov >= dbsub_parameters.hmmcov:
                self.Genes.append(dbSub_record(lines))
    
    def __iter__(self):
        return iter(self.Genes)
    
    def GeneID2gene(self):
        ### each gene may contain more than two dbsub records
        geneid2gene = {}
        for gene in self:
            geneid2gene.setdefault(gene.GeneID,[]).append(gene)
        self.geneid2gene = geneid2gene
        return geneid2gene


class dbSub_record(object):
    '''
    design for dbCAN_sub output, each line
    '''
    def __init__(self,lines):
        self.dbcan_sub_subfam = lines[0]
        self.Subfam_Composition = lines[1]
        self.Subfam_EC = lines[2]
        self.Substrate = lines[3] if lines[3]!= "-" else ""
        self.hmm_Length = lines[4]
        self.GeneID = lines[5]
        self.GeneLen = lines[6]
        self.E_Value = lines[7]
        self.hmm_Start = lines[8]
        self.hmm_End = lines[9]
        self.Gene_Start = lines[10]
        self.Gene_End = lines[11]
        self.Cov = lines[12]
    
    def __repr__(self):
        return "\t".join([self.__dict__[name] for name in self.__dict__])
    def __str__(self):
        return "\t".join([self.__dict__[name] for name in self.__dict__])

def replace_black_underline(df):
    col_names = df.columns.tolist()
    for index,value in enumerate(col_names):
        col_names[index]= value.replace(" ","_")
    df.columns=col_names

class dbCAN_substrate_prediction(object):
    '''
    design for running substrate prediciton process with the output of dbCAN3. Run the following step.
    Step 1: run dbCAN-PUL searcing substrate prediction. -> dbCAN_PUL_substrate_predict
    Step 2: run dbCAN_sub-family substrate prediciton. -> dbcan_sub_subfamily_substrate_prediction
    Step 3: combine eCAM-family and dbCAN-PUL substrate prediciton.
    '''
    
    def __init__(self,args):
        '''
        prepare all the input file and outfile names required by substrate prediction
        '''
        ### input
        self.input_folder = args.input if args.input.endswith("/") else args.input +"/"
        self.cgc_out = self.input_folder + "cgc.out"
        self.cgc_standard_out = self.input_folder + "cgc_standard.out"
        self.dbsub_out = self.input_folder +"dbcan-sub.hmm.out"
        self.overview_txt = self.input_folder +"overview.txt" 
        self.protein_db     = self.input_folder +"uniInput"

        ### output
        self.out = args.out
        self.random_str = uuid.uuid4().hex ### tmp folder to save some tmp files
        #self.random_str = "59f220bd5fc4422187679301976a3d76" ## for debug
        #self.tmp_folder = "/dev/shm/" + self.random_str ### tmp folder to save some tmp files
        self.tmp_folder = self.input_folder
        self.tmp_blastp_out = self.tmp_folder + "PUL_blast.out"
        self.tmp_CAZyme_pep = f"{self.tmp_folder}CGC.faa"

        ### parameters
        self.PULdb  = f"{ROOT_FOLDR}PUL.faa"
        self.pul_excel_filename = f"{ROOT_FOLDR}dbCAN-PUL.xlsx"
        self.homologous_parameters  = HitParamter(args)
        self.dbsub_parameters  = dbcan_sub_parameter(args)
        
        ### output parameters, intermediate results
        self.odbcan_sub = args.odbcan_sub
        self.odbcanpul = args.odbcanpul
        self.dbcanpul_tmp = "dbcanpul.tmp.txt"
        self.dbcan_sub_tmp = "dbcansubpul.tmp.txt"

        ### Method to predict substrate 
        self.run_dbCAN_sub = True
        self.run_dbCAN_PUL = True
        
        ### 
        self.dbcan_sub_CGC2substrates = {}
        self.queryCGC2hit = {}
        self.dbcan_sub_CGC2maxscore = {}
    
    def check_input(self):
        '''
        check input files 
        '''
        ### input_folder
        if not os.path.exists(self.input_folder):
            print(f"input folder: {self.input_folder} dose not exist!",file=sys.stderr);exit()
        if not os.path.exists(self.cgc_out):
            print(f"dbCAN3 cgc finder out file: {self.cgc_out} dose not exist!",file=sys.stderr)
        if not os.path.exists(self.cgc_standard_out):
            print(f"dbCAN3 cgc finder standard out file: {self.cgc_out} dose not exist!",file=sys.stderr);exit()
        if not os.path.exists(self.dbsub_out):
            print(f"dbCAN3 dbsub output file: {self.dbsub_out} dose not exist!",file=sys.stderr)
            print(f"Substrate prediciton based on major voting will not applied!",file=sys.stderr)
            self.run_dbCAN_sub = False
        if not os.path.exists(self.overview_txt):
            print(f"dbCAN3 overview file: {self.overview_txt} dose not exist!",file=sys.stderr)
        if not os.path.exists(self.protein_db):
            print(f"dbCAN3 protein sequences: {self.protein_db} dose not exist!",file=sys.stderr);exit()
        if not os.path.exists(self.PULdb):
            print(f"dbCAN-PUL database: {self.PULdb} dose not exist!",file=sys.stderr);exit()
        if not os.path.exists(self.pul_excel_filename):
            print(f"dbCAN-PUL substrate excel file: {self.pul_excel_filename} dose not exist!",file=sys.stderr);exit()

    def __repr__(self):
        '''
        return the all the string variables of class
        '''
        return "\n".join([self.__dict__[name] for name in self.__dict__ if isinstance(self.__dict__[name],str)])
    
    def extract_seq_in_CGC(self):
        '''read protein sequence and save them in self.seqs
        '''
        print("Start extracting sequence!")
        self.dbCAN_hits = dbCAN_Out(self.cgc_standard_out)
        ### 
        #for gene in self.dbCAN_hits:
        #    print(gene)
        #exit()
        
        self.seqid2seq = SeqIO.to_dict(SeqIO.parse(self.protein_db,'fasta'))
        self.protid2gene = self.dbCAN_hits.ProteinID2genes()
        self.seqs = []
        for seqid in self.seqid2seq:
            if seqid in self.protid2gene:
                #print(seqid,self.seqid2seq[seqid].id,self.protid2gene[seqid].Get_CAzyID())
                self.seqid2seq[seqid].id = self.protid2gene[seqid].Get_CAzyID()
                self.seqs.append(self.seqid2seq[seqid])

    def do_blastp_against_dbCANPUL(self):
        ''' sequences save to dish, and then do blastp.
        to accelerate, the sequence can be save to /dev/shm/
        '''
        print(f"Start blastp CAZyme sequences against to database: {self.PULdb}",file=sys.stderr)
        _ = self.extract_seq_in_CGC() if not self.seqs else 1 ### to run
        os.makedirs(self.tmp_folder, exist_ok=True)
        ###
        SeqIO.write(self.seqs,self.tmp_CAZyme_pep,'fasta')
        
        outfmt = '"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
        self.blastp_command = f"blastp -max_hsps 1 -query {self.tmp_CAZyme_pep} -db {self.PULdb} -outfmt {outfmt} -evalue 0.01 -out {self.tmp_blastp_out} -num_threads 32 "
        print(self.blastp_command)
        print("[whether PUL db exists]", os.path.exists(self.PULdb))
        
        ### checking the blastp out
        
        if not os.path.exists(self.tmp_blastp_out):
            status = os.system(self.blastp_command)
        else:
            status = 0 

        #status = os.system(self.blastp_command) if not os.path.exists(self.tmp_blastp_out) else 0
        if status != 0:
            print(f'comand line: "{self.blastp_command}" runs return error.', file=sys.stderr)
    
    def read_dbCAN_PUL(self):
        self.Puls = pd.read_excel(self.pul_excel_filename)
        replace_black_underline(self.Puls)
    
    def read_blastp_result(self,filename):
        '''
        design for reading blastp 6 fmtout, return a dictionary with key as query id, and value is array of blastp out, one element one line.
        I don't want to the function has many parameters, so design it as class function. 
        '''
        print(f"Reading blastp result {filename}")
        querydict = {}
        for line in open(filename):
            qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,qlen,slen = line.split() 
            ## SeqID|CGCID|signature|family
            queryids = qseqid.split("|") ### so the query seqid can not named with "|"
            if len(queryids) < 2:
                print(f"{queryid} can not be splited with |.",file=sys.stderr)
                exit()
            queryid = queryids[0] + "|" + queryids[1]
            ### add filter conditions from user defined
            if float(evalue) > self.homologous_parameters.evalue_cutoff:
                continue
            if float(pident) < self.homologous_parameters.identity_cutoff:
                continue
            if (float(qend) - float(qstart) + 1)/float(qlen) < self.homologous_parameters.coverage_cutoff:
                continue
            querydict.setdefault(queryid,[]).append(blastp_hit(line.split()))
        self.hitdict = querydict
    
    def Uniq_blastp_hit(self,blast_list):
        uniqs = []; genes = []; uniqss = []; 
        ### hit types, CAZy
        homologous_pairs = []
        for tmp in blast_list:
            if tmp.qseqid not in uniqs:### the seqid in 
                genes.append(tmp) ### gene should be the homologous pairs
                uniqs.append(tmp.sseqid)
            uniqss.append(tmp.qseqid)
            #print(tmp)
            hit_type   = tmp.sseqid.split(":")[-2]
            query_type = tmp.qseqid.split("|")[3] ## this requires the input fasta seqid can not be named with "|"
            #print(hit_type,query_type)
            #exit()
            #print (tmp.sseqid)
            #print (tmp.qseqid)
            if hit_type == query_type:
                homologous_pairs.append(query_type + "-" + hit_type)
            CAZyme_pairs_num = homologous_pairs.count("CAZyme-CAZyme")
        if len(uniqs) >= self.homologous_parameters.uqcgn and len(set(uniqss))>= self.homologous_parameters.upghn and CAZyme_pairs_num >= self.homologous_parameters.cpn and len(homologous_pairs) >= self.homologous_parameters.tpn:
            if not self.homologous_parameters.ept: ### extra signature pairs?
                score = Sum_bitscore(genes) ### sum score
                if score/len(uniqs) >= self.homologous_parameters.bitscore_cutoff: ### >= average bitscore_cutoff
                    ### for blastp hit
                    ### self.queryCGC_CGChits_genes_blastp_hit
                    return score,homologous_pairs
                else:
                    return -1,homologous_pairs
            else: ### require extra signature pairs
                signature_pairs = self.homologous_parameters.ept
                signature_pairs_num = self.homologous_parameters.eptn
                pair_signal = 0
                for i,signature_pair in enumerate(signature_pairs):
                    if homologous_pairs.count(signature_pair) >= int(signature_pairs_num[i]):
                        pair_signal += 1
                if pair_signal == len(signature_pairs): ### satifies all signature_pair
                    score = Sum_bitscore(genes) ### sum score
                    if score/len(uniqs) >= self.homologous_parameters.bitscore_cutoff: ### >= average bitscore_cutoff
                        return score,homologous_pairs
                    else:
                        return -1,homologous_pairs
                else:
                    return -1,homologous_pairs
        else:
            return -1,homologous_pairs
    
    def dbcan_sub_read_cgc(self):
        if not self.cgcid2cgc:
            ### need to read cgc
            self.dbCAN_hits = dbCAN_Out(self.cgc_standard_out)
            self.cgcs = CGC_hub(self.dbCAN_hits)
            self.cgcid2cgc = self.cgcs.CGCID2CGC()
    
    def analyze_blastp_out(self):
        print("Start analyzing blastp result",file=sys.stderr)
        self.read_blastp_result(self.tmp_blastp_out)
        self.cgcs = CGC_hub(self.dbCAN_hits)  ### dbCAN_hits has been generated by running extract_seq
        self.cgcid2cgc = self.cgcs.CGCID2CGC() ### cgcid 2 cgc
        
        #self.bitscore_cutoff = 50; self.uniq_query_cgc_gene_num = 3; self.uniq_hit_pul_gene_num = 3
        #self.CAZyme_pairs_num_cutoff = 1;self.signature_pair = 2
        
        queryCGC2scores = {}; queryCGC2pulids = {}; queryCGC2Mapped_types = {}
        self.queryCGC2blastp_hits = {} ### to save the process results
        for hit in self.hitdict:
            tmp_dict = {}
            for hit2 in self.hitdict[hit]: ### for each query 
                pulid = hit2.sseqid.split("_")[0]  ### why? PUL0350_2 -> PUL0350
                tmp_dict.setdefault(pulid,[]).append(hit2)
            scores = []; pulids = []; maped_types = []
            ### for each query
            for pulid in tmp_dict:
                score,homologous_pairs = self.Uniq_blastp_hit(tmp_dict[pulid])
                if score > 0: ### if score == -1,  now we got an homologous hits, here
                    scores.append(score)
                    pulids.append(pulid)
                    maped_types.append(homologous_pairs)
                    ### maybe here, we need to track the proces, so that we can debug in the feature, just kept the blastp hits here.
                    self.queryCGC2blastp_hits.setdefault(hit,[]).append(tmp_dict[pulid])
            queryCGC2scores[hit] = scores ; queryCGC2pulids[hit] = pulids ; queryCGC2Mapped_types[hit] = maped_types
        ### save the query CGC hits results
        self.queryCGC2scores = queryCGC2scores ; self.queryCGC2pulids = queryCGC2pulids ; self.queryCGCmapedtypes = queryCGC2Mapped_types
    
    def get_best_pul_hit(self):
        self.read_dbCAN_PUL() ### get self.Puls
        queryCGC2hit = {}
        for queryCGC in self.queryCGC2scores:
            scores = self.queryCGC2scores[queryCGC]
            if len(scores) == 0: ### not hit
                continue
            ### sorted the best values
            score_orders = np.argsort(scores)
            max_score_index = score_orders[-1]
            score = scores[max_score_index]
            bestpulid = self.queryCGC2pulids[queryCGC][max_score_index]
            substrate = self.Puls[self.Puls["ID"] == bestpulid].substrate_final.values[0]
            mapped_types = self.queryCGCmapedtypes[queryCGC][max_score_index]
            queryCGC2hit[queryCGC] = PULhit(score,bestpulid,substrate,mapped_types)
        self.queryCGC2hit = queryCGC2hit
    
    def get_best_pul_hit_and_blastphit(self):
        print("Get best pul hit ",file=sys.stderr)
        self.read_dbCAN_PUL() ### get self.Puls
        queryCGC2hit = {}; self.queryCGC_best_genes_blastp_hit = {}; self.queryCGC_CGChits_genes_blastp_hit = {}
        #if self.odbcanpul: ### output
        #    f = open(self.dbcanpul_tmp,'w')
        for queryCGC in self.queryCGC2scores:
            scores = self.queryCGC2scores[queryCGC]
            if len(scores) == 0: ### not hit
                continue
            ### sorted the best values
            score_orders = np.argsort(scores)
            max_score_index = score_orders[-1] ### get best hits
            best_score = scores[max_score_index]
            bestpulid = self.queryCGC2pulids[queryCGC][max_score_index]
            substrate = self.Puls[self.Puls["ID"] == bestpulid].substrate_final.values[0]
            mapped_types = self.queryCGCmapedtypes[queryCGC][max_score_index]
            queryCGC2hit[queryCGC] = PULhit(best_score,bestpulid,substrate,mapped_types)
            self.queryCGC_best_genes_blastp_hit[queryCGC] = self.queryCGC2blastp_hits[queryCGC][max_score_index] ### save the best PUL, genes blastp records
            ###
            ### save all the blastp results 

            for idx in score_orders:
                self.queryCGC_CGChits_genes_blastp_hit.setdefault(queryCGC,[]).append(self.queryCGC2blastp_hits[queryCGC][idx]) 

            #if self.odbcanpul:
            #    for sorted_idx in score_orders[::-1]:
            #        score = scores[sorted_idx]
            #        pulid = self.queryCGC2pulids[queryCGC][sorted_idx]
            #        substrate = self.Puls[self.Puls["ID"] == pulid].substrate_final.values[0]
            #        mapped_types = ";".join(self.queryCGCmapedtypes[queryCGC][sorted_idx])
            #        f.write(f"{queryCGC}\t{pulid}\t{substrate}\t{score}\t{mapped_types}\n") ### save the all the potential hits for each CGC
        self.queryCGC2hit = queryCGC2hit
        #f.close()
    
    def print_best_result_and_blastp(self):
        '''
        output substrate prediction result and the genes homologous 
        '''
        for queryCGC in self.queryCGC2hit:
            pulhit = self.queryCGC2hit[queryCGC]
            print(queryCGC,pulhit)
            for blastp in self.queryCGC_best_genes_blastp_hit[queryCGC]:
                print(blastp)
    
    def print_result_and_blastp(self):
        with open("dbcanpul.hit.blastp.txt",'w') as f:
            for queryCGC in self.queryCGC_CGChits_genes_blastp_hit:
                for hit_PULID_blastps in self.queryCGC_CGChits_genes_blastp_hit[queryCGC]:
                    for blastp in hit_PULID_blastps:
                        f.write(blastp.format_str()+"\n")

    def dbCAN_PUL_substrate_predict(self):
        '''
        substrate prediction based on dbCAN-PUL searching.
        Step 1: extract cgc protein sequences from dbCAN3 output. -> extract_seq_in_CGC
        Step 2: blastp cgc protein sequences against dbCAN-PUL. -> do_blastp_against_dbCANPUL
        Step 3: deal with the blastp result. This step includes customed parameters: average bitscore, evalue, signature pairs. and then run dbCAN-PUL homologous substrate prediction. -> analyze_blastp_out
        Step 4: get the best hit and predict the substrate -> get_best_pul_hit or get_best_pul_hit_and_blastphit
        '''
        self.extract_seq_in_CGC()
        self.do_blastp_against_dbCANPUL()
        self.analyze_blastp_out()
        #self.get_best_pul_hit()
        self.get_best_pul_hit_and_blastphit()
        if self.odbcanpul:
            self.print_result_and_blastp()
    
    def Read_CAZyme_substrate(self):
        '''
        read dbsub.out into RAM
        '''
        print(f"Reading dbsub outfile:{self.dbsub_out}")

        self.CAZyme2substrate = dbSub(self.dbsub_out,self.dbsub_parameters)  ### need to filter dbsub results with conditions 
        self.geneid2dbsub = self.CAZyme2substrate.GeneID2gene() ### one geneID map to a list, the element is one line in dbsub.out
        cgcid2sub = {}
        for cgcid in self.cgcid2cgc: ## NC_000913.3|CGC26
            for gene in self.cgcid2cgc[cgcid]: ### 
                if gene.Protein_ID in self.geneid2dbsub:
                    cgcid2sub.setdefault(cgcid,[]).append(self.geneid2dbsub[gene.Protein_ID])
                    #print(cgcid,gene.Protein_ID,self.geneid2dbsub[gene.Protein_ID])


    def dbcan_sub_subfamily_substrate_prediction(self):
        '''
        substrate prediction based on dbCAN-sub subfamily.
        this part includes 4 steps.
        Step 1: reading dbsub prediciton result -> Read_CAZyme_substrate
        Step 2: reading cgc result -> dbcan_sub_read_cgc
        Step 3: combine dbsub result and cgc result. -> CGC2substrate_dbcan_sub
        Step 4: scoring the CGC substrate predicted by dbCAN-sub, and get the best substrate -> substrate_scoring_dbcan_sub
        '''
        
        self.Read_CAZyme_substrate()
        self.dbcan_sub_read_cgc()
        self.CGC2substrate_dbcan_sub()
        self.substrate_scoring_dbcan_sub()

    def substrate_scoring_dbcan_sub(self):
        print("Start dbCAN-sub subfamily substrate scoring")
        finalsub = {}; finalscores = {}; finalranks = {}; finalmaxscore = {}
        cgcid2sub = self.cgcid2substrate_dbcan_sub 
        ### cgcid map to dbsub records, two dimensions list. 
        ### one dimension is from CAZyme to substrate due to one CAZyme can have more than two substrate or two sub-family
        ### another dimension is from cgc, one cgc includes several CAZymes
        ### a complex logicl as following
        for cgcid in cgcid2sub:
            if self.cgcid2CAZyme_substrate_num[cgcid] < self.dbsub_parameters.num_of_domains_substrate_cutoff: ## domains
                ### only has one CAZyme substrate prediction,so discard this part
                continue
            if self.cgcid2CAZyme_domain_substrate_num[cgcid] < self.dbsub_parameters.num_of_protein_shared_substrate_cutoff: ## seqeunces
                ### sequences
                continue
            scores = {};ranks = []
            for subs_list in cgcid2sub[cgcid]: ### dealing with substrates, two dimensions
                for subs in subs_list:
                    subs = subs.Substrate
                    subs = subs.replace("and",",") ### based some naming habit by someone
                    subs = subs.replace(" ","") ### remove some blanks
                    subss = subs.split(",")
                    subss = set(subss)
                    #print(cgcid,subss)
                    for tmp_sub in subss: ### loop for substrate
                        if not tmp_sub: ### exclude "" substrate come from "-"
                            continue
                        tmp_sub = tmp_sub.strip(" ") ### remove blank in the two ends of substrate, 
                        tmp_subs = tmp_sub.split("|") ### some substrates combined by "|"
                        #tmp_subs = list(set(tmp_subs)) ### unique
                        #print(cgcid,tmp_subs)
                        for i in range(len(tmp_subs)): ### 
                            scores.setdefault(tmp_subs[i],[]).append(math.pow(2,-i)) 
                            ### "|" indicates the priority of the substate for CAZyme, so different score was assigned
            for sub in scores: ### loop for substrate
                scores[sub] = sum(scores[sub])
                ranks.append(f"{sub}:{scores[sub]}")
            finalscores[cgcid] = scores
            max_score = max(scores.values())
            if max_score < self.dbsub_parameters.dbcan_substrate_scors: ### subsrate score less than cutoff
                continue ## next cgc
            finalmaxscore[cgcid] = max_score
            final_subs = []
            for sub,score in scores.items():
                if score == max_score: ### each substrate may includes more than two substrate with identity score
                    final_subs.append(sub)
            finalsub[cgcid] = ",".join(final_subs)
            finalranks[cgcid] = ranks
        
        self.dbcan_sub_CGC2substrates = finalsub  ### save the cgc substrate
        self.dbcan_sub_CGC2scores     = finalscores ### almost the same as dbcan_sub_substrate_score
        self.dbcan_sub_substrate_score= finalranks  ###
        self.dbcan_sub_CGC2maxscore   = finalmaxscore ### max score
       

    def dbcan_sub_intermediate_file(self):
        f = open("dbCAN-sub.tmp.txt",'w')
        geneids_uniq = []
        for cgcid in self.cgcid2cgc:
            if cgcid in self.dbcan_sub_CGC2substrates:
                for gene in self.cgcid2cgc[cgcid].genes:
                    geneid = gene.Protein_ID
                    if geneid not in geneids_uniq:
                        geneids_uniq.append(geneid)
                        ECs,subs,esubfam = self.from_geneid_get_dbcan_sub(geneid)
                        f.write(cgcid+"\t"+geneid+"\t"+",".join(ECs)+"\t"+",".join(subs)+"\t"+",".join(esubfam)+"\n")
    
    def from_geneid_get_dbcan_sub(self,geneid):
        genes = self.geneid2dbsub.get(geneid,"")
        if not genes:
            return [],[],[]
        ECs = [] ; subs = [] ; esubfam = []
        for i,gene in enumerate(genes):
            ECs.extend(clean_EC(gene))
            subs.extend(clean_sub(gene))
            esubfam.append(gene.dbcan_sub_subfam)
        ECs = set(ECs)
        subs = set(subs)
        return ECs,subs,esubfam

        #print(geneid,ECs,subs,esubfam)
        #ECs = set([gene.Subfam_EC.split(":")[0] for gene in genes])
        #subs = set([gene.Substrate.strip(" ") for gene in genes])
        #esubfam = [gene.dbcan_sub_subfam for gene in genes]
        #f.write(geneid+"\t"+",".join(ECs)+"\t"+",".join(subs)+"\t"+",".join(esubfam)+"\n")
        #f.close()

    def dbcan_sub_sub_print_result(self):
        for cgc in self.dbcan_sub_CGC2substrates:
            print(cgc,self.dbcan_sub_CGC2substrates[cgc]) ### cgcid and substrates
            print("\t".join(self.dbcan_sub_substrate_score[cgc]))  ### ranking
            #print(self.dbcan_sub_CGC2scores[cgc])  ### ranking
            tmp_lines = ""
            for dbsubs in self.cgcid2substrate_dbcan_sub[cgc]:
                for dbsub in dbsubs:
                    print(dbsub)
            print("-"*20)
    
    def CGC2substrate_dbcan_sub(self):
        cgcid2sub = {} ; cgcid2substrate_CAZyme_num = {}
        for cgcid in self.cgcid2cgc: ## NC_000913.3|CGC26
            for gene in self.cgcid2cgc[cgcid]: ### 
                if gene.Protein_ID in self.geneid2dbsub:
                    cgcid2sub.setdefault(cgcid,[]).append(self.geneid2dbsub[gene.Protein_ID])
                    #print(cgcid,gene.Protein_ID,self.geneid2dbsub[gene.Protein_ID])
        self.cgcid2substrate_dbcan_sub = cgcid2sub
        
        self.cgcid2CAZyme_domain_substrate_num = {}
        ### count how many sequences in CAZyme has a substrate, and then calcuate the cgc potential substrate number
        for cgcid in cgcid2sub:
            cgcid_uniq_sequences = []
            for dbsub_records in cgcid2sub[cgcid]: ## CGC loop 
                for dbsub_record in dbsub_records: ## CAZyme loop
                    if dbsub_record.Substrate: ### substrate for CAZyme
                        cgcid_uniq_sequences.append(dbsub_record.GeneID)
            self.cgcid2CAZyme_domain_substrate_num[cgcid]  = len(set(cgcid_uniq_sequences))
        ### count how many domains in CAZyme has a substrate, and then calcuate the cgc potential substrate number
        for cgcid in cgcid2sub:
            cgcid_CAZyme_sub_num = 0
            for dbsub_records in cgcid2sub[cgcid]: ## CGC loop 
                for dbsub_record in dbsub_records: ## CAZyme loop
                    if dbsub_record.Substrate: ### substrate for CAZyme
                        cgcid_CAZyme_sub_num += 1
            cgcid2substrate_CAZyme_num[cgcid] = cgcid_CAZyme_sub_num
        self.cgcid2CAZyme_substrate_num = cgcid2substrate_CAZyme_num

    def substrate_predict(self):
        if self.run_dbCAN_PUL:
            self.dbCAN_PUL_substrate_predict()
        if self.run_dbCAN_sub:
            self.dbcan_sub_subfamily_substrate_prediction()

    def __del__(self):
        ''' remove tmp folder
        shutil.rmtree()
        '''

        ### no need to remove the tmp_folder, because is the working folder for job
        ## print(f"Rmoving tmp file:{self.tmp_folder}")
        ## shutil.rmtree(self.tmp_folder)
    
    def integrate_dbCANPUL_dbcan_sub(self): ### maybe need, in the future. 
        '''
        combine two methods
        '''
        pass

    def print_result(self):
        ### self.queryCGC2hit ### dbCAN-PUL hit
        ### self.dbcan_sub_CGC2substrates ### dbcan_sub substrate
        shared_cgcids = self.queryCGC2hit.keys() | self.dbcan_sub_CGC2substrates.keys()
        print("#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignarture pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score")
        for cgcid in shared_cgcids:
            dbcan_pul_part = self.queryCGC2hit.get(cgcid,"")
            dbcan_sub_substate = self.dbcan_sub_CGC2substrates.get(cgcid,"")
            PULID = dbcan_pul_part.pulid if dbcan_pul_part else ""
            dbcan_pul_sub = dbcan_pul_part.substrate if dbcan_pul_part else ""
            bitscore = dbcan_pul_part.score if dbcan_pul_part else ""
            sig_pairs = ";".join(dbcan_pul_part.maped_types) if dbcan_pul_part else ""
            dbcan_sub_maxscore = self.dbcan_sub_CGC2maxscore.get(cgcid,"")
            print(f"{cgcid}\t{PULID}\t{dbcan_pul_sub}\t{bitscore}\t{sig_pairs}\t{dbcan_sub_substate}\t{dbcan_sub_maxscore}")
    
    def result_print_to_file(self):
        ### self.queryCGC2hit ### dbCAN-PUL hit
        ### self.dbcan_sub_CGC2substrates ### dbCAN-sub substrate
        shared_cgcids = self.queryCGC2hit.keys() | self.dbcan_sub_CGC2substrates.keys()
        print (f"Writing substrate prediction result to file:{self.input_folder+self.out}")
        with open(self.input_folder+self.out,'w') as f:
            f.write("#cgcid\tPULID\tdbCAN-PUL substrate\tbitscore\tsignature pairs\tdbCAN-sub substrate\tdbCAN-sub substrate score\n")
            for cgcid in shared_cgcids:
                dbcan_pul_part = self.queryCGC2hit.get(cgcid,"")
                dbcan_sub_substate = self.dbcan_sub_CGC2substrates.get(cgcid,"")
                PULID = dbcan_pul_part.pulid if dbcan_pul_part else ""
                dbcan_pul_sub = dbcan_pul_part.substrate if dbcan_pul_part else ""
                bitscore = dbcan_pul_part.score if dbcan_pul_part else ""
                sig_pairs = ";".join(dbcan_pul_part.maped_types) if dbcan_pul_part else ""
                dbcan_sub_maxscore = self.dbcan_sub_CGC2maxscore.get(cgcid,"")
                f.write(f"{cgcid}\t{PULID}\t{dbcan_pul_sub}\t{bitscore}\t{sig_pairs}\t{dbcan_sub_substate}\t{dbcan_sub_maxscore}\n")

class PULhit(object):
    '''
    design for cgc search against PUL dababase, save bestscore, substrate, mapped_types
    '''
    def __init__(self,score,pulid,substrate,mapped_types):
        self.score = score
        self.substrate = substrate
        self.maped_types = mapped_types
        self.pulid = pulid
    def __repr__(self):
        return "\t".join([self.pulid,str(self.score),self.substrate,",".join(self.maped_types)])


def Modified_json(parameter_file):
    paramters_sub_dict =  json.load(open(parameter_file))
    paramters_sub_dict["Sub_Pred"] = "sub_pred"
    paramters_sub_dict["dbCAN_Sub"] = "dbCAN_Sub"
    paramters_sub_dict["dbCAN_PUL"] = "dbCAN_PUL"
    ### no need to adjust the dbcansub due to it require to run dbCAN-sub first
    print(f"Write parameters to file: {parameter_file}",file=sys.stderr);
    with open(parameter_file,'w') as f:
        json.dump(paramters_sub_dict,f)

def cgc_prediction_webserver(args,sub_pred):
    ### change wo working directory
    args.workdir = args.workdir if args.workdir.endswith("/") else args.workdir +"/"
    os.chdir(args.workdir) ### change the working directory to webserver working
    script_folder = sys.path[0]
    db_folder = script_folder  ### the script and database are in the same folder
    
    ### update parameters
    sub_pred.tmp_folder = args.workdir ### tmp folder to save some tmp files
    sub_pred.tmp_blastp_out = sub_pred.tmp_folder + "PUL_blast.out"
    sub_pred.tmp_CAZyme_pep = sub_pred.tmp_folder + "CAZyme.faa"

    sub_pred.PULdb  = f"{script_folder}/PUL.faa"
    sub_pred.pul_excel_filename = f"{script_folder}/dbCAN-PUL.xlsx"
    
    ### loading parameters file from php blast.php named with parameters.json
    parameter_file = "parameters.json"
    
    if args.rerun:
        Modified_json(parameter_file)

    ### adjust the paremeter to run or not?
    paramters_sub_dict =  json.load(open(parameter_file))
    
    if paramters_sub_dict["dbCAN_Sub"] and paramters_sub_dict["Sub_Pred"] and paramters_sub_dict["dbcansub"]:
        sub_pred.run_dbCAN_sub = True
    else:
        sub_pred.run_dbCAN_sub = False
    
    if paramters_sub_dict["Sub_Pred"] and paramters_sub_dict["dbCAN_PUL"]:
        sub_pred.run_dbCAN_PUL = True
    else:
        sub_pred.run_dbCAN_PUL = False
    
    ## {jobid:20221122142425,dbcansub:dbcansub,Sub_Pred:sub_pred,dbCAN_PUL:dbCAN_PUL,dbCAN_Sub:dbCAN_Sub}
    ## or {"jobid":"20221122161713","dbcansub":"","Sub_Pred":"sub_pred","dbCAN_PUL":"dbCAN_PUL","dbCAN_Sub":""}

def dbCAN3_paramters_prepare(args):
    args.input = args.out_dir if args.out_dir.endswith("/") else args.out_dir+"/"
    args.workdir = args.out_dir if args.out_dir.endswith("/") else args.out_dir+"/"
    ##
    global ROOT_FOLDR 
    ROOT_FOLDR = args.db_dir if args.db_dir.endswith("/") else args.db_dir+"/"

def cgc_substrate_prediction(args):
    ### add parameters for dealing with substrates
    
    if args.cgc_substrate: ### this means dbCAN3 call cgc substrate prediction
        dbCAN3_paramters_prepare(args)
    
    sub_pred = dbCAN_substrate_prediction(args)
    if args.env == "server":
        cgc_prediction_webserver(args,sub_pred) ### dealing with webserver

    sub_pred.check_input()
    time_start = time.time()
    sub_pred.substrate_predict()
    time_end = time.time()
    print(f"Substrate prediciton done! {(time_end-time_start)}s")
    sub_pred.result_print_to_file()
    
    if sub_pred.odbcan_sub:
        sub_pred.dbcan_sub_intermediate_file()
    time_end = time.time()
    
    ### plot the syntenic block 
    if args.cgc_substrate:
        os.chdir(args.workdir)
    if args.db_dir.startswith("/"):
        plot_command = f"syntenic_plot syntenic_plot -b PUL_blast.out --cgc cgc_standard.out -i {args.out} --db {args.db_dir}"
    else:
        plot_command = f"syntenic_plot syntenic_plot -b PUL_blast.out --cgc cgc_standard.out -i {args.out} --db ../{args.db_dir}"
    #print command
    print(plot_command)
    os.system(plot_command)
    print(f"All done! {(time_end-time_start)}s")

class HitParamter(object):
    '''
    design for parameters, how to identify real homologous hit for genes in GCG
    '''
    def __init__(self,args):
        self.upghn = args.uniq_pul_gene_hit_num
        self.uqcgn = args.uniq_query_cgc_gene_num
        self.cpn = args.CAZyme_pair_num
        self.tpn = args.total_pair_num
        self.ept = args.extra_pair_type.split(",") if args.extra_pair_type else None
        self.eptn = args.extra_pair_type_num.split(",") if self.ept else 0
        self.identity_cutoff = args.identity_cutoff
        self.coverage_cutoff  = args.coverage_cutoff
        self.bitscore_cutoff = args.bitscore_cutoff
        self.evalue_cutoff = args.evalue_cutoff
        
        ### check the additional requires for the other signature pairs
        if self.ept and len(self.ept) != len(self.eptn):
            print(f"The optional chocices of {self.ept} is not equal to {self.eptn}.",file=sys.stderr)
            exit()

    def __repr__(self):
        return "\n".join([name + ": " +str(self.__dict__[name]) for name in self.__dict__])

class dbcan_sub_parameter(object):
    def __init__(self,args):
        self.hmmevalue = args.hmmevalue
        self.hmmcov    = args.hmmcov
        self.num_of_protein_shared_substrate_cutoff = args.num_of_protein_substrate_cutoff
        self.num_of_domains_substrate_cutoff = args.num_of_domains_substrate_cutoff
        self.dbcan_substrate_scors =  args.substrate_scors
    
    def __repr__(self):
        return "\n".join([name + ": " +str(self.__dict__[name]) for name in self.__dict__])

def clean_sub(sub):
    subs = sub.Substrate
    subs = subs.replace("and",",") ### based some naming habit by someone
    subss = subs.split(",")
    subss = set(subss)
    tmp_subs = []
    for tmp_sub in subss: ### loop for substrate
        if not tmp_sub: ### exclude "" substrate come from "-"
            continue
        tmp_sub = tmp_sub.strip(" ") ### remove blank in the two ends of substrate, 
        tmp_subs.extend(tmp_sub.split("|")) ### some substrates combined by "|"
    return list(set(tmp_subs))

def clean_EC(sub):
    subs = sub.Subfam_EC
    tmp_subs = []
    for tmp_sub in subs.split("|"): ### loop for substrate
        if not tmp_sub or tmp_sub =="-": ### exclude "" substrate come from "-"
            continue
        tmp_sub = tmp_sub.strip(" ") ### remove blank in the two ends of substrate, 
        tmp_subs.append(tmp_sub.split(":")[0])
    return list(set(tmp_subs))

def parse_argv():
    parser = argparse.ArgumentParser(description='run_dbCAN substrate prediction.')
    parser.add_argument('function', help='what function will be used to analyze.')
    group = parser.add_argument_group('general optional arguments')
    group.add_argument('-i','--input',help="input file: dbCAN3 output folder")
    group.add_argument('--cgc')
    group.add_argument('--pul',help="dbCAN-PUL PUL.faa")
    group.add_argument('-f','--fasta')
    group.add_argument('-b','--blastp')
    group.add_argument('-o','--out',default="substrate.out")
    group.add_argument('-w','--workdir',type=str,default=".")
    group.add_argument('-rerun','--rerun',type=bool,default=False)
    group.add_argument('-env','--env',type=str,default="local")
    group.add_argument('-odbcan_sub','--odbcan_sub', help="output dbcan_sub prediction intermediate result?")
    group.add_argument('-odbcanpul','--odbcanpul',type=bool,default=True,help="output dbCAN-PUL prediction intermediate result?")
    parser.add_argument('--db_dir', default="db", help='Database directory')
    
    ### paramters to identify a homologous PUL
    ### including blastp evalue,number of CAZyme pair, number of pairs, extra pair, bitscore_cutoff, uniq query cgc gene hits.
    ### uniq PUL gene hits. identity cutoff. query coverage cutoff.
    group1 = parser.add_argument_group('dbCAN-PUL homologous conditons', 'how to define homologous gene hits and PUL hits')
    group1.add_argument('-upghn','--uniq_pul_gene_hit_num',default = 2,type=int)
    group1.add_argument('-uqcgn','--uniq_query_cgc_gene_num',default = 2,type=int)
    group1.add_argument('-cpn','--CAZyme_pair_num',default = 1,type=int)
    group1.add_argument('-tpn','--total_pair_num',default = 2,type=int)
    group1.add_argument('-ept','--extra_pair_type',default = None,type=str,help="None[TC-TC,STP-STP]. Some like sigunature hits")
    group1.add_argument('-eptn','--extra_pair_type_num',default ="0",type=str,help="specify signature pair cutoff.1,2")
    group1.add_argument('-iden','--identity_cutoff',default = 0.,type=float,help="identity to identify a homologous hit")
    group1.add_argument('-cov','--coverage_cutoff',default = 0.,type=float,help="query coverage cutoff to identify a homologous hit")
    group1.add_argument('-bsc','--bitscore_cutoff',default = 50,type=float,help="bitscore cutoff to identify a homologous hit")
    group1.add_argument('-evalue','--evalue_cutoff',default = 0.01,type=float,help="evalue cutoff to identify a homologous hit")

    group2 = parser.add_argument_group('dbCAN-sub conditons', 'how to define dbsub hits and dbCAN-sub subfamily substrate')
    group2.add_argument('-hmmcov','--hmmcov',default = 0.,type=float)
    group2.add_argument('-hmmevalue','--hmmevalue',default = 0.01,type=float)
    group2.add_argument('-ndsc','--num_of_domains_substrate_cutoff',default = 2,type=int,help="define how many domains share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-npsc','--num_of_protein_substrate_cutoff',default = 2,type=int,help="define how many sequences share substrates in a CGC, one protein may include several subfamily domains.")
    group2.add_argument('-subs','--substrate_scors',default = 2,type=int,help="each cgc contains with substrate must more than this value")
    
    args = parser.parse_args()
    return args

def main():
    args = parse_argv()
    cgc_substrate_prediction(args)

if __name__=="__main__":
    args = parse_argv()
    if args.function == "cgc_substrate_prediction":
        # python3 cgc_substrate_prediction.py cgc_substrate_prediction -i output -cpn 0 -upghn 1 -uqcgn 1 -bsc 20 
        # python3 cgc_substrate_prediction.py cgc_substrate_prediction -i output -cpn 2 -cov [0-1] -ept "TC-TC" -eptn 2
        # python3 /array1/www/dbCAN3/ty/cgc_substrate_prediction.py cgc_substrate_prediction -i /array1/www/dbCAN3/data/blast/20221121163906 -w /array1/www/dbCAN3/data/blast/20221121163906 -env server 
        cgc_substrate_prediction(args)
