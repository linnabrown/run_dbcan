import os
import sys

def HLError(mess):
    return f"\033[1;31;40m{mess}:\033[0m"
### paf record sample
'''
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)
13      attribute

'''
### diamond 

'''
1.  qseqid      query or source (gene) sequence id
2.  sseqid      subject or target (reference genome) sequence id
3.  pident      percentage of identical positions
4.  length      alignment length (sequence overlap)
5.  mismatch    number of mismatches
6.  gapopen     number of gap openings
7.  qstart      start of alignment in query
8.  qend        end of alignment in query
9.  sstart      start of alignment in subject
10.  send        end of alignment in subject
11.  evalue      expect value
12.  bitscore    bit score
13.  qlen        query length
14.  slen       subject length
'''

def CAZy_filter(cazy):
    return set([aa for aa in cazy])
    #return set([aa.split("_")[0] for aa in cazy])

### need to convert to blastp 6
class PafRecord(object):
    def __init__(self,lines):
        self.Qsn = lines[0]
        self.Qsl = lines[12]
        self.Qs  = int(lines[6]) -1
        self.Qe  = lines[7]
        self.Strand = lines[4]
        self.Tsn = lines[1]
        self.Tsl = lines[13]
        self.Ts  = int(lines[8]) -1
        self.Te  = lines[9]
        self.Nrm = lines[11]
        self.Abl = lines[3]
        self.Mq  = lines[10] ### if the paf was converted from sam, Mq here stands for the MAPQ 
        ### deal information
        self.SeqID = self.Tsn.split('|')[0]
        self.CAZys = CAZy_filter(self.Tsn.strip("|").split("|")[1:]) ### seqid|cazy1|cazy2|...| ## not subfamily
        self.UniReadId = lines[0].split(".")[0]
    def __str__(self):
        return "\t".join([str(getattr(self, value)) for value in vars(self) if value != "CAZys"])

class Paf(object):
    def __init__(self,filename):
        self.records = [PafRecord(line.split()) for line in open(filename)]
    def __iter__(self):
        return iter(self.records)
    ### get reads id
    def GetReadId(self):
        return [record.Qsn for record in self]
    ### get protein id
    def GetSeqId(self): 
        return [record.SeqID for record in self]
    ### get protein id: protein length dictory
    def GetSeqLen(self):
        return {record.SeqID:record.Tsl for record in self}
    ### get CAZy family id 2 protein id: one-many
    def CAZy2SeqID(self,CazySeqId):
        for record in self:
            for cazy in record.CAZys:
                CazySeqId.setdefault(cazy,[]).append(record.SeqID)
    ## get protein id 2 read is: one-many
    def SeqID2ReadID(self,aa):
        for record in self:
            aa.setdefault(record.SeqID,[]).append(record.Qsn)
    def ReadID2Record(self):
        return {record.Qsn:record for record in self}
    def Output(self):
        [print (record) for record in self]
    ## the CAZy information for megahit are not Qsn instead of they are in the 
    def Assign_CAZy_megahit(self):
        for cazy in self:
            cazy.CAZys = CAZy_filter(cazy.Qsn.strip("|").split("|")[1:])
    def Assign_subfam(self,CAZyID2subfam):
        for hit in self:
            hit.subfams = CAZyID2subfam.get(hit.Tsn,"")
    def Get_subfam2SeqID(self,subfam2SeqID):
        for record in self:
            for cazy in record.subfams:
                subfam2SeqID.setdefault(cazy,[]).append(record.SeqID)

def CAZyReadCount(cazyid,cazy2seqid,readtable):
    tmp_sum = 0
    for seqid in cazy2seqid[cazyid]:
        tmp_sum += readtable[seqid]
    return tmp_sum

def FPKMToCsv(args,tool,cazyfpkm,readtable,cazy2seqid):
    outfilename = args.output 
    with open(outfilename,'w') as f:
        f.write(f"Family\tAbundance\tSeqNum\tReadCount\n")
        for cazyid in cazyfpkm:
            seqnum = len(cazy2seqid[cazyid])
            readcount = CAZyReadCount(cazyid,cazy2seqid,readtable)
            fpkm = cazyfpkm[cazyid]
            if not cazyid[0].isdigit():
                f.write(f"{cazyid}\t{fpkm}\t{seqnum}\t{readcount}\n")

def check_read_type(filename):
    if filename.endswith("fq") or filename.endswith("fq.gz"):
        return "fq"
    elif filename.endswith("fa") or filename.endswith("fa.gz"):
        return "fa"
    else:
        sys.stderr.write(HLError("Error") + " File type not supported, please provide .fa(fa.gz) or (fq)fq.gz reads file.\n")
        exit(1)

def get_count_reads(file):
    if file.endswith("fq.gz"):
        r = os.popen("zcat " + file + " | echo $((`wc -l`/4))")
    elif filename.endswith(".fq"):
        r = os.popen("cat " + file + " | echo $((`wc -l`/4))")
    elif file.endswith("fa.gz"):
        r = os.popen("zcat " + file + " | grep '>' " + " | wc -l")
    elif filename.endswith(".fa"):
        r = os.popen("grep '>' " + file + " | wc -l")
    text = r.read()
    r.close()
    return float(text)

### need to modify codes for single-end sequencing

def diamond_unassemble_data(args):
    check_read_type(args.raw_reads)
    paf1 = Paf(args.paf1)
    if args.paf2:
        paf2 = Paf(args.paf2)
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber)*2
    ### FPKM or TPM is based on args.normalized
    cazyfpkm,readtable,cazy2seqid = Cal_FPKM(paf1,paf2,totalreadnumber,args.normalized)
    FPKMToCsv(args,"Diamond",cazyfpkm,readtable,cazy2seqid)

def diamond_filter(args):
    print_seqids = {}
    for line in open(args.paf1):
        lines = line.split()
        if lines[0] not in print_seqids:
            print(line.rstrip("\n"))
            print_seqids[lines[0]] = 1

def getSeqlen(paf1,paf2):
    x = paf1.GetSeqLen()
    y = paf2.GetSeqLen()
    return merge_two_dicts(x,y)

def getCazySeqId(paf1,paf2):
    cazy2seqid = {}
    paf1.CAZy2SeqID(cazy2seqid)
    paf2.CAZy2SeqID(cazy2seqid)
    for cazy in cazy2seqid:
        cazy2seqid[cazy] = set(cazy2seqid[cazy])
    return cazy2seqid

def get_subfam2seqid(paf1,paf2):
    subfam2seqid = {}
    paf1.Get_subfam2SeqID(subfam2seqid)
    paf2.Get_subfam2SeqID(subfam2seqid)
    for subfam in subfam2seqid:
        subfam2seqid[subfam] = set(subfam2seqid[subfam])
    return subfam2seqid

def getSeqReadID(paf1,paf2):
    seqid2readid = {}
    paf1.SeqID2ReadID(seqid2readid)
    paf2.SeqID2ReadID(seqid2readid)
    return seqid2readid

def SeqReadCount(seqid2readid):
    ## 0.5 two reads should be one because of the input is pair end
    return{seqid:len(seqid2readid[seqid]) for seqid in seqid2readid}

def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def SequenceFPKM(readtable,seq2len,totalreadnumber):
    seqfpkm = {}
    for seqid in readtable:
        tmp_total_read = float(totalreadnumber)/pow(10,6)
        tmp_trans_len  = float(seq2len[seqid])/1000
        read_count = float(readtable[seqid])
        tmp_fpkm = read_count/tmp_total_read/tmp_trans_len
        #print(seqid,totalreadnumber,seq2len[seqid],read_count)
        seqfpkm[seqid] = tmp_fpkm 
    return seqfpkm

###           Ni/Li*10^6
###  TPM =    ------------------------------
###           sum(N1/L1+N2/L2 + ... + Nn/Ln)

### but the TPM cann't not applied in assembly-free, because it will increase the value
### very limit reads can align to CAZyme
def SequenceTPM(readtable,seq2len,totalreadnumber):
    seqtpm = {}
    normalized_tpm = 0.
    ### calculate normalized_tpm
    for seqid in readtable:
        read_count = float(readtable[seqid])
        seqlen = float(seq2len[seqid])
        normalized_tpm += read_count/seqlen
    ### calculate tpm
    for seqid in readtable:
        read_count = float(readtable[seqid])
        seqlen = float(seq2len[seqid])
        normalized_reads_counts = read_count/seqlen*pow(10,6)
        tmp_seqtpm = normalized_reads_counts/normalized_tpm
        seqtpm[seqid] = tmp_seqtpm
    return seqtpm

###         read count *10^6
###  RPM = ------------------------
###         total read count(mapped)

def SequenceRPM(readtable,seq2len,totalreadnumber):
    seqrpm = {}
    for seqid in readtable:
        read_count = float(readtable[seqid])
        seqlen = float(seq2len[seqid])
        rpm    = read_count*pow(10,6)/totalreadnumber
        seqrpm[seqid] = rpm
    return seqrpm

def CAZyFPKM(seqfpkm,cazy2seqid): ## named as FPKM, but can apply in tpm
    cazyfpkm = {}
    for cazy in cazy2seqid:
        tmp_fpkm = 0.
        for seqid in cazy2seqid[cazy]:
            tmp_fpkm += float(seqfpkm[seqid])
        cazyfpkm[cazy] = tmp_fpkm
    return cazyfpkm

def Cal_FPKM(paf1,paf2,totalreadnumber,normalized):
    ## get sequence length from paf
    seq2len = getSeqlen(paf1,paf2)
    
    # get CAZy family to seq mapping table: CAZy ID 2 protein ID
    cazy2seqid = getCazySeqId(paf1,paf2)
    # outdict_list(cazy2seqid)

    ## get SeqID2ReadID to generate mapping table: protein ID 2 read ID
    seqid2readid = getSeqReadID(paf1,paf2)
    ### read table: protein ID 2 read count 
    readtable = SeqReadCount(seqid2readid)
    ## outdict(readtable)
    ## calculate fpkm for each protein seq
    if normalized == "FPKM":
        seqfpkm = SequenceFPKM(readtable,seq2len,totalreadnumber)
    elif normalized == "RPM":
        seqfpkm = SequenceRPM(readtable,seq2len,totalreadnumber)
    else:
        seqfpkm = SequenceTPM(readtable,seq2len,totalreadnumber)
    ## outdict(seqfpkm)
    cazyfpkm = CAZyFPKM(seqfpkm,cazy2seqid)
    #outdict(cazyfpkm)
    return cazyfpkm,readtable,cazy2seqid

import argparse

## dict: str -> []

def read_EC2substrate_table(args):
    famEC2substrate = {}
    map_table = f"{args.db}fam-substrate-mapping.tsv"
    map_table_lines = open(map_table).readlines()
    for line in map_table_lines[1:]:
        lines = line.rstrip("\n").split("\t")
        substrates = [sub_tmp.strip(" ") for sub_tmp in lines[0].strip().replace("and","").split(',')]
        #famEC2substrate.setdefault(lines[2],[]).extend(substrates)
        famEC2substrate.setdefault(lines[-1],[]).extend(substrates)
        #famEC2substrate[lines[-1]] = lines[0]
    for fam in famEC2substrate:
        famEC2substrate[fam] = list(set(famEC2substrate[fam]))
    return famEC2substrate

### each protein may includes more than 1 eCAMI subfam
def read_CAZyID2subfam_table(args):
    CAZyID2subfam = {}
    map_table = f"{args.db}CAZyID_subfam_mapping.tsv"
    map_table_lines = open(map_table).readlines()
    for line in map_table_lines:
        lines = line.rstrip("\n").split("\t")
        CAZyID2subfam.setdefault(lines[-1],[]).append(lines[0])
    return CAZyID2subfam

def read_subfam2ECosub_table(args):
    subfam2EC = {};subfam2subtrate = {}
    map_table = f"{args.db}subfam_EC_mapping.tsv"
    map_table_lines = open(map_table).readlines()
    for line in map_table_lines:
        lines = line.rstrip("\n").split("\t")
        if lines[-1] != "-":
            substrates = [sub.strip() for sub in lines[-1].strip().replace("and","").split(",")]
            subfam2subtrate.setdefault(lines[0],[]).extend(substrates)
        if lines[1] != "-":
            subfam2EC.setdefault(lines[0],[]).append(lines[1])
    
    for subfam in subfam2EC:
        subfam2EC[subfam] = list(set(subfam2EC[subfam]))
    for subfam in subfam2subtrate:
        subfam2subtrate[subfam] = list(set(subfam2subtrate[subfam]))
    
    ### dict, sub -> []
    return subfam2EC,subfam2subtrate


def diamond_EC_abund(args):
    if not args.db.endswith("/"):
        args.db += "/"
    subfam2EC,subfam2subtrate = read_subfam2ECosub_table(args)
    
    EC2Abund = {} ; EC2subfam = {}
    for line in open(args.input):
        subfam,FPKM,ReadCount,SeqNum = line.rstrip("\n").split("\t")
        if subfam in subfam2EC:
            ECs = subfam2EC[subfam]
            for EC in ECs:
                subfams = EC2subfam.get(EC,[])
                if subfam not in subfams:
                    EC2subfam.setdefault(EC,[]).append(subfam)
                    EC2Abund.setdefault(EC,[]).append(float(FPKM))
    
    outfilename = args.output 
    with open(outfilename,'w') as f:
        f.write("EC\tAbundance\tsubfam\n")
        for sub in EC2Abund:
            f.write(sub+"\t"+str(sum(EC2Abund[sub]))+"\t"+";".join(EC2subfam[sub])+"\n")

def CAZyme_substrate(args):
    if not args.db.endswith("/"):
        args.db += "/"
    
    EC2substrate = read_EC2substrate_table(args)
    subfam2EC,subfam2subtrate = read_subfam2ECosub_table(args)
    
    Sub2Abund = {}; Sub2subfam = {}
    for line in open(args.input):
        #subfam,FPKM,ReadCount,SeqNum = line.rstrip("\n").split("\t")
        #Subfamily	Abundance	SeqNum	ReadCount
        subfam,FPKM,SeqNum,ReadCount = line.rstrip("\n").split("\t")
        ### route1, subfam->EC->substrate
        if subfam in subfam2EC:
            ECs = subfam2EC[subfam]
            if ECs:
                for EC in ECs:
                    substrates = EC2substrate.get(EC,"")
                    if substrates:
                        for sub in substrates:
                            subfams = Sub2subfam.get(sub,[])
                            if subfam not in subfams:
                                Sub2Abund.setdefault(sub,[]).append(float(FPKM))
                                Sub2subfam.setdefault(sub,[]).append(subfam)
        ### route2, subfam -> substrate
        substrates = subfam2subtrate.get(subfam,"")
        if substrates:
            for sub in substrates:
                subfams = Sub2subfam.get(sub,[])
                if subfam not in subfams:
                    Sub2Abund.setdefault(sub,[]).append(float(FPKM))
                    Sub2subfam.setdefault(sub,[]).append(subfam)
    
    outfilename = args.output 
    with open(outfilename,'w') as f:
        f.write("Substrate\tAbundance\tsubfam\n")
        for sub in Sub2Abund:
            f.write(sub+"\t"+str(sum(Sub2Abund[sub]))+"\t"+";".join(Sub2subfam[sub])+"\n")

def Cal_subfam_FPKM(paf1,paf2,totalreadnumber,normalized):
    ## get sequence length from paf
    seq2len = getSeqlen(paf1,paf2)
    
    # get CAZy family to seq mapping table: CAZy ID 2 protein ID
    #cazy2seqid = getCazySeqId(paf1,paf2)
    subfam2seqid = get_subfam2seqid(paf1,paf2)
    # outdict_list(cazy2seqid)

    ## get SeqID2ReadID to generate mapping table: protein ID 2 read ID
    seqid2readid = getSeqReadID(paf1,paf2)
    ### read table: protein ID 2 read count 
    readtable = SeqReadCount(seqid2readid)
    ## outdict(readtable)
    ## calculate fpkm for each protein seq
    if normalized == "FPKM":
        seqfpkm = SequenceFPKM(readtable,seq2len,totalreadnumber)
    elif normalized == "RPM":
        seqfpkm = SequenceRPM(readtable,seq2len,totalreadnumber)
    else:
        seqfpkm = SequenceTPM(readtable,seq2len,totalreadnumber)
    ## outdict(seqfpkm)
    cazyfpkm = CAZyFPKM(seqfpkm,subfam2seqid)
    #outdict(cazyfpkm)
    return cazyfpkm,readtable,subfam2seqid


def diamond_subfam_abund(args):
    if not args.db.endswith("/"):
        args.db += "/"
    check_read_type(args.raw_reads)
    ### FPKM or TPM is based on args.normalized
    CAZyID2subfam = read_CAZyID2subfam_table(args)
    paf1 = Paf(args.paf1)
    if args.paf2:
        paf2 = Paf(args.paf2)
    paf1.Assign_subfam(CAZyID2subfam)
    paf2.Assign_subfam(CAZyID2subfam)
    totalreadnumber = get_count_reads(args.raw_reads)
    if args.paf2:
        totalreadnumber = float(totalreadnumber)*2
    
    subfamfpkm,readtable,subfam2seqid = Cal_subfam_FPKM(paf1,paf2,totalreadnumber,args.normalized)
    FPKMToCsv(args,"Diamond",subfamfpkm,readtable,subfam2seqid)

def arg_parse():
    parser = argparse.ArgumentParser(description='diamond assembly-free method')
    parser.add_argument('function', help='which function will be used.',choices=["diamond_fam_abund","diamond_substrate_abund","diamond_subfam_abund","diamond_EC_abund"])
    parser.add_argument('-paf1', type=str, help='R1 reads diamond blastx',default="")
    parser.add_argument('-paf2', type=str, help='R2 reads diamond blastx',default="")
    parser.add_argument('-i','--input', type=str)
    parser.add_argument('-d','--db', type=str,default="./db")
    parser.add_argument('-o','--output', type=str ,default="asmfree_fam_abund")
    parser.add_argument('--raw_reads', type=str ,default=" ",help="compress or uncompress fq/fa type of reads.")
    parser.add_argument('-n','--normalized', type=str ,help="FPKM, TPM, RPM",default = "TPM",choices=['FPKM', 'RPM', 'TPM'])
    return parser.parse_args()

def main():
    args = arg_parse()
    if args.function == "diamond_fam_abund":
        ### dbcan_asmfree diamond_fam_abund -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx --raw_reads Dry2014_1_val_1.fq.gz -n FPKM -o Dry2014_fam_abund
        ### dbcan_asmfree diamond_fam_abund -paf1 Wet2014_1.blastx -paf2 Wet2014_2.blastx --raw_reads Wet2014_1_val_1.fq.gz -n FPKM -o Wet2014_fam_abund
        diamond_unassemble_data(args)
    if args.function == "diamond_subfam_abund":
        ### dbcan_asmfree diamond_subfam_abund -paf1 Dry2014_1.blastx -paf2 Dry2014_2.blastx --raw_reads Dry2014_1_val_1.fq.gz -o Dry2014_subfam_abund -n FPKM
        ### dbcan_asmfree diamond_subfam_abund -paf1 Wet2014_1.blastx -paf2 Wet2014_2.blastx --raw_reads Wet2014_1_val_1.fq.gz -o Wet2014_subfam_abund -n FPKM
        diamond_subfam_abund(args)
    if args.function == "diamond_EC_abund":
        ### dbcan_asmfree diamond_EC_abund -i Dry2014_subfam_abund -o Dry2014_EC_abund
        ### dbcan_asmfree diamond_EC_abund -i Wet2014_subfam_abund -o Wet2014_EC_abund
        diamond_EC_abund(args)
    if args.function == "diamond_substrate_abund":
        ### dbcan_asmfree diamond_substrate_abund -i Dry2014_subfam_abund -o Dry2014_substrate_abund
        ### dbcan_asmfree diamond_substrate_abund -i Wet2014_subfam_abund -o Wet2014_substrate_abund
        CAZyme_substrate(args)
if __name__== "__main__":
    main()
