##################################################
# dbCAN2 Driver Script (Stand Alone Version)
#
# Written by Tanner Yohe in the Yin Lab at NIU
# Revised by Huang Le in the Zhang Lab at NKU
#
# Last updated 4/23/18
#
# Accepts user input
# Predicts genes if needed
# Runs input against HMMER, DIAMOND, and Hotpep
# Optionally predicts CGCs with CGCFinder
#
####################################################
from subprocess import Popen, call, check_output
import os
import argparse
import sys


parser = argparse.ArgumentParser(description='dbCAN2 Driver Script')

parser.add_argument('inputFile', help='User input file. Must be in FASTA format.')
parser.add_argument('inputType', choices=['protein', 'prok', 'meta'], help='Type of sequence input. protein=proteome; prok=prokaryote; meta=metagenome') #protein=proteome, prok=prokaryote nucleotide, meta=metagenome nucleotide
parser.add_argument('--cluster', '-c', help='Predict CGCs via CGCFinder. This argument requires an auxillary locations file if a protein input is being used')
parser.add_argument('--dia_eval', default=1e-121,type=float, help='DIAMOND E Value')
parser.add_argument('--dia_cpu', default=5, type=int, help='Number of CPU cores that DIAMOND is allowed to use')
parser.add_argument('--hmm_eval', default=1e-15, type=float, help='HMMER E Value')
parser.add_argument('--hmm_cov', default=0.35, type=float, help='HMMER Coverage value')
parser.add_argument('--hmm_cpu', default=1, type=int, help='Number of CPU cores that HMMER is allowed to use')
parser.add_argument('--hot_hits', default=4, type=int, help='Hotpep Hit value')
parser.add_argument('--hot_freq', default=2.0, type=float, help='Hotpep Frequency value')
parser.add_argument('--hot_cpu', default=4, type=int, help='Number of CPU cores that Hotpep is allowed to use')
parser.add_argument('--out_pre', default="", help='Output files prefix')
parser.add_argument('--out_dir', default="", help='Output directory')
parser.add_argument('--db_dir', default="db/", help='Database directory')
parser.add_argument('--cgc_dis', default=2, help='CGCFinder Distance value')
parser.add_argument('--cgc_sig_genes', default='tp', choices=['tp', 'tf','all'], help='CGCFinder Signature Genes value')
parser.add_argument('--tools', '-t', nargs='+', choices=['hmmer', 'diamond', 'hotpep', 'all'], default='all', help='Choose a combination of tools to run')

args = parser.parse_args()

####
#python run_dbcan.py [inputFile] [inputType]
####

##########################
# Begin Setup and Input Checks

dbDir = args.db_dir
prefix = args.out_pre
outDir = args.out_dir
auxFile = ""
input = args.inputFile
inputType = args.inputType
find_clusters = False
if args.cluster != None:
	find_clusters = True
	auxFile = args.cluster
if not dbDir.endswith("/") and len(dbDir) > 0:
	dbDir += "/"
if not os.path.isdir(dbDir):
	print("ERROR: The database directory does not exist")
	exit()
if not os.path.isfile(dbDir+'CAZy.dmnd'):
	print("ERROR: No CAZy DIAMOND database found. Please make sure that your CAZy DIAMOND databased is named 'CAZy.dmnd' and is located in your database directory")
	exit()
if not os.path.isfile(dbDir+'dbCAN.txt'):
	print("ERROR: No dbCAN HMM database found. Please make sure that your dbCAN HMM database is named 'dbCAN.txt', has been through hmmpress, and is located in your database directory")
	exit()
if not outDir.endswith("/") and len(outDir) > 0:
	outDir += "/"
if not os.path.isdir(outDir):
	call(['mkdir', outDir])
if find_clusters and inputType == "protein":
	if len(auxFile) > 0:
		if not os.path.isfile(auxFile):
				print("ERROR: It seems that the auxillary filename that you provided does not exist, or is not a file")
				exit()
	else:
		print("ERROR: Please provide an auxillary input file with the position of each gene. This file can either be in BED or GFF format")
		exit()
tools = [True, True, True] #DIAMOND, HMMER, Hotpep
if args.tools != 'all':
	if 'diamond' not in args.tools:
		tools[0] = False
	if 'hmmer' not in args.tools:
		tools[1] = False
	if 'hotpep' not in args.tools:
		tools[2] = False

# End Setup and Input Checks
#########################
# Begin Gene Prediction Tools

if inputType == 'prok':
    call(['prodigal', '-i', input, '-a', outDir+prefix+'uniInput', '-o', outDir+prefix+'prodigal.gff', '-f', 'gff', '-q'])
if inputType == 'meta':
    call(['FragGeneScan1.30/run_FragGeneScan.pl', '-genome='+input, '-out='+outDir+prefix+'fragGeneScan', '-complete=1', '-train=complete', '-thread=10'])

#Frag Gene Scan
if inputType == 'meta':
    call(['cp', outDir+prefix+'fragGeneScan.faa', outDir+prefix+'uniInput'])

#Proteome
if inputType == 'protein':
    call(['cp', input, outDir+prefix+'uniInput'])

# End Gene Prediction Tools
#######################
# Begin SignalP

signalpos = Popen('signalp -t gram+ '+outDir+prefix+'uniInput > '+outDir+prefix+'signalp.neg', shell=True)
signalpneg = Popen('signalp -t gram- '+outDir+prefix+'uniInput > '+outDir+prefix+'signalp.pos', shell=True)

# End SignalP
#######################
# Begin Core Tools

if tools[0]:
	diamond = Popen(['diamond', 'blastp', '-d', dbDir+'CAZy.dmnd', '-e', str(args.dia_eval), '-q', outDir+prefix+'uniInput', '-k', '1', '-p', str(args.dia_cpu), '-o', outDir+prefix+'diamond.out', '-f', '6'])

if tools[1]:
	hmmer = Popen(['hmmscan', '--domtblout', outDir+prefix+'h.out', '--cpu', str(args.hmm_cpu), '-o', '/dev/null', dbDir+'dbCAN.txt', outDir+prefix+'uniInput'])

if tools[2]:
	count = int(check_output("tr -cd '>' < "+outDir+prefix+"uniInput | wc -c", shell=True))    #number of genes in input file
	numThreads = args.hot_cpu	    														#number of cores for Hotpep to use
	count_per_file = count/numThreads														#number of genes per core
	directory = input.split('.')[0]
	call(['mkdir','-m','777','Hotpep/'+directory])
	num_files = 1
	num_genes = 0
	out = open("Hotpep/"+directory+"/orfs"+str(num_files)+".txt", 'w')
	with open(outDir+prefix+'uniInput', 'r') as f:
		for line in f:
			if line.startswith(">"):
				num_genes += 1
				if num_genes > count_per_file and num_files != numThreads:
					out.close()
					num_files += 1
					num_genes = 0
					out = open("Hotpep/"+directory+"/orfs"+str(num_files)+".txt", 'w')
			out.write(line)
				
	os.chdir('Hotpep/')
	hotpep = Popen(['python', 'train_many_organisms_many_families.py', directory, str(numThreads), str(args.hot_hits), str(args.hot_freq)])
	os.chdir('../')
	hotpep.wait()

	hotpepDir = 'Hotpep/'+directory
	call(['mv', hotpepDir+'/Results/output.txt', outDir+prefix+'Hotpep.out'])

if tools[0]:
	print("Waiting on DIAMOND")
	diamond.wait()
	print("DIAMOND complete")
if tools[1]:
	print("Waiting on HMMER")
	hmmer.wait()
	print("HMMER complete")

	call('python hmmscan-parser.py '+outDir+prefix+'h.out '+str(args.hmm_eval)+' '+str(args.hmm_cov)+' > '+outDir+prefix+'hmmer.out', shell=True)
	call(['rm', outDir+prefix+'h.out'])

# End Core Tools
########################
# Begin Adding Column Headers

if tools[2]:
	with open(outDir+prefix+'Hotpep.out') as f:
		with open(outDir+prefix+'temp', 'w') as out:
			out.write('CAZy Family\tPPR Subfamily\tGene ID\tFrequency\tHits\tSignature Peptides\n')
			for line in f:
				out.write(line)
	call(['mv', outDir+prefix+'temp', outDir+prefix+'Hotpep.out'])
if tools[1]:
	with open(outDir+prefix+'hmmer.out') as f:
		with open(outDir+prefix+'temp', 'w') as out:
			out.write('HMM Profile\tProfile Length\tGene ID\tGene Length\tE Value\tProfile Start\tProfile End\tGene Start\tGene End\tCoverage\n')
			for line in f:
				out.write(line)
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

########################
# Begin TF and TP prediction

	call(['diamond', 'blastp', '-d', dbDir+'tf.dmnd', '-e', '1e-10', '-q', outDir+prefix+'uniInput', '-k', '1', '-p', '1', '-o', outDir+prefix+'tf.out', '-f', '6'])
	call(['diamond', 'blastp', '-d', dbDir+'tcdb.dmnd', '-e', '1e-10', '-q', outDir+prefix+'uniInput', '-k', '1', '-p', '1', '-o', outDir+prefix+'tp.out', '-f', '6'])
	tp = set()
	tf = set()
	tp_genes = {}
	tf_genes = {}
	with open(outDir+prefix+'tf.out') as f:
		for line in f:
			row = line.rstrip().split('\t')
			tf.add(row[0])
			if not row[0] in tf_genes:
				tf_genes[row[0]] = row[1]
			else:
				tf_genes[row[0]] += ','+row[1]
	with open(outDir+prefix+'tp.out') as f:
		for line in f:
			row = line.rstrip().split('\t')
			tp.add(row[0])
			if not row[0] in tp_genes:
				tp_genes[row[0]] = row[1]
			else:
				tp_genes[row[0]] += ','+row[1]
# End TF and TP prediction
##########################
# Begine CAZyme Extraction
	cazyme_genes = {}
	dia = set()
	hot = set()
	hmm = set()
	if tools[0]:
		with open(outDir+prefix+'diamond.out') as f:
			next(f)
			for line in f:
				row = line.rstrip().split('\t')
				dia.add(row[0])
				if row[0] not in cazyme_genes:
					cazyme_genes[row[0]] = set()
				cazyme_genes[row[0]].add(row[1].split('|')[1])
	if tools[1]:
		with open(outDir+prefix+'hmmer.out') as f:
			next(f)
			for line in f:
				row = line.rstrip().split('\t')
				hmm.add(row[2])
				if row[2] not in cazyme_genes:
					cazyme_genes[row[2]] = set()
				cazyme_genes[row[2]].add(row[0].split('.')[0])
	if tools[2]:
		with open(outDir+prefix+'Hotpep.out') as f:
			next(f)
			for line in f:
				row = line.rstrip().split('\t')
				hot.add(row[2])
				if row[2] not in cazyme_genes:
					cazyme_genes[row[2]] = set()
				cazyme_genes[row[2]].add(row[0])
	if tools.count(True) > 1:
		temp1 = hmm.intersection(hot)
		temp2 = hmm.intersection(dia)
		temp3 = dia.intersection(hot)
		cazyme = temp1.union(temp2, temp3)
	else:
		cazyme = hmm.union(dia, hot)
# End CAZyme Extraction
######################
# Begin GFF preperation

	if inputType == "prok":   #use Prodigal GFF output
		with open(outDir+prefix+'prodigal.gff') as f:
			with open(outDir+prefix+'cgc.gff', 'w') as out:
				for line in f:
					if not line.startswith("#"):
						row = line.rstrip().split('\t')
						num = row[-1].split(";")[0].split('_')[-1]
						gene = row[0] + '_' + num
						row[8] = ""
						if gene in tf:
							row[2] = "TF"
							row[8] = "DB="+tf_genes[gene]
						elif gene in tp:
							row[2] = "TC"
							row[8] = "DB="+tp_genes[gene]
						elif gene in cazyme:
							row[2] = "CAZyme"
							row[8] = "DB="+'|'.join(cazyme_genes[gene])
						row[8] += ";ID="+gene
						out.write('\t'.join(row)+'\n')
						
						
	elif inputType == "meta":  #use FragGeneScan GFF output
		with open(outDir+prefix+'fragGeneScan.gff') as f:
			with open(outDir+prefix+'cgc.gff', 'w') as out:
				for line in f:
					if not line.startswith("#"):
						row = line.rstrip().split('\t')
						gene = row[-1].split(";")[0].split("=")[1]
						if gene in tf:
							row[2] = "TF"
							row.insert(8, "DB="+tf_genes[gene])
						elif gene in tp:
							row[2] = "TC"
							row.insert(8, "DB="+tp_genes[gene])
						elif gene in cazyme:
							row[2] = "CAZyme"
							row.insert(8, "DB="+'|'.join(cazyme_genes[gene]))
						else:
							row.insert(8, "")
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
								note = row[8].split(";")
								gene = ""
								notes = {}
								for x in note:
									temp = x.split('=')
									notes[temp[0]] = temp[1]
								if "Name" in notes:
									gene = notes["Name"]
								elif "ID" in notes:
									gene = notes["ID"]
								else:
									continue
								if gene in tf:
									row[2] = "TF"
									row[8] = "DB="+tf_genes[gene]
								elif gene in tp:
									row[2] = "TC"
									row[8] = "DB="+tp_genes[gene]
								elif gene in cazyme:
									row[2] = "CAZyme"
									row[8] = "DB="+'|'.join(cazyme_genes[gene])
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
						row = line.rstrip().split('\t')
						outrow = ['.','.','.','.','.','.','.','.','']
						gene = row[1]
						if gene in tf:
							outrow[2] = 'TF'
							outrow[8] =  "DB="+tf_genes[gene]
						elif gene in tp:
							outrow[2] = 'TC'
							outrow[8] = "DB="+tp_genes[gene]
						elif gene in cazyme:
							outrow[2] = 'CAZyme'
							outrow[8] = "DB="+'|'.join(cazyme_genes[gene])
						else:
							outrow[2] = 'CDS'
						outrow[0] = row[0]
						outrow[3] = row[2]
						outrow[4] = row[3]
						outrow[6] = row[4]
						outrow[8] += ";ID="+gene
						out.write('\t'.join(outrow)+'\n')
# End GFF preperation
####################
# Begin CGCFinder call

	call(['python', 'CGCFinder.py', outDir+prefix+'cgc.gff', '-o', outDir+prefix+'cgc.out', '-s', args.cgc_sig_genes, '-d', str(args.cgc_dis)])

# End CGCFinder call
# End CGCFinder
####################
# Begin SignalP combination

print("Waiting on signalP")
signalpos.wait()
signalpneg.wait()
print("SignalP complete")
with open(outDir+prefix+'temp', 'w') as out:
	with open(outDir+prefix+'signalp.pos') as f:
		for line in f:
			if not line.startswith('#'):
				row = line.split(' ')
				row = [x for x in row if x != '']
				if row[9] == 'Y':
					out.write(line)
	with open(outDir+prefix+'signalp.neg') as f:
		for line in f:
			if not line.startswith('#'):
				row = line.split(' ')
				row = [x for x in row if x != '']
				if row[9] == 'Y':
					out.write(line)
call('sort -u '+outDir+prefix+'temp > '+outDir+prefix+'signalp.out', shell=True)
call(['rm', outDir+prefix+'temp', outDir+prefix+'signalp.pos', outDir+prefix+'signalp.neg'])

# End SignalP combination
#######################
# End script
