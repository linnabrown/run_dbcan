#Runs parallel_group_many_proteins_many_patterns.rb for many organisms and types of proteins
from subprocess import call
import os
import sys

###
#python train_many_organisms_many_families.py [inputFolder] [threads] [hits] [freq]
###

organism_array = ["Chaetomium_globosum_cbs_148_51"]
#Chaetomium_globosum_cbs_148_51 fungus_fungus

cazyme_array = ["CE", "GH", "AA", "PL", "GT", "CBM"]
threads = 1
peptide_length = 6 #length of conserved peptides
hit_cut_off = 3 #number of conserved peptides necessary to classify a protein
freq_cut_off = 1.0 #minimum sum of frequencies necessary to classify a protein
list_multidomain_enzymes = "no"

if len(sys.argv) > 1:
	organism_array = [sys.argv[1]]
if len(sys.argv) > 2:
	threads = int(sys.argv[2])
if len(sys.argv) > 4:
	hit_cut_off = int(sys.argv[3])
	freq_cut_off = float(sys.argv[4])
try:
	call(['rm', organism_array[0]+'/Results/output.txt'])
except:
	pass
for protein_dir_name in organism_array:
	print("Screening "+protein_dir_name+" for")
	for cazy_class in cazyme_array:
		print(cazy_class)
		peptide_dir_name = "CAZY_PPR_patterns/"+cazy_class
		variables =  [threads, protein_dir_name, peptide_dir_name, peptide_length, hit_cut_off, freq_cut_off]
		call("python parallel_group_many_proteins_many_patterns_noDNA.py "+" ".join(str(x) for x in variables), shell=True)
		#call(["python", "add_functions_orf.py", protein_dir_name, peptide_dir_name])
		var1 = 1
		while var1 <= threads:
			try:
				os.remove(protein_dir_name+"/thread"+str(var1)+".txt")
			except:
				pass
			var1 += 1
			
	if list_multidomain_enzymes == "yes":
		call("python list_multidomain_proteins.py "+protein_dir_name+" "+"_".join(cazyme_array), shell=True)
print("\nScreened\n"+"\n".join(organism_array))
print("for proteins of the types\n"+", ".join(cazyme_array))
