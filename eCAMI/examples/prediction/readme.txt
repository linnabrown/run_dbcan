The input is a fasta file with sequences to be annotated for an EC or CAZyme family.

The output has one single file mimicing a fasta format (multi-family proteins can occur multiple times in multiple lines in this file):
>protein_name    fam_name:group_number    subfam_name_of_the_group:subfam_name_count
matching kmers(start position in the query sequence)

The protein_name is the protein ID provided in the input fasta. 
The fam_name:group_number is the CAZyme family or EC number assigned by eCAMI by comparing the query protein against the k-mer peptide library of CAZymes or ECs (fam_name is the CAZyme or EC cluster from the library, and group_number is the cluster/subfam number/ID of that CAZyme or EC cluster). 
The subfam_name_of_the_group:subfam_name_count is extracted from the k-mer peptide library of CAZymes or ECs (subfam_name_of_the_group is the CAZy subfam or EC at the 4th level, and subfam_name_count is the count).
