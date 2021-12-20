The input is a fasta file to be clustered. The known label of each sequence can be added after the sequence ID and separated with a '|', e.g., the EC numbers or CAZy subfamilies known for this protein. 

Two files are output from CAMI clustering for each cluster: (i) fasta sequence, and (ii) characteristic k-mer peptides. 

For (ii), the output file has its first line being the cluster number/ID (starting from 0), 
and the second line has two tab separated fields: 1 (e.g., GH5:41) is the CAZyme famID or EC at the 3rd level, colon, and then the count of proteins; 2 (e.g., GH5_2:40|CBM28:31|CBM17:25|3.2.1.4:15|GH5:1) is a string separated by |, and each part has the known label (e.g., the labels will be the EC numbers or CAZy subfamilies provided in the input fasta file), colon, the counts of proteins; if there are multiple lables, they are separated |, 
and the third line includes each k-mer peptide followed by their frequency.
