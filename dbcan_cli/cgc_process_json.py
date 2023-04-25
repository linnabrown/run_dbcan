##########################
# to generate json file for all cgc_stardard.out file from run_dbcan
# use: python cgc_process_json.py -i cgc_standard.out -o cgc_standard.out.json
# written by Roland Madadjim in Cui's lab at Soc, UNL
# last updated: 12/09/2022
##########################

#from __future__ import print_function
import os
import json
import time
import argparse
import pandas as pd
import numpy as np


class PrePro:

    def __init__(self, data):
        self.df = data


    def extract_gs(self, dataList):
        i = 0
        geneL = []
        gene = list(map(lambda e : "{}".format(e['Gene_Type']),dataList))
        pfam = list(map(lambda e : "{}".format(e['Protein_Family']),dataList))
        while i < len(dataList):  
            if (gene[i] == 'CAZyme'):
                s = pfam[i]
                geneL.append(s)
            elif (gene[i] == 'TC'):
                s = pfam[i]
                geneL.append(s)
            elif (gene[i] == 'TF'):
                s = pfam[i]
                geneL.append(s)
            elif (gene[i] == 'STP'):
                s = pfam[i]
                geneL.append(s)
            elif (gene[i] == 'Null_Type'):
                s = 'NA'
                geneL.append(s)
            i=i+1
        gene_st = '-'.join(geneL)
        return gene_st


    def pul_section(self):
        for (cgc_id), df_pul_grouped in self.df.groupby("cgc_id"):
            datalist = list(self.cluster_section(df_pul_grouped))
            gene_str = self.extract_gs(datalist)
            yield {cgc_id: {
                "changelog": [],
                "Cluster_ID": cgc_id,
                "Gene_String":gene_str,
                "Contig_ID": self.df.loc[self.df['cgc_id'] == cgc_id, 'contig_id'].iloc[0],
                "ncbi_species_tax_id": [],
                "organism_name": [],
                "publication": [],
                "Protein": list(self.cluster_section(df_pul_grouped))
                #"dbCan_Pul_accession": ID, # as string or integer?
                #"publication": df.loc[df['ID'] == ID, 'PMID'].iloc[0],
                }
            }

            
    def cluster_section(self, df_pul_grouped):
        for (contig_id,gene_type,contig_id,protein_id,gene_start,gene_stop,direction,protein_family), df_puls in df_pul_grouped.groupby(
                ["contig_id","gene_type","contig_id","protein_id","gene_start","gene_stop","direction","protein_family"]
        ):
            yield {
                "protein_id": protein_id,
                "Gene_Type": gene_type,
                "Gene_Start": gene_start,
                "Gene_Stop": gene_stop,
                "Strand": direction,
                "Protein_Family": protein_family,
            }
            
    def run_dbCan_section(self, df_puls):
        for row in df_puls.itertuples():
            yield {
                #"Gene_Type": row.gene_type,
                #"Gene_Start": row.gene_start,
                #"Gene_Stop": row.gene_stop,
                #"Strand": row.direction,
            # "Protein_Family": row.protein_family
            }

def file_ext(choices,fname):
    ext = os.path.splitext(fname)[1][1:]
    if ext not in choices:
        parser.error("File needs to be a .out or .csv")
    return fname

class CustomEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return str(obj)
        return super().default(obj)

def main():
    parser = argparse.ArgumentParser(description='Compiling Json from cgc_standard.out')
    parser.add_argument('-i',required=True,help='path to output file (cgc_standard.out) file', type=lambda s:file_ext(("out","csv"),s))
    parser.add_argument('-o','--output')
    args = parser.parse_args()

    with open(args.i) as file:  ### input files
        data = pd.read_csv(file,  sep ='\t')
        data.rename(columns = {'CGC#':'cgc_id','Gene Type':'gene_type','Contig ID':'contig_id','Protein ID':'protein_id','Gene Start':'gene_start','Gene Stop':'gene_stop','Direction':'direction','Protein Family':'protein_family'}, inplace = True)
        data['gene_type'].fillna('Null_Type', inplace=True)
        data['protein_family'].fillna('0', inplace=True)
        p = PrePro(data)
    
    pul_list = list(p.pul_section())
    pul_dict = {}
    for sub_dict in pul_list:
        pul_dict.update(sub_dict)
    jsonPuls = json.dumps(pul_dict, indent=4, cls=CustomEncoder)
    
    with open(args.output,"w") as outfile:
    #with open("Json"+time.strftime("%Y%m%d%H%M%S")+".json","w") as outfile:
        outfile.write(jsonPuls)

if __name__ == "__main__":
    main()
