#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 20:08:33 2019
Revised on Nov 11, 2021 by Le Huang
@author: joy1314bubian
"""

import numpy as np

import os
import datetime
import re
#import psutil
import argparse,sys
from multiprocessing import Pool 
# add by Le Nov 14, 2021

# from eCAMI_data import eCAMI_data_path
#add by Le Nov 14, 2021 end
def read_file(database_dir,input_fasta_file):
    f = open(database_dir +input_fasta_file, 'r')
    text=f.read()
    f.close()
    text=text.split('\n')
    while '' in text:
        text.remove('')
    proteinname=[]
    new_text={}

    for i in range(len(text)):
        if '>' in text[i]:
            if ' ' in text[i]:
                proteinname.append(text[i].split(' ')[0])
            else:
                proteinname.append(text[i])
        else:
            if '>' not in text[i] and '>' in text[i-1]:
                new_text[proteinname[len(proteinname)-1]]=text[i]
            else:
                new_text[proteinname[len(proteinname)-1]]+=text[i]
    protein_string=[]
#    print(len(proteinname))
    for each_name in proteinname:
        protein_string.append(new_text[each_name])
    return [proteinname,protein_string]


def find_all(sub,s):
	index_list = []
	index = s.find(sub)
	while index != -1:
		index_list.append(index+1)
		index = s.find(sub,index+1)
	
	if len(index_list) > 0:
		return index_list
	else:
		return -1


def get_kmer_dict(n_mer_dir_path):
    kmer_dir_list=os.listdir(n_mer_dir_path)
    fam_kmer_dict={}
    for k in range(len(kmer_dir_list)):
        dir_name=n_mer_dir_path
        kmer_label=[]
        kmer_message=[]
        all_cluster_kmer=[]
        kmer_file_list=os.listdir(dir_name+'/'+kmer_dir_list[k]+'/kmer_for_each_cluster')
        for each_file in kmer_file_list:
            f=open(dir_name+'/'+kmer_dir_list[k]+'/kmer_for_each_cluster'+'/'+each_file,'r')
            text=f.read()
            text=text.split('\n')
            f.close()
            kmer_label.append(kmer_dir_list[k]+':'+text[0].split(':')[1])
            temp_name=((text[1].split('\t')[0]).split(':')[0]).split('_')[0]
            kmer_message.append(text[1].split('\t')[1])
            each_cluster_kmer=text[2]
            each_cluster_kmer=each_cluster_kmer.split('\t')
            while '' in each_cluster_kmer:
                each_cluster_kmer.remove('')
            cluster_kmer={}
            for each_kmer in each_cluster_kmer:
                temp_mess=each_kmer.split(':')
                cluster_kmer[temp_mess[0]]=float(temp_mess[1])
            all_cluster_kmer.append(cluster_kmer)

        temp_name=temp_name.split('.')
        while '-' in temp_name:
            temp_name.remove('-')
        temp_name='.'.join(temp_name)

        fam_kmer_dict[temp_name]=[kmer_label,kmer_message,all_cluster_kmer]
    return fam_kmer_dict

def write_line(selected_fam_message,temp_protein_string):
    string_line=''
    string_line+=str(selected_fam_message[4])+'\t'+str(selected_fam_message[3])+'\n'
    temp_same_kmer=selected_fam_message[2]
    kmer_index_dict={}
    for each_kmer in temp_same_kmer:
        index=find_all(each_kmer,temp_protein_string)
        kmer_index_dict[np.min(index)]=[each_kmer,index]
    sorted_kmer_index_dict=sorted(kmer_index_dict.items(),key=lambda x:x[0])
    for each_member in sorted_kmer_index_dict:
        string_line+=each_member[1][0]
        string_line+='('
        index=each_member[1][1]
        for m in range(len(index)):
            if m==len(index)-1:
                string_line+=str(index[m])
            else:
                string_line+=str(index[m])+','
        string_line+=')'+','
    string_line+='\n'
    return string_line

def get_cluster_number(file_name,fam_kmer_dict,output_dir,n_mer,piece_number,protein_name,protein_string,important_n_mer_number,beta):


    fw=open(output_dir+file_name,'w')
    for i in piece_number:
        kmer=[]
        for j in range(n_mer):
            text=protein_string[i][j:len(protein_string[i])]
            kmer+=re.findall('.{'+str(n_mer)+'}', text)
        kmer=list(set(kmer))
        selected_fam_dict={}
        sort_fam=[]
        for each_fam in fam_kmer_dict.keys():
            all_cluster_kmer=fam_kmer_dict[each_fam][2]
            kmer_message=fam_kmer_dict[each_fam][1]
            kmer_label=fam_kmer_dict[each_fam][0]
            each_cluster_score=[]

            for j in range(len(all_cluster_kmer)):
                score=0
                number_score=0
                each_cluster_kmer=all_cluster_kmer[j]
                same_kmer=list(set(kmer).intersection(set(each_cluster_kmer.keys())))
                temp_same_kmer=same_kmer[:]
                for k in temp_same_kmer:
                    number_score+=each_cluster_kmer[k]
                score+=len(same_kmer)
                if score:
                    each_cluster_score.append([j,score,number_score,same_kmer])

            if len(each_cluster_score)==0:
                continue
            each_cluster_score=sorted(each_cluster_score,key=(lambda x:x[2]),reverse=True)
            each_cluster_score=sorted(each_cluster_score,key=(lambda x:x[1]),reverse=True)
            if each_cluster_score[0][2]>=beta and each_cluster_score[0][1]>=important_n_mer_number:
                sort_fam.append([each_cluster_score[0][1],each_cluster_score[0][2],each_fam])
                selected_fam_dict[each_fam]=[each_cluster_score[0][1]],\
                                  each_cluster_score[0][2],each_cluster_score[0][3],\
                                  kmer_message[each_cluster_score[0][0]],kmer_label[each_cluster_score[0][0]]
        if len(sort_fam)==0:
            continue
        sort_fam=sorted(sort_fam,key=(lambda x:x[1]),reverse=True)
        sort_fam=sorted(sort_fam,key=(lambda x:x[0]),reverse=True)


        first_fam=sort_fam[0][2]
        fam_key=list(selected_fam_dict.keys())
        fam_key.remove(first_fam)
        for each_fam_key in fam_key:
            if len(set(list(''.join(selected_fam_dict[each_fam_key][2]))))<5:
                del selected_fam_dict[each_fam_key]
        flag=len(first_fam.split('.'))
        fw.write(protein_name[i]+'\t'+write_line(selected_fam_dict[first_fam],protein_string[i]))
        used_fam=[first_fam]
        other_fams=selected_fam_dict[first_fam][3]
        other_fams=other_fams.split('|')
        for each_fam in other_fams:
            temp_fam=each_fam.split(':')
            temp_name=temp_fam[0].split('_')[0]
            if flag>1:
                temp_name=temp_name.split('.')
                temp_name='.'.join(temp_name[0:flag])
            if len(temp_fam)==2 and temp_name not in used_fam:
                used_fam.append(temp_fam[0].split('_')[0])
                if int(temp_fam[1])>1 and temp_name in selected_fam_dict.keys():
                    fw.write(protein_name[i]+'\t'+write_line(selected_fam_dict[temp_name],protein_string[i]))
                    

    fw.close()






def get_validation_results(input_fasta_file,database_dir,output_dir,output_file_name,n_mer,important_n_mer_number,beta,fam_kmer_dict,jobs):
    if database_dir:
        database_dir=database_dir+'/'
    [protein_name,protein_string] = read_file(database_dir,input_fasta_file)

    current_path=os.getcwd()

    all_output_dir_name=output_dir.split('/')
    while '' in all_output_dir_name:
        all_output_dir_name.remove('')
    for i in range(len(all_output_dir_name)):
        if i==0:
            temp_dir_list=os.listdir(current_path)
            if all_output_dir_name[i] not in temp_dir_list:
                os.mkdir(all_output_dir_name[i])
        else:
            temp_dir_name='/'.join(all_output_dir_name[0:i])
            temp_dir_list=os.listdir(temp_dir_name)
            if all_output_dir_name[i] not in temp_dir_list:
                os.mkdir(temp_dir_name+'/'+all_output_dir_name[i])
    temp_output_dir='/'.join(all_output_dir_name)  
    pred_file_name=output_file_name
    if temp_output_dir:
        temp_output_dir=output_dir+'/'

    pool = Pool(processes=jobs)
    piece_number=np.array_split(list(range(len(protein_name))),jobs)

    for i in range(jobs):
        file_name=input_fasta_file+'_thread_'+str(i)+'.txt'
#        get_cluster_number(file_name,fam_kmer_dict,temp_output_dir,n_mer,piece_number[i],protein_name,protein_string,important_n_mer_number,beta)

        pool.apply_async(get_cluster_number, (file_name,fam_kmer_dict,temp_output_dir,n_mer,piece_number[i],protein_name,protein_string,important_n_mer_number,beta,))
    pool.close()
    pool.join()
    fw=open(temp_output_dir+pred_file_name,'w')
    for i in range(jobs):
        f=open(temp_output_dir+input_fasta_file+'_thread_'+str(i)+'.txt')
        text=f.read()
        f.close()
        fw.write(text)
        os.remove(temp_output_dir+input_fasta_file+'_thread_'+str(i)+'.txt')
    fw.close()
        






# def arg_parser(args): 
#     '''
#     Process command-line user arguments and prepare for further directories
#     '''

#     ####################################
#     #Folders and files:
#     parser = argparse.ArgumentParser(description='The input parameters for the identification technology')

#     parser.add_argument('-kmer_db',     default="CAZyme",type=str,        help="Change n_mer directories path for prediction")
#     parser.add_argument('-output',     default="examples/prediction/output/test_pred_cluster_labels.txt",type=str,        help="file name for prediction results saving")
#     parser.add_argument('-input',     default="examples/prediction/input/test.faa",type=str,        help="Define the fasta file name")
#     parser.add_argument('-k_mer',           default=8,         type=int,                 help="Peptide length for prediction")
#     parser.add_argument('-jobs',            default=8,         type=int,                 help='Number of processor for use for prediction')
#     parser.add_argument('-important_k_mer_number',        default=5,      type=int,               help="Minimum number of n_mer for prediction")
#     parser.add_argument('-beta',        default=2,      type=float,               help="Minimum sum of percentage of frequency of n_mer for prediction")


#     args = parser.parse_args(args)

#     return (args)


class eCAMI_config(object):
    def __init__(
        self,
        db_type = 'Cazyme',
        output = 'examples/prediction/output/test_pred_cluster_labels.txt',
        input = 'examples/prediction/input/test.faa',
        k_mer = 8,
        jobs = 8,
        important_k_mer_number = 5,
        beta = 2.0
    ) -> None:
        if db_type == 'CAZyme':
            self.kmer_db = f'{os.path.dirname(__file__)}/CAZyme'
        else:
            # print(__file__)
            self.kmer_db = f'{os.path.dirname(__file__)}/EC'
        # print(self.kmer_db)
        self.output = output
        self.input = input
        self.k_mer = k_mer
        self.jobs = jobs
        self.important_k_mer_number = important_k_mer_number
        self.beta = beta


def eCAMI_main(args):
    starttime = datetime.datetime.now()
    # args = arg_parser(sys.argv[1:])
    
    database_dir=os.path.dirname(args.input)
    input_fasta_file = os.path.basename(args.input)
    
    output_dir=os.path.dirname(args.output)
    output_file_name = os.path.basename(args.output)
    
    n_mer=args.k_mer
    n_mer_dir_path=args.kmer_db


    important_n_mer_number=args.important_k_mer_number
    jobs=args.jobs
    beta=args.beta

    # print(input_fasta_file)
    fam_kmer_dict=get_kmer_dict(n_mer_dir_path)
    get_validation_results(input_fasta_file,database_dir,output_dir,output_file_name,n_mer,important_n_mer_number,beta,fam_kmer_dict,jobs)

    endtime = datetime.datetime.now()
    print('total time:'+str((endtime - starttime).total_seconds())+'s')
