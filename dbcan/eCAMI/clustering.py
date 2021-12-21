#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 20:57:51 2019

@author: joy1314bubian
"""


import psutil,argparse,sys

import numpy as np
import linecache
import os
import datetime
from multiprocessing import Pool 
import re
import gc
from collections import Counter


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




def clusterkmercount(all_kmer,selected_label):
    selected_cluster_kmer_number=[]
    for i in range(len(selected_label)):
        each_selected_cluster_kmer_number={}
        for k in selected_label[i]:
            for j in all_kmer[k]:
                if j not in each_selected_cluster_kmer_number.keys():
                    each_selected_cluster_kmer_number[j]=1
                else:
                    each_selected_cluster_kmer_number[j]+=1
        selected_cluster_kmer_number.append(each_selected_cluster_kmer_number)
    return selected_cluster_kmer_number

    

                    
def clusteradd(selected_cluster_kmer_number,all_kmer,new_selected_label,new_unselected_label,beta):
    selected_label=new_selected_label[:]
    unselected_label=new_unselected_label[:]
    all_unselected_score=[]
    all_unselected_number_score=[]
    for i in unselected_label:
        each_cluster_score=[]
        each_cluster_number_score=[]
        for j in range(len(selected_cluster_kmer_number)):
            score=0
            number_score=0
            each_cluster_kmer=selected_cluster_kmer_number[j]
            same_kmer=list(set(all_kmer[i]).intersection(set(each_cluster_kmer.keys())))
            for k in same_kmer:
                each_kmer_percent=each_cluster_kmer[k]/len(selected_label[j])
                if each_kmer_percent>=0.2:
#                    score+=1
                    score+=1
                    number_score+=each_kmer_percent
           
            each_cluster_score.append(score)
            each_cluster_number_score.append(number_score)
        all_unselected_score.append(each_cluster_score)
        all_unselected_number_score.append(each_cluster_number_score)
    max_similary_cluster_label=np.argmax(all_unselected_score,1).tolist()
    remove_label=[]
    for i in range(len(max_similary_cluster_label)):
        if all_unselected_number_score[i][max_similary_cluster_label[i]]>=beta:

            selected_label[max_similary_cluster_label[i]].append(unselected_label[i])
            remove_label.append(unselected_label[i])

    for i in remove_label:
        unselected_label.remove(i)
    
    return [selected_label,unselected_label]


def clusterunion(selected_cluster_kmer_number,all_kmer,newselected_label,alpha):
    selected_label=newselected_label[:]

    between_cluster_score=[]
    kmer_more_than_number=[]
    for i in range(len(selected_cluster_kmer_number)):
        each_cluster_score=[]
        each_kmer_more_than_number=[]

        for j in range(len(selected_cluster_kmer_number)):
            score=0

            more_than_number=0
            if i!=j:
                
                same_key=dict.fromkeys([x for x in selected_cluster_kmer_number[i] if x in selected_cluster_kmer_number[j]]).keys()
                for each_key in same_key:
                    first_selected_cluster_kmer_percent=selected_cluster_kmer_number[i][each_key]/len(selected_label[i])
                    second_selected_cluster_kmer_percent=selected_cluster_kmer_number[j][each_key]/len(selected_label[j])
                    if first_selected_cluster_kmer_percent>=0.2 and second_selected_cluster_kmer_percent>=0.2:
                      
                        score+=1
                        more_than_number+=(selected_cluster_kmer_number[i][each_key]+selected_cluster_kmer_number[j][each_key])/(len(selected_label[i])+len(selected_label[j]))
            each_kmer_more_than_number.append(more_than_number)            
            each_cluster_score.append(score)
        between_cluster_score.append(each_cluster_score)
        kmer_more_than_number.append(each_kmer_more_than_number)
    max_similary_cluster_label=np.argmax(between_cluster_score,1).tolist()
    labels_cluster={}
    k=0
    for i in range(len(max_similary_cluster_label)):
        if kmer_more_than_number[i][max_similary_cluster_label[i]]>=alpha:
            labels_cluster[k]=[]
            labels_cluster[k].append(i)
            labels_cluster[k].append(max_similary_cluster_label[i])
            k=k+1
    
    
    for i in labels_cluster.keys():
        for j in labels_cluster.keys():
            if i!=j and list(set(labels_cluster[j]).intersection(set(labels_cluster[i])))!=[]:
                labels_cluster[i]=list(set(labels_cluster[j]).union(set(labels_cluster[i])))
                labels_cluster[j]=[]

    for i in labels_cluster.keys():
        if labels_cluster[i]!=[]:
            for j in range(len(labels_cluster[i])):
                if j!=0:
                    selected_label[labels_cluster[i][0]]+=selected_label[labels_cluster[i][j]]
                    selected_label[labels_cluster[i][j]]=[]
                    selected_label[labels_cluster[i][0]]=list(set(selected_label[labels_cluster[i][0]]))
    for item in selected_label[:]:
        if item == []:
            selected_label.remove(item)
    
    return selected_label





def correctrate(union_cluster,familyname,proteinname):
    ErrorNumber=0
    allProteinNumber=0
    f=open(familyname+'/'+'result.txt','w')
    subfamfilename=familyname+'/'+familyname+'.subfam.list'
    count=len(open(subfamfilename).readlines())
    pro_label={}
    for i in range(count):
        pro=linecache.getline(subfamfilename,i+1)
        proname=pro.split('\t')[0]
        pro_label['>'+proname]=int(pro.split('\t')[1])
    cluster_number=0
    for i in range(len(union_cluster)):
        t=[]
        if len(union_cluster[i])>=5:
            cluster_number+=1
            for j in union_cluster[i]:
                
                t.append(pro_label[proteinname[j]])
        if t!=[]:
            count = np.bincount(t)
            max_number_label=np.argmax(count)
            ErrorNumber=ErrorNumber+sum(t!=max_number_label)
            f.write('have a fam label number:'+str(len(t))+'\n')
            allProteinNumber=len(t)+allProteinNumber

    f.write('error number:'+str(ErrorNumber)+'\n')
    f.write('have fam label sum number:'+str(allProteinNumber)+'\n')
    f.write('correct rate:'+str(1-ErrorNumber/allProteinNumber)+'\n')
    f.close()
    print('cluster number:'+str(cluster_number))
    print('error number:'+str(ErrorNumber))
    print('have fam label sum number:'+str(allProteinNumber))
    print('correct rate:'+str(1-ErrorNumber/allProteinNumber))





def chooseFeature(all_kmer,sorted_kmer_all_number,kmer_all_number,kmer_id,new_selected_label,new_unselected_label,important_n_mer_number,del_percent,):
#    print ("### RAM usage during clustering: %s" % (ram()))
    selected_label=new_selected_label[:]
    unselected_label=new_unselected_label[:]
    i=0
    delete_feature=[]
    selected_feature=[]
    
    
    inter_cluster_sample_label=[]
    temp_cluster_sample_label=[]
    inter_cluster_label=[]
    selected_cluster_feature=[]

    while sorted_kmer_all_number!=[]:
        sample_label=[]
        for j in range(len(all_kmer)):
            if sorted_kmer_all_number[0] in all_kmer[j]:
                sample_label.append(j)
        
        temp_selected_feature=set(all_kmer[sample_label[0]])
        for j in range(1,len(sample_label)):
            temp_selected_feature=temp_selected_feature.intersection(set(all_kmer[sample_label[j]]))
        temp_selected_feature=temp_selected_feature.intersection(set(sorted_kmer_all_number))
        temp_selected_feature=list(temp_selected_feature)
        sorted_kmer_all_number.remove(sorted_kmer_all_number[0])
        if len(temp_selected_feature)>=important_n_mer_number:
            temp = temp_selected_feature[:]
            temp_delete_feature=[]

            for j in temp:
                if len(sample_label)/kmer_all_number[kmer_id[j]][1]<del_percent:
                    temp_delete_feature.append(j)
                    temp_selected_feature.remove(j)

            if len(temp_selected_feature)>=important_n_mer_number:
                flag=True
                
                
                for cluster_number in range(len(selected_label)):
                    if len(list(set(sample_label).difference(set(selected_label[cluster_number]))))<2 \
                    or len(list(set(selected_label[cluster_number]).intersection(set(sample_label))))>0.9*len(sample_label):
                        flag=False
                        break
                    

                

                if flag:
                    
                    for cluster_number in range(len(selected_label)):
                        if list(set(selected_label[cluster_number]).intersection(set(sample_label)))!=[]:
                            inter_cluster_label.append([cluster_number,len(selected_label)])
                            temp_cluster_sample_label.append(list(set(selected_label[cluster_number]).intersection(set(sample_label))))
                            for temp_label in list(set(selected_label[cluster_number]).intersection(set(sample_label))):
                                inter_cluster_sample_label.append([len(inter_cluster_label)-1,temp_label])
#                    print(flag)
                    selected_label.append(sample_label)
                    for j in sample_label:
                        if j in unselected_label:
                            unselected_label.remove(j)

                    delete_feature=list(set(delete_feature).union(set(temp_delete_feature)))

                    selected_feature=list(set(selected_feature).union(set(temp_selected_feature)))
            
                    selected_cluster_feature.append(temp_selected_feature) 
                    

        
        
        
    for i in delete_feature:
        if i in selected_feature:
            delete_feature.remove(i)

    return [delete_feature,selected_label]         




def write_cluster_kmer_file(cluster_kmer_number,union_cluster,familyname,protein_name,output_dir):
    all_cluster_message=[]
    for i in range(len(union_cluster)):
        cluster_message=[]
        for each_cluster_label in union_cluster[i]:
            each_protein_message=protein_name[each_cluster_label]
            each_protein_message=each_protein_message.split('|')
            each_protein_message=each_protein_message[1:len(each_protein_message)]
            cluster_message+=each_protein_message
        cluster_message=Counter(cluster_message)
        cluster_message=sorted(cluster_message.items(),key = lambda x:x[1],reverse = True)
        all_cluster_message.append(cluster_message)
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
    Dir_list=os.listdir(temp_output_dir)
    if 'kmer_for_each_cluster' not in Dir_list:
        os.mkdir(temp_output_dir+'/kmer_for_each_cluster')
    temp_output_dir=temp_output_dir+'/kmer_for_each_cluster'
    del_list=os.listdir(temp_output_dir)
    for del_file in del_list:
        os.remove(temp_output_dir+'/'+del_file)
    for i in range(len(union_cluster)):
        f=open(temp_output_dir+'/'+'cluster_'+str(i)+'.txt','w')
        f.write('cluster number:'+str(i)+'\n')
        f.write(familyname+':'+str(len(union_cluster[i]))+'\t')
        for each_message in all_cluster_message[i]:
            if each_message==all_cluster_message[i][-1]:
                f.write(each_message[0]+':'+str(each_message[1]))
            else:
                f.write(each_message[0]+':'+str(each_message[1])+'|')
        f.write('\n')
        for each_kmer in cluster_kmer_number[i].keys():
            percent=cluster_kmer_number[i][each_kmer]/len(union_cluster[i])
            if percent>=0.2:
                f.write(each_kmer+':'+str(percent)+'\t')
        f.close()




def writeclusterfile(union_cluster,uncluster_label,protein_name,protein_string,output_dir,fasta_extension):
#    filename=database_dir+'/'+familyname+fasta_extension
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
    Dir_list=os.listdir(temp_output_dir)
    if 'fasta_for_each_cluster' not in Dir_list:
        os.mkdir(temp_output_dir+'/fasta_for_each_cluster')
    temp_output_dir=temp_output_dir+'/fasta_for_each_cluster'
    delList = os.listdir(temp_output_dir)
    for delfile in delList:
        if 'cluster' in delfile:
            os.remove(temp_output_dir + os.sep + delfile)
    for i in range(len(union_cluster)):
        f=open(temp_output_dir+'/cluster_'+str(i)+'.txt','w')
        for j in union_cluster[i]:
            f.write(protein_name[j]+'\n')
            f.write(protein_string[j]+'\n')
        f.close()
    if uncluster_label:
        f=open(temp_output_dir+'/unclustered.txt','w')
        for j in uncluster_label:
            f.write(protein_name[j]+'\n')
            f.write(protein_string[j]+'\n')
        f.close()
def getNewkmer(kmer,delete_kmer):
    return list(set(kmer).difference(delete_kmer))







def arg_parser(args): 
    '''
    Process command-line user arguments and prepare for further directories
    '''

    ####################################
    #Folders and files:
    parser = argparse.ArgumentParser(description='The input parameters for the classification technology')
    parser.add_argument('-output_dir',     default="examples/clustering/output/GH5",type=str,        help="Change output directory for both clustering results and n_mer message for prediction")
    parser.add_argument('-input', default='examples/clustering/input/GH5.faa',     type=str,                 help="Define the fasta file name")

   
    ####################################
    #CUPPclustering parameters:
    parser.add_argument('-piece_number',           default=8,         type=int,                 help="Number of spliting sorted n_mer features for making groups")

    parser.add_argument('-k_mer',           default=8,         type=int,                 help="Peptide length for clustering")

    parser.add_argument('-minimum_group_size',default=5,       type=int,               help="Ignore groups of this number of protein members under minimum group size")
    parser.add_argument('-jobs',            default=8,         type=int,                 help='Number of processor for use for clustering')
    parser.add_argument('-alpha',        default=7,      type=float,               help="Minimum frequency of unioning current groups")
    parser.add_argument('-beta',        default=1,      type=float,               help="Minimum frequency of adding members for current groups")
    parser.add_argument('-important_k_mer_number',        default=10,      type=int,               help="Minimum number of n_mer for making groups in the first step")
    parser.add_argument('-del_percent',        default=0.2,      type=float,               help="Minimum percentage of deleting unimportant n_mer for reducing running time")
    parser.add_argument('-selected_k_mer_number_cut',           default=2,         type=int,                 help="Minimum number of deleting n_mer of making group for reducing running time")
    args = parser.parse_args(args)

    return (args)




def ram():
    '''
    Determine the current RAM usage by the process
    '''
    py = psutil.Process(os.getpid())
    return round(py.memory_info()[0]/2.**30,2)










if __name__ == "__main__":
    args = arg_parser(sys.argv[1:])
    inputs=args.input
    inputs=inputs.split('/')
    database_dir='/'.join(inputs[0:len(inputs)-1])
    temp_name=inputs[-1].split('.')
    allfamilyname=['.'.join(temp_name[0:len(temp_name)-1])]
    fasta_extension='.'+temp_name[-1]

    output_dir=args.output_dir
    minimum_group_size=args.minimum_group_size
    jobs=args.jobs
    n_mer=args.k_mer
    alpha=args.alpha
    beta=args.beta
    important_n_mer_number=args.important_k_mer_number
    selected_kmer_number_cut=args.selected_k_mer_number_cut
    del_percent=args.del_percent

    pool_number=args.piece_number

    for familyname in allfamilyname:
        
        starttime = datetime.datetime.now()
        print(familyname)
        pool_new = Pool(processes=jobs)
        if database_dir:
            database_dir=database_dir+'/'
#        file_object = open(database_dir + familyname + fasta_extension, 'r')
        [protein_name,protein_string] = read_file(database_dir,familyname+fasta_extension)
        

        

        all_kmer = []
        kmer_all_number=[]
        kmer_id={}

        
        print('protein number:'+str(len(protein_name)))

        for i in range(len(protein_name)):
            kmer=[]
            for j in range(n_mer):
                text=protein_string[i][j:len(protein_string[i])]
                kmer+=re.findall('.{'+str(n_mer)+'}', text)
            temp_kmer=kmer[:]
            kmer=list(set(kmer))
            kmer.sort(key=temp_kmer.index)
            all_kmer.append(kmer)
            kmerString=','.join(kmer)

            for each_protein_kmer in kmer:
                if each_protein_kmer not in kmer_id.keys():
                    kmer_id[each_protein_kmer]=len(kmer_all_number)
                    kmer_all_number.append([each_protein_kmer,1])
                else:
                    kmer_all_number[kmer_id[each_protein_kmer]][1]+=1
#        f.close()
        sorted_kmer_all_number=sorted(kmer_all_number,key=lambda x:x[1], reverse=True)
        v = [x[1] for x in sorted_kmer_all_number] 
        v=np.array(v)
        remain_line=np.where(v<selected_kmer_number_cut)[0]
        if len(remain_line)>0:
            last_remain_line=remain_line[0]
        else:
            last_remain_line=len(sorted_kmer_all_number)
        sorted_kmer_all_number=sorted_kmer_all_number[0:last_remain_line]
        sorted_kmer_all_number=[x[0] for x in sorted_kmer_all_number] 

        
        selected_label=locals()
        unselected_label=locals()
        each_sorted_kmer_all_number=locals()

        if len(sorted_kmer_all_number)>1000:
            each_sorted_kmer_all_number_size=len(sorted_kmer_all_number)//pool_number
            for pool_no in range(pool_number-1):
                selected_label['selected_label_'+str(pool_no)]=[]
                unselected_label['unselected_label_'+str(pool_no)]=list(range(len(all_kmer)))
                each_sorted_kmer_all_number['each_sorted_kmer_all_number_'+str(pool_no)]=sorted_kmer_all_number[each_sorted_kmer_all_number_size*pool_no:each_sorted_kmer_all_number_size*(pool_no+1)]
            each_sorted_kmer_all_number['each_sorted_kmer_all_number_'+str(pool_number-1)]=sorted_kmer_all_number[each_sorted_kmer_all_number_size*(pool_number-1):len(sorted_kmer_all_number)]
            selected_label['selected_label_'+str(pool_number-1)]=[]
            unselected_label['unselected_label_'+str(pool_number-1)]=list(range(len(all_kmer)))
        else:
            each_sorted_kmer_all_number_size=0
            for pool_no in range(pool_number-1):
                selected_label['selected_label_'+str(pool_no)]=[]
                unselected_label['unselected_label_'+str(pool_no)]=list(range(len(all_kmer)))
                each_sorted_kmer_all_number['each_sorted_kmer_all_number_'+str(pool_no)]=[]
            each_sorted_kmer_all_number['each_sorted_kmer_all_number_'+str(pool_number-1)]=sorted_kmer_all_number[:]
            selected_label['selected_label_'+str(pool_number-1)]=[]
            unselected_label['unselected_label_'+str(pool_number-1)]=list(range(len(all_kmer)))   
        
        
        
        
        del sorted_kmer_all_number
        del v
        del last_remain_line
        del remain_line
#        del middletime
        del each_sorted_kmer_all_number_size
        gc.collect()

        results=[]
#        print ("### RAM usage during clustering: %s" % (ram()))
        for pool_no in range(pool_number):
#            print ("### RAM usage during clustering: %s" % (ram()))
            results.append(pool_new.apply_async(chooseFeature,(all_kmer,\
                                              each_sorted_kmer_all_number['each_sorted_kmer_all_number_'+str(pool_no)],\
                                              kmer_all_number,kmer_id,selected_label['selected_label_'+str(pool_no)],\
                                              unselected_label['unselected_label_'+str(pool_no)],important_n_mer_number,del_percent,)))
        pool_new.close()
        pool_new.join()
        all_selected_label=[]
        
        new_kmer=[]
        all_delete_feature=[]
        all_delete_feature=set(all_delete_feature)
        for res in results:
            each_feature_cluster_result=res.get()

            all_delete_feature=all_delete_feature.union(each_feature_cluster_result[0])
            all_selected_label+=each_feature_cluster_result[1]

        new_all_kmer=[]
        for each_kmer in all_kmer:
            new_all_kmer.append(list(set(each_kmer).difference(all_delete_feature)))

     
        del all_kmer
        gc.collect()
     
        unique_selected_label=[]
        ununique_labels=[]
        for each_selected_label in all_selected_label:
            ununique_labels+=list(set(unique_selected_label).intersection(set(each_selected_label)))
            unique_selected_label=list(set(unique_selected_label).union(set(each_selected_label)))
        all_unselected_label=list(set(list(range(len(protein_name)))).difference(set(unique_selected_label)))
        selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)   
        if all_unselected_label:
            [all_selected_label,all_unselected_label]=clusteradd(selected_cluster_kmer_number,new_all_kmer,all_selected_label,all_unselected_label,beta)
            
            
        temp_selected_label=[]
        while temp_selected_label!=all_selected_label and len(all_selected_label)>1:
            temp_selected_label=all_selected_label[:]
            selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)   
            all_selected_label=clusterunion(selected_cluster_kmer_number,new_all_kmer,all_selected_label,alpha)
            unique_selected_label=[]
            ununique_labels=[]
            for each_selected_label in all_selected_label:
                ununique_labels+=list(set(unique_selected_label).intersection(set(each_selected_label)))
                unique_selected_label=list(set(unique_selected_label).union(set(each_selected_label)))
            all_unselected_label=list(set(list(range(len(protein_name)))).difference(set(unique_selected_label)))
            selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)   
            if all_unselected_label:
                [all_selected_label,all_unselected_label]=clusteradd(selected_cluster_kmer_number,new_all_kmer,all_selected_label,all_unselected_label,beta)
#            break
            if len(ununique_labels)<5:
                break
#        print('ununique label length:'+str(len(ununique_labels)))


        
        for i in range(len(all_selected_label)):
            for j in range(len(all_selected_label)):
                if i!=j and all_selected_label[i]!=[] and all_selected_label[j]!=[]:
                    same_label=list(set(all_selected_label[i]).intersection(set(all_selected_label[j])))
                    if same_label!=[]:
                        kmer_key_1=new_all_kmer[all_selected_label[i][0]]
                        for each_label in all_selected_label[i]:
                            kmer_key_1=list(set(kmer_key_1).intersection(set(new_all_kmer[each_label])))
                        kmer_key_2=new_all_kmer[all_selected_label[j][0]]
                        for each_label in all_selected_label[j]:
                            kmer_key_2=list(set(kmer_key_2).intersection(set(new_all_kmer[each_label])))
                        if len(kmer_key_1)>len(kmer_key_2) or \
                        (len(kmer_key_1)==len(kmer_key_2) and len(all_selected_label[i])>=len(all_selected_label[j])):
                            for each_same_label in same_label:
                                all_selected_label[j].remove(each_same_label)
                        else:
                            for each_same_label in same_label:
                                all_selected_label[i].remove(each_same_label)
        

        
        
        proteinname=[]
        for i in range(len(protein_name)):
            proteinname.append(protein_name[i].split('|')[0])
        

        while [] in all_selected_label:
            all_selected_label.remove([])
        temp_selected_label=all_selected_label[:]
        for each_temp_selected_label in temp_selected_label:
            if len(each_temp_selected_label)==1:
                all_selected_label.remove(each_temp_selected_label)
        all_unselected_label=list(range(len(new_all_kmer)))
        for each_temp_selected_label in all_selected_label:
            all_unselected_label=list(set(all_unselected_label).difference(set(each_temp_selected_label)))

        cluster_number=0
        for each_selected_label in all_selected_label:
            cluster_number+=len(each_selected_label)
#        print('all protein cluster number:'+str(cluster_number))

        for round in range(2,minimum_group_size):
            temp_selected_label=[]
            temp_unselected_label=all_unselected_label[:]
            for each_selected_label in all_selected_label:
                if len(each_selected_label)>round:
                    temp_selected_label.append(each_selected_label)
                else:
                    for each_label in each_selected_label:
                        temp_unselected_label.append(each_label)
            all_selected_label=temp_selected_label[:]
            add_unselected_label=list(set(temp_unselected_label).difference(set(all_unselected_label)))
            all_unselected_label=temp_unselected_label[:]
            selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)
            temp_selected_label=all_selected_label[:]
            if add_unselected_label!=[]:
                [all_selected_label,all_unselected_label]=clusteradd(selected_cluster_kmer_number,new_all_kmer,all_selected_label,all_unselected_label,beta)
            
            if all_selected_label!=temp_selected_label:
                selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)
                all_selected_label=clusterunion(selected_cluster_kmer_number,new_all_kmer,all_selected_label,alpha)

        temp_selected_label=[]
        new_temp_selected_label=all_selected_label[:]
        while temp_selected_label!=all_selected_label and len(all_selected_label)>2:
            temp_selected_label=all_selected_label[:]
            selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)
            all_selected_label=clusterunion(selected_cluster_kmer_number,new_all_kmer,all_selected_label,alpha)
        
        if new_temp_selected_label!=all_selected_label and all_unselected_label:
            selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)
            [all_selected_label,all_unselected_label]=clusteradd(selected_cluster_kmer_number,new_all_kmer,all_selected_label,all_unselected_label,beta)
        temp_selected_label=all_selected_label[:]
        for each_selected_label in temp_selected_label:
            if len(each_selected_label)<minimum_group_size:
                all_selected_label.remove(each_selected_label)
                all_unselected_label+=each_selected_label
        all_unselected_label=list(set(all_unselected_label))
        cluster_number=0
        for each_selected_label in all_selected_label:
            cluster_number+=len(each_selected_label)
        print('all protein cluster number:'+str(cluster_number))
        writeclusterfile(all_selected_label,all_unselected_label,protein_name,protein_string,output_dir,fasta_extension)
        selected_cluster_kmer_number=clusterkmercount(new_all_kmer,all_selected_label)
        write_cluster_kmer_file(selected_cluster_kmer_number,all_selected_label,familyname,protein_name,output_dir)

##        
        endtime = datetime.datetime.now()

        print('total time:'+str((endtime - starttime).total_seconds())+'s')
#        get_results(familyname)
        print('#####################################################')
