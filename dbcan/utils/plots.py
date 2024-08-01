import time
from subprocess import Popen, call, check_output
import argparse,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import seaborn as sns
from matplotlib import pyplot
from matplotlib.patches import Patch
#plt.style.use('seaborn')
from dbcan.utils.utils import cgc_standard_line
from dbcan.cli.syntenic_plot import syntenic_plot,read_blast_result_cgc,read_UHGG_CGC_stanrdard_out,read_PUL_cgcgff
from dbcan.cli.syntenic_plot import Get_parameters_for_plot,plot_Polygon_homologous,plot_syntenic_block
from dbcan.cli.syntenic_plot import Get_Position as synGet_Position
from dbcan.cli.syntenic_plot import plot_genome_line as synplot_genome_line
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.colors as colors

class CGC_Standard_Out(object):
    def __init__(self,filename):
        hits = open(filename).readlines()[1:]
        self.genes = []
        for line in hits:
            if line.startswith("CGC#"):
                continue
            lines = line.split()
            self.genes.append(cgc_standard_line(lines))
    def __iter__(self):
        return iter(self.genes)

    def CGCID2genes(self):
        cgcdict = {}
        for gene in self:
            cgcdict.setdefault(gene.cgcid,[]).append(gene)
        return cgcdict

class CGC(object):
    def __init__(self,genes):
        self.genes = genes
        self.ID = genes[0].cgcid ### get cgc id
        self.start = min([gene.gene_start for gene in genes])
        self.end = max([gene.gene_end for gene in genes])
        self.gene_num = len(genes)

    def __iter__(self):
        return iter(self.genes)
    def __repr__(self):
        return "\t".join([self.ID,str(self.start),str(self.end),str(self.gene_num)])
    def __len__(self):
        return len(self.genes)
    
    def get_positions(self):
        starts = [] ; ends = [] ; strands = []
        for gene in self:
            starts.append(gene.gene_start)
            ends.append(gene.gene_end)
            strands.append(gene.strand)
        return starts,ends,strands
    
    def get_proteinID(self):
        return [gene.seqid for gene in self]
    def get_cgc_CAZyme(self):
        return [gene.gene_type for gene in self]

class CGC_standard_out_2CGC(object):
    def __init__(self,dbcan):
        self.CGCs = []
        cgcdict = dbcan.CGCID2genes()
        for cgc in cgcdict:
            self.CGCs.append(CGC(cgcdict[cgc]))
    def __iter__(self):
        return iter(self.CGCs)
    
    def cgcid2CGC(self):
        return {cgc.ID:cgc for cgc in self}

def CGC_plot(args):
    paras = plot_parameters(args)
    dbCAN_standard_out  = CGC_Standard_Out(paras.PUL_annotation)
    cgcs = CGC_standard_out_2CGC(dbCAN_standard_out)
    cgcid2cgc = cgcs.cgcid2CGC()
    cgc = cgcid2cgc[args.cgcid]
    genetypes = cgc.get_proteinID()
    ## get gene starts
    starts,ends,strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    print(f"{args.cgcid.split('|')[0]}:{min(starts)}-{max(ends)}")
    cgc_fig_plot(starts,ends,strands,types,genetypes)

def read_location_reads_count(filename):
    xs2ys = {}
    for line in open(filename):
        lines = line.split()
        xs2ys[int(lines[1])] = int(lines[2])
    return xs2ys

def CGC_plot_reads_count(args):
    paras = plot_parameters(args)
    dbCAN_standard_out  = CGC_Standard_Out(paras.PUL_annotation)
    cgcs = CGC_standard_out_2CGC(dbCAN_standard_out)
    cgcid2cgc = cgcs.cgcid2CGC()
    cgc = cgcid2cgc[args.cgcid]
    genetypes = cgc.get_proteinID()
    ## get gene starts
    starts,ends,strands = cgc.get_positions()
    types = cgc.get_cgc_CAZyme()
    cgc_fig_plot_abund(starts,ends,strands,types,genetypes,paras)

def Get_Position(starts,ends,strands,labels,yshift=0):
    Width = 1000 ; Height = 160; 
    poly_heigth = 10
    Triangle_length = 4
    plot_start_x, plot_start_y = [0,Height/2 - poly_heigth-yshift]
    shfit_pos = starts[0]
    maxbp = max(ends) - min(starts)
    pixeachbp =  Width / maxbp  
    for i in range(len(starts)):
        starts[i] = starts[i] - shfit_pos
        ends[i]   = ends[i] - shfit_pos
    #maxbp = max(ends) - min(starts)  
    ###  5     4
    ###            3         
    ###  1     2
    lines = [] ;polygens = [];texts = []
    for i in range(len(starts)):
        if strands[i] == "+":
            positions_str = str( starts[i] * pixeachbp) + " " + str(plot_start_y) + " " ## first point x,y
            positions_str += str( ends[i] * pixeachbp - Triangle_length) + " " + str(plot_start_y) + " "## second point
            positions_str += str( ends[i] * pixeachbp) + " " + str(plot_start_y + poly_heigth) + " " ## 3
            positions_str += str( ends[i] * pixeachbp - Triangle_length) + " " + str( plot_start_y + 2*poly_heigth) + " " ### 4
            positions_str += str( starts[i] * pixeachbp )+ " " + str(plot_start_y + 2*poly_heigth)
        if strands[i] == "-":
            positions_str = str( starts[i] * pixeachbp ) + " " + str(plot_start_y + poly_heigth) + " "
            positions_str += str( starts[i] * pixeachbp + Triangle_length) + " " + str(plot_start_y) + " "
            positions_str += str(ends[i] * pixeachbp) + " " + str(plot_start_y) + " "
            positions_str += str( ends[i] *pixeachbp ) + " " + str(plot_start_y + 2* poly_heigth) + " "
            positions_str += str( starts[i]* pixeachbp +Triangle_length) + " " + str(plot_start_y + 2* poly_heigth)
        polygens.append(positions_str)
        ### for genome line
        if i < len(starts) -1:
            positions_str = str( ends[i] *pixeachbp) + " " + str(plot_start_y + poly_heigth)  + " "
            positions_str += str( starts[i+1]*pixeachbp) + " " + str(plot_start_y + poly_heigth)
            lines.append(positions_str)
        texts.append(labels[i].split('.')[0])
    
    scale_number = 10
    each_scale_bp = maxbp / scale_number
    each_scale_pix = each_scale_bp * pixeachbp
    
    plot_start_y -= 50
    scale_positions = []; scale_positions_texts = [] ; scale_text = []
    scale_positions.append("0 " + str(plot_start_y + 3*poly_heigth) + " " + str(10*each_scale_pix) + " " + str(plot_start_y + 3*poly_heigth))
    plot_start_y -= 1
    for i in range(scale_number+1):
        positions_str = str(i*each_scale_pix) + " "
        positions_str += str(plot_start_y + 3* poly_heigth) + " "
        positions_str += str(i*each_scale_pix) + " "
        positions_str += str(plot_start_y + 3*poly_heigth + 0.6* poly_heigth)
        scale_positions.append(positions_str)
        positions_str = str(i*each_scale_pix) + " " + str(plot_start_y + 3*poly_heigth + 0.6* poly_heigth)
        scale_positions_texts.append(positions_str)
        scale_text.append(str(int(each_scale_bp*i)+ shfit_pos))

    return polygens,lines,texts,scale_positions,scale_text

def plot_Polygon(polygens1,types1,ax):
    colors_map = {"CAZyme":"#FF0000","null":"#808080","other":"#808080",
    "TC":"#9400D3","CDS":"#00FFFF","STP":"#0000FF","TF":"#1E90FF"}
    for j in range(len(polygens1)):
        polygen = polygens1[j].split()
        points = []
        color  = colors_map[types1[j]]
        for i in range(int(len(polygen)/2)):
            points.append([float(polygen[2*i]),float(polygen[2*i+1])])
        ax.add_patch(Polygon(points, color=color, alpha=0.5,edgecolor=None,facecolor=None,lw=0))

def plot_genome_line(lines,ax):
    for line in lines:
        x1,y1,x2,y2 = points2(line)
        ax.add_patch(Polygon([(x1,y1),(x2,y2)], color="gray",lw=1,edgecolor=None))

def plot_scale_line(lines,label,ax):
    for i,line in enumerate(lines):
        x1,y1,x2,y2 = points2(line)
        ax.add_patch(Polygon([(x1,y1),(x2,y2)], color="gray",lw=1,edgecolor=None))
        if i>=1:
            ax.text(float(x1),float(y1)-20,label[i-1],va='bottom', ha='center')

def points2(coord):
    x1,y1,x2,y2 = coord.split()
    return x1,y1,x2,y2

def cgc_fig_plot(starts,ends,strands,types,labels):
    custom_lines = [Line2D([0], [0], color="red", lw=4,alpha=0.5),
        Line2D([0], [0], color="blue", lw=4,alpha=0.5),
        Line2D([0], [0], color="green", lw=4,alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4,alpha=0.5),
        Line2D([0], [0], color="gray", lw=4,alpha=0.5)]

    labelcolor=["red","blue","green","cyan","gray"]

    genecustom_lines = [Patch(color="#FF0000",alpha=0.5),
        Patch(color="#808080", alpha=0.5),
        Patch(color="#9400D3", alpha=0.5),
        Patch(color="#0000FF", alpha=0.5),
        Patch(color="#1E90FF", alpha=0.5)]
    genelabelcolor=["#FF0000","#808080","#9400D3","#0000FF","#1E90FF"]
    geneslabels    = ["CAZyme","Other","TC","STP","TF"]
    ### for legends
    px = 1/plt.rcParams['figure.dpi'] ## px
    Width = 1400 ; Height = 100
    fig = plt.figure(figsize=(Width*px*1.2,Height*px*2))
    ax  = fig.add_subplot(111)
    maxbp = max(ends) - min(starts)
    polygens,lines,texts,scale_positions,scale_text = Get_Position(starts,ends,strands,labels)
    #print (texts,scale_positions,scale_text)
    plot_Polygon(polygens,types,ax)
    plot_genome_line(lines,ax)
    plot_scale_line(scale_positions,scale_text,ax)
    ax.plot()
    legend = pyplot.legend(genecustom_lines,geneslabels,frameon=False,labelcolor=genelabelcolor,loc='best',title_fontsize="x-large")
    ax.add_artist(legend)
    plt.ylim(0,150)
    plt.xlim(-50,1100)
    plt.tight_layout(pad=0.1)
    plt.axis('off')
    #plt.show()
    file_name = "cgc.pdf"
    print(f"Save figure to file {file_name}!")
    plt.savefig(f"{file_name}")
    plt.close()

def cgc_fig_plot_abund(starts,ends,strands,types,labels,parameters):
    ori_starts = starts.copy(); ori_ends = ends.copy() #### starts will be shift, kept them
    custom_lines = [Line2D([0], [0], color="red", lw=4,alpha=0.5),
        Line2D([0], [0], color="blue", lw=4,alpha=0.5),
        Line2D([0], [0], color="green", lw=4,alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4,alpha=0.5),
        Line2D([0], [0], color="gray", lw=4,alpha=0.5)]

    labelcolor=["red","blue","green","cyan","gray"]

    genecustom_lines = [Patch(color="#FF0000",alpha=0.5,lw=0),
        Patch(color="#808080", alpha=0.5,lw=0),
        Patch(color="#9400D3", alpha=0.5,lw=0),
        Patch(color="#0000FF", alpha=0.5,lw=0),
        Patch(color="#1E90FF", alpha=0.5,lw=0)]
    genelabelcolor=["#FF0000","#808080","#9400D3","#0000FF","#1E90FF"]
    geneslabels    = ["CAZyme","Other","TC","STP","TF"]
    ### for legends
    px = 1/plt.rcParams['figure.dpi'] ## px
    Width = 1400 ; Height = 100
    fig = plt.figure(figsize=(Width*px*1.2,Height*px*4))

    #plt.subplots_adjust(bottom=-0.5)
    ax  = fig.add_subplot(212)
    maxbp = max(ends) - min(starts)
    polygens,lines,texts,scale_positions,scale_text = Get_Position(starts,ends,strands,labels)
    #print (texts,scale_positions,scale_text)
    plot_Polygon(polygens,types,ax)
    plot_genome_line(lines,ax)
    plot_scale_line(scale_positions,scale_text,ax)
    ax.plot()
    legend = pyplot.legend(genecustom_lines,geneslabels,frameon=False,labelcolor=genelabelcolor,loc='best',title_fontsize="x-large")
    ax.add_artist(legend)
    plt.ylim(0,150)
    xlim_x1,xlim_x2 = (-10,1100)
    plt.xlim(xlim_x1,xlim_x2)
    #plt.tight_layout(pad=0.1)
    plt.axis('off')
    
    ### here we need to plot the reads_count of each position
    ### layout 2
    xs2ys = read_location_reads_count(parameters.reads_count)
    max_y = max(xs2ys.values())
    add_readcount_layout(fig,ori_starts,ori_ends,xs2ys,max_y,-3,max_y+10,xlim_x1,xlim_x2,maxbp)
    #plt.show()
    file_name = "cgc-coverage.pdf"
    print(f"Save figure to file {file_name}!")
    plt.savefig(f"{file_name}")
    #plt.close()

def add_readcount_layout(fig,starts,ends,xs2ys,max_y,ylim_y1,ylim_y2,xlim_x1,xlim_x2,syn_maxbp):
    maxbp = max(ends) -min(starts)
    Width = 1000
    pixeachbp =  Width / syn_maxbp ### here the maxbp should from the whole max
    ax  = fig.add_subplot(211)
    #ax = plt.axes([0, 0.4, 1, 0.3])
    plt.ylim(ylim_y1,ylim_y2)
    plt.xlim(xlim_x1,xlim_x2)
    plt.tight_layout(pad=0.1)
    plt.plot((0,1000),(0,0),color='gray',lw=1)
    all_xs = []; all_ys = []
    start = min(starts)
    for i in range(1,maxbp+1): ### 
        all_xs.append(pixeachbp*i)
        if i+start in xs2ys:
            all_ys.append(xs2ys[i+start])
        else:
            all_ys.append(0)
    #print(starts,all_xs[0],all_xs[-1])
    plt.plot(all_xs,all_ys,'-',alpha=0.5,color='red',lw=1)
    #ax.fill_between(all_xs,0*len(all_xs),all_ys,facecolor='red', alpha=0.3,edgecolor="white")
    ax.fill_between(all_xs,all_ys,0,facecolor='red', alpha=0.3,edgecolor="white")
    for pos in ['top', 'right', 'bottom']:
        ax.spines[pos].set_visible(False)
    ax.tick_params(bottom=False, top=False, left=True, right=False)
    ax.set_xticks([])

class plot_parameters():
    def __init__(self,args):
        self.input = args.input if args.input.endswith("/") else args.input +"/"
        #self.R1 = args.R1
        #self.R2 = args.R2
        self.bedtools = args.bedtools
        self.reads_count = args.readscount
        self.output = args.function + "_" + args.output
        self.CAZyme_annotation  = self.input + "overview.txt"
        self.dbCANsub_substrate_annotation  = self.input + "dbcan-sub.hmm.out"
        self.PUL_substrate_annotation  = self.input + "substrate.out"
        self.PUL_annotation  = self.input + "cgc_standard.out"
        self.function = args.function
        self.parameters_check()
    
    def parameters_check(self):
        if self.function == "CGC_plot":
            print("You are plotting the CGC regarding abundance!")
            if not os.path.exists(self.PUL_annotation):
                print(f"PUL annotation file {self.PUL_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
        if self.function == "CGC_coverage_plot":
            print("You are plotting the CGC with reads count!")
            if not os.path.exists(self.PUL_annotation):
                print(f"PUL annotation file {self.PUL_annotation} dose not exit, please check run_dbcan finish status!")
                exit(1)
            if not os.path.exists(self.reads_count):
                print(f"Reads count file {self.reads_count} dose not exit, please run samtools depth first!")
                exit(1)
        
        if self.function == "CGC_synteny_plot":
            self.blastp = self.input + "PUL_blast.out"
        if self.function == "CGC_synteny_coverage_plot":
            self.blastp = self.input + "PUL_blast.out"

def generate_syntenic_block(cgcpul,cgcpul_blastp,genes1,genes2):
    blocks = []
    for record in cgcpul_blastp[cgcpul]: ### generate block information
        query = record.qseqid
        hit   = record.sseqid
        cgc_proteinid = query.split("|")[2]
        pul_proteinid = hit.split(":")[3]
        if not pul_proteinid:
            pul_proteinid = hit.split(":")[2]
        try:
            index1 = genes1.index(cgc_proteinid)
            index2 = genes2.index(pul_proteinid)
            blocks.append(f"{index1}-{index2}-{record.pident}")
        except:
            print (cgcpul,query,hit,cgc_proteinid,pul_proteinid,genes1,genes2)
            continue
    return blocks
        #print (cgc_proteinid2gene[cgc_proteinid],pul_proteinid2gene[pul_proteinid])

def CGC_syntenic_with_PUL(args):
    paras = plot_parameters(args)
    cgcid2pulid = {line.rstrip().split("\t")[0]:line.rstrip().split("\t")[1] for line in open(paras.PUL_substrate_annotation).readlines()[1:]}
    cgc = args.cgcid
    pul = cgcid2pulid.get(cgc,"")
    if pul:
        cgcpul_blastp = read_blast_result_cgc(paras.blastp)
        cgc_proteinid2gene,cgcid2gene,cgcid2geneid = read_UHGG_CGC_stanrdard_out(paras.PUL_annotation)
        PULid_proteinid2gene,PULid2gene,PULid2geneid = read_PUL_cgcgff(args) ### read db_dir
        cgcpul = cgc+":"+pul
        bed_cgc = cgcid2gene[cgc]
        bed_pul = PULid2gene[pul]
        starts1,ends1,strands1,types1 = Get_parameters_for_plot(bed_cgc)
        starts2,ends2,strands2,types2 = Get_parameters_for_plot(bed_pul)
        genes1 = cgcid2geneid[cgc]
        genes2 = PULid2geneid[pul]
        #print (cgc,pul)
        #print (genes1)
        #print (genes2)
        blocks = generate_syntenic_block(cgcpul,cgcpul_blastp,genes1,genes2)
        syntenic_plot(starts1,starts2,ends1,ends2,strands1,strands2,types1,types2,blocks,cgc,pul)
    else:
        print(f"Does not find homolog PUL for CGC: {cgc}!")
        exit(1)

def CGC_syntenic_with_PUL_abund(args):
    paras = plot_parameters(args)
    cgcid2pulid = {line.rstrip().split("\t")[0]:line.rstrip().split("\t")[1] for line in open(paras.PUL_substrate_annotation).readlines()[1:]}
    cgc = args.cgcid
    pul = cgcid2pulid.get(cgc,"")
    if pul:
        cgcpul_blastp = read_blast_result_cgc(paras.blastp)
        cgc_proteinid2gene,cgcid2gene,cgcid2geneid = read_UHGG_CGC_stanrdard_out(paras.PUL_annotation)
        PULid_proteinid2gene,PULid2gene,PULid2geneid = read_PUL_cgcgff(args) ### read db_dir
        cgcpul = cgc+":"+pul
        bed_cgc = cgcid2gene[cgc]
        bed_pul = PULid2gene[pul]
        starts1,ends1,strands1,types1 = Get_parameters_for_plot(bed_cgc)
        starts2,ends2,strands2,types2 = Get_parameters_for_plot(bed_pul)
        genes1 = cgcid2geneid[cgc]
        genes2 = PULid2geneid[pul]
        #print (cgc,pul)
        #print (genes1)
        #print (genes2)
        blocks = generate_syntenic_block(cgcpul,cgcpul_blastp,genes1,genes2)
        syntenic_plot_with_abund(starts1,starts2,ends1,ends2,strands1,strands2,types1,types2,blocks,cgc,pul,paras)
    else:
        print(f"Does not find homolog PUL for CGC: {cgc}!")
        exit(1)

def syntenic_plot_with_abund(starts,starts1,ends,ends1,strands,strands1,Types,Types1,blocks,cgcid,pulid,paras):
    ### for legends
    custom_lines = [Line2D([0], [0], color="red", lw=4,alpha=0.5),
        Line2D([0], [0], color="blue", lw=4,alpha=0.5),
        Line2D([0], [0], color="green", lw=4,alpha=0.5),
        Line2D([0], [0], color="cyan", lw=4,alpha=0.5),
        Line2D([0], [0], color="gray", lw=4,alpha=0.5)]
    #print(starts,starts1,ends,ends1,strands,strands1,Types,Types1,blocks,cgcid,pulid)
    ### syntenic block colors 
    labelcolor=["red","blue","green","cyan"]
    labels    = ["80-100","60-80","40-60","20-40"]
    genecustom_lines = [Patch(color="#FF0000",alpha=0.5),
        Patch(color="#808080", alpha=0.5),
        Patch(color="#9400D3", alpha=0.5),
        Patch(color="#0000FF", alpha=0.5),
        Patch(color="#1E90FF", alpha=0.5)]
    genelabelcolor=["#FF0000","#808080","#9400D3","#0000FF","#1E90FF"]
    geneslabels    = ["CAZyme","Other","TC","STP","TF"]

    ### for legends

    px = 1/plt.rcParams['figure.dpi'] ## px
    Width = 1600 ; Height = 320*2

    fig = plt.figure(figsize=(Width*px,Height*px*2/2.5))
    ax  = fig.add_subplot(212)
    ### decide which 
    maxbp = max([max(ends) - min(starts),max(ends1) - min(starts1)])

    ori_starts = starts.copy(); ori_ends = ends.copy() #### starts will be shift, kept them

    ### get postion for all elements of CGC 
    polygens,blocks_coor,lines_coor,scale_positions,scale_text = synGet_Position(starts,ends,strands,maxbp,yshift=40,up=2) ## CGC
    plot_scale_line(scale_positions,scale_text,ax)
    
    polygens1,blocks1_coor,lines_coor1,_,_ = synGet_Position(starts1,ends1,strands1,maxbp,yshift=0,up=1) ### PUL
    ### 
    plot_Polygon_homologous(polygens,polygens1,Types,Types1,2,ax)
    ###
    plot_syntenic_block(blocks,blocks_coor,blocks1_coor,ax)
    synplot_genome_line(lines_coor,lines_coor1,ax)
    ### need to add the genome postion scale 
    ### legend1
    legend1 = pyplot.legend(custom_lines,labels,frameon=False,labelcolor=labelcolor,
        loc='upper right',title="Identity")
    ax.add_artist(legend1)
    
    legend2 = pyplot.legend(genecustom_lines,geneslabels,frameon=False,
        labelcolor=genelabelcolor,loc='lower right',title="Gene")
    ax.add_artist(legend2)

    plt.text(500,10,cgcid,fontsize=10,horizontalalignment='center')
    plt.text(500,90,pulid,fontsize=10,horizontalalignment='center')
    xlim_x1,xlim_x2 = (-10,1100)
    ylim_y1,ylim_y2 = (0,100)
    plt.ylim(ylim_y1,ylim_y2)
    plt.xlim(xlim_x1,xlim_x2)
    plt.axis('off')
    ax.plot()
    #plt.tight_layout(pad=0.1)
    cgcid = cgcid.replace("|","_") ### need to replace "|" to "_", because | is a special chara for system
    ### for local 
    xs2ys = read_location_reads_count(paras.reads_count)
    max_y = max(xs2ys.values())
    add_readcount_layout(fig,ori_starts,ori_ends,xs2ys,max_y,ylim_y1,max_y,xlim_x1,xlim_x2,maxbp)
    #plt.show();exit()
    print(f"Saving figure to file {cgcid}-syntenic-cov.pdf!")
    plt.savefig(f"{cgcid}-syntenic-cov.pdf")
    plt.close()

def combined_datafram_based_on_first_col(pd_lists,samples):
    if len(pd_lists) <= 1:
        return pd_lists[0]
    else:
        col_name = pd_lists[0].columns
        on_merge_col = col_name[0] ###
        merged_table = pd.merge(pd_lists[0],pd_lists[1],on=[on_merge_col],how="outer")
        
        for i in range(len(pd_lists)):
            ori_names = pd_lists[i].columns
            mod_names = [ori_names[0]] + [ori_names[j] +"_"+ samples[i] for j in range(1,len(ori_names))]
            pd_lists[i].columns = mod_names
        
        for i in range(2,len(pd_lists)):
            merged_table = pd.merge(merged_table,pd_lists[i],on=[on_merge_col],how="outer")
    #print(merged_table.columns)
    abundance_col = col_name[1]
    merged_table.fillna(0,inplace=True)
    merged_table["diff_abs"] = np.abs(merged_table[abundance_col+"_x"] - merged_table[abundance_col+"_y"])
    merged_table["diff"] = merged_table[abundance_col+"_x"] - merged_table[abundance_col+"_y"]
    
    ### rename columns names
    merged_columns = merged_table.columns
    rename_columns = []
    abund_index = 0 
    for column in merged_columns:
        if abundance_col in column:
            rename_columns.append(samples[abund_index])
            abund_index += 1
        else:
            rename_columns.append(column)
    #merged_table.rename(columns={abundance_col+"_x": samples[0], abundance_col+"_y": samples[1]},inplace=True)
    merged_table.columns = rename_columns
    merged_table.sort_values("diff_abs",inplace=True,ascending=False)
    return merged_table,on_merge_col

def filter_out_enzyme_number(table):
    bools = []
    for i in table.iloc[:,0]: ### first col
        if i[0].isdigit():
            bools.append(True)
        elif i in ["PL0","GH0","GT0","CBM0","AA0","CE0"]:
            bools.append(False)
        else:
            bools.append(True)
    table = table[bools]
    #print (table)
    return table

import re
def add_column_type(table): ### Like 
    cols = []
    for i in table["CAZy"]:
        fam = re.sub(r'[0-9]+', '', i)
        cols.append(fam)
        #print (i,fam)
    table["fam"] = cols
    return table

def heatmap_plot(args):
    pds = [filter_out_enzyme_number(pd.read_csv(filename,sep="\t")) for filename in args.input.split(",")]
    samples = args.samples.split(",")
    plt.style.use(args.plot_style)
    if len(pds) != len(samples):
        print("The number of samples is not eaqul the abundance!")
        exit(1)
    for i in range(len(pds)):
        pds[i]["sample"] = samples[i]
    data,x = combined_datafram_based_on_first_col(pds,samples)
    
    if not args.col:
        data = data.iloc[0:int(args.top),:]
    else:
        if args.value:
            data = data.loc[data[args.col].isin(args.value.split(","))]
        else:
            data = data.iloc[0:int(args.top),:]
    data = data.set_index(data.iloc[:,0])
    data = data[samples]
    sns.set_style("whitegrid")
    sns.set_context("paper")
    ### user defined color map 

    ### default color 
    if args.palette:
        cmap = args.palette
    else:
        mycolor=['aliceblue','skyblue','deepskyblue','orange','tomato','red']
        cmap = colors.LinearSegmentedColormap.from_list('my_list', mycolor)
        cmap.set_under('white')
    
    if args.cluster_map:
        sns.clustermap(data, cmap=cmap,cbar=True,vmin=0.1,dendrogram_ratio=0.03,cbar_pos=(0.1, 1, 0.1, 0.1),
                       col_cluster=False,cbar_kws={"shrink": 0.3},
                       figsize=(len(data.columns)*1.2,len(data.index)/3))
        #plt.tight_layout(pad=0.1)
        if args.show_fig:
            plt.show()
        else:
            plt.savefig("heatmap_cluster.pdf")
    else:
        plt.figure(figsize=(len(data.columns)*1.2,len(data.index)/4))
        if args.show_abund:
            ax = sns.heatmap(data, cmap=cmap,yticklabels=True,annot=True,fmt=".0f",linewidths=.5,cbar=True,vmin=0.1,
                         cbar_kws={"shrink": 0.3,"anchor":(0, 0.0)})
        else:
            ax = sns.heatmap(data, cmap=cmap,yticklabels=True,annot=args.show_abund,fmt=".0f",linewidths=0,cbar=True,vmin=0.1,
                         cbar_kws={"shrink": 0.3,"anchor":(0, 0.0)})
        ax.collections[0].colorbar.ax.tick_params(labelsize=6)
    
        plt.xticks(rotation=30)
        plt.tight_layout(pad=0.1)
        if args.show_fig:
            plt.show()
        else:
            plt.savefig("heatmap.pdf")

def bar_plot(args):
    ### --input CAZyme_abund_output,CAZyme_abund_output
    ### --samples fefifo_8002_1,fefifo_8002_7
    ### show top10? or most different?
    pds = [filter_out_enzyme_number(pd.read_csv(filename,sep="\t")) for filename in args.input.split(",")]
    samples = args.samples.split(",")
    plt.style.use(args.plot_style)
    if len(pds) != len(samples):
        print("The number of samples is not eaqul the abundance!")
        exit(1)
    for i in range(len(pds)):
        pds[i]["sample"] = samples[i]
    data,x = combined_datafram_based_on_first_col(pds,samples)
    ### top different abundance --value 'host glycan,starch' --col Substrate
    if not args.col:
        data = data.iloc[0:int(args.top),:]
    else:
        if args.value:
            data = data.loc[data[args.col].isin(args.value.split(","))]
        else:
            data = data.iloc[0:int(args.top),:]
    ### normal
    if args.vertical_bar:
        ax = data.plot.barh(x=x, y=samples)
        plt.ylabel("")
        plt.xlabel("Abundance")
    else:
        ax = data.plot(x=x, y=samples, kind="bar")
        plt.xticks(rotation=90)
        plt.xlabel("")
        plt.ylabel("Abundance")

    #plt.gca().invert_yaxis()
    #plt.gca().invert_xaxis()
    ### switch

    #leg = plt.legend(frameon=False,handletextpad=-2.0, handlelength=0)
    #for item in leg.legendHandles:
    #    item.set_visible(False)
    #axins = inset_axes(ax,  "30%", "30%" ,loc="upper center", borderpad=1)
    
    plt.title(f"The most top{args.top} different families")
    #plt.show()

    if not args.pdf.endswith(".pdf"):
        args.pdf += args.pdf+".pdf"

    plt.savefig(f"{args.pdf}")
    print(f"Saving plot to file: {args.pdf}")

def parse_argv():
    usage = '''
    %(prog)s [positional arguments] [options]
    -----------------------------------------
    Plot CGC
    dbcan_plot CGC_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1'
    -----------------------------------------
    Plot CGC with abundance
    dbcan_plot CGC_coverage_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1' --readscount cgc.depth.txt
    -----------------------------------------
    Plot syntenic blocks between CGC and dbCAN-PUL
    dbcan_plot CGC_synteny_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1'
    -----------------------------------------
    Plot syntenic blocks between CGC and dbCAN-PUL with abundance
    dbcan_plot CGC_synteny_coverage_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1' --readscount cgc.depth.txt
    -----------------------------------------
    barplot of abundance for CAZyme,substate across samples
    dbcan_plot bar_plot -i ../fefifo_8002_1.CAZyme_abund,../fefifo_8002_7.CAZyme_abund --samples fefifo_8002_1,fefifo_8002_7 --top 40 --vertical_bar
    -----------------------------------------
    Heatmap of abundance for CAZyme,substate across samples
    dbcan_plot heatmap_plot -i ../fefifo_8002_1.CAZyme_abund,../fefifo_8002_7.CAZyme_abund --samples fefifo_8002_1,fefifo_8002_7 --show_abund
    -----------------------------------------
    '''
    parser = argparse.ArgumentParser(description='dbCAN plot ulities.',usage=usage,prog='dbcan_plot')
    parser.add_argument('function', help='What function will be used to analyze?')
    parser.add_argument('-i','--input',help='dbCAN CAZyme annotation oupput folder.',default="output",required=True)
    parser.add_argument('-bt','--bedtools',help='bedtools gene reads count results.')
    parser.add_argument('-o','--output',help='output files',default="output")
    parser.add_argument('--db_dir', default="db", help='Database directory')
    parser.add_argument('--cgcid', help='CGC id, consists of contig_ID|cgc_order',type=str,default=None)
    parser.add_argument('--readscount', help='Read counts file generated by samtools depth!')
    parser.add_argument('--samples', help='Samples seperated by ",".')
    parser.add_argument('--top', help='Plot top num, the family or substrate was sorted by abs different.',default=20,type=int)
    parser.add_argument('--plot_style', help='Style for barplot and heatmap.',choices=["ggplot",'seaborn','seaborn-poster'],default="ggplot")
    parser.add_argument('--vertical_bar', help='Vertical bar',action='store_true')
    parser.add_argument('--show_fig', help='Show plot figures or not',action='store_true')
    parser.add_argument('--show_abund', help='Show abundance in heatmap?',action='store_true')
    parser.add_argument('--palette', help='palettes or colormaps defined in matplotlib',default=None)
    parser.add_argument('--cluster_map', action='store_true')
    parser.add_argument('--col',default=None)
    parser.add_argument('--value',default=None)
    parser.add_argument('--pdf',default="bar_plot.pdf",help="bar plot output pdf file name")
    args = parser.parse_args()
    return args


def main():
    args = parse_argv()
    if args.function == "CGC_plot":
        ### dbcan_plot CGC_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1'
        CGC_plot(args)
    if args.function == "CGC_coverage_plot":
        ### dbcan_plot CGC_coverage_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1' --readscount cgc.depth.txt
        CGC_plot_reads_count(args)
    if args.function == "CGC_synteny_plot":
        ### dbcan_plot CGC_synteny_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1'
        CGC_syntenic_with_PUL(args)
    if args.function == "CGC_synteny_coverage_plot":
        ### dbcan_plot CGC_synteny_coverage_plot -i ../fefifo_8002_1.dbCAN --cgcid 'k141_145331|CGC1' --readscount cgc.depth.txt
        CGC_syntenic_with_PUL_abund(args)
    if args.function == "bar_plot":
        # dbcan_plot bar_plot -i fefifo_8002_1.CAZymeSub_abund,fefifo_8002_7.CAZymeSub_abund --samples fefifo_8002_1,fefifo_8002_7 --vertical_bar --col Substrate --value 'host glycan,starch,pectin,xylan,glycogen,cellulose,sucrose,mucin,arabinoxylan'
        bar_plot(args)
    if args.function == "heatmap_plot":
        heatmap_plot(args)
if __name__ =="__main__": ### for test
    #split_uniInput_dbcansub("uniInput",32,"a","a",10,0.5)
    pass
