import os
import json
import shutil
import numpy as np
import pysam
import fusion_utils as fu 

#takes in the multi-chrom fasta file and rads it by seq/chrom
#reads each sequence and then writes that portion to the chrom_fasta_dir
#using the chrom_base='' seq_name.fa
#write_fasta_by_chrom(path,path[0:path.rfind('/')],'')

def get_reverse_complement(seq):
    m = {'A':'T','a':'t','T':'A','t':'a','C':'G','c':'g','G':'C','g':'c','N':'N'}
    return ''.join([m[k] for k in seq[::-1]])

def read_fasta_substring(fasta_path,chrom,pos,end):
    ss = ''
    with pysam.FastaFile(fasta_path) as f:
        ss = f.fetch(chrom)
    return ss[pos:end]
    
def read_fasta_chrom(fasta_path,chrom):
    ss = ''
    with pysam.FastaFile(fasta_path) as f:
        ss = f.fetch(chrom)
    return ss

def read_fasta(fasta_path):
    ss = {}
    with pysam.FastaFile(fasta_path) as f:
        names = f.references
        for chrom in names:
            ss[chrom] = f.fetch(chrom)
    return ss

def get_fasta_seq_names(fasta_path):
    ss = []
    with pysam.FastaFile(fasta_path) as f:
        ss = list(f.references)
    return ss

def get_fasta_seq_lens(fasta_path):
    ss = []
    with pysam.FastaFile(fasta_path) as f:
        ss = list(f.lengths)
    return ss

def get_fasta_seq_names_lens(fasta_path):
    ss = {}
    with pysam.FastaFile(fasta_path) as f:
        names = f.references
        lens  = f.lengths
    for i in range(len(names)):
        ss[names[i]]=lens[i]
    return ss   
    
def write_fasta(seqs,fasta_path,index=True):
    with open(fasta_path,'w') as fasta:
        for k in seqs:
            fasta.write('\n'.join(['>%s'%k]+[seqs[k][i:(i+80)] for i in range(0,len(seqs[k]),80)]+['\n']))
    if index: pysam.faidx(fasta_path) #reindex
    return True    
    
#ss is a HTSeq Sequence list
def write_fasta_by_chrom(seqs, chrom_fasta_dir, chrom_base=''):
    names = []
    for k in seqs:
        name = chrom_fasta_dir+'/'+chrom_base+k+'.fa'
        names += [name]
        with open(name, 'w') as fasta:
            fasta.write('\n'.join(['>%s'%k]+[seqs[k][i:(i+80)] for i in range(0,len(seqs[k]),80)]+['\n']))
    return names

def write_fasta_mask(M,json_path):
    with open(json_path,'w') as f:
        json.dump(M,f)
    return True

#compute an expectation given randomly distributed short reads for the RD windows (hist bins)
def expected_window(depth=20,length=100,target=100):
    return int(round((1.0*target/depth)*2.0*length,0))           

def get_local_path(local_path):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    return path+local_path
   
def get_coordinate_offsets(json_path):
    info,O = {},{}
    with open(json_path,'r') as f:
        info = json.load(f)
    for i in info:
        O[str(i)] = info[i]
    return O

def get_stage_map(json_name):
    raw = {}
    with open(json_name,'r') as f:
        raw = json.load(f)
    stage_map = {}
    for r in raw:
        try:
            k = int(r)
            stage_map[k] = str(raw[r]) #byte string in python 2.7.10+
        except ValueError:
            print(r)
            pass
    return stage_map

def write_coordinate_offsets(fasta_path,json_path,trim_chr):
    #read in the fasta lengths and then sort them by length into a json offset map
    L = get_fasta_seq_names_lens(fasta_path)
    l_i = list(np.argsort([L[k] for k in L]))[::-1] #sort by max length
    O,offset = {},0 #starts at zero
    for i in l_i:
        if trim_chr:
            chrom = L.keys()[i]
            chrom_tag = chrom.split('CHROM')[-1]
            chrom_tag = chrom_tag.split('CHR')[-1]
            chrom_tag = chrom_tag.split('chrom')[-1]
            chrom_tag = chrom_tag.split('chr')[-1]
            O[chrom_tag] = offset
        else:
            O[L.keys()[i]] = offset
        offset = offset+L[L.keys()[i]] #old + new + 1
    with open(json_path,'w') as f:
        json.dump(O,f)
    return True

def get_chrom_dict(json_name):
    path = os.path.dirname(os.path.abspath(__file__))+'/data/'
    info,I = {},{}
    with open(path+json_name,'r') as f:
        info = json.load(f)
    for i in info:
        I[str(i)] = int(info[i])
    return I
    
#input is json data store for svmask regions which are the chrom and start,stop pairs
#in the 0-chrom_length form and a offsetmap object O that will update coordnates
#output is a sorted by ref start pos list of list to filter on
def get_mask_regions(json_path,O,complement=False):
    M = {}
    with open(json_path,'r') as f:
        M = json.load(f) #load the svmask per sequence
    N = []    
    for k in M:
        for i in range(len(M[k])):
            N += [[np.uint32(O[k])+np.uint32(M[k][i][0]),np.uint32(O[k])+np.uint32(M[k][i][1])]] #apply offsets
    N = sorted(N,key=lambda x: x[0])
    if complement: N = complement_mask_regions(N,O)
    return N

#get the opposite part of the mask in the offset mask format
def complement_mask_regions(R,O):
    #make a dummy SVU and use LR
    D = [np.uint32(0),[],np.uint32(0),np.uint32(0),{}]
    T,S = [],[[np.uint32(0),np.uint32(sorted([O[k] for k in O])[-1])]+D]
    for i in range(len(R)):
        T += [R[i]+D]
    I,U,D1,D2 = fu.LR(T,S) #D2 is the complement now
    T = []
    for i in range(len(D2)):
        T += [D2[i][0:2]]
    return T

#flatten the input mask regions in 2D space
def flatten_mask_regions(R,O,complement=False):
    C = sorted(R,key=lambda x: x[0])
    M = []
    n = len(C)
    if n > 0:
        i = 0
        while i < n-1:
            j = i+1
            b = C[i][1]
            while b+1 >= C[j][0] and j < n-1:
                if b < C[j][1]: b = C[j][1]
                j += 1
            M += [[C[i][0],b]]
            i = j                              #----------span of i to j here-------------
        if len(M)>0 and M[-1][1]+1>=C[i][0]:   #----------potential span of i-1 to i here-
            if M[-1][1]<C[i][1]:
                M[-1][1] = C[i][1]
        else:                                  #------------only i is here----------------
            M += [[C[i][0],C[i][1]]]
    if complement: M = complement_mask_regions(M,O)
    return M

def write_mask_regions(json_name):
    return True

def copy_json_mask(in_path,out_path):
    shutil.copy2(in_path,out_path)
    return True

def bed_mask_to_json_mask(bed_path,json_path):
    bed_data = []
    with open(bed_path, 'r') as f:
        bed_data = [i.split('\t') for i in f.read().split('\n')]
    data = {}
    for row in bed_data: #chrom,start,stop
        if len(row)>=3:
            if data.has_key(row[0]): data[row[0]] += [[int(row[1]),int(row[2])-1]]
            else:                    data[row[0]]  = [[int(row[1]),int(row[2])-1]]
    for k in data:
        data[k] = sorted(data[k], key = lambda x: x[0])
    #no coordinate sort per contig
    with open(json_path,'w') as f:
        json.dump(data,f)
    return True
    
def get_offsets(chroms,order):
    x,m = 0,{k:0 for k in chroms}
    for k in order:
        m[k] = x
        x += chroms[k]
    return m

def ucsc_clip_chr(chrom):
    if chrom.startswith('chr'):
        return chrom.replace('chr','')
    else:
        return chrom

def write_ucsc_knownGene(tsv_path,json_path,offset_map):
    return False


#[1]
#chrom chromStart chromEnd ... no header
def write_bed_json(tsv_path,json_path,offset_map):
    S,raw,err = {k:[] for k in offset_map},[],[]
    with open(tsv_path, 'r') as f:
        raw = f.readlines()
        for i in range(len(raw)): #skip header line if skip==1
            try:
                row = raw[i].split('\t')
                ctg = ucsc_clip_chr(row[0])
                if S.has_key(ctg):
                    S[ctg] += [[int(row[1]),int(row[2])]]
            except Exception as e:
                err += [e]
                pass
        if len(err)>0: print(err)
    with open(json_path,'w') as f:
        json.dump(S,f)       
        return True
    return False

#[1]
#bin	 chrom chromStart chromEnd forms: SimpleRepeatx,Gap,SegmentalDups,MicroSatellite
def write_ucsc_json(tsv_path,json_path,offset_map):
    S,raw,err = {k:[] for k in offset_map},[],[]
    with open(tsv_path, 'r') as f:
        raw = f.readlines()
        for i in range(1,len(raw)): #skip header line if skip==1
            try:
                row = raw[i].split('\t')
                ctg = ucsc_clip_chr(row[1])
                if S.has_key(ctg):
                    S[ctg] += [[int(row[2]),int(row[3])]]
            except Exception as e:
                err += [e]
                pass
        if len(err)>0: print(err)
    with open(json_path,'w') as f:
        json.dump(S,f)       
        return True
    return False
    

#[2]
#bin swScore milliDiv milliDel milliIns genoName genoStart genoEnd genoLeft strand 
#repName repClass repFamily repStart repEnd repLeft id
def write_ucsc_repeat_masker(tsv_path,json_path,offset_map,
                             repClasses=['SIMPLE_REPEAT','LOW_COMPLEXITY','SATELLITE']):
    S,raw,err = {k:[] for k in offset_map},[],[]
    rep_class_counts = {k:0 for k in repClasses}
    with open(tsv_path, 'r') as f:
        raw = f.readlines()
        for i in range(1,len(raw)): #skip header line
            try:
                row = raw[i].split('\t')
                ctg = ucsc_clip_chr(row[5])
                cls = row[11].upper()
                if rep_class_counts.has_key(cls) and S.has_key(ctg):
                    rep_class_counts[cls] += 1
                    S[ctg] += [[int(row[6]),int(row[7])]]
            except Exception as e:
                err += [e]
                pass
        if len(err)>0: print(err)
        for k in rep_class_counts:
            print('parsed %s rows with repClass value %s'%(rep_class_counts[k],k))
    with open(json_path,'w') as f:
        json.dump(S,f)       
        return True
    return False

#M = read_fasta_mask('/home/tbecker/data/human_g1k_v37_decoy/human_g1k_v37_decoy_S15.fa.svmask.fasta')
#write_fasta_mask(M,'/home/tbecker/software/fusionSVU/data/human_g1k_v37_decoy_svmask.json')
   
   