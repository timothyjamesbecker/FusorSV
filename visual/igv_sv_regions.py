#start IGV server (that has socket serving on local host)
#and then send requests for auto generating SV region pics

import glob
import os
import socket
import sys
import time
import numpy as np
relink = os.path.dirname(os.path.abspath(__file__))+'/../'
sys.path.append(relink) #go up one in the modules
import svu_utils as su

class IGV:
    def __init__(self,ip='127.0.0.1',port=60151): #default IGV port
        self.con = (ip,port)
    
    def __enter__(self):
        return self
    
    def __exit__(self,type,value,traceback):
        return 0
    
    def execute(self,command,buff_size=4096):
        try:
            client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            client.connect(self.con)
            client.send(command+'\n')
            response = client.recv(buff_size)
            print(response)
            client.close()
            return True
        except Exception:
            return False
    
    #reset the IGV session for loading in a new sample
    def new(self):
        return self.execute('new')
        
    #load the needed bam files and vcfs
    def load(self,uri):
        return self.execute('load %s'%uri)
    
    #scroll through the VCF file and snapshot each area
    def goto(self,pos):
        return self.execute('goto %s'%pos)

    def region(self,chrom,start,end):
        return self.execute('region %s %s %s'%(chrom,start,end))
    
    def snapshot(self,image_path):
        return self.execute('snapshot %s'%image_path)
    
    def maxPanelHeight(self,height):
        return self.execute('maxPanelHeight %s'%height)
    
    def collapse(self):
        return self.execute('collapse')
        
    #assume fusorSV VCF4.2 row is a list:
    #[0]CHROM [1]POS [2]ID [3]REF [4]ALT [5]QUAL [6]FILTER [7]INFO [8]FORMAT
    def parse_chrom_pos_range(self,vcf_row,tigra_key=38,hs=3,
                              he={'DEL':0.3,'DUP':0.1,'INV':0.2},flanking=0.5):
        chrom,start,end,svtype,f_id = '','','','',''
        try:
            chrom  = vcf_row[0]
            start  = int(vcf_row[1])
            end    = int(vcf_row[7].rsplit('END=')[-1].rsplit(';')[0])
            svlen  = abs(end-start)
            start  = max(0,start-int(svlen*flanking))
            end    = end+int(svlen*flanking)
            svtype = vcf_row[4].replace('<','').replace('>','')
            f_id   = vcf_row[2].replace('fusorSV_','')
            target = 1
            try:              target = int(vcf_row[7].split('TARGET=')[-1].split(';')[0])
            except Exception: pass
            if target == 0: #novel to target
                idx = su.info_to_idx(row[7])
                tigra   = idx.has_key(tigra_key)
                #check for high confidence
                exp     = float(row[7].split('SVEX=')[-1].split(';')[0].split(',')[0]) >= he[svtype]
                #check for more than sc supporting callers
                support = len(set(idx.keys()).difference(set([tigra_key]))) >= hs
                s = ''
                if support: s += '_HS'
                if exp:     s += '_HE'
                if tigra:   s += '_TA'
                f_id += s
        except Exception:
            pass
        return chrom,start,end,svtype,f_id
    
    #read a fusorSV VCF file
    def read_vcf(self,vcf_file):
        sname,header,data = '',[],[]
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header += [line]
                else:
                    data += [line.split('\t')]
            sname = header[-1].split('\t')[-1] #sname should be here
        return sname,header,data
        
if __name__ == '__main__':   
    igv = IGV(ip='127.0.0.1',port=60151)
#    igv = IGV(ip='127.0.0.1',port=60152)
#    snames = ['HG00096', 'HG00268', 'HG00419', 'HG00759', 
#              'HG01051', 'HG01112', 'HG01500', 'HG01565', 
#              'HG01583', 'HG01595', 'HG01879', 'HG02568', 
#              'HG02922', 'HG03006', 'HG03052', 'HG03642', 
#              'HG03742', 'NA12878', 'NA18525', 'NA18939', 
#              'NA19017', 'NA19238', 'NA19239', 'NA19625', 
#              'NA19648', 'NA20502', 'NA20845']
    snames = ['NA19017','NA12878','HG00419']
#    snames = ['NA19238','NA19239','NA19625','NA18525']
    #VCF files---------------------------------------------------------------------
    p3v = glob.glob('/home/tbecker/data/meta_caller_R4/*/*_S0.vcf')
    fsv = glob.glob('/home/tbecker/data/meta_caller_final/vcf/*_S-1.vcf')
    ssv = glob.glob('/home/tbecker/data/meta_caller_final/vcf/*.support.vcf')
             
    #bam files to look at in relation to VCF----------------------------------------
    tars = glob.glob('/home/tbecker/data/tigra_Y/*/tigra.sorted.bam')
    bams = glob.glob('/home/tbecker/data/g1k_high_low/*high_cov*.bam')
    out_dir = '/home/tbecker/data/meta_caller_all_vcf_tracks_igv/'
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    for sname in snames:
        print("trying to process sample %s"%sname)
        if not os.path.exists(out_dir+'/'+sname): os.makedirs(out_dir+'/'+sname)
        igv.new()
        #load the bam file if available
        for bam in bams:
            if bam.find(sname)>=0:
                print('loading track %s'%bam)
                igv.load(bam)
        for tar in tars:
            if tar.find(sname)>=0:
                print('loading track %s'%tar)
                igv.load(tar)
        #vcf files now-----------------
        ps,exclude = [],1
        for p3 in p3v:
            if p3.find(sname)>=0:
                print('loading track %s'%p3)
                igv.load(p3)
        for fs in fsv:
            if fs.find(sname)>=0:
                print('loading track %s'%fs)
                igv.load(fs)
        for ss in ssv:
            if ss.find(sname)>=0:
                print('loading track %s'%ss)
                igv.load(ss)
        igv.collapse()
        print('data tracks are loaded for %s'%sname)
        #for each row in HE,HS,TA view the region and save as a file: sname_fid_type_HE|HS|TA.jpg
        calls = {} #get the unique rows for all agregate filtered files
        for fs in fsv:
            if fs.find(sname)>=0:
                sample,header,data = igv.read_vcf(fs)
                for row in data:
                    calls[tuple(row)] = 1       
        #now works on multiple filtered vcf file inputs
        for row in [list(k) for k in sorted(calls.keys())]:
            chrom,start,end,svtype,f_id = igv.parse_chrom_pos_range(row)
            print('sample=%s\tchrom=%s\tstart=%s\tend=%s\ttype=%s\tf_id=%s'%(sname,chrom,start,end,svtype,f_id))
            if not os.path.exists(out_dir+'/'+sname+'/'+sname+'_'+f_id+'_'+svtype+'.jpg'):
                if f_id.find('_HE')>=0 or f_id.find('HS')>=0 or f_id.find('TA')>=0:
                    pos = str(chrom)+':'+str(start)+'-'+str(end)
                    igv.goto(pos)
                    igv.snapshot(out_dir+'/'+sname+'/'+sname+'_'+f_id+'_'+svtype+'.jpg')
        
