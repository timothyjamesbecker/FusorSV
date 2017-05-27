#merged genotype IGV visualization sockets script

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
    #type dependant flanking regions settings
    def parse_chrom_pos_range(self,vcf_row,flanking={'DEL':0.5,'DUP':0.7,'INV':0.25}):
        chrom,start,end,svtype,c_id = '','','','',''
        try:
            c_id    = vcf_row[2]
            chrom   = vcf_row[0]
            start   = int(vcf_row[1])
            end     = int(vcf_row[7].rsplit('END=')[-1].rsplit(';')[0])
            svlen   = abs(end-start)
            svtype  = vcf_row[4].replace('<','').replace('>','')
            start   = max(0,start-int(svlen*flanking[svtype]))
            end     = end+int(svlen*flanking[svtype])
            snames  = vcf_row[7].rsplit('TARGET=')[-1].split(',')
            snames  = [sname.replace('0_','') for sname in snames]
        except Exception:
            pass
        return chrom,start,end,svtype,c_id,snames
    
    #read a fusorSV VCF file
    def read_vcf(self,vcf_file):
        snames,header,data = [],[],[]
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header += [line]
                else:
                    data += [line.split('\t')]
            snames = header[-1].split('\t')[9:] #sname should be here
        return snames,header,data
        
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
    snames = ['NA19017','NA12878','HG00419','NA19238','NA19239','NA19625','NA18525']
    #VCF files---------------------------------------------------------------------
    p3 = '/home/tbecker/data/validation_analysis_vcfs/g1k.P3.vcf.gz' #clusters now
    fs = '/home/tbecker/data/validation_analysis_vcfs/all_samples_merged.vcf.forQihui.vcf'         
    #bam files to look at in relation to VCF----------------------------------------
    tars = glob.glob('/home/tbecker/data/tigra_Y/*/tigra.sorted.bam')
    bams = glob.glob('/home/tbecker/data/g1k_high_low/*high_cov*.bam')
    out_dir = '/home/tbecker/data/validation_analysis_vcfs/igv/' 
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    vcf_samples,header,data = igv.read_vcf(fs)
    for row in sorted(data,key=lambda x: (x[0],x[1])):
        igv.new()
        igv.load(p3)
        igv.load(fs)
        chrom,start,end,svtype,c_id,row_samples = igv.parse_chrom_pos_range(row)
        print('snames=%s\tchrom=%s\tstart=%s\tend=%s\ttype=%s\tf_id=%s'%(row_samples,chrom,start,end,svtype,c_id))
        one,test,control = True,[],'/'
        for sname in row_samples:  #look at the snames found in the vcf
            for bam in bams:
                if bam.find(sname)>=0:
                    test += [bam]
                elif one and not bam.find('NA19240')>=0:     #pick a control to display against
                    control = bam
                    one = False
        print('loading control track %s'%control)
        igv.load(control)
        for t in test:
            print('loading test track %s'%t)
            igv.load(t)
        igv.collapse()
        #jump to the svtype dependant coordinate now and plot
        pos = str(chrom)+':'+str(start)+'-'+str(end)
        igv.goto(pos)
        igv.snapshot(out_dir+'/'+c_id+'_'+svtype+'_'+'_'.join(row_samples)+'.jpg')
        
