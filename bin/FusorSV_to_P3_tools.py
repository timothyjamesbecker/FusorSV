#!/usr/bin/env python
#Timothy James Becker, PhD candidate, UCONN 12/09/2017-12/09/2017
#pure python 2.7.10+ implementation of VCF overlap, intersection,difference
#design to query FusorSV merged VCF files with g1k.P3.vcf file
import argparse
import gzip
import time

des = """
FusorSV all_samples.vcf to g1k.P3.vcf intersection, difference, in silico tool
Timothy James Becker, PhD candidate, UCONN 12/09/2017-12/09/2017"""
parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i', '--fusorsv_vcf',type=str,help='all_samples.vcf file path\t[None]')
parser.add_argument('-t', '--g1k_P3_vcf',type=str,help='g1k.P3.vcf.gz file path\t[None]')
parser.add_argument('-o', '--out_dir',type=str,help='vcf output file path\t[None]')
parser.add_argument('-r', '--overlap',type=float,help='reciprocal overlap\t[0.5]')
args = parser.parse_args()

if args.fusorsv_vcf is not None:
    vcf_path = args.fusorsv_vcf
else:
    print('fusorsv_vcf path not given!')
    raise IOError
if args.g1k_P3_vcf is not None:
    p3_vcf_path = args.g1k_P3_vcf
else:
    print('g1.P3.vcf path not given!')
    raise IOError
if args.out_dir is not None:
    out_dir = args.out_dir
else:
    print('out_dir not given!')
    raise IOError
if args.overlap is not None:
    r = args.overlap
else:
    r = 0.5

def get_info_item(row,key,INFO=7):
    return row[INFO].rsplit(key)[-1].rsplit(';')[0]

def get_methods(row):
    method = get_info_item(row,'METHOD=').rsplit('|')
    T = {get_info_item(row,'SVTYPE='):[m.rsplit(':')[0] for m in method]}
    return T

def read_vcf(path):
    header, data, raw, last = [], [], [], ""
    if not path.endswith('.gz'):
        with open(path, 'r') as f:
            raw = f.readlines()
            for line in raw:
                if line.startswith('#'):
                    header += [line.replace('\n', '').rsplit('\t')]
                else:
                    if last != line:
                        data += [line.replace('\n', '').rsplit('\t')]
                    last = line
    else:
        with gzip.GzipFile(path,'rb') as f:
            raw = f.readlines()
            for line in raw:
                if line.startswith('#'):
                    header += [line.replace('\n', '').rsplit('\t')]
                else:
                    if last != line:
                        data += [line.replace('\n', '').rsplit('\t')]
                    last = line
    return header,data

def get_in_silico(data,in_silico_key='38'):
    F = {}
    for row in data:
        m = get_methods(row)
        if F.has_key(m.keys()[0]):
            F[m.keys()[0]] += [m[m.keys()[0]]]
        else:
            F[m.keys()[0]] = [m[m.keys()[0]]]
    # count the in silico for each type
    S = {}
    for t in sorted(F.keys()):
        S[t] = [0,len(F[t])]
        for i in range(len(F[t])):
            if in_silico_key in F[t][i]:
                S[t][0] += 1
        print('in silico support for SVTYPE=%s is %s / %s'%(t,S[t][0],S[t][1]))
    return S

def get_samples(data):
    S = set([])
    for row in data:
        S.add(get_info_item(row,'SAMPLE='))
    return sorted(list(S))

#amount of reciprocal overlap between two ranges:
#c1=[chr1,x1,x2],c2=[chr1,x1,x2]=>1.0
def overlap(c1, c2):
    if c1[0]==c2[0] and c1[2]>=c1[1] and c2[2]>=c2[1]:
        l = (c1[2]-c2[1]+1)+(c2[2]-c1[1]+1)
        u = float(min(l,max(c1[2],c2[2])-min(c1[1],c2[1])+1))
        i = 1.0*float(abs(c1[1]-c2[1])+abs(c1[2]-c2[2]))
        x = max(0.0,u-i)/u
    else:
        x = 0.0
    return x

def intersect(C1,C2,r=0.5):
    ts = sorted(list(set(C1.keys()+C2.keys())))
    I = {t:[] for t in ts}
    for t in ts:
        if C1.has_key(t) and C2.has_key(t):
            for i in range(len(C1[t])):
                for j in range(len(C2[t])):
                    if overlap(C1[t][i],C2[t][j]) >= r:
                        I[t] += [C1[t][i]]
                        break
    return I

def difference(C1,C2,r=0.5):
    ts = sorted(list(set(C1.keys() + C2.keys())))
    D = {t: [] for t in ts}
    for t in ts:
        if C1.has_key(t) and C2.has_key(t):
            for i in range(len(C1[t])):
                unique = True
                for j in range(len(C2[t])):
                    if overlap(C1[t][i], C2[t][j]) >= r:
                        unique = False
                        break
                if unique: D[t] += [C1[t][i]]
    return D

def get_t_ranges(data,unique=True):
    if unique:
        D,R = {},{}
        for row in data:
            t = get_info_item(row,'SVTYPE=')
            chrom,x1,x2 = row[0].replace('chr',''),int(row[1]),int(get_info_item(row,'END='))
            if D.has_key(t):
                if (chrom,x1,x2) not in D[t]:
                    D[t].add((chrom,x1,x2))
            else:
                D[t]  = set([(chrom,x1,x2)])
        for t in D:
            R[t] = sorted([[d[0],d[1],d[2]] for d in D[t]],key=lambda x: (x[0],x[1],x[2]))
    else:
        R = {}
        for row in data:
            t = get_info_item(row,'SVTYPE=')
            chrom,x1,x2 = row[0].replace('chr',''),int(row[1]),int(get_info_item(row,'END='))
            if R.has_key(t):
                R[t] += [[chrom,x1,x2,'\t'.join(row)]]
            else:
                R[t]  = [[chrom,x1,x2,'\t'.join(row)]]
    return R

start = time.time()
print('reading %s'%vcf_path)
h1,d1 = read_vcf(vcf_path)
print('getting ranges for %s'%vcf_path)
C1 = get_t_ranges(d1,unique=False)

print('gathering in silico support statistics for %s'%vcf_path)
#process the calls by type and get each method
S = get_in_silico(d1)
print('finished gathering in silico support statistics for %s'%vcf_path)

print('reading %s'%p3_vcf_path)
h2,d2 = read_vcf(p3_vcf_path)
print('getting ranges for %s'%p3_vcf_path)
C2 = get_t_ranges(d2)

print('gathering intersection and forward difference from %s to %s'%(vcf_path,p3_vcf_path))
print('counts are computed with reciprocal overlap of %s'%r)
I = intersect(C1,C2,r=r)
D = difference(C1,C2,r=r)
print('SVTYPE, intersection, difference, total are as follows:')
for t in sorted(C1.keys()):
    print('%s:I=%s,D=%s,T=%s'%(t,len(I[t]),len(D[t]),len(C1[t])))
stop = time.time()
print('total execution time was %s sec'%(round(stop-start,2)))

