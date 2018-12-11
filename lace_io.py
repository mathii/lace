#############################################################################
#
#   Copyright 2018 Iain Mathieson
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
#############################################################################

# Input/output functions for the nearest_neighbour script.

from __future__ import division
import sys, getopt, gzip
from math import exp, log, fsum
import numpy as np

##########################################################################################################

def parse_individual(arg):
    """
    Try to open arg as a file and return list of each line, otherwise assume it is a comma separated 
    list and use split to return an actual list
    """

    inds=[]
    try:
        file=open(arg, "r")
        inds=file.readlines()
        inds=[i[:-1] for i in inds] # remove newlines
    except IOError:
        inds=arg.split(",")

    return inds

##########################################################################################################

def load_gen_data(gen_file, sample_file):
    """
    Load genotype data from an IMPUTE format .gen and .smaple file. 
    Will pick the most likely genotype call - even if the probability 
    is very low. Genotype qc should be implented elsewhere.
    """

    gen_data=open(gen_file, "r")
    sample_data=open(sample_file, "r")

    sample_names=[]

    # First two lines are headers
    sample_data.next()
    sample_data.next()

    for line in sample_data:
        sample_names.append(line.split(" ")[0])
        
    snp_names=[]
    snp_pos=[]
    genotype_data=[]

    for line in gen_data:
        data=line[:-1]          # Remove \n from the end of the line
        data=data.rstrip().split(" ")
        snp_names.append(data[1])
        snp_pos.append(int(data[2]))
        gt = data[5:]
        if not len(gt)==3*len(sample_names): 
            print gt
            raise Exception("Bad data line: %d samples and %d entries" % (len(sample_names), len(gt)) )
        gt = [x.index(max(x)) for x in zip(gt[1::3],gt[2::3], gt[3::3])]
        genotype_data.append(gt)

    return {"sample_names":sample_names, "snp_names":snp_names, "snp_pos":snp_pos, "genotype_data":genotype_data}

##########################################################################################################

def load_minimal_data(test_file):
    """
    Load gentotype data from a minimal test format - just a matrix of 1's and 0's
    with postitions as row names and sample names as column names. 
    """
    
    genotype_data=[]
    
    test_data=None
    if(test_file[-3:]==".gz"):
        test_data=gzip.open(test_file, "r")
    else:
        test_data=open(test_file, "r")
    
    snp_pos=[]
    i=0
    for line in test_data:
       if i==0: 
           sample_names=line.split()
       else:
           if not all([x in ["0", "1", "2", "."] for x in line.split()[1:]]):
               raise Exception("Could not read line: " + lines)  
           genotype_data.append([3 if x=="." else int(x) for x in line.split()[1:]])
           snp_pos.append(int(line.split()[0]))
       i+=1

    snp_names=["SNP"+str(x) for x in snp_pos]
    test_data.close()

    return {"sample_names":sample_names, "snp_names":snp_names, "snp_pos":snp_pos, "genotype_data":genotype_data}

##########################################################################################################

def load_eigenstrat_data(file_root):
    """
    Load gentotype data from an [unpacked] eigenstrat file. file_root{.snp,.ind,.geno}
    Missing data coded as 9. 
    """

    ind_file=open(file_root+".ind", "r")
    snp_file=open(file_root+".snp", "r")
    gen_file=open(file_root+".geno", "r")
    
    sample_names=ind_file.readlines()
    sample_names=[x.strip() for x in sample_names]
    sample_names=[x.split()[0] for x in sample_names]
    ind_file.close()
    
    snp_data=snp_file.readlines()
    snp_data=[x.strip() for x in snp_data]
    snp_names=[x.split()[0] for x in snp_data]
    snp_pos=[int(x.split()[3]) for x in snp_data]
    snp_file.close()

    genotype_data=np.genfromtxt(file_root+".geno", dtype=np.int, delimiter=1)
    genotype_data[genotype_data==9]=3
    return {"sample_names":sample_names, "snp_names":snp_names, "snp_pos":snp_pos, "genotype_data":genotype_data}

    
##########################################################################################################

def load_vcf_data(vcf_file):
    """
    Load gentotype data from VCF - We're not parsing the vcf properly, so can't 
    guarantee it's not buggy. Missing data ("./.") coded as 3 
    """
    
    if(vcf_file[-3:]==".gz"):
        vcf_data=gzip.open(vcf_file, "r")
    else:
        vcf_data=open(vcf_file, "r")
        
    snp_names=[]
    snp_pos=[]
    genotype_data=[]

    missing=0
    
    for line in vcf_data:

        if line[0:2] == '##':
            continue
        elif line[0:1] == '#':
            data=line[1:-1]
            data=data.split("\t")
            if data[0:9]==["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                sample_names=data[9:]
            else:
                print data[0:9]
                raise Exception("Bad vcf header line")
        else:
            data=line[:-1]
            data=data.split("\t")

            if len(data[4].split(","))>1: 
                print "Warning: ignoring multi alleleic site at " + data[0]+":"+data[1] 
                continue # multi-allelic sites. 

            if data[2] != ".":
                snp_names.append(data[2])
            else:
                snp_names.append(data[0]+":"+data[1])

            snp_pos.append(int(data[1]))

            if not all([(x[0]=="." and x[2]==".") or (x[0] in ["0", "1"] and x[2] in ["0", "1"]) for x in data[9:]]):
                raise Exception("Could not read line: " + line)  
            
            genotype_data.append([ 3 if x[0]=="." and x[2]=="." else int(x[0])+int(x[2]) for x in data[9:] ])

    return {"sample_names":sample_names, "snp_names":snp_names, "snp_pos":snp_pos, "genotype_data":genotype_data}

##########################################################################################################

def print_traceback_summary(summary, sample_name):
    
    ordering = sorted(summary, key=summary.__getitem__, reverse = True)
    print "Relatedness proportions for " + sample_name
    print
    print "Sample\tDistance"
    for s in ordering:
        print "%s\t%1.3f" % (s,summary[s])

##########################################################################################################

def parents_to_string(parent_tuple):
    """
    Convert parent strings to a slash string
    """
    return str(parent_tuple[0])+" "+str(parent_tuple[1])

##########################################################################################################

def output_phased_data(phasing, sample_names, snp_names, options):
    """
    Output phased data - and quality scores. 
    """
    things_to_output=[]
    things_to_output.append( ("la", "local_ancestry", parents_to_string))
    if options.get("best_parents", None): things_to_output.append( ("bp", "best_parents", parents_to_string) )
   
    # Output phased data
    for suffix, tag, format_func in things_to_output:

        if(options.get("gzip", None)):
            file_name = options["out"]+"."+suffix+".txt.gz"
            out_file = gzip.open(file_name, "w")
        else:
            file_name = options["out"]+"."+suffix+".txt"
            out_file = open(file_name, "w")
        
        #out_file.write( "\t".join(["POS"]+sample_names) + "\n" )
        for i in range(len(phasing[sample_names[0]][tag])):
            #out_file.write( "\t".join([snp_names[i]]+[format_func(phasing[s][tag][i]) for s in sample_names] ) + "\n")
            out_file.write( " ".join([format_func(phasing[s][tag][i]) for s in sample_names] ) + "\n")

        out_file.close()


##########################################################################################################

