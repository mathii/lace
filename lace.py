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

from __future__ import division
import sys, getopt, ancestry
import lace_io as io
import recombination as rec   
import numpy as np  
import preclustering as pre
from numpy import array
from collections import defaultdict
from multiprocessing import Pool

##########################################################################################################

# algorithm name to module map:
# viterbi - full 2 parent viterbi algorithm, implemented in cython
algo_defs = { "viterbi":"c_viterbi3" }

##########################################################################################################

def help():
    print "Run the lace algorithm on genotype data"
    print "Usage: python lace.py {options}"
    print "Options with asterisks require arguments"
    print
    print "Input/Output:"
    print "-m*   [minimal] input file"
    print "-v*   [vcf] input file"
    print "-e*   [eigenstrat] input root"
    print "-r*   [recombination] map file or constant cm/mb"
    print "-o*   [out]put file root"
    print "-t*   [population] labels - one per line"
    print "-b    output [best_parents]"
    print "-z    output [gzip]ped files"
    print
    print "Options:"
    print "-s    Input data is [pseudo_haploid]"
    print "-w    [window] number of snps to smooth"
    print "-i*   calculate for these [individual]s - comma sep list or file with one per line"
    print "-n*   use a [panel] of only these individuals - as -i option"
    print "-u*   [multi_process]ing: use this many processes"
    print "-x*   Only consider the first [max_snps] snps"
    print "-c*   Select only this many [closest] samples to query for each individual"
    print
    print "Other settings"
    print "--Ne*  Change Ne. Presumably you know what you're doing"
    print "--tbk* Number of steps to check traceback chunks - for viterbi: a memory/speed tradeoff"    
    print "--mtp* Mutation probability - probability of imperfect copying. Default 0.01"
    print "--thw* Triple heterozgote weight - use to downweight the trple het probability. Default 0.01"
    print "--npt* Number of traceback paths to use for ancestry - the more you use, the more you phase"
    print "--smo Smooth output"


##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "Ne": 14000, "out":"pace.out", "algorithm":"viterbi", "traceback_lookback_k":100, "recombination_map":"1", "mutation_probability":0.01, "pseudo_haploid":False, "populations":None,  "triple_het_weight":0.01, "n_traceback_paths":9, "window":1, "smooth_output":False}

    try:
        opts, args = getopt.getopt(sys.argv[1:], "m:v:e:r:o:p:bzsi:n:u:x:c:w:", ["help", "eigenstrat=",  "minimal=", "vcf=", "recombination=", "max_snps=", "out=", "best_parents", "pseudo_haploid", "gzip", "phase", "individual=", "multi_process=", "closest=", "Ne=", "tbk=","mtp=",  "panel=", "populations=", "npt=", "window=", "smo"])
    except Exception as err:
        print str(err)
        help()
        sys.exit()

    if not opts:
        help()
        sys.exit()        
        
    for o, a in opts:
        if o in ["-h","--help"]:
            help()
            sys.exit()
        elif o in ["-m","--minimal"]:        options["test_file"] = a
        elif o in ["-v","--vcf"]:            options["vcf_file"] = a
        elif o in ["-e","--eig"]:            options["eigenstrat_root"] = a
        elif o in ["-r","--recombination"]:  options["recombination_map"] = a
        elif o in ["-w","--window"]:         options["window"] = int(a)
        elif o in ["-o","--out"]:            options["out"] = a      
        elif o in ["-q","--quality"]:        options["quality"] = True      
        elif o in ["-b","--best_parents"]:   options["best_parents"] = True      
        elif o in ["-z","--gzip"]:           options["gzip"] = True      
        elif o in ["-s","--pseudo_haploid"]: options["pseudo_haploid"] = True
        elif o in ["-i","--individual"]:     options["individual"] = io.parse_individual(a)
        elif o in ["-n","--panel"]:          options["panel"] = io.parse_individual(a)
        elif o in ["-p","--populations"]:    options["populations"] = io.parse_individual(a)
        elif o in ["-x","--max_snps"]:       options["max_snps"] = int(a)      
        elif o in ["-u","--multi_process"]:  options["multi_process"] = int(a)      
        elif o in ["-c","--closest"]:        options["closest"] = int(a)      
        elif o in ["--Ne"]:                  options["Ne"] = int(a)      
        elif o in ["--tbk"]:                 options["traceback_lookback_k"] = int(a)      
        elif o in ["--mtp"]:                 options["mutation_probability"] = float(a)      
        elif o in ["--thw"]:                 options["triple_het_weight"] = float(a)      
        elif o in ["--npt"]:                 options["n_traceback_paths"] = int(a)      
        elif o in ["--smo"]:                 options["smooth_output"] = True      

    # Check we entered some sensible data
    validate_options(options)
    print
    print "Running with the following options:"
    for o, a in options.items():
        print o + "\t"*(4- int(1+len(o)/7.5)) + str(a)
    print

    return options

##########################################################################################################

def validate_options(options):
    """
    Check that the options we entered are sensible
    """
    if ("test_file" in options) + ("vcf_file" in options) + ("eigenstrat_root" in options) !=1:
        raise Exception("Must specify exactly one data source")
    if not options.get("recombination_map"):
        raise Exception("Must specify recombination map")
    if not options["populations"]:
        raise Exception("Must specify population labels (-p/--populations) to call local ancestry")
    
##########################################################################################################

def safe_run_for_one_sample(args):
    """
    Safe wrapper for nn_for_one_sample. If something goes wrong, just return an empty summary
    This might give wierd results, but at least it won't just blow up. 
    """
    try:
        out=run_for_one_sample(args)
    except Exception as ex:
        print "Caught an exception in nn_for_one_sample for sample " + args[0]
        print "Safe mode: returning result for " + args[0]
        print "Details: " + str(ex)
        print
        out = {}

    return out

##########################################################################################################

def run_for_one_sample(args):
    """
    Estimate the nearest neighbour in the sample at each snp by using the 
    supplied algorithm. Returns a list of length equal to the number
    of snps with the nearest neighbour in each position. 
    what, should be either nn or phase. Sample indices maps used sample index
    back to actual sample index, which is from a superset. 
    """
    (sample_name, data, recombinator, options, summary_function) = args

    # Load the right (cython/python) module 
    algorithm=__import__(algo_defs[options["algorithm"]])

    i = data["sample_names"].index(sample_name) # This is the index of the sample to be queried


    include=np.ones(len(data["sample_names"]), dtype=np.bool)
    if "panel" in options:
        include=np.in1d(np.array(data["sample_names"]), np.array(options["panel"]))
    # exclude current snp
    include[i]=False
    
    # Data with the current snp excluded
    observations=data["genotype_data"][:,i]
    used_sample_names=[x for x,i in zip(data["sample_names"],include) if i]
    used_genotype_data=data["genotype_data"][:,include]
    
    used_sample_indices=np.where(include)[0]
    
    N_samples = sum(include)
    N_snps = len(data["snp_pos"])
    
    # Subsampling
    if "closest" in options:
        used_genotype_data, used_sample_names = pre.closest_n( used_genotype_data, used_sample_names, observations, options["closest"] )

    used_genotype_frequency=None
    used_options=options.copy()
    used_genotype_data_na=used_genotype_data.astype(float)
    used_genotype_data_na[used_genotype_data_na>2.0]=np.nan
    used_genotype_frequency=np.nanmean(used_genotype_data_na,axis=1)/2
    used_options["used_genotype_frequency"]=used_genotype_frequency

    trans=algorithm.transition( N_samples, options["Ne"], recombinator, data["snp_pos"])
    emiss=algorithm.emission(N_samples, options)
    if options["pseudo_haploid"]:
        emiss=algorithm.pseudohaploid_emission(N_samples, options)
    vit=algorithm.calculator(used_genotype_data, trans, emiss, observations, used_options)

    vit.calculate()
    out = summary_function(vit, used_sample_indices, data["snp_pos"], options, used_genotype_data, observations, i ) 

    if "multi_process" in options: # This is a bit of a hack to get some output from the multiprocess. 
        print "\033[1ACompleted: "+str(data["sample_names"].index(sample_name))+"/"+str(len(data["sample_names"]))

    return out

##########################################################################################################

def calculate_full_matrix(data, samples_to_run, recomb, options, summary_function):
    """
    Calculate the full relatedness matrix for all the samples.
    Can use multiple processes
    """
    
    if( options.get("safe", False) ):
        nn_function=safe_run_for_one_sample
    else:
        nn_function=run_for_one_sample

    info=summary_function.__doc__.strip()

    if options.get("multi_process",0)>1:
        mp = options["multi_process"]
        print info +" using %d processes\n" %(mp)
        pool = Pool(mp)
        args = [ (s, data, recomb, options, summary_function) for s in samples_to_run ]
        pool_results=pool.map(nn_function, args, mp )
        results=dict(zip(samples_to_run, pool_results))
    else:
        print info + ":\n"
        results = {}
        for i, sample in enumerate(samples_to_run):
            print "\033[1A"+sample+" ["+str(i+1)+"/"+str(len(samples_to_run))+"]"
            results[sample] = nn_function((sample, data, recomb, options, summary_function))

    return results    

##########################################################################################################


def main(options):

    if options.get("test_file"):
       data = io.load_minimal_data(options["test_file"])
    elif options.get("vcf_file"):
        data = io.load_vcf_data(options["vcf_file"])
    elif options.get("eigenstrat_root"):
        data = io.load_eigenstrat_data(options["eigenstrat_root"])
    else:
        raise Exception("No input file specified")
    
    recomb = rec.get_recombinator(options["recombination_map"])

    if options["pseudo_haploid"] and np.any(np.equal(data["genotype_data"], 1)):
        raise Exception("Cannot use pseudohaploid algorithm on data with hetozygote sites. "+
                        "Pseudohaploids should be coded as 0 and 2.")

    if not options["pseudo_haploid"] and not np.any(np.equal(data["genotype_data"], 1)):
        raise Exception("All your data is 0 or 2. Are you sure you don't want the pseudohaploid "+
                        "algorithm (-s)?")
    
    # Cut down data if specified, and turn into an array
    max_snps = options.get("max_snps",None)
    if max_snps:
        data["snp_names"]=data["snp_names"][0:max_snps]
        data["snp_pos"]=data["snp_pos"][0:max_snps]
        data["genotype_data"]=data["genotype_data"][0:max_snps]
    data["genotype_data"] = array(data["genotype_data"])
    options["missing_probability"]=np.mean(data["genotype_data"]==3)
    if options["missing_probability"]>0:
        print "Found "+str(int(np.round(options["missing_probability"]*100))) + "% missing genotypes"
    
    samples_to_run=data["sample_names"]
    if options.get("individual", None):
        samples_to_run=options["individual"]

    # Phasing
    results=calculate_full_matrix(data, samples_to_run, recomb, options, ancestry.ancestry_n_tracebacks)
    io.output_phased_data(results, samples_to_run, data["snp_names"], options)
        
##########################################################################################################

if __name__ == "__main__" :
    options = parse_options()
    main(options)
