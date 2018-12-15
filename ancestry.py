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

# Contains the code which actually does the phasing from tracebacks
# At some point to be replaced with cython code for speed. 

from __future__ import division
import random
import numpy as np
from hmmlearn import hmm

########################################################################################################## 

def ancestry_n_tracebacks(viterbi_object, sample_indices, snp_pos, options, genotypes, observations, sample_index):
    """
    Phasing
    """

    #Phase a bunch of tracebacks and then combine them.
    phases=[]
    phase=None
    parents=None
    ancestry=None
    tbs=viterbi_object.traceback(n_paths= options["n_traceback_paths"], use_everything=False )
                
    parents=[[(sample_indices[p1],sample_indices[p2]) for p1,p2 in order_parents(this_tb)] for this_tb in tbs]
    ancestries=[[(options["populations"][p1],options["populations"][p2]) for p1,p2 in pars] for pars in parents]
    ancestry=combine_ancestry(ancestries)
    ancestry=smooth_ancestry(ancestry, options, snp_pos)

    return {"best_parents":parents[0], "local_ancestry":ancestry}

########################################################################################################## 

def smooth_ancestry(ancestry, options, snp_pos):
    """
    Smooth phasing with a window.
    """
    window=options["window"]
    if window < 2: # no smoothing 
        return ancestry
    
    half_window=int(window/2)
    
    n_sites=len(ancestry)
    smoothed_ancestry=[(None,None)]*n_sites

    for i in range(n_sites):
        start=max(0,i-half_window)
        end=min(n_sites, i+half_window)
        votes=ancestry[start:end]
        unique=list(set(votes))
        scores=[votes.count(u) for u in unique]
        winner=unique[scores.index(max(scores))]
        smoothed_ancestry[i]=winner
    #pdb.set_trace()
    
    return smoothed_ancestry

########################################################################################################## 

def smooth_ancestry_hmm(ancestry, options, snp_pos):
    """
    Smooth phasing with an hmm.
    """
    transition_prob=options["generations"]/2/len(snp_pos)
    entries=list(set(ancestry))
    entries.sort()
    observations=np.array([entries.index(a) for a in ancestry])
    observations=observations.reshape(-1,1)
    
    n_components=len(entries)
    startprob_prior=np.zeros(n_components)+1/n_components
    transmat_prior=np.zeros((n_components,n_components))+transition_prob
    np.fill_diagonal(transmat_prior, 1-(n_components-1)*transition_prob)

    model=hmm.MultinomialHMM(n_components=n_components, startprob_prior=startprob_prior,
                                 transmat_prior=transmat_prior, params='e', init_params='es')

    model.transmat_=transmat_prior
    
    model.fit(observations)

    smoothed_ancestry=[entries[a] for a in model.decode(observations, algorithm="viterbi")[1]]

    pdb.set_trace()

    return smoothed_ancestry
        
########################################################################################################## 

def combine_ancestry(ancestries):
    """
    Combine phasing - takes an array of arrays of phase tuples (and None) and combines them to get a consesus, with 
    all the paths voting equally.
    """
    n_paths=len(ancestries)
    n_sites=len(ancestries[0])

    if(n_paths==1): 
        ancestry=ancestries[0]
        return ancestry

    ancestry=[(None,None)]*n_sites

    for i in range(n_sites):
        votes=[(min(a[i]),max(a[i])) for a in ancestries]
        unique=list(set(votes))
        scores=[votes.count(u) for u in unique]
        winner=votes[scores.index(max(scores))]
        ancestry[i]=winner

    return ancestry
        
########################################################################################################## 

def order_parents(traceback):
    """
    Make sure that the order of the parents in the traceback is consistent - the traceback
    returns each pair ordered by number so we need to make sure that [(1,2), (2,3)] actually
    shows up as [(1,2),(3,2)].
    """
    
    ln=len(traceback)
    for i in range(1,ln):
        last_pair=traceback[i-1]
        this_pair=traceback[i]            

        if((this_pair[0]!=last_pair[0]) and (this_pair[1]!=last_pair[1])): # Need to flip
            traceback[i:]=phase_flip(traceback[i:])
                                  
    return traceback

########################################################################################################## 

def phase_flip(phase):
    """
    Flip phasing
    """
    return [(y,x) for x,y in phase]

########################################################################################################## 
