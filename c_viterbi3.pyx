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

# Viterbi algorithm classes for 2-parent viterbi algorithm
# This is the cython version of veiterbi2.py. It should be much faster

from __future__ import division
from scipy import interpolate, sparse
from math import exp, log, fsum
from collections import defaultdict
from viterbi_2d_helpers import transition, emission, pseudohaploid_emission
import numpy as np
cimport numpy as np

import cython
cimport cython

max_num_sparse_elems = 1e6

##########################################################################################################

class calculator(object):
    """
    implement the viterbi algorithm using the supplied transition and emission probabilities
    This finds the two closest related individuals, i.e. people with genotypes compatable with
    being your parents. 
    """

    def __init__(self, data, transition, emission, observed, options={}):
        self.data=data
        self.transition=transition

        self.observed=np.array(observed)
        self.emission=emission
        self.Nx = len(data)
        self.Ny = len(data[0])
        self.options=options
        self.stored_ordered_states=None
        self.frequency=options["used_genotype_frequency"]

        Ns = (self.Ny*(self.Ny-1))//2      # Number of states ( samples^2 but state[0] >= state[1] )

        self.Ns = Ns
        self.initial_p= [1/Ns]*Ns

        # Build a map of (2D) states to indices and reversed. 
        # this: [ (1,0),(2,0),(2,1),..., (N-1,N-2) ]
        self.states = []
        for i in range(self.Ny):
            for j in range(i):
                self.states.append((i,j))

    def convert_emission_to_matrix(self, f):
        """
        Convert the dictionary based lookup from the emission object to an array based lookup
        That is, a 3D array indexed by {0,1,2} - genotypes, containing probabilities.
        """
        cdef np.ndarray[np.float64_t, ndim=3] em=np.zeros((4,4,4), dtype=np.float64)
        prob_dict = self.emission.probabilities
        
        for key,value in prob_dict.items():
            prob=self.emission.emission_probability(key[0:2], key[2], f)
            if callable(prob):
                 em[key]=prob(f)
            else:
                 em[key]=prob

        return(em)


    @cython.boundscheck(False)
    @cython.nonecheck(False)
    @cython.wraparound(False)
    def calculate(self):
        """
        Calculate the viterbi matrix and the traceback matrix
        For the traceback matrix, when there is more than one possible
        choice, we take the first one. We only allow one state
        to change at each state, so you can go from (0,0)->(1,0) but 
        not from (0,0)->(1,1). This is for speed. 

        To save memory, we don't store the whole viterbi matrix, and we 
        only store the traceback elements where something changes. 
        """
        
        cdef bint everything = self.options.get("everything", False) # Do we have to try and phase everything?

        cdef int Nx = self.Nx     # Number of markers
        cdef int Ny = self.Ny  # Number of samples
        cdef int Ns = self.Ns     # Number of states ( samples^2 but state[0] > state[1] )
        cdef int nse = 0        # Counter for number of elements in traceback
        cdef int i,j,k,s0,s1, best_idx, idx, next_idx, tb_k, best_move_idx, gt_sum

        cdef int max_nse =  max_num_sparse_elems
        cdef int max_tb_k = self.options["traceback_lookback_k"]
        cdef np.float64_t best, best_value, tp0, tp1, best_last_V, best_last_V_tp1, best_i, best_j, best_move, thisVal

        # Type all the members we need locally as numpy arrays. 
        cdef np.ndarray[np.int_t, ndim=2] data=self.data
        cdef np.ndarray[np.int_t, ndim=2] states=np.array(self.states)
        cdef np.ndarray[np.int_t, ndim=1] observed=self.observed

        #These are used for building the sparse traceback matrix
        cdef np.ndarray[np.int_t, ndim=1] I=np.zeros(max_nse, dtype=int)
        cdef np.ndarray[np.int_t, ndim=1] J=np.zeros(max_nse, dtype=int)
        cdef np.ndarray[np.int_t, ndim=1] V=np.zeros(max_nse, dtype=int)
        cdef np.ndarray[np.int_t, ndim=2] tb_arr=np.zeros((max_tb_k,Ns), dtype=int)

        # used for bits of the Viterbi algorithm
        cdef np.ndarray[double, ndim=1] initial_p=np.array(self.initial_p)
        cdef np.ndarray[np.float64_t, ndim=1] thisV=np.zeros(Ns, dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=1] lastV=np.zeros(Ns, dtype=np.float64)
        cdef np.ndarray[np.float64_t, ndim=1] best_values=np.zeros(Ny, dtype=np.float64)
        cdef np.ndarray[np.int_t, ndim=1] best_idxes=np.zeros(Ny, dtype=int)
        cdef np.ndarray[np.float64_t, ndim=3] em=np.zeros((4,4,4), dtype=np.float64)    

        # Cache the indices to look up. 
        cdef np.ndarray[np.int_t, ndim=2] state_indices=np.zeros((Ny,Ny), dtype=int)

        # Sparse matrix - not a cdef. Would be nice if we had one.
        t = sparse.lil_matrix((Nx,Ns), dtype=int).tocoo()

        # Fill in the indices to look up
        for i from 1 <= i < Ny:
            for j from 0 <= j < i:
                idx = i*(i-1)//2 + j
                state_indices[i,j]=idx
                state_indices[j,i]=idx
                
        # Setup first row. 
        best_last_V=-1.0
        tb_k=0

        for j from 0 <= j < Ns:
            lastV[j]=initial_p[j]*self.emission.emission_probability((data[0,states[j,0]], data[0,states[j,1]]), observed[0], self.frequency[0])
            tb_arr[0,j]=j
            if best_last_V < lastV[j]:
                best_last_V = lastV[j]
        
        if best_last_V > 0:
            for j from 0 <= j < Ns:
                lastV[j] = lastV[j]/best_last_V 
        else:
            for j from 0 <= j < Ns:
                lastV[j] = 1/Ns 
            
        # For each snp
        for i from 1 <= i < Nx:

            # Get the transition probabilities
            tp = self.transition.single_transition_probability(i)
            tp0 = tp[0]
            tp1 = tp[1]
            tp2 = tp[2]         # Only used if we are phasing everything. 
            best_last_V_tp1=best_last_V*tp1

	    # Get the emission matrix
            em=self.convert_emission_to_matrix(self.frequency[i])
            #Calculate the best transitions for each i, j                
            for j from 0<=j<Ny:
                best = -1.0
                for k from 0<=k<Ny:
                    idx = state_indices[j,k]
                    
                    if lastV[idx] > best:
                        best_values[j]=lastV[idx]*tp1
                        best_idxes[j]=idx
                        best=lastV[idx]

            #traceback if we're on a multiple of chunk size, or at the end
            traceback_this_iteration=(tb_k+1==max_tb_k) or (i+1==Nx)

            # For each state see what the most likely previous state was. 
            for j from 0 <= j < Ns:
                best = lastV[j]*tp0
                best_idx = j
                s0 = states[j,0]
                s1 = states[j,1]

                if best < best_last_V_tp1: # If it might be better to move

                    best_i = best_values[s0] # best value if we let first index vary
                    best_j = best_values[s1] # best value if we let second index vary
                
                    if best_i < best_j:
                        best_move = best_j
                        best_move_idx = best_idxes[s1]
                    else:
                        best_move = best_i
                        best_move_idx = best_idxes[s0]
                    
                    if best < best_move:
                        best=best_move
                        best_idx=best_move_idx

                thisV[j]=best*em[data[i,s0], data[i,s1], observed[i]]
                tb_arr[tb_k,j]=best_idx
                              
                # If we have demanded that we phase *everything* and site we're going to is not phasable 
                # then try and find the best informative state. If none of them are informative give up
                # (but we might phase randomly later)
                if everything:
                    if observed[i-1]==1 and data[i-1,states[best_idx,0]]==data[i-1,states[best_idx,1]]: 
                        best=0
                        best_idx=-1
                        for k from 0<=k<Ns:
                            if data[i-1,states[k,0]]!=data[i-1,states[k,1]]:
                                thisVal=lastV[k]
                                if states[k,0]!=s0 and states[k,0]!=s1 and states[k,1]!=s0 and states[k,1]!=s1:
                                    thisVal *= tp2
                                elif states[k,0]==s0 and states[k,1]==s1:
                                    thisVal *= tp0
                                else:
                                    thisVal *= tp1
                                
                                if thisVal>best:
                                    best=thisVal
                                    best_idx=k

                        if best_idx>-1: # If we found something - flip traceback
                            tb_arr[tb_k,j]=-best_idx-2 # this is negative so that we can tell that we only came here for one site. 

                if traceback_this_iteration:
                    idx = j
                    for k from 0 <= k < tb_k:
                        if idx == -1:
                            break
                        next_idx = tb_arr[tb_k-k,idx]                            

                        tb_arr[tb_k-k,idx]=-1

                        if next_idx!=idx:
                            I[nse]=i-k
                            J[nse]=idx
                            V[nse]=next_idx+1 # Sparse matrix returns 0 for empty, so shift everything by 1
                            nse+=1
                            if nse > max_nse:
                                t=t+sparse.coo_matrix((V,(I,J)),shape=(Nx,Ns),dtype=int)
                                nse=0            
                            if next_idx>=0:
                                idx=next_idx

            # Back to outer loop ( i over Nx ) 
            # Move traceback on, wrapping round if required
            tb_k = (tb_k + 1) % max_tb_k

            best_last_V=-1.0
            for j from 0 <= j < Ns:
                if best_last_V < thisV[j]:
                    best_last_V = thisV[j]

            for j from 0 <= j < Ns: 
                thisV[j] = thisV[j]/best_last_V
                lastV[j] = thisV[j]

            best_last_V=1.0 # We are normalising everything so that the best element==1
                
        # Finally save everything
        self.viterbi = thisV
        t=t+sparse.coo_matrix((V[:nse],(I[:nse],J[:nse])),shape=(Nx,Ns),dtype=int)
        self.traceback_matrix = t.tocsr()
        self.stored_ordered_states=None

    def ordered_viterbi_states(self):
        """
        Get the states in order of the best Viterbi score
        """
        if self.stored_ordered_states==None:
            V = self.viterbi
            to_sort = zip([-v for v in V], range(len(V)))
            to_sort.sort()
            self.stored_ordered_states=[s[1] for s in to_sort]
    
        return self.stored_ordered_states
        
    @cython.boundscheck(False)
    @cython.nonecheck(False)
    @cython.wraparound(False)
    def traceback(self, n_paths=1, use_everything=True):
        """
        Get the traceback of one of the most likely paths
        which path=i gets the i+1th best path.  
        """
        cdef int max_tb_k = self.options["traceback_lookback_k"]
        cdef np.ndarray[np.int_t, ndim=2] t=self.traceback_matrix[(self.Nx-max_tb_k):self.Nx,:].todense()
        cdef int Nx = self.Nx     # Number of markers
        cdef int i,j, t_i, s_i
        cdef np.ndarray[np.int_t, ndim=1] index=np.zeros(n_paths, dtype=int)

        i=0
        j=0

        tb = [[None]*len(self.observed) for k in range(n_paths)]
        one_step_index=[None]*n_paths

        ordered_elems=self.ordered_viterbi_states()
        index=np.array([ordered_elems[s] for s in range(n_paths)])
       
        for i from 0<=i<Nx:
            t_i=i % max_tb_k
            s_i=max_tb_k-t_i-1

            if t_i==0 and i>0:
                start=max(0,self.Nx-i-max_tb_k)-(self.Nx-i-max_tb_k)
                t[start:max_tb_k,]=self.traceback_matrix[max(0,self.Nx-i-max_tb_k):(self.Nx-i),:].todense() # copy the traceback matrix into a c array in chunks for fast indexing
            for j from 0<=j<n_paths:
                if one_step_index[j]:
                    tb[j][i]=self.states[one_step_index[j]]
                    one_step_index[j]=None
                else: 
                    tb[j][i]=self.states[index[j]]

                    back_trace=t[s_i,index[j]]-1
                    if back_trace>=0:
                        index[j]=back_trace
                    if back_trace<=-2 and use_everything:   # jump off the optimal path for one step to avoid unphasable site. 
                        one_step_index[j]=-back_trace-2
       
        [tr.reverse() for tr in tb]
        return tb
        


##########################################################################################################
