#############################################################################
#
#   Copyright 2011 Iain Mathieson
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

# Contains functions to precluster the data. Generally when we have many individuals, 
# we try to cut down the N^2 dimensionality by selecting only some subset of the individuals
# to run on, 

from __future__ import division
import numpy as np
from numpy import array

##########################################################################################################

def closest_n( data, sample_names, observations, n, method="incompatable" ):
    """
    Select the n individuals from the dataset which are closest to observation, 
    using a given method to rank them. Returns the reduced dataset and the 
    selected sample names
    """

    dists = distance_methods[method](data, observations)

    dists = zip(dists, sample_names)
    dists.sort()
    top_n_samples=[x[1] for x in dists[:n]]

    new_samples = [s for s in sample_names if s in top_n_samples]
    include = array([(s in top_n_samples) for s in sample_names])
    new_data = data[:,include]

    return new_data, new_samples

##########################################################################################################
def incompatable_distance(data, observations):
    """
    count the number of incompatable sites (0 vs 2), weight by
    allele frequency, rank
    """
    
    frequency = data.mean(axis=1)
    frequency = np.choose(frequency>0, (-1,frequency))
    weights = np.choose(frequency>0, (0, 1/frequency))

    def  obs_mult(x): return (abs(x-observations)==2) * weights
    
    incomp_array = np.apply_along_axis(obs_mult, 0, data)
    scores = incomp_array.sum(axis=0)

    return scores

##########################################################################################################

distance_methods = { "incompatable": incompatable_distance }

##########################################################################################################
