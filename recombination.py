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

# Dealing with recombination rates and the genetic map

from scipy import interpolate
import gzip

##########################################################################################################

def get_recombinator(recombination_data):
    """
    If recombinationd data looks like float, then return a constant recombinator. If it looks like 
    it might be a file, then try get a map based recombinator
    """
    try:
        rec_rate=float(recombination_data)
        print "Using constant recombination rate " + str(rec_rate) + " cm/Mb"
        return constant_recombinator(rec_rate)
    except ValueError:
        print "Loading recombination map from " + recombination_data
        return recombinator(recombination_data)

##########################################################################################################

class recombinator(object):
    """
    Class to interpolate the genetic distance between two positions, given am IMPUTE stype recombination map.
    This is c_j d_j from Li and Stephens. 
    """
    
    def __init__(self, recombination_map):
        
        recombination_file=None
        if recombination_map[-3:]==".gz":
            recombination_file=gzip.open(recombination_map, "r")
        else:
            recombination_file=open(recombination_map, "r")
        
        self.position=[]
        self.rate=[]
        self.dist=[]

        header=recombination_file.next()
        header_map=self.parse_header(header)
        
        for line in recombination_file:
            data=line[:-1].split()
            self.position.append(int(data[header_map["POS"]]))
            self.rate.append(float(data[header_map["RATE"]]))
            self.dist.append(float(data[header_map["MAP"]]))
        
        recombination_file.close()

        self.max_pos = max(self.position)
        self.min_pos = min(self.position)

        self.fitter = interpolate.UnivariateSpline(self.position, self.dist, k=1, s=0) # Linear interpolation

    def parse_header(self, header_line):
        """
        Try and parse the header line and figure out what are the position, rate and dist columns
        """
        bits=header_line.split()
        
        header_map={"POS":None, "RATE":None, "MAP":None}
        
        for what in header_map.keys():
            where=[y.upper()[:len(what)]==what for y in bits]
            if sum(where)!=1:
                raise Exception("Check that map file has headers \"Pos\", \"Rate\" and \"Map\" (in CM)")
            else:
                header_map[what]=[i for i in range(len(where)) if where[i]][0]

        return header_map
                
    def distance(self, position_1, position_2):
        """
        Return the genetic distance in cm between two points
        """
        map_pos = self.fitter([position_1, position_2])
        return map_pos[1] - map_pos[0]

##########################################################################################################

class constant_recombinator(object):
    """
    Class to give the genetic distance between populations, given a constant recombination rate (cm/mb) 
    """
    
    def __init__(self, recombination_rate):
        self.rate=recombination_rate

    def distance(self, position_1, position_2):
        """
        Return the genetic distance between two points
        in cM. The rec rate is given in cm/mb
        """
        return self.rate*(position_2-position_1)*10e-6 

##########################################################################################################

