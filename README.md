##Requirements:

LACE requires Python version 2.6 or 2.7. In addition it requires the standard packages
numpy, scipy and  cython. It also requires the pyEigenstrat package available from https://github.com/mathii/pyEigenstrat.

##Install:

python setup.py install

##Test
We include some test data to check that LACE has installed correctly

#Diploid mode:
python lace.py -m testdata.gt.txt -o testout -p testpops.txt

#Pseudohaploid mode: 
python lace.py -m testdata.phgt.txt -o testout.ph -p testpops.txt -s

##Options

To see available options run:

python lace.py -h

Key options include: 
-e specify input data in Eigenstrat format (i.e. "-e data" implies data.geno, data.snp and data.ind as input files)
-n list of individuals to use as source populations (file, with one ID per line)
-i list of individuals to find ancestry (file, with one ID per line, or comma separated list, or single individual)

