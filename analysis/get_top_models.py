#!/usr/bin/env python

#Jared Adolf-Bryfogle
#Basic script to get output of top models for Rosetta MPI.
#No error checking, just read the parser options and pass at least --file
import os
import re
import sys
import numpy
import heapq

from optparse import OptionParser, IndentedHelpFormatter
from collections import defaultdict



def get_score_indexes(score_terms, score_term_line):
    """
    Gets the index of the score_terms listed from the second line of the .sc file
    Returns as map score_term:index
    """
    headerSP = score_term_line.strip().split()
    delta_indexes = defaultdict();
    for score_term in score_terms:
        for delta_index in range(0, len(headerSP)):
            if headerSP[delta_index] == score_term:
                delta_indexes[score_term] = delta_index
                break;
            
    return delta_indexes

def get_model_of_energy(data, e_term, energy):
    """
    Get the model of a particular score term of particular energy.
    """
    return data[e_term]['model'][data[e_term]['e'].index(energy)]
    
    
def get_lowest_energies(data, e_term, n):
    """
    Get a list of the lowest n energies of a particular score term
    """
    #Heapq returned an error
    energies = sorted(data[e_term]['e'])
    result = [energies[i] for i in range(0, n)]
    repr(result)

    return result
    
    
if __name__== '__main__':
    
    parser = OptionParser()
    args = sys.argv
    
    parser.add_option("--file", "-f",

        help = "Input Score file")
    
    parser.add_option("--top", "-t",
        type = 'int',
        default=10,
        help = "Number of structures top structures")
    
    parser.add_option("--result", "-r",
        default="results.txt",
        help = "File to write results to")
    
    parser.add_option("--score_terms", "-s",
        default='total_score',
        help = "Score terms to get top models.  For multiple score terms, split by, Ex: total_score,interface_delta_x")
    
(options, args) = parser.parse_args(args=args[1:])

FILE = open(options.file, 'r')
FILE.readline()
header = FILE.readline()

score_terms = options.score_terms.split(',')
delta_indexes = get_score_indexes(score_terms, header)



data = defaultdict(dict)

for line in FILE:
    line = line.strip()
    lineSP = line.split()
    file = lineSP[-1]
    
    for score_term in sorted(delta_indexes):
        if not data.has_key(score_term):
            data[score_term]['e'] = []
            data[score_term]['model'] = []
        delta = float(lineSP[delta_indexes[score_term]])
        data[score_term]['e'].append(delta)
        data[score_term]['model'].append(file)
    
    
FOUT = open(options.result, 'w')

for e_term in score_terms:
    FOUT.write("#e_term min max avg std min_model\n")
    delta_min = min(data[e_term]['e'])
    delta_max = max(data[e_term]['e'])
    delta_mean = numpy.mean(data[e_term]['e'])
    delta_min_model = get_model_of_energy(data, e_term, delta_min)
    delta_sd = numpy.std(data[e_term]['e'])
    line = e_term+" %.3f"%delta_min+" %.3f"%delta_max + " %.3f"%delta_mean+" %.3f"%delta_sd+" "+repr(delta_min_model)+"\n"

    #print line
    FOUT.write(line+"\n")
    FOUT.write("#e_term top_num model energy\n")
    energies = get_lowest_energies(data, e_term, options.top)
    models = []
    i = 1
    for e in energies:
        model = get_model_of_energy(data, e_term, e)
        models.append(model+".pdb")
        
        FOUT.write(e_term+" "+" "+repr(i)+" "+model+"   "+"%.3f"%e+"\n")
        i+=1
    for_pymol = " ".join(models)
    print e_term
    print "pymol "+for_pymol
    FOUT.write("\n")
    FOUT.write("#"+e_term+" "+"PYMOL COMMAND:\t")
    FOUT.write("pymol "+for_pymol+"\n\n")
FOUT.close()
FILE.close()
