#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   calculate_orientational_distances.py
## @brief  Calculates orientational distance between all antibodies in database
## @author Nick Marze

import os, sys, numpy

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main(args):
    lines = [line.strip() for line in open('angles.sc')]
    description = []
    distance = []
    opening = []
    opposite = []
    packing = []
    for i in range(1, len(lines)):
        split_line = lines[i].split(",")
        description.append(split_line[0])
        distance.append(float(split_line[1]))
        opening.append(float(split_line[2]))
        opposite.append(float(split_line[3]))
        packing.append(float(split_line[4]))
    matrix = []
    stdev = [numpy.std(distance), numpy.std(opening), numpy.std(opposite), numpy.std(packing)]
    for j in range(0, len(description)):
        distance_column = []
        for k in range(0, len(description)):
            a = (distance[j]-distance[k])/stdev[0]
            b = (opening[j]-opening[k])/stdev[1]
            c = (opposite[j]-opposite[k])/stdev[2]
            d = (packing[j]-packing[k])/stdev[3]
            entry = (a * a) + (b * b) + (c * c) + (d * d)
            distance_column.append(entry)
        matrix.append(distance_column)
    with open('comparisons.txt', 'w') as f:
        f.write("pdb_code ")
        f.write("%s\n" % " ".join([str(n) for n in description]))
        for x in range(0,len(matrix)):
            f.write("%s " % description[x])
            f.write("%s\n" % " ".join([str(n) for n in matrix[x]]))

if __name__ == "__main__": main(sys.argv)
