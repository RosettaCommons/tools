# FILENAME
# compute_pnear.py
#
# AUTHOR
# Vikram K. Mulligan, Flatiron Institute (vmulligan@flatironinstitute.org).
#
# DATE
# Created 10 October 2019.
#
# DESCRIPTION
# A Python script to convert fold propensity (PNear) from a table of scores and RMSD values.  PNear is a value that
# ranges from 0 to 1, with 0 representing no propensity to favour the desired conformation and 1 representing 100%
# propensity to favour that conformation.  The PNear computation takes threeinputs: a file containing scores and RMSDs
# in columns, a value "lambda", and a value "kbt".  The lambda value,  in Angstroms, indicates the extent to which
# a structure may deviate from the native state and still be considered "native-like"; typical values are 1.5 to 2.0 A
# for peptides, and perhaps 2.0 to 4.0 A for proteins.  The kbt value, in kcal/mol, determines how large the energy
# gap between native and nonnative must be in order for a design to be considered a good folder.  Typically, 0.62 is
# used (a value corresponding to physiological temperature).
#
# USAGE
# -- First, prepare a file containing two columns, the first with score values and the second with RMSD values.  If you
# have a binary silent file from ab initio structure prediction, for example, you can use:
#     grep SCORE myfile.silent | awk '{print $2, $22}` > score_vs_rmsd.txt
# Note that in the above, the column indices ($2 and $22) may need to be updated.
# -- Next, run this script, providing the input file, lambda, and kbt values as input:
#     python3 compute_pnear.py score_vs_rmsd.txt 1.5 0.62 > pnear_calculation.log
# -- The PNear value will appear at the aned of pnear_calculation.log

# Includes:
from math import exp
import sys

# Given a vector of scores, a matching vector of rmsds, and values for lambda and kbt,
# compute the PNear value.
def calculate_pnear( scores, rmsds, lambda_val=1.5, kbt=0.62 ) :
    nscores = len(scores)
    assert nscores == len(rmsds), "Error in calculate_pnear(): The scores and rmsds lists must be of the same length."
    assert nscores > 0, "Error in calculate_pnear(): At least one score/rmsd pair must be provided."
    assert kbt > 1e-15, "Error in calculate_pnear(): kbt must be greater than zero!"
    assert lambda_val > 1e-15, "Error in calculate_pnear(): lambda must be greater than zero!"
    minscore = min( scores )
    weighted_sum = 0.0
    Z = 0.0
    lambdasq = lambda_val * lambda_val
    for i in range( nscores ) :
        val1 = exp( -( rmsds[i] * rmsds[i] ) / lambdasq )
        val2 = exp( -( scores[i] - minscore ) / kbt )
        weighted_sum += val1*val2
        Z += val2
    assert Z > 1e-15, "Math error in calculate_pnear()!  This shouldn't happen."
    return weighted_sum/Z


# Given a file, a lambda value (in Angstroms) , and a kbt value (in kcal/mol), compute
# the PNear value.
def calculate_pnear_given_file( filename, lambda_val=1.5, kbt=0.62 ) :
    print( "Reading file \"" + filename + "\".")
    with open( filename, 'r' ) as file:
        lines = file.readlines()
    print( "Read file \"" + filename + "\".  Beginning parse...")

    first = True
    scores = []
    rmsds = []
    
    for line in lines :
        linestripped = line.strip()
        linesplit = linestripped.split()
        assert len(linesplit) == 2, "Error in calculate_pnear_given_file(): Could not parse line \"" + linestripped + "\" in file \"" + filename + "\"."
        try:
            curscore = float( linesplit[0] )
            currmsd = float( linesplit[1] )
        except ValueError:
            print( "Skipping line \"" + linestripped + "\" in file \"" + filename + "\"." )
            continue
        if( first == True ):
            first = False
        scores.append( curscore )
        rmsds.append( currmsd )
    assert first == False, "No valid score-vs-rmsd lines were found in file \"" + filename + "!"

    return calculate_pnear( scores, rmsds, lambda_val=lambda_val, kbt=kbt )

assert len( sys.argv ) == 4, "Error in compute_pnear.py!  Three arguments are expected: the score-vs-rmsd file name, the lambda value, and the kbt value."
pnear = calculate_pnear_given_file( sys.argv[1], lambda_val=float(sys.argv[2]), kbt=float(sys.argv[3]) )
print( "PNear:\t" + str(pnear) )
