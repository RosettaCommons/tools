#!/usr/bin/env python
#coding: utf8 

note = '	author: Chris King 4/2011 (dr.chris.king@gmail.com)\n\
	\n\
	WHAT IS THIS?\n\
	This script is for post-processing the output of Rosetta\'s pepspec application\n\
	It will sort peptides by Rosetta score, and generate a position-specific frequency matrix in a *.pwm file\n\
	see: Structure-based prediction of protein–peptide specificity in Rosetta, Proteins 2010, 78, 3437 \n\
\n\
	HOW DO I USE IT? \n\
	gen_pepspec_pwm.py pepspec_output n_nterm_res score_frac_cutoff score_name background_pwm\n\
	Essential:\n\
		Pepspec output file.\n\
		Additionally, enter the number of residues n-terminal (to the left) to the anchor residue.\n\
	Recommended:\n\
		Additionally, enter the fraction of low scoring sequences to use\n\
		(default 0.1, or 10%), what score term to use (default binding score minus protein score),\n\
		and the location of the background PWM\n\
		The default backgound PWM can be found in <path_to_rosetta_database>/pepspec_background.binding-prot-0.1.pwm\n\
\n\
	WHAT DOES IT DO?\n\
	This script sorts output peptide sequences by given score type, takes the lowest-energy N% peptide sequences,\n\
	and generates a frequency-based position-weight matrix (PWM). Optionally, a background reference PWM can be\n\
	supplied by the user for normalization. For normalization, the peptide p+0 sequence position is ignored.\n\
'


import sys
import os
import math

#quicksort stolen/modified from http://code.activestate.com/recipes/66473-just-for-fun-quicksort-in-3-lines/
#sort list L of tuples by tuple element n, e.g. qsort( L, 1 ) will sort L on 2nd column
#this sort dies for large N! too many recursions...
#def qsort( L, n ):
#	if len( L ) <= 1: return L
#	return qsort( [ lt for lt in L[ 1: ] if lt[ n ] < L[ 0 ][ n ] ], n ) + [ L[ 0 ] ]  + qsort( [ ge for ge in L[ 1: ] if ge[ n ] >= L[ 0 ][ n ] ], n )

#single-letter AA types
restypes = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

#read command line
if( len( sys.argv ) < 2 ):
	print 'args:\tpepspec_output\tn_nterm_res\t<score_frac_cutoff(0.1)>\t<score_name(binding-prot_score)>\t<background_pwm>\n\"gen_pepspec_pwm.py help\" for more info'
	sys.exit()
sys.argv.pop( 0 )
specfilename = sys.argv.pop( 0 )
if specfilename == 'help':
	print note
	sys.exit()
assert len( sys.argv ) > 0, 'automatic nterm detection not yet implemented!'
nterm = -1 * int( sys.argv.pop( 0 ) )
if len( sys.argv ) > 0: scorecut = float( sys.argv.pop( 0 ) )
else: scorecut = 0.1
if len( sys.argv ) > 0: scorename = sys.argv.pop( 0 )
else: scorename = "binding-prot_score"
if len( sys.argv ) > 0: bgpwmname = sys.argv.pop( 0 )
else: bgpwmname = None

nametag = '.'.join( specfilename.split( '.' )[ :-1 ] )
nametag = nametag + '.' + scorename + '.' + str( scorecut )

#load pepspec output
specfile = open( specfilename, 'r' ).readlines()
#get n lines and scorecut fraction
n_seqs = len( specfile )
n_seqs_cut = int( round( scorecut * n_seqs, 0 ) )

##			WARNING			##
##!!sequence column hardcoded for now!!##
seqcol = 1

#get nres in peptide
nres = len( specfile[ 0 ].split()[ seqcol ] )
cterm = nres + nterm - 1

#get score col
scorecol = None
parse = specfile[ 0 ].split() 
for itag in range( len( parse ) ): 
	tag = parse[ itag ].rstrip( ':' )
	#actual value is next column
	if tag == scorename: scorecol = int( itag + 1 )
if scorecol is None:
	print 'score_name', scorename, 'was not found on first line of', specfilename, '!!'
	sys.exit()

#load seqs and scores as list of string, float 2ples
seqscores = []
for line in specfile:
	parse = line.split()
	seq = parse[ seqcol ]
	score = float( parse[ scorecol ] )
	seqscores.append( ( seq, score ) )


#sort seqscore data
seqscores = sorted( seqscores, key=lambda seqscore: seqscore[1] )

#filter seqscore data
seqscores = seqscores[ :n_seqs_cut ]
seqname = nametag + '.seq'
seqfile = open( seqname, 'w' )

#fill seq counts, nres x 20aa
#replicating a list with * doesn't create copies, it only creates references to the existing objects!
freqmat = [ None ] * nres
for ires in range( nres ):
	freqmat[ ires ] = [ 0.0 ] * len( restypes )
	
#increment matrix from sequences
for seqscore in seqscores:
	seq = seqscore[ 0 ]
	for iseq in range( len( seq ) ):
		res = seq[ iseq ]
		aaidx = restypes.index( res )
		freqmat[ iseq ][ aaidx ] += ( 1.00000 / n_seqs_cut )
	#and write sequence to fasta file
	seqfile.write( seq + '\n' )
seqfile.close()

#print pwm
pwmname = nametag + '.pwm'
pwmfile = open( pwmname, 'w' )
#col header ( seqpos )
for iseq in range( nres ):
	seqpos = str( iseq + nterm )
	pwmfile.write( '\t' + seqpos )
pwmfile.write( '\n' )
#print data lines
for iaa in range( len( restypes ) ):
	aa = restypes[ iaa ]
	pwmfile.write( aa )
	for iseq in range( nres ):
		#pwmfile.write( '\t' + str( freqmat[ iseq ][ iaa ] ) )
		pwmfile.write( '\t' + '%.3f' % freqmat[ iseq ][ iaa ] )
	pwmfile.write( '\n' )

#normalize by bg pwm?
if bgpwmname is None: sys.exit()

#load bgpwm
bgpwm = open( bgpwmname, 'r' ).readlines()
header = bgpwm.pop( 0 ).split()
bgpwm.pop( 0 ) #get rid of newline
bgnterm = int( header[ 0 ] )
bgcterm = int( header[ -1 ] )
bgnres = bgcterm - bgnterm + 1

#init bg freq matrix
bgfreqmat = [ None ] * bgnres
for ires in range( bgnres ):
	bgfreqmat[ ires ] = [ 0.0 ] * len( restypes )

#load bg freqs
assert len( bgpwm ) == len( restypes ), 'background pwm has wrong number of rows!'
for ires in range( len( bgpwm ) ):
	line = bgpwm[ ires ]
	bgfreqs = line.split()
	bgfreqs.pop( 0 )
	for ifreq in range( len( bgfreqs ) ):
		bgfreqmat[ ifreq ][ ires ] = float( bgfreqs[ ifreq ] )

#calc offset, BG PWM MUST SPAN INPUT PWM!
bgoffset = nterm - bgnterm
assert bgnterm <= nterm and bgcterm >= cterm, 'background pwm does not have enough columns!' 

#subtract matrices
colsum = [ 0.0 ] * nres
for iseq in range( nres ):
	#if anchor position, leave it alone
	if iseq + nterm == 0:
		colsum[ iseq ] = 1
		continue
	for ires in range( len( restypes ) ):
		#if bg starts before pwm just skip offset in bg
		diff = freqmat[ iseq ][ ires ] - bgfreqmat[ iseq + bgoffset ][ ires ] 
		if diff < 0: diff = 0
		freqmat[ iseq ][ ires ] = diff
		colsum[ iseq ] += diff

#print bg pwm
normpwmname = nametag + '.norm.pwm'
normpwmfile = open( normpwmname, 'w' )

for iseq in range( nres ):
	seqpos = str( iseq + nterm )
	normpwmfile.write( '\t' + seqpos )
normpwmfile.write( '\n' )
#renormalize and print
for iaa in range( len( restypes ) ):
	aa = restypes[ iaa ]
	normpwmfile.write( aa )
	for iseq in range( nres ):
		#if colsum[ iseq ] == 0, then pwm matched background, so zero information, so equal probabilities!
		if colsum[ iseq ] == 0.0: normpwmfile.write( '\t' +  '%.3f' % ( 1.000 / len( restypes ) ) )
		else: normpwmfile.write( '\t' + '%.3f' % ( freqmat[ iseq ][ iaa ] / colsum[ iseq ] ) )
	normpwmfile.write( '\n' )

