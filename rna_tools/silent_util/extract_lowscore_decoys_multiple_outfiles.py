#!/usr/bin/env python

from sys import argv,exit,stderr,stdout
from os import popen, system
from os.path import basename
import string,subprocess

def Help():
    print
    print 'Usage: '+argv[0]+' <silent out file 1> < silent file 2> ... <N> '
    print '  Will extract N models with lowest score from each silent file.'
    print '  If you want to select based on another column, say 12 (Rg), the'
    print '    last arguments should be -12 <N>  (for lowest Rg) or +12 <N>'
    print '    (for highest Rg).'
    print '  If N is a float, it will be treated as a score cutoff, rather than'
    print '    desired number of models.'
    print
    if ( subprocess.call( ['/bin/bash','-i','-c','alias exo']) == 1 ):
        print ' You might consider using an alias "exo" for this function. Add '
        print '   alias exo="extract_lowscore_decoys_multiple_outfiles.py'
        print ' to your ~/.bashrc script.'

    exit()

if len(argv)<2:
    Help()

output_to_file = True

SCORE_CUTOFF = 0

try:
    NSTRUCT = float(argv[-1])

    if argv[-1].count( '.' ) > 0: # user specified a float...
        SCORE_CUTOFF = NSTRUCT
        NSTRUCT = 0
    else:
        NSTRUCT = int( NSTRUCT )

    del(argv[-1])
except:
    NSTRUCT = 2

scorecol_defined = 0
try:
    scorecol = int(argv[-1])
    del(argv[-1])
    scorecol_defined = 1
except:
    scorecol = -1

REVERSE = ''
if scorecol > 0:
    REVERSE = ' --reverse '

#Another possibility... user supplies -rms or +rms
scorecol_name_defined = 0
if not scorecol_defined:
    scorecol_name = argv[-1]
    if scorecol_name[0] == '-':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        del( argv[-1] )
        REVERSE = ''
    if scorecol_name[0] == '+':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        REVERSE = '-r'
        del( argv[-1] )

infiles = argv[1:]
IS_OUTFILE = 1


for infile in infiles:

    score_plus_lines = []
    firstlines = []

    if len( firstlines ) == 0:
        firstlines = popen('head -n 5 '+infile).readlines()
        scoretags = string.split( firstlines[1] )
        if firstlines[0].find( "SEQUENCE" ) < 0  and   firstlines[0].find("SCORE:") >= 0:
            IS_OUTFILE = 0
            scoretags = string.split( firstlines[0] )

    scoretag=''
    if scorecol_defined:
        scoretag = scoretags[ abs(scorecol) ]

    if scorecol_name_defined:
        scorecol_names = string.split( scorecol_name,',' )
        scorecols = []
        for s in scorecol_names:
            assert( scoretags.count( s ))
            scorecol = scoretags.index( s )
            scorecols.append( scorecol )
        scoretag = scorecol_name
    else:
        scorecols  = [scorecol]

    #####################################
    # Make the list of decoys to extract
    #####################################
    command = 'grep SCORE '+infile+' | grep -v NATIVE'
    lines = popen( command ).readlines()

    for line in lines:
        cols = string.split( line )
        score = 0.0

        try:
            for scorecol in scorecols: score += float( cols[ abs(scorecol) ] )
        except:
            continue

        if ( NSTRUCT == 0 ):
            if REVERSE:
                if (score < SCORE_CUTOFF ): continue
            else:
                if (score > SCORE_CUTOFF ): continue

        if REVERSE: score *= -1
        score_plus_lines.append( ( score, line, infile ))

    score_plus_lines.sort()

    tags_for_infile = {}
    tags_for_infile[ infile ] = []

    count = 0
    for score_plus_line in score_plus_lines:
        line = score_plus_line[1]
        cols = string.split(line)
        tag = cols[-1]

        infile_for_tag = score_plus_line[2]
        tags_for_infile[ infile_for_tag ].append( tag )

        count += 1
        if ( NSTRUCT > 0 and count >= NSTRUCT ): break

    fid = stdout
    if output_to_file:
        newfile = infile
        if scorecol_name_defined: newfile = newfile.replace( ".out", ".%s.out" % scorecol_name )
        newfile = newfile.replace( ".out", ".top%d.out" % NSTRUCT  )
        fid = open( newfile , 'w' )

    # Report to command-line
    stderr.write( 'Extracting %d models from %s to %s\n' % (count, infile, newfile ) )

    #####################################
    # Output header
    #####################################
    if not IS_OUTFILE:
        command = 'head -n 1 '+infile
        lines = popen(command).readlines()
    elif (firstlines[2][:6] == 'REMARK' ):
        command = 'head -n 3 '+infile
        lines = popen(command).readlines()
    elif (firstlines[4][:6] == 'REMARK' ):
        command = 'head -n 5 '+infile
        lines = popen(command).readlines()
        del( lines[2] )
        del( lines[2] )
    else:
        command = 'head -n 2 '+infile
        lines = popen(command).readlines()

    for line in lines: fid.write( line )

    ok_tags = tags_for_infile[ infile ]
    if len( ok_tags ) == 0: continue

    data = open(infile,'r')
    line = data.readline() # Skip first two lines
    line = data.readline()

    ###################################################
    # output silent file lines with correct decoy tag
    ###################################################
    writeout = 0
    while line:
        line = data.readline()[:-1]
        if len( line ) < 2: continue
        if line[:9] == 'SEQUENCE:': continue # Should not be any more sequence lines!
        cols = string.split( line )
        tag = cols[-1]

        if cols[0][:5] == "SCORE":
            if ok_tags.count( tag ) > 0: writeout = 1
            else: writeout = 0

        if not writeout: continue

        # following was for concatenating amongst several outfiles
        #if (n > 0):
        #    tagcols = string.split( tag, '_')
        #    try:
        #        tagnum = int( tagcols[-1] )
        #       tagcols[-1] = '%04d_%06d' %  ( n, tagnum )
        #        newtag = string.join( tagcols,'_')
        #        description_index = line.find( tag )
        #        line = line[:description_index] + newtag
        #    except:
        #        continue

        fid.write(  line+"\n" )

    data.close()

    if output_to_file:
        fid.close()
        #stderr.write( "Outputted models to: %s\n" % newfile )
