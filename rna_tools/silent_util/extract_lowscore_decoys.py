#!/usr/bin/env python

from __future__ import print_function
from sys import argv,exit
from os import popen, system
from os.path import basename,exists,expanduser,expandvars
import string,subprocess
from glob import glob

def Help():
    print()
    print('Usage: '+argv[0]+' <silent out file 1> < silent file 2> ... <N> ')
    print('  Will extract N decoys with lowest score from each silent file.')
    print('  If you want to select based on another column, say 12 (Rg), the')
    print('    last arguments should be -12 <N>  (for lowest Rg) or +12 <N>')
    print('    (for highest Rg).')
    print()
    if ( subprocess.call( ['/bin/bash','-i','-c','alias ex']) == 1 ):
        print(' You might consider using an alias "ex" for this function. Add ')
        print('  alias ex="extract_lowscore_decoys.py" ')
        print(' to your ~/.bashrc script.')
    print()
    exit()


if len(argv)<2:
    Help()

replace_names = 1
if argv.count('-no_replace_names'):
    pos = argv.index('-no_replace_names')
    del( argv[pos] )
    replace_names = 0

extract_first_chain = 0
if argv.count('-extract_first_chain'):
    pos = argv.index('-extract_first_chain')
    del( argv[pos] )
    extract_first_chain = 1

start_at_zero = 0
if argv.count('-start_at_zero'):
    pos = argv.index('-start_at_zero')
    del( argv[pos] )
    start_at_zero = 1

use_start_pdb = 0
if argv.count('-start_pdb'):
    pos = argv.index('-start_pdb')
    del( argv[pos] )
    start_pdb_file = argv[ pos ]
    del( argv[pos] )
    use_start_pdb = 1

output_virtual = 0
if argv.count('-output_virtual'):
    pos = argv.index('-output_virtual')
    del( argv[pos] )
    output_virtual = 1

extra_res_fa = []
if argv.count('-extra_res_fa'):
    pos = argv.index('-extra_res_fa')
    del( argv[pos] )
    while pos < len( argv ) and argv[ pos ][0] != '-':
        extra_res_fa.append( argv[ pos ] )
        del( argv[ pos ] )

rosetta_folder = ""
if argv.count('-rosetta_folder'):
	pos = argv.index('-rosetta_folder')
	rosetta_folder = argv[ pos+1 ]
	del( argv[pos+1] )
	del( argv[pos] )
NSTRUCT_defined = 0
try:
    NSTRUCT = int(argv[-1])
    del(argv[-1])
    NSTRUCT_defined = 1
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
match_score_defined = 0
if not scorecol_defined:
    scorecol_name = argv[-1]
    if scorecol_name[0] == '-':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        REVERSE = ''
    if scorecol_name[0] == '+':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        REVERSE = '-r'
    if '=' in scorecol_name:
        scorecol_name_defined = 1
        scorecol_name_split = scorecol_name.split('=')
        scorecol_name = scorecol_name_split[0]
        try:
            match_score = float(scorecol_name_split[-1])
            match_score_defined = 1
            if not REVERSE:
                REVERSE = ''
            if not NSTRUCT_defined:
                NSTRUCT = None
        except:
            match_score_defined = 0
    if scorecol_name_defined:
        del( argv[-1] )

infiles = argv[1:]

if rosetta_folder == "":
	rosetta_folder = expandvars("$ROSETTA")
MINI_DIR = rosetta_folder + '/main/source/bin/'
DB = rosetta_folder + '/main/database/'

# deprecated -- set up for old CASP runs, but probably
#  will never be used again...
HOMEDIR = expanduser('~')

for infile in infiles:
    tags = []

    scoretags = popen('head -n 2 '+infile).readlines()[1].split()
    scoretag=''
    if scorecol_defined:
        scoretag = scoretags[ abs(scorecol) ]

    if scorecol_name_defined:
        scorecol_names = scorecol_name.split(',')
        scorecols = []
        for s in scorecol_names:
            assert( scoretags.count( s ))
            scorecol = scoretags.index( s )
            scorecols.append( scorecol )
        scoretag = scorecol_name
        if match_score_defined:
            scoretag += "_%d" % int(match_score)
    else:
        scorecols  = [scorecol]


    binary_silentfile = 0
    remark_lines = popen('head -n 7 '+infile).readlines()
    for line in remark_lines:
        if ( len( line ) > 6 and line[:6] == "REMARK" ):
            remark_tags = line.split()
            if remark_tags.count('BINARY_SILENTFILE'):
                binary_silentfile = 1
            if remark_tags.count('BINARY'):
                binary_silentfile = 1

    coarse = 0
    if exists( 'remark_tags') and remark_tags.count('COARSE'):
        coarse = 1

    assert(infile[-3:] == 'out')
#    lines = popen('grep SCORE '+infile+' |  sort -k %d -n %s | head -n %d' % (abs(SCORECOL)+1, REVERSE, NSTRUCT+1) ).readlines()


    # Check if this run appeared to use -termini
    terminiflag = ''
    fid = open( infile, 'r')
    line = 'ATOM'
    while (line.count('ATOM') or line.count('SCORE') or
           line.count('SEQU') or line.count('JUMP') or line.count('FOLD')):
        line = fid.readline()
    if line.count('AAV'):
        terminiflag = ' -termini '


    # Make the list of decoys to extract
    lines = popen( 'grep SCORE '+infile+' | grep -v NATIVE').readlines()

    score_plus_lines = []
    for line in lines:
        cols = line.split()
        score = 0.0
        try:
            for scorecol in scorecols: score += float( cols[ abs(scorecol) ] )
        except:
            continue
        if REVERSE: score *= -1
        score_plus_lines.append( ( score, line ))

    score_plus_lines.sort()
    if match_score_defined:
        score_plus_lines = filter( lambda x: float(x[0]) == match_score, score_plus_lines)
    lines = map( lambda x:x[-1], score_plus_lines[:NSTRUCT] )

    templist_name = 'temp.%s.list'% basename(infile)

    fid = open(templist_name,'w')
    count = 0
    for line in lines:
        cols = line.split()
        tag = cols[-1]
        if tag.find('desc') < 0:
            fid.write(tag+'\n')
            tags.append(tag)
            count = count+1
        if NSTRUCT and count >= NSTRUCT:
            break
    outfilename = infile

    fid.close()

    startpdbflag = ''
    if (use_start_pdb): startpdbflag = '-start_pdb '+start_pdb_file

    extract_first_chain_tag = ''
    if (extract_first_chain): extract_first_chain_tag = ' -extract_first_chain '

    #Set up bonds file?
    softlink_bonds_file = 0
    wanted_bonds_file = infile+'.bonds'
    wanted_rot_templates_file = infile+'.rot_templates'
    bonds_files = glob( '*.bonds')
    if len( bonds_files ) > 0:
        if not exists( wanted_bonds_file ):
            softlink_bonds_file = 1
            system( 'ln -fs '+bonds_files[0]+' '+wanted_bonds_file )
            system( 'ln -fs '+bonds_files[0].replace('.bonds','.rot_templates') \
                    +' '+wanted_rot_templates_file )


    # Centroid readout?
    MINI_EXE = MINI_DIR+'extract_pdbs'
    if not exists( MINI_EXE):
        MINI_EXE = MINI_DIR+'extract_pdbs.linuxgccrelease'
    if not exists( MINI_EXE):
        MINI_EXE = MINI_DIR+'/extract_pdbs.macosgccrelease'
        assert( exists( MINI_EXE ) )

    command = '%s -load_PDB_components -in:file:silent  %s   -in:file:tags %s -database %s -out:file:residue_type_set centroid ' % \
                  ( MINI_EXE, outfilename, " ".join(tags), DB )

    # This section (rosetta++) is really old and probably could be removed
    old_rosetta = 0
    scorelabels = popen( 'head -n 2 '+outfilename ).readlines()[-1].split()
    if "SCORE" in scorelabels:
        EXE = HOMEDIR+'/src/rosetta++/rosetta.gcc'
        if not exists( EXE ):
            EXE = 'rm boinc* ros*txt; '+HOMEDIR+'/src/rosetta++/rosetta.mactelboincgraphics '
        assert( exists( EXE ) )
        command = '%s -extract -l %s -paths %s/src/rosetta++/paths.txt -s %s %s %s '% (EXE, templist_name, HOMEDIR,outfilename, terminiflag, startpdbflag+extract_first_chain_tag)
        old_rosetta = 1
        print("OLD_ROSETTA", old_rosetta)

    # Check if this is an RNA run.
    fid = open( infile, 'r')
    line = fid.readline(); # Should be the sequence.
    print(line)
    rna = 0
    sequence = line.split()[-1]
    rna = 1
    for c in sequence:
        if not ( c == 'a' or c == 'c' or c == 'u' or c == 'g'):
            rna = 0
            break
    if rna:     command  += ' -enable_dna -enable_rna '


    #        command = command.replace('rosetta++','rosetta_rna')
    #print "RNA? ", rna

    # Check if this is full atom.
    lines = popen('head -n 8 '+outfilename).readlines()
    if len(lines[6].split()) > 10:
        command += ' -fa_input'

    # Hey this could be a new mini RNA file
    if rna and not old_rosetta:

        command = '%s -load_PDB_components -database %s -in::file::silent %s -tags %s  -extract' % \
                  ( MINI_EXE, DB, outfilename, " ".join(tags))

        if binary_silentfile:
            silent_struct_type = 'binary_rna'
        else:
            silent_struct_type = 'rna'

        command = '%s -load_PDB_components -database %s -in:file:silent %s -in:file:tags %s -in:file:silent_struct_type %s  ' % \
                  ( MINI_EXE, DB,outfilename, " ".join(tags), silent_struct_type )

        if coarse:
            command += " -out:file:residue_type_set coarse_rna "
        else:
            #command += " -out:file:residue_type_set rna "
            # will default to fa_standard, which holds rna residue types now.
            #command += " -chemical:patch_selectors VIRTUAL_RNA_RESIDUE VIRTUAL_PHOSPHATE VIRTUAL_RIBOSE "
            pass

        if output_virtual: command += " -output_virtual "

    elif ( binary_silentfile ):

        command = '%s -load_PDB_components -in:file:silent  %s  -in:file:silent_struct_type binary -in:file:fullatom -in:file:tags %s -database %s  ' % \
                  ( MINI_EXE, outfilename, " ".join(tags), DB )
        if output_virtual: command += " -output_virtual "
        if len( extra_res_fa ) > 0:
            command += " -extra_res_fa "
            for m in extra_res_fa: command += " "+m

        if (scoretags.count('vdw')): command += ' -out:file:residue_type_set centroid '


    print(command)
    system(command)


    # deprecated -- set up for old CASP runs, but probably
    #  will never be used again...
    if outfilename.find('t343')>0:
        command = HOMEDIR+'/python/extract_t343.py %s %s' % (outfilename,
                                                                 " ".join(tags))
        print(command)
        system(command)


    count = 1
    if start_at_zero: count = 0

    if replace_names:
        for tag in tags:
            if scorecol_defined or scorecol_name_defined:
                command = 'mv %s.pdb %s.%s.%d.pdb' % (tag,basename(infile),scoretag,count)
            else:
                command = 'mv %s.pdb %s.%d.pdb' % (tag,basename(infile),count)
            print(command)
            system(command)
            count += 1

    command = 'rm '+templist_name
    print(command)
    system(command)

    if (softlink_bonds_file):
        #system( 'rm '+wanted_bonds_file+' '+wanted_rot_templates_file )
        print(' WARNING! WARNING')
        print(' Found a .bonds and .rot_templates file and used it!')

