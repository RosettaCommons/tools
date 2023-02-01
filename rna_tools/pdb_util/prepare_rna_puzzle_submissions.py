#!env python

from sys import argv
from os import system,remove
from os.path import exists
import renumber_pdb_in_place
import parse_tag

def Help():
    print()
    print(argv[0], ' <Puzzle Number> [-renumber <Res Num>, -unvirt <Res Num>]')
    print()
    print('  pdb file names should be specified in order in submit_order.list')
    print()
    exit( 0 )
if len( argv ) == 1:
    Help()

lines = open( 'submit_order.list' ).readlines()
pdbs = [x.replace( '\n', '') for x in lines]

problem_number = argv[1]

new_numbers, chains = [], []
if '-renumber' in argv:
    new_numbers, chains = parse_tag.parse_tag(argv[argv.index('-renumber')+1])
    print(new_numbers, chains)

unvirt_res = [1]
if '-unvirt' in argv:
    unvirt_res += parse_tag.parse_tag(argv[argv.index('-unvirt')+1])[0]
unvirt_res = ' '.join(list(set(map(str, unvirt_res))))


final_submit_pdb = 'DASLAB_Problem%s.pdb' % (problem_number)
if exists( final_submit_pdb): remove(final_submit_pdb)

count = 0
for pdb in pdbs:
    if len(pdb) == 0: continue

    count = count+1

    if len(new_numbers) and not pdb.startswith('default') :
        new_pdb = pdb.replace('.pdb','.renumbered.pdb')
        cmd = ' '.join(['cp', pdb, new_pdb])
        print(cmd)
        system(cmd)
        pdb = new_pdb

        renumber_pdb_in_place.renumber_pdb([pdb], new_numbers, chains = chains)

        new_pdb = pdb.replace('.pdb','.REORDER.pdb')
        cmd = 'reorder_pdb.py '+pdb
        print(cmd)
        system(cmd)
        pdb = new_pdb

    new_pdb = pdb.replace('.pdb','.unvirtualized.pdb')
    command = 'rna_graft -s %s %s -o %s -unvirtualize_phosphate_res %s -graft_backbone_only' % (
        pdb, pdb, new_pdb, unvirt_res
    )
    print(command)
    print(system( command ))
    pdb = new_pdb

    new_pdb = 'DASLAB_Problem%s_Rank%d.pdb' % (problem_number, count )
    command = 'reorder_to_standard_pdb.py %s > %s' % ( pdb, new_pdb )
    print(command)
    system( command )

    command = 'echo MODEL %d >> %s' % (count,final_submit_pdb)
    print(command)
    system( command )

    command = 'cat %s >> %s' % (new_pdb,final_submit_pdb)
    print(command)
    system( command )

    command = 'echo TER >> %s' % (final_submit_pdb)
    print(command)
    system( command )

    command = 'echo ENDMDL >> %s' % (final_submit_pdb)
    print(command)
    system( command )


#command = 'tar cvfz DASLAB_Problem%s.tgz DASLAB_Problem%s_Rank*.pdb' % (problem_number, problem_number)
#print(command)
#system( command )

print( '\nSubmit this combined pdb: ', final_submit_pdb )
