#!/usr/bin/env python

from optparse import OptionParser
from sys import exit, stderr, stdout
from os import popen, system, path
from os.path import exists
from glob import glob

from amino_acids import longer_names
from renumber_structure_profile_checkpoint_back_to_seqres import renumber_structure_profile_checkpoint_back_to_seqres
from parse_dssp_results import parse_dssp_results
from jump_over_missing_density import jump_over_missing_density
from numbering_back_to_pdbseqres import numbering_back_to_pdbseqres
import ConfigParser

#parser = ArgumentParser()
parser = OptionParser()
parser.add_option("-p", dest="pdb_fn", help="five-letter code + chain pdb id  eg. 2oxgZ")
parser.add_option("--no_structure_profile", action="store_true", default=False, help="don't get structure_profile check point, use sequence profile instead.")
parser.add_option("-d", dest="debug", action="store_true", help="debug: print out some information")
parser.add_option("-n", dest="dry_run", action="store_true", help="dry run: don't run the system( cmd )")
(options,args) = parser.parse_args()
#options = parser.parse_args()

if not options.pdb_fn:
    parser.print_help()
    exit()

if len(options.pdb_fn) != 5 :
    stderr.write("ERROR: five-letter pdb name needed: pdbid + chain\n")
    parser.print_help()
    exit()

options.pdb_fn += ".pdb"

if options.dry_run:
    print "-----------------------------------"
    print "             dry_run"
    print "-----------------------------------"

## get PDB2VALL_PATH
PDB2VALL_PATH = path.abspath(path.dirname(__file__)) + "/"

## read config
config = ConfigParser.RawConfigParser(allow_no_value=True)
config.read(PDB2VALL_PATH + "pdb2vall.cfg")

## EXTERNAL PROGRAMS PATHES
if not PDB2VALL_PATH:
    stderr.write("ERROR: you should specify the path where your packages are first.\n"); exit()

script_make_sequence_profile_checkpoint = PDB2VALL_PATH + "sequence_profile_scripts/run_psiblast_filtnr_tight.pl"
script_get_structure_profile_checkpoint = PDB2VALL_PATH + "structure_profile_scripts/get_structure_profile_checkpoint.py"
relax_sequence_file = PDB2VALL_PATH + "relax_seqeuence_file.txt"

## ROSETTA APP
ROSETTA_PATH = PDB2VALL_PATH + "../../../";
if config.get('pdb2vall', 'rosetta_path'):
    ROSETTA_PATH = config.get('pdb2vall', 'rosetta_path')
if not ROSETTA_PATH.endswith(path.sep):
    ROSETTA_PATH += path.sep
idealization_app = ROSETTA_PATH + "main/source/bin/idealize_jd2.linuxgccrelease"
relax_app = ROSETTA_PATH + "main/source/bin/relax.linuxgccrelease"
rosetta_database = ROSETTA_PATH + "main/database/"

####################

pdb_id           = options.pdb_fn[:4]
pdb_chain        = options.pdb_fn[4]
pdb              = pdb_id.lower() + pdb_chain.upper()
fasta_output_fn  = pdb_id.lower() + pdb_chain.upper()+".fasta"


## GET FASTA FROM pdbseqres.txt
print "pdb2vall.py: getting fasta from pdb_seqres.txt ..."
tmp_results = numbering_back_to_pdbseqres( pdb )

if not tmp_results:
    stderr.write("ERROR: can't get results from numbering_back_to_pdbseqres().\n"); exit()

missing_den_list = tmp_results[ "missing_density_rsn_list" ]


## Dictionary used
fasta_dict                 = tmp_results[ "seq_from_pdbseqres_Dict" ]
seq_checkpoint_dict        = {}  ##
strpro_checkpoint_dict     = {}  ## from renumber_structure_profile_checkpoint_back_to_seqres
secstr_dict                = {}  ## from idealized rosetta output "xxxx_0001_0001.pdb"; the numbering is according to the missing density regions, recorded in disorder_dict.
idealized_pdb_xyz_dict     = {}  ## from idealized rosetta output "xxxx_0001_0001.pdb"; the numbering is according to the missing density regions, recorded in disorder_dict.
idealized_pdb_torsion_dict = {}  ## from idealized rosetta output "xxxx_0001_0001.pdb"; the numbering is according to the missing density regions, recorded in disorder_dict.


## GET SEQUENCE CHECKPOINT FILE
cmd = script_make_sequence_profile_checkpoint + ' ' +  fasta_output_fn
print "pdb2vall.py: running run_psiblast_filtnr_tight.pl to make sequence profile checkpoint..."
if not options.dry_run:
    system( cmd )

if exists( fasta_output_fn[:5] + ".checkpoint" ):
    print "pdb2vall.py: checking for each position whether rsd in seq_checkpoint file is the same as that in pdbseqres_fasta."
    checkpoint_count = 1
    checkpoint_file = open( fasta_output_fn[:5] + ".checkpoint", "r").readlines()
    for line in checkpoint_file:
        if len(line) > 10:
            #print checkpoint_count
            if checkpoint_count in missing_den_list:
                #print "missing", checkpoint_count, missing_den_list
                checkpoint_count += 1
                continue
            else:
                rsn = checkpoint_count
                line_edit = line.strip().split()[1:]
                if fasta_dict[ rsn ] == line.strip().split()[0]:
                    if options.debug:
                        print checkpoint_count, line_edit
                    seq_checkpoint_dict[ rsn ] = " ".join( x for x in line_edit)
                    checkpoint_count += 1
                else:
                    stderr.write("ERROR:pdb2vall.py: sequence checkpoint rsd doesn't match the fasta from seqres\n"); exit()
else:
    stderr.write("ERROR:pdb2vall.py: no ", fasta_output_fn[:5] + ".checkpoint file has been made!\n"); exit()



## GET STRUCTURE PROFILE CHECKPOINT
if not options.no_structure_profile:
    if exists("%s.50.9mers.ali.fasta.new.blast.checkpoint" % pdb):
        print "checkpoint file", glob("*9mer*checkpoint")[0], "has existed!"
        print "-"*75
    else:
        cmd = script_get_structure_profile_checkpoint + " " + options.pdb_fn
        print "pdb2vall.py get_structure_profile_checkpoint.py: making structure profile checkpoint..."
        print "cmd", cmd, "\n"
        print "-"*75
        if not options.dry_run:
            system( cmd )

    # RENUMBERED IT AND GET A STRUCTURE PROFILE CHECKPOINT_DICT BACK
    if exists("%s.50.9mers.ali.fasta.new.blast.checkpoint" % pdb):
        print "pdb2vall.py: renumber_structure_profile_checkpoint_back_to_seqres(): renumbering structure profile checkpoint back to pdbseqres..."
        strpro_checkpoint_dict = renumber_structure_profile_checkpoint_back_to_seqres( fasta_output_fn, "%s.50.9mers.ali.fasta.new.blast.checkpoint" %pdb, fasta_dict )
    else:
        stderr.write("ERROR:pdb2vall.py: %s.50.9mers.ali.fasta.new.blast.checkpoint doesn't exist!\n" % pdb )
        stderr.write("ERROR:pdb2vall.py: structure_profile_checkpoint generation fails\n"); exit()

# IDEALIZE PDB AND GET ITS TORSION ANGLE
if exists( pdb + "_0001.pdb"):
    print "pdb2vall.py: detected the idealized pdb - %s_0001.pdb" % pdb
else:
    print "pdb2vall.py: trying to run idealize app."
    cmd = idealization_app + " -in:file:s " + pdb + ".pdb -out:file:output_torsions -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm -chemical:exclude_patches ActeylatedProteinNterm -chainbreaks -overwrite -database %s" % rosetta_database
    print cmd
    if not options.dry_run:
        system( cmd )
    if not exists( pdb + "_0001.pdb"):
        stderr.write("ERROR:pdb2vall.py: idealized app failed at making %s_0001.pdb\n" % pdb ); exit()

if len( fasta_dict ) > 800:
		system("cp  %s_0001.pdb %s_0001_0001.pdb" %(pdb, pdb))

# RELAX WITH -relax::constrain_relax_to_start_coords
if exists( pdb + "_0001_0001.pdb"):
    print "pdb2vall.py: detected the relaxed pdb - %s_0001_0001.pdb" % pdb
else:
    print "pdb2vall.py: trying to run relax app."
    cmd = relax_app + " -in:file:s " + pdb + "_0001.pdb  -relax::default_repeats 1 -relax:fast -in:file:fullatom -ignore_unrecognized_res -use_input_sc -relax::sequence_file " + relax_sequence_file + " -out:file:output_torsions -chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm -chemical:exclude_patches ActeylatedProteinNterm -relax::constrain_relax_to_start_coords -overwrite -database %s" % rosetta_database
    print cmd
    if not options.dry_run:
        system( cmd )
    if not exists( pdb + "_0001_0001.pdb"):
        stderr.write("ERROR:pdb2vall.py: relax app failed at making %s_0001_0001.pdb\n" % pdb ); exit()

# RUN DSSP AND PARSE THE RESULTS
print "pdb2vall.py: parse_dssp_results(): using the UNidealized pdb to run dssp to get  1.solvent assessbilty  2.dssp_phi  3.dssp_psi ..."
dssp_Dict = parse_dssp_results( options.pdb_fn, fasta_dict )

# RENUMBER IDEALIZED PDB BACK TO PDB_SEQRES
if not exists( pdb + "_0001_0001.pdb"):
    stderr.write("ERROR:pdb2vall.py: cannot detected the idealized and relaxed pdb - %s_0001_0001.pdb\n" % pdb ); exit()
else:
    ## XYZ COORDINATES FROM IDEALIZED PDB
    idealized_pdb_xyz_lines = jump_over_missing_density( fasta_output_fn, fasta_output_fn[:5], "xyz")
    for line in idealized_pdb_xyz_lines:
        rsd   = line[17:20].strip()
        rsn   = int( line[22:26].strip())

        if rsn > len( fasta_dict.keys()): break

        if options.debug:
            print rsn, "fasta_dict", fasta_dict[ rsn ], longer_names[ rsd ], "rsd from pdb"

        if fasta_dict[ rsn ] == longer_names[ rsd ]:
            xcord = line[30:38].strip()
            ycord = line[38:46].strip()
            zcord = line[46:54].strip()
            idealized_pdb_xyz_dict[ rsn ] = [ xcord, ycord, zcord ]
        else:
            stderr.write("ERROR: pdb2vall.py: for position %s of %s.pdb, xyz rsd from idealized pdb doesn't match the fasta from seqres\n" %( rsn, pdb )); exit()

    if options.debug:
        for index in idealized_pdb_xyz_dict.keys():
            print index, idealized_pdb_xyz_dict[ index ]

    ## TORSIONS FROM IDEALIZED PDB
    idealized_pdb_torsion_lines = jump_over_missing_density( fasta_output_fn, fasta_output_fn[:5], "torsion")
    for line in idealized_pdb_torsion_lines:
        rsn = int( line.split()[0] )
        rsd = line.split()[1]
        secstr = line.split()[2]

        if rsn > len( fasta_dict.keys()): break
        if rsd == 'X': break

        if fasta_dict[ rsn ] == rsd :
            phi   = line.split()[3]
            psi   = line.split()[4]
            omega = line.split()[5]

            idealized_pdb_torsion_dict[ rsn ] = [ phi, psi, omega ]
            secstr_dict[ rsn ]                = secstr

            ## CHECK WHETHER THE IDEALIZED_PHI, PSI IS SIMILAR TO DSSP_PHI, PSI - should I?

            #print rsn, rsd, secstr, phi, psi, omega
        else:
            stderr.write("ERROR:pdb2vall.py: for position %s of %s.pdb, torsions of rsd from idealized pdb doesn't match the fasta from seqres %s\n" %( rsn, pdb, rsd )); exit()

    if options.debug:
        for index in idealized_pdb_xyz_dict.keys():
            print index, idealized_pdb_xyz_dict[ index ]

print "-"*75
print "pdb2vall.py is trying to combine %s.vall from:" % pdb
print "     1. fasta_dict               (from pdb_seqres.txt)"
print "     2. sectr_dict               (from dssp in rosetta idealizer)"
print "     3. idealized_pdb_xyz_dict"
print "     4. idealized_torsion_dict"
print "     5. dssp_Dict"
print "     6. seq_pro_checkpoint_dict"
print "     7. str_pro_checkpoint_dict / seq_pro_checkpoint_dict"
print "-"*75

bfactor = 0
CB_x = 0
CB_y = 0
CB_z = 0
CEN_x = 0
CEN_y = 0
CEN_z = 0

vall_out = open( pdb + ".vall", "w")
for rsn in fasta_dict.keys():
    if rsn > len( secstr_dict.keys()): continue
    if rsn > len( fasta_dict.keys()): continue
    if rsn in missing_den_list:
        continue
    else:
        vall_out.write("%5s %1s %1s %5d %6.2f " %( pdb, fasta_dict[ rsn ], secstr_dict[ rsn ], rsn, bfactor ))
        vall_out.write("%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f " %( float(idealized_pdb_xyz_dict[ rsn ][0]), float(idealized_pdb_xyz_dict[ rsn ][1]), float(idealized_pdb_xyz_dict[ rsn ][2]), CB_x, CB_y, CB_z, CEN_x, CEN_y, CEN_z ))
        vall_out.write("%8.3f %8.3f %8.3f " %( float(idealized_pdb_torsion_dict[ rsn ][0]), float(idealized_pdb_torsion_dict[ rsn ][1]), float(idealized_pdb_torsion_dict[ rsn ][2]) ))
        vall_out.write("%8.3f %8.3f %3d %3d " %( float(dssp_Dict[ rsn ][1]), float(dssp_Dict[ rsn ][2]), int(dssp_Dict[ rsn ][0]), 0.0))
        vall_out.write("%s "  %( seq_checkpoint_dict[ rsn ]))
        if not options.no_structure_profile:
            vall_out.write("%s\n" %( strpro_checkpoint_dict[ rsn ]))
        else:
            vall_out.write("%s\n"  %( seq_checkpoint_dict[ rsn ]))
vall_out.close()
print "pdb2vall.py is done making %s.vall" % pdb

