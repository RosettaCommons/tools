#!/usr/bin/env python

from sys import exit, stderr
from os import popen, system, path
from os.path import exists
from optparse import OptionParser
from get_missing_den_from_disss import get_fulllength_fasta
import ConfigParser

'''
the difference between this script and the previous version is introducing
sequence alignment program bl2seq to find the missing density regions.

the reference is taken from pdb_seqres.txt, and trying to match the pdb from get_pdb_new.py
'''

def numbering_back_to_pdbseqres( pdb_target_name ):
    PDB2VALL_PATH = path.abspath(path.dirname(__file__)) + "/"
    if not PDB2VALL_PATH:
            stderr.write("ERROR: you should specify the path where your packages are first.\n")
            return 0

    ## read config
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(PDB2VALL_PATH + "pdb2vall.cfg")

    ## EXTERNAL PROGRAMS
    script_get_single_chain_pdb = PDB2VALL_PATH + "pdb_scripts/get_pdb_new.py"
    script_pdb2fasta            = PDB2VALL_PATH + "pdb_scripts/pdb2fasta.py"
    bl2seq = ""
    if not exists( bl2seq ):
        bl2seq = PDB2VALL_PATH + "../../../../../src/blast/bin/bl2seq" # Robetta location
    if not exists( bl2seq ):
        bl2seq = PDB2VALL_PATH + "../blast/bin/bl2seq" # fragment_tools location
    if not exists( bl2seq ):
        bl2seq = PDB2VALL_PATH + "pdb_scripts/bl2seq" # pdb2val location
    if config.get('pdb2vall', 'bl2seq'):
        bl2seq = config.get('pdb2vall', 'bl2seq')

    print bl2seq

    ## FILES
    pdb_seqres_fn = PDB2VALL_PATH + "database/rcsb_data/derived_data/pdb_seqres.txt"

    if len( pdb_target_name.split(".")[0] ) == 5:
        pdb_id    = pdb_target_name[:4]
        pdb_chain = pdb_target_name[4]
    else:
        stderr.write("ERROR: you should use only five letter name pdb. eg: pdbid + chain.pdb \n")
        return 0


    ## MAKING TWO THINGS IN THIS WHILE LOOP:
    ## 1. fasta_out_fn based on pdb_seqres.txt ( w/o missing_density_rsn )
    ## 2. seq_from_pdbseqres_Dict based on pdb_seqres.txt
    seq_from_pdbseqres_Dict = {}
    fasta_out_fn            = open( pdb_id.lower() + pdb_chain.upper() + ".fasta", "w")

    file  = open( pdb_seqres_fn, "r")
    line  = file.readline()
    tag   = pdb_id.lower() + "_" + pdb_chain.upper()
    count = 1
    seq_from_pdbseqres = ""
    while line:
        if line.startswith(">") and tag in line:
            #print line.strip()
            header    = line
            seq_from_pdbseqres = file.readline()
            for rsd in seq_from_pdbseqres.strip():
                seq_from_pdbseqres_Dict[ count ] = rsd
                count += 1

            fasta_out_fn.write( line + seq_from_pdbseqres )
            fasta_out_fn.close()
            break
        line = file.readline()

    if not seq_from_pdbseqres:
        stderr.write("ERROR: numbering_back_to_pdbseqres(): can't find %s in pdb_seqres.txt\n\n" % ( pdb_id.lower()+pdb_chain.upper()))
        return 0

    ## RENUMBER PDB FROM get_pdb_new.py
    ## make seq_from_getpdbpy_Dict
    seq_from_getpdbpy_Dict = {}

    cmd = script_get_single_chain_pdb + " " + pdb_id.lower() + " " + pdb_chain.upper()
    system( cmd )

    pdb_fn = pdb_id.lower() + pdb_chain.upper() + ".pdb"
    count  = 1
    if exists( pdb_fn ):
        cmd = script_pdb2fasta + " " + pdb_fn + " | tee " + pdb_fn + ".fasta"
        fasta_from_pdb_lines = popen( cmd ).readlines()
        for line in fasta_from_pdb_lines:
            if line.startswith(">"):
                continue
            fasta_from_getpdbpy = line.strip()
            for rsd in fasta_from_getpdbpy:
                seq_from_getpdbpy_Dict[ count ] = rsd
                count += 1
    else:
        stderr.write("ERROR:numbering_back_to_pdbseqres(): fails to fetch pdb via get_pdb_new.py\n")
        return 0
        exit()

    ## RUN bl2seq
    cmd = bl2seq + " -i " + pdb_fn + ".fasta -j " + pdb_fn[:5] + ".fasta -F F -p blastp"
    bl2seq_results = popen( cmd ).readlines()
    if not bl2seq_results:
        stderr.write("ERROR: not able to get results from bl2seq\n")
        return 0

    tmp_fasta_pdbseqres_bl2seq_list     = []
    tmp_fasta_from_getpdbpy_bl2seq_list = []
    fasta_pdbseqres_bl2seq              = ""
    fasta_from_getpdbpy_bl2seq          = ""

    for line in bl2seq_results:
        if "Sbjct" in line:
            fasta_pdbseqres_bl2seq += line.strip().split()[2]
            tmp_fasta_pdbseqres_bl2seq_list.append( line.strip() )

        if "Query" in line and not pdb_fn in line:
            fasta_from_getpdbpy_bl2seq += line.strip().split()[2]
            tmp_fasta_from_getpdbpy_bl2seq_list.append( line.strip() )

    '''### debug
    print  "fasta_pdbseqres_bl2seq",     fasta_pdbseqres_bl2seq
    print  "fasta_from_getpdbpy_bl2seq", fasta_from_getpdbpy_bl2seq
    print tmp_fasta_pdbseqres_bl2seq_list #'''

    numbering_ref_Dict     = {}
    numbering_ref_rev_Dict = {}

    sbjct_starting_number = int( tmp_fasta_pdbseqres_bl2seq_list[0].split()[1] )
    query_starting_number = int( tmp_fasta_from_getpdbpy_bl2seq_list[0].split()[1] )
    newrsdnum       = sbjct_starting_number
    count           = 1
    missing_density_rsn_list = []
    for rsd in fasta_from_getpdbpy_bl2seq:
        if rsd == "-":
            missing_density_rsn_list.append( newrsdnum )
            newrsdnum += 1
            continue
        numbering_ref_Dict[ count ]     = newrsdnum
        numbering_ref_rev_Dict[ newrsdnum ] = count
        count += 1
        newrsdnum += 1
    #print numbering_ref_rev_Dict

    ## ADDING HEAD AND TAIL THAT THE BL2SEQ DIDN'T COVER
    '''debug
    print "sbjct_starting_numbers", sbjct_starting_number
    print "query_starting_number", query_starting_number#'''

    if newrsdnum - 1 < len( seq_from_pdbseqres.strip()):
        print ( newrsdnum - 1 , len( seq_from_pdbseqres ))
        missing_density_rsn_list += range( newrsdnum , len( seq_from_pdbseqres.strip()) + 1 )
    if sbjct_starting_number != query_starting_number:
        if ( sbjct_starting_number > 1 ) and ( query_starting_number == 1 ):
            missing_density_rsn_list += range( 1, sbjct_starting_number )
        elif ( sbjct_starting_number > 1 ) and ( query_starting_number > 1 ):
            missing_density_rsn_list += range( 1, sbjct_starting_number - query_starting_number + 1 )
            #print "something funnny happens again"
            #return 0
            #exit()

    disss_missing_den_list = get_fulllength_fasta( pdb_target_name, "missing_den_list" )
    '''debug
    print disss_missing_den_list
    print missing_density_rsn_list #'''

    for missing_rsn_from_dissstxt in disss_missing_den_list:
        if missing_rsn_from_dissstxt not in missing_density_rsn_list:
            #print  missing_rsn_from_dissstxt
            ## check b factor believable or not
            pdblines = open( pdb_fn, "r").readlines()
            for pdbline in pdblines:
                if not "TER" in pdbline:
                    resnum   = int( pdbline[22:26] )
                    atomname = pdbline[12:16].strip()
                    #print atomname
                    #print resnum
                    if resnum == numbering_ref_rev_Dict[ missing_rsn_from_dissstxt ]:
                        if atomname in ["CA", "N", "C" ]:
                            print pdbline.strip()
                            #print resnum, seq_from_getpdbpy_Dict[ resnum ], numbering_ref_rev_Dict[ missing_rsn_from_dissstxt ], seq_from_getpdbpy_Dict[ numbering_ref_rev_Dict[ missing_rsn_from_dissstxt ]]
                            bfactor = float( pdbline[60:66] )
                            #print bfactor
                            if bfactor >= 95.0:
                                missing_density_rsn_list.append( missing_rsn_from_dissstxt )
                                print "rsn %s has a bfactor over 100, adding that into missing_density_list"
                            else:
                                #print "!#$! that"
                                #stderr.write("ERROR:numbering_back_to_pdbseqres(): dis_ss.txt messed that up.\n")
                                stderr.write("WARNING:numbering_back_to_pdbseqres(): dis_ss.txt messed that up.\n")
                                stderr.write( str(missing_rsn_from_dissstxt) + " " + str(numbering_ref_rev_Dict[ missing_rsn_from_dissstxt ]) + " at %s should be included as missing density\n" % pdb_fn)

                                #return 0
                                #exit()





    ## WRITE OUT AGAIN THE FASTA_FROM_PDBSEQRES
    missing_density_rsn_list.sort()
    string = " ".join( str(x) for x in missing_density_rsn_list ) + "\n"
    fasta_out_fn              = open( pdb_id.lower() + pdb_chain.upper() + ".fasta", "w")
    fasta_from_pdbseqres = header.strip() + "  missing_density_rsn: " +  string + seq_from_pdbseqres
    fasta_out_fn.write( fasta_from_pdbseqres )
    fasta_out_fn.close()



    results_Dict = { "numbering_ref_Dict"      :numbering_ref_Dict,
                     "seq_from_pdbseqres_Dict" :seq_from_pdbseqres_Dict,
                     "seq_from_getpdbpy_Dict"  :seq_from_getpdbpy_Dict,
                     "missing_density_rsn_list":missing_density_rsn_list,
                     "fasta_lines"             :fasta_from_pdbseqres}

    return results_Dict



if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-p", dest="pdb_id", help="five letter pdb_id eg. 2oxgZ.pdb")
    parser.add_option("-r", dest="request", help="request for returning, numbering_ref_Dict / seq_from_pdbseqres_Dict / seq_from_getpdbpy_Dict / missing_density_rsn_list / fasta_lines", default='fasta_lines')
    (options,args) = parser.parse_args()

    if not options.pdb_id:
        print "you missed something, please look at info below"
        parser.print_help()
        exit()

    results = numbering_back_to_pdbseqres( options.pdb_id )
    if results:
        results = results[ options.request ]
        #print numbering_back_to_pdbseqres( options.pdb_id )
        #print Dict
        if "Dict" in options.request:
            for rsn in results.keys():
                print rsn, results[ rsn ]
        else:
            print results
