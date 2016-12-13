#!/usr/bin/env python

from sys import exit, stderr, stdout
from os import popen, system, path
from optparse import OptionParser
import ConfigParser

PDB2VALL_PATH = path.abspath(path.dirname(__file__)) + "/"

## read config
config = ConfigParser.RawConfigParser()
config.read(PDB2VALL_PATH + "pdb2vall.cfg")


def get_fulllength_fasta( pdb_target_name, request ):
    if not PDB2VALL_PATH:
        stderr.write("ERROR: you should specify the path where your packages are first.\n")
        return 0

    stdout.write("get_fulllength_fasta(): updating database/rcsb_data/\n")

    cmd = PDB2VALL_PATH + "database/rcsb_data/update_rcsb_data.pl"

    if config.has_option('pdb2vall', 'wwpdb_host'):
      netpdbhost = config.get('pdb2vall', 'wwpdb_host')
      if netpdbhost:
		    cmd = 'ssh %s %s' % ( netpdbhost, cmd )

    # update the database first
    print cmd
    system( cmd )

    ss_dis_fn              = PDB2VALL_PATH + 'database/rcsb_data/ss_dis.txt'
    pdb_seqres_fn          = PDB2VALL_PATH + 'database/rcsb_data/derived_data/pdb_seqres.txt'

    fasta_dict                 = {}  ## from dis_ss.txt
    disorder_dict              = {}  ## from dis_ss.txt

    pdb_id           = pdb_target_name[:4].upper()
    pdb_chain        = pdb_target_name[4].upper()
    new_pdb_id_chain = pdb_id +':'+ pdb_chain

    file = open( ss_dis_fn, "r")
    line = file.readline()
    fasta_seq  = ''
    #secstr     = ''
    disorder   = ''

    while line:
        if new_pdb_id_chain+':sequence' in line:
            # SEQUENCE HEADER EG. >103l:a:SEQUENCE
            #print line.strip()
            line = file.readline()
            while line:
                if pdb_id in line:
                    break
                # SEQUENCE INFO
                #print line.strip()
                fasta_seq = fasta_seq + line.strip()
                line = file.readline()

        '''NOTE: THERE IS SOME DISCREPANCY IN BETWEEN THE NEW VERSION DENOTATION FOR SECSTR AND THE OLD VERSION.
                 I TOOK THE DSSP ROSETTA CALCULATES FOR IDEALIZED MODEL.
        if new_pdb_id_chain+':secstr' in line:
            # sequence header eg. >103L:A:sequence
            #print line.strip()
            line = file.readline()
            while line:
                if pdb_id in line:
                    break
                # sequence info
                #print line.strip()
                secstr = secstr + line[0:75]
                line = file.readline() '''

        if new_pdb_id_chain+':disorder' in line:
            #print line.strip()
            line = file.readline()
            while line:
                if ">" and 'sequence' in line:
                    line = 'die'
                    break
                #print line.strip()
                disorder = disorder + line.strip()
                line = file.readline()

        if line == 'die':
            break

        line = file.readline()

    cmd = 'grep %s_%s %s' %( pdb_id.lower(), pdb_chain, pdb_seqres_fn )
    fasta_header = popen( cmd ).readline().strip()
    length_from_pdb_seqres = fasta_header.split()[2].split(":")[1]

    #secstr = secstr.replace("\n", "")
    '''if options.debug:
        print "information from dis_ss.txt"
        print fasta_seq
        #print secstr
        print disorder'''

    ###################################################################################################
    ## GET THE MISSING RSDS
    ##
    missing_den_list = []
    for index in range(0, len( fasta_seq )):
        rsn = index + 1
        fasta_dict[rsn]    = fasta_seq[index]
        disorder_dict[rsn] = disorder[index]
        #secstr_dict[rsn]   = ss3state( secstr[index] )

        if disorder[index] == 'X':
            missing_den_list.append( rsn )
        '''if options.debug:
            print '%3s %s %s' %( rsn, fasta_seq[index], disorder[index] )'''

    '''if options.debug:
        for index in range(0, len( fasta_seq )):
            rsn = index + 1
            print "%3s %s %s" %( rsn, fasta_dict[rsn ], disorder_dict[rsn ])
            print "%3s %s %s %s" %( rsn+1, fasta_dict[rsn + 1], secstr_dict[rsn + 1], disorder_dict[rsn + 1]'''

    ## WRITE OUT A FASTA FILE WITH A HEADER CONTAINING MISSING DENSITY RESIDUES
    missing_den_line = " ".join( str(x) for x in missing_den_list )
    '''fasta_output_fn = pdb_id.lower() + pdb_chain.upper()+".fasta"
    if len( fasta_seq ) == int( length_from_pdb_seqres ):
        fasta_output = open( fasta_output_fn, "w")
        fasta_output.write( fasta_header + " missing_density_rsn: " +  missing_den_line + "\n")
        fasta_output.write( fasta_seq + "\n")
        fasta_output.close()
        print fasta_output_fn, "done!"
    else:
        print "fasta is wrong!!!"
        return 0
        exit()'''
    fasta_output = ( fasta_header + " missing_density_rsn: " +  missing_den_line + "\n") + fasta_seq + "\n"

    if   request == 'fasta_lines':
        return fasta_output

    elif request == 'fasta_dict':
        return fasta_dict

    elif request == 'missing_den_list':
        return missing_den_list

    elif request == 'disorder_dict':
        return disorder_dict

    else:
        return 0

if __name__ == '__main__':
    parser = OptionParser()

    parser.add_option("-p", dest="pdb_id", help="five letter pdb_id eg. 2oxgZ.pdb")
    parser.add_option("-d", action="store_true", dest="debug", help="print out some information")
    parser.add_option("-r", dest="request", help="request for returning, fasta_lines / fasta_dict / missing_den_list / disorder_dict", default='fasta_lines')
    (options,args) = parser.parse_args()

    if not options.pdb_id:
        print "you miss some files, please look at info below"
        parser.print_help()
        exit()

    if options.debug:
        print get_fulllength_fasta_debug( options.pdb_id, options.request )
    else:
        print get_fulllength_fasta( options.pdb_id, options.request )


