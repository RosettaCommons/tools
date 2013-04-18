#!/usr/bin/env python


## @file   antibody_repertoire.py
## @brief  Processes a fasta file of an antibody repertoire to create a dB of structures
## @author Jeffrey Gray
## @usage  antibody_repertoire.py repertoire.fasta


#Class antibody:
#    name = ""
#    Hseq = ""
#    Lseq = ""


import os, sys


def main(args):
    '''Script to process a fasta file of an antibody repertoire and create a file hierarchy for structure prediction
Usage: antibody_repertoire.py repertoire.fasta'''

    try:
        args[1]
    except:
        print main.__doc__
        sys.exit()

    repfilename = args[1]
    if repfilename[-6:] == '.fasta':
        repname = repfilename[:-6]
    else:
        repname = repfilename
        repfilename += '.fasta'

    repertoire_file_unpack(repname,repfilename)
    #create_grafted_templates()  ## replaced by bash script


def repertoire_file_unpack(repname,repfilename):

    print 'repertoire',repname
    print 'opening repertoire file',repfilename

    # input from repertoirefile #####################
    repertoirefile = open(repfilename)

    H = {}
    L = {}
    N_complete_Abs = 0
    N_H_only = 0
    N_L_only = 0
    ATTnum = 0

    for l in repertoirefile:

        if l[0] == '>':  # new sequence
            abchainname = l[1:].split()[0]
            abname = abchainname[:-1]
            chain = abchainname[-1]
            if abname=="V": # for Daisuke's Antigen Type Test
                if chain == "L":
                    ATTnum +=1
                abname = "%s%02i" % (repname,ATTnum)
            print 'Found ab', abname, '\t chain ', chain

        else:
            seq = l  #.rstrip()
            if chain == 'H':
                if H.has_key(abname): H[abname] += seq
                else: H[abname] = seq

            elif chain == 'L':
                if L.has_key(abname): L[abname] += seq
                else: L[abname] = seq

    # output to individual directories #####################
    try:
        os.stat(repname)
    except:
        os.mkdir(repname)

    for ab,Hseq in H.iteritems():
        try:
            Lseq = L[ab]
        except:
            print 'Warning: found chain',ab,'H but not L...skipping'
            N_H_only += 1
            break
        abdir = repname + '/' + ab
        try:
            os.stat(abdir)
        except:
            os.mkdir(abdir)

        hf=open(abdir + '/' + ab + 'H.fasta','w')
        hf.write('> '+ab+'H from '+repname+'\n')
        hf.write(Hseq)
        hf.close()

        lf=open(abdir + '/' + ab + 'L.fasta','w')
        lf.write('> '+ab+'L from '+repname+'\n')
        lf.write(Lseq)
        lf.close()

        hlf=open(abdir + '/' + ab + 'HL.fasta','w')
        hlf.write('> '+ab+'H from '+repname+'\n')
        hlf.write(Hseq)
        hlf.write('> '+ab+'L from '+repname+'\n')
        hlf.write(Lseq)
        hlf.close()
        print 'Completed Ab',ab
        N_complete_Abs += 1

    for ab,Lseq in L.iteritems():
        try:
            Hseq = H[ab]
        except:
            print 'Warning: found chain',ab,'L but not H...skipping'
            N_L_only += 1

    # run complete #####################
    print '******************************'
    print '** Repertoire process complete'
    print '** Complete antibodies:',N_complete_Abs
    print '** H only:',N_H_only
    print '** L only:',N_L_only



if __name__ == "__main__": main(sys.argv)
