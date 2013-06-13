#!/usr/bin/python
import loadPDB
import calccontacts
import string
import sys
import os

#Set parameters in calccontacts.py

def check_if_already_run():
    fileList = os.listdir('.')
    if '.disulfOK' in fileList:
        print 'makefixdisulf has already run in this directory'
        print 'exiting without running again'
        sys.exit(0)

def getcys(loadedPDB):
    """
    filters out the CYS residue from PDBrespdb (the output of loadPDB.loadProt)
    """
    cyslist = []
    cysindex = 0
    for res in loadedPDB:
        if res.restype == 'CYS':
            cysindex = cysindex + 1
            setattr(res, 'cysnum', cysindex)
            cyslist.append(res)
    return cyslist

def main(pdb):
    """
    To calculate disulfides:
    1) load the PDB using loadPDB module
    2) Filter out the cysteines -> cysPDB
    3) Use calccontacts.py with disulfide mode to calculate disulfides directly from
    cysPDB
    4) Translate the disulfide contact array into a fixdisulf file compatible with
    rosetta
    """
    #Don't run if I have run before in this directory.
    check_if_already_run()
    loadedPDB = loadPDB.loadProt(pdb)
    cysPDB = getcys(loadedPDB)
    Ncys = len(cysPDB)
    if Ncys == 0:
        print pdb, 'has no cysteines'
        print 'No fixdisulf file will be made'
        sys.exit(0)
    
    disulfObject = calccontacts.contactProtein('disulf')
    disulfObject.contactPDB = cysPDB
    disulfObject.find_contacts('all')
    disulfList = disulfObject.contactMap
    Ndisulf = len(disulfList)

    #Make an empty file '.disulfOK'
    #This will let the script know that makefixdisulf.py has been run
    #once and does not need to be run again.
    disulf_record = open('.disulfOK', 'w')
    disulf_record.write('makefixdisulf.py has already been run\n')
    
    if Ndisulf == 0:
        print pdb, 'has no disulfides'
        print 'No fixdisulf file will be made'
        sys.exit(0)        
    
    print pdb, 'has', Ncys, 'cysteines'
    print pdb, 'has', Ndisulf, 'disulfides:'
    disulfObject.outputMap('print', pdb)
    
    #Create the .fixdisulf file
    pdbid = string.split(pdb, '.')[0]
    disulfname = pdbid + '.fixdisulf'
    print 'making fixdisulf file', disulfname
    fixdisulf = open(disulfname, 'w')
    
    #Write out the disulfides
    for disulf in disulfList:
        outstr = string.ljust(str(disulf.res1.cysnum), 3) + \
                 string.ljust(str(disulf.res2.cysnum), 3) + '\n'
        fixdisulf.write(outstr)



if len(sys.argv) < 2:
    print 'Usage:  makefixdisulf.py <pdb>'
    sys.exit(0)

pdb = sys.argv[1]
main(pdb)
