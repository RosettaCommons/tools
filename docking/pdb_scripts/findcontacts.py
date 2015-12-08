#!/usr/bin/python
import loadPDB_contacts
import prottorsions
import sys
import string


if len(sys.argv) < 4:
    print 'Usage: findcontacts.py <pdbfile> <mode = C | I | CI> <radius>'
    print 'C = contacts only, I = interface only, CI = contacts and interface'
    print 'Not enough arguments detected ... exiting'
    sys.exit(1)

def dist_sq(atom1, atom2):
    return (atom2[0]-atom1[0])**2 + (atom2[1]-atom1[1])**2 + (atom2[2]-atom1[2])**2

def num(resstring):
    #For dealing with noncanonical residue numbers (e.g. '134A')
    nums = ''
    for item in resstring:
        if item >= '0' and item <= '9':
            nums = nums + item
    return int(nums)

def maxlength(restype):
    """
    input: restype
    output: maxlength
    includes the distance between CA and the farthest atom in the fully
    trans conformation plus a C-H bond (1.5A).
    """
    if restype in ['LEU', 'ILE', 'ASP', 'ASN']:
        return 4.0
    elif restype in ['ALA', 'GLY']:
        return 2.4
    elif restype in ['PRO', 'SER', 'THR', 'VAL', 'CYS']:
        return 2.8
    elif restype in ['GLU', 'GLN', 'PHE']:
        return 5.4
    elif restype == 'LYS':
        return 6.4
    elif restype == 'ARG':
        return 7.3
    elif restype in ['TYR', 'TRP']:
        return 6.8
    elif restype == 'HIS':
        return 4.8
    elif restype == 'MET':
        return 5.8
  
def makeMaxLengthList(loadedPDB):
    maxLengthList = []
    for item in loadedPDB:
        maxLengthList.append(maxlength(item.restype))
    return maxLengthList

def extractPartners(pdb):
    """
    Input:  a pdb (with docking partners split by a 'TER')
    output:  [partner1, partner2] Where partner1 and partner2
    are still in the PDB format (have not been loaded yet)
    """
    partner1 = []
    partner2 = []
    endoffirstchain = 0
    for item in pdb:
        if not endoffirstchain:
            if len(item) >= 4 and item[0:4] == 'ATOM':
                partner1.append(item)
            elif len(item) >= 3 and item[0:3] == 'TER':
                endoffirstchain = 1
        else:
            if len(item) >= 4 and item[0:4] == 'ATOM':
                partner2.append(item)
    return [partner1, partner2]

def pair_close_residues(loadedPDB1, loadedPDB2, torsionPDB1, torsionPDB2, contactfilter):
    """
    This is the time-saver function.  It compares every CA in partner 1 against
    every CA in partner 2 and saves ONLY the residue pairs with CA closer than 20A.
    Input:  The loadedPDB for partner 1 and partner2 (see module loadPDB_contacts.py), torsionPDBs
    for both partners (see prottorsions.py for details), and the filter
    Output:  A list of residue pairs with at least one contact (two atoms 4A apart)
    [res1, res2] where res1 and res2 are loadedPDB residues
    """
    len1 = len(torsionPDB1)
    len2 = len(torsionPDB2)
    maxLengthList1 = makeMaxLengthList(torsionPDB1)
    maxLengthList2 = makeMaxLengthList(torsionPDB2)
#    sys.stderr.write(str(len1*len2) + ' residue comparisons possible\n')
    paired_residues = []
    #Make len1*len2 comparisons
    for i in range(len1):
        for j in range(len2):
            filter = maxLengthList1[i] + maxLengthList2[j] + contactfilter + 1.0
            filter_sq = filter*filter
            #(Use torsionPDB; has CA available)
            ca_dist_sq = dist_sq(torsionPDB1[i].CA.coords, torsionPDB2[j].CA.coords)
            #save loadedPDB residue, which is not missing some atoms like the torsionPDBresidue
            if ca_dist_sq <= filter_sq:
                paired_residues.append([loadedPDB1[i], loadedPDB2[j]])
    return paired_residues
    
def count_contacts(residuepairlist, filter):
    """
    Input:  paired residue list (output of pair_close_residues)
    Output:  list of residue pairs and number of contacts for each pair that has
    more than one contact, i.e.
    [res1 (loadedPDBres), res2, # of contacts]
    """
    contactlist = []
#    sys.stderr.write('contact filter is ' + str(filter) + ' A\n')
    filter_sq = filter*filter
    for item in residuepairlist:
        #atomlist is the list of atoms in a loadedPDBres
        res1atomlist = item[0].atomlist
        res2atomlist = item[1].atomlist
        Natoms1 = len(res1atomlist)
        Natoms2 = len(res2atomlist)
        contactcounter = 0
        hbondcounter = 0
        #make Natoms1*Natoms2 comparisons
        for i in range(Natoms1):
            for j in range(Natoms2):
                atom_dist_sq = dist_sq(res1atomlist[i].coords, res2atomlist[j].coords)
                if atom_dist_sq < filter_sq:
                    contactcounter = contactcounter + 1
                    if [res1atomlist[i].atomtype[0], res2atomlist[j].atomtype[0]] in [['N', 'O'], ['O', 'N']]:
                        hbondcounter = hbondcounter + 1
        #Only save residue pairs that have contacts
        if contactcounter > 0:
            contactlist.append(item + [contactcounter] + [hbondcounter])
    return contactlist

def sortintlist(intlist):
    #sorts an intlist based on resid and chain
    newintlist = intlist
    sortedlist = []
    length = len(newintlist)
    for i in range(0, length):
        smallestindex = i
        for j in range(i+1, length):
            smallestchain = intlist[smallestindex][1]
            smallestresnum = num(intlist[smallestindex][2])
            currentchain = intlist[j][1]
            currentresnum = num(intlist[j][2])
            if currentchain < smallestchain or (currentchain == smallestchain and currentresnum < smallestresnum):
                smallestindex = j
        smallest = newintlist.pop(smallestindex)
        newintlist.insert(i, smallest)
    return newintlist


def getintlist(contactlist):
    """
    Input : a contactlist (output of count_contacts)
    output: a list of all UNIQUE residues in the contact list,
    sorted by partner and within partner, by chain and residue number
    """
    partner1list = [[contactlist[0][0].restype, contactlist[0][0].chain, contactlist[0][0].resid]]
    partner2list = [[contactlist[0][1].restype, contactlist[0][1].chain, contactlist[0][1].resid]]
    for i in range(0, len(contactlist)):
        partner1res = [contactlist[i][0].restype, contactlist[i][0].chain, contactlist[i][0].resid]
        partner2res = [contactlist[i][1].restype, contactlist[i][1].chain, contactlist[i][1].resid]
        #Add partner1res to partner1list if it is not already in partner1list
        if not (partner1res in partner1list):
            partner1list.append(partner1res)
        #Add partner2res to partner2list if it is not already in partner2list
        if not (partner2res in partner2list):
            partner2list.append(partner2res)
    """
    Partner1list can be returned as is, because partner 1 is the first
    dimension of iteration in function countcontacts, the source of the
    contactlist above.  Partner 2, however, is out of order because it is
    the second dimension of iteration in countcontacts.  Thus, I had to
    write sort routine sortintlist to sort partner 2.
    """
    sortpartner2list = sortintlist(partner2list)
    return [partner1list, sortpartner2list]
    
def makeContactlist(pdbfile, contactfilter, mode):
    """
    The master function
    inputs a pdb file and outputs a list of contacts.
    """
    pdb = open(pdbfile, 'r').readlines()
    if not (mode in ['C', 'I', 'CI', 'IC']):
        print 'Invalid mode', mode
        print 'Valid modes are C (contacts), I (interface), CI (contacts and interface)'
        sys.exit(1)



    length = len(pdbfile)
    pdbid = pdbfile[0:length-4]
#    pdbid = string.split(pdbfile, '.')[0]
    #Split the PDB into docking partners
    partnerlist = extractPartners(pdb)
    partner1pdb = partnerlist[0]
    partner2pdb = partnerlist[1]
    #load the PDB into PDBres class (loadPDB.py)
    partner1_loaded = loadPDB_contacts.loadProt(partner1pdb, 0)
    partner2_loaded = loadPDB_contacts.loadProt(partner2pdb, 0)
    #load the PDB into torsionRes class (prottorsions.py)
    partner1_torsionPDB = prottorsions.makeTorsionPDB(partner1_loaded, 0)
    partner2_torsionPDB = prottorsions.makeTorsionPDB(partner2_loaded, 0)
    #Filter out residues without close CAs
    residuepairlist = pair_close_residues(partner1_loaded, partner2_loaded, partner1_torsionPDB, partner2_torsionPDB, contactfilter)
#    sys.stderr.write('making ' + str(len(residuepairlist)) + ' residue comparisons\n')
    contactlist = count_contacts(residuepairlist, contactfilter)
#    sys.stderr.write(str(len(contactlist)) +  ' contacts found\n')
    #Exit if no contacts found
    if len(contactlist) == 0:
        sys.stderr.write('No contacts found; not making a contact or interface file\n')
        sys.exit(0)
    #Make a contact file if user requested it
    if mode in ['C', 'CI','IC']:
        contactfname = pdbid + '.contacts'
#        print 'Making contact file', contactfname
        contactoutput = open(contactfname, 'w')
        for item in contactlist:
            outstring = item[0].restype + ' ' + item[0].chain + ' ' + string.rjust(item[0].resid, 5) + ' -- ' + item[1].restype + ' ' + item[1].chain + ' ' + string.rjust(item[1].resid, 5) + ' ' + string.rjust(str(item[2]), 2) + ' contacts'
            if item[3] > 0:
                outstring = outstring + ' ' + str(item[3]) + ' Hbonds'
            outstring = outstring + '\n'
            contactoutput.write(outstring)
    #Make an interface file if user requested it
    if mode in ['I', 'CI', 'IC']:
        intfname = pdbid + '.int'
#        print 'Making interface file', intfname
        intoutput = open(intfname, 'w')
        intlist = getintlist(contactlist)
        intoutput.write('Partner 1:\n')
        partner1=intlist[0]
        for item in partner1:
            intoutput.write(item[0] + ' ' + item[1] + ' ' + string.rjust(item[2], 5) + '\n')
        intoutput.write('Partner 2:\n')
        partner2 = intlist[1]
        for item in partner2:
            intoutput.write(item[0] + ' ' + item[1] + ' ' + string.rjust(item[2], 5) + '\n')
      
pdbfile = sys.argv[1]
mode = string.upper(sys.argv[2])
radius=float(sys.argv[3])
makeContactlist(pdbfile, radius, mode)

