##############################################
#Module loadPDB.py
#created 10/13/03, Mike Daily
#last modified 3/22/04, Mike Daily
##############################################
import string
import coordlib
import sys
import copy
import os

#*****************************************************************************
#FUNCTIONS FOR LOADING PDB INTO PDB CLASS

output_clean_pdb = 0

def corr_int(resid):
    resid_num = ''
    for char in resid:
        if char in string.digits:
            resid_num = resid_num + char
    return int(resid_num)

class PDBatom:
    """
    reduce a pdb atom line (columns of text) to a set of objects
    describing that atom
    """
    def __init__(self, PDBline):
        recordType=PDBline[0:4]
        if recordType == 'ATOM':
            self.groupType='std'
        elif recordType == 'HETA':
            self.groupType='het'
        self.atomtype = string.strip(PDBline[13:16])
        self.restype = PDBline[17:20]
        self.chain = PDBline[21]
        #Add in a dummy chain if there is none from the pdb file
        #(Makes later data processing easier)
        if self.chain == ' ':
            self.chain = '0'
        self.resid = string.strip(PDBline[22:27])
        #Please note resid is a STRING not an int (works better
        #if resid is something like '14A')
        self.coords = [float(PDBline[30:38]), float(PDBline[38:46]),
                       float(PDBline[46:54])]

class new_PDBatom:
    #For creation of PDBatoms from coordinates
    def __init__(self, groupType, atomtype, restype, chain, resid, coords):
        self.groupType = groupType
        self.atomtype = atomtype
        self.restype = restype
        self.chain = chain
        self.resid = resid
        self.coords = coords

#============================================================================

class PDBres:
    """
    take a resatomlist (output of getResidues) and reduces it to a
    residue class, consisting of a type, chain, resid, and a list of
    atoms.  Each atom has been reduced from a pdbatom (output of class
    PDBatom) to a reduced atom (output of reducedAtom).  The pdbatom
    information is no longer needed, since it is the same for every residue.

    All the methods below are part of this class, down to the next line of
    '=' characters.
    """
    def __init__(self, resatomlist):
        self.groupType=resatomlist[0].groupType
        self.restype = resatomlist[0].restype
        try:
            self.resaa = restypedict[self.restype]
        except KeyError:
            self.resaa = 'NA'
        self.chain = resatomlist[0].chain
        self.resid = resatomlist[0].resid
        self.resnum = corr_int(self.resid)
        self.printid = self.restype + ' ' + self.chain + ' ' + self.resid
        self.initializeAtoms(resatomlist)

    def initializeAtoms(self, resatomList):
        """
        Initializes several definitions of the atoms in resatomlist:

        1) direct-callable:  calling on the name of the atom records
        coordinates of that atom
        2) list of atoms:  list of coordinates; for blind iteration over
        atom coordinates (e.g. calculating centroid)
        3) list of types:  which atom types does the residue contain
        """
        self.coordlist = []
        self.atomlist = []
        for atom in resatomList:
            #don't load in multiple copies of the same atom (e.g. CB A, CB A)
            if not(atom.atomtype in self.atomlist):
                #initializes objects for direct calling of atoms (e.g. res.N = ?)
                coords = atom.coords
                setattr(self, atom.atomtype, coords)
                #makes coordlist, atomtypelist, namedAtomlist
                self.coordlist.append(coords)
                self.atomlist.append(atom.atomtype)
    #I have made getitem equivalent to getattr so that I can easily reference
    #atoms using the [] overloading convention, e.g. res['N'] if I want to
    #loop over a list of atomtypes without using getattr(res, atomtype)

    def copy(self):
        tmp_copy = copy.copy(self)
        tmp_copy.coordlist = self.coordlist[:]
        tmp_copy.atomlist = self.atomlist[:]
        return tmp_copy
    
    def __getitem__(self, name):
        return getattr(self, name)

    def __setitem__(self, name, value):
        return setattr(self, name, value)

    global bbAtomList
    bbAtomList = ['N', 'CA', 'C', 'O']
    
    def getScAtoms(self):
        scAtomList = []
        for atomtype in self.atomlist:
            if not (atomtype in bbAtomList):
                scAtomList.append(atomtype)
        self.scatomlist = scAtomList
        
    def addCentroid(self):
        self.centroid = coordlib.vavg(self.coordlist)
        
    def fetchAtoms(self, desiredAtomList):
        """
        input:  list of desired atom names from this residue
        output:  list of coordinates in the order of the list of
        desired atom names.
        If one of the atoms is not found, an error message is printed and
        an empty list is returned.
        The program which called this function can then interpret the empty
        list,which could mean a few things:
        1) coding error:  searched for atom types inappropriate for this
        residue type
        2) incomplete sidechain or backbone
        """
        coordlist = []
        for atomname in desiredAtomList:
            try:
                #get the coordinates
                coordlist.append(self[atomname])
            #atom not found in this residue!
            except AttributeError:
                sys.stderr.write('atom ' + atomname + ' not found in \
                residue ' + self.printid + '\n' )
                return []
        return coordlist

    def fetchBB(self):
        return self.fetchAtoms(bbAtomList)

    def fetchSC(self):
        if not(hasattr(self, 'scatomlist')):
            self.getScAtoms()
        return self.fetchAtoms(self.scatomlist)

    def addscCentroid(self):
        scatoms = self.fetchSC()
        #For glycines, set the sidechain centroid to the whole residue
        #centroid
        if scatoms == []:
            self.addCentroid()
            self.sccentroid = self.centroid
        else:
            self.sccentroid = coordlib.vavg(scatoms)
#==========================================================================

restypedict = {
    'ALA':  'A',
    'CYS':  'C',
    'ASP':  'D',
    'GLU':  'E',
    'PHE':  'F',
    'GLY':  'G',
    'HIS':  'H',
    'ILE':  'I',
    'LYS':  'K',
    'LEU':  'L',
    'MET':  'M',
    'ASN':  'N',
    'PRO':  'P',
    'GLN':  'Q',
    'ARG':  'R',
    'SER':  'S',
    'THR':  'T',
    'VAL':  'V',
    'TRP':  'W',
    'TYR':  'Y',
    }

def is_bb_complete(PDBres):
    #Crashes the script if a residue in the PDB has a missing backbone
    atomlist=PDBres.atomlist
    bb_complete = (('N' in atomlist) * ('CA' in atomlist) *
                   ('C' in atomlist) * ('O' in atomlist))
    if not bb_complete:
        print 'loadPDB.py has found a missing backbone in PDB residue', \
              PDBres.printid
        print 'generating a new pdb without the missing backbone residue'
        print '(same name as original pdb)'
        print 'The old pdb has been moved to <pdb>.badbb.pdb'
    return bb_complete

def loadAtoms(pdb):
    """
    takes opened pdb file 'pdb' and reduces each ATOM entry to a
    PDBatom as defined
    in the PDBatom class above.
    the output is a list of PDBatoms
    At this step, waters and hydrogens are removed
    """
    atomlist = []
    for item in pdb:
        if item[0:4] in ['ATOM', 'HETA']:
            newatom = PDBatom(item)
            #remove waters and hydrogens!
            if newatom.restype != 'HOH' and newatom.atomtype[0] != 'H':
                atomlist.append(newatom)
    return atomlist

def parseAtomList(pdbatoms):
    """
    takes a list of pdbatoms (output of loadAtoms above) and splits
    it into sublists, one for each residue
    The output looks like the following:
    [[atom1, atom2, ... , atom-Natoms(res1)] , [atom1, atom2, ... ,
    atom-Natoms(res2)] , ... , Nres]
    The format of the pdbatom entries does not change, only the list
    structure.
    """
    outputlist = []
    resatomlist = [pdbatoms[0]]
    for i in range (1, len(pdbatoms)):
        current = pdbatoms[i]
        previous = pdbatoms[i-1]
        if [current.chain, current.resid] == [previous.chain,
                                              previous.resid]:
            resatomlist.append(current)
        else:
            outputlist.append(resatomlist)
            resatomlist = [current]
    outputlist.append(resatomlist)
    return outputlist

def convert_to_loadedPDB(resAtomList):
    """
    Input:  a list of PDBatoms parsed into residues
    Output: a loadedPDB (each residue is converted to a PDBres class instance)

    Residues with missing backbone atoms are skipped, and the global flag
    output_clean_pdb is turned on so that a new pdb will be created without the
    missing backbone
    """
    global output_clean_pdb #Set to 0 at beginning of file
    PDBresList = []
    for res in resAtomList:
        newRes = PDBres(res)
        if newRes.groupType == 'std' and not(is_bb_complete(newRes)):
            output_clean_pdb = 1
        else:
            PDBresList.append(newRes)
    return PDBresList

def loadProt(pdbfile):
    """
    combines above functions
    Input is a pdb file
    Output is a the output of function reduceReslist above.
    Function proceeds in several steps:
    1) Extract ATOM lines and reduce them to pdbatom classes (function
    loadAtoms)
    2) Parse pdbatom list (result of (1)) into residues (function
    getResidues) and convert to PDBres classes
    """
    pdb=open(pdbfile, 'r').readlines()
    pdbid = string.split(pdbfile, '.')[0]
    
    pdbatomlist = loadAtoms(pdb)
    resAtomList = parseAtomList(pdbatomlist)
    loadedPDB = convert_to_loadedPDB(resAtomList)
    addIndices(loadedPDB)
    if output_clean_pdb:
        bad_pdb_filename = pdbid + '.badbb.pdb'
        os.system('mv ' + pdbfile + ' ' + bad_pdb_filename)
        writeSinglePDB(loadedPDB, pdbfile)
    return loadedPDB

#*******************************************************************************
#FUNCTIONS FOR OUTPUTTING LOADEDPDB TO A FILE

def outputCoords(coords):
    #Converts coordinates from list [x,y,z] to string (PDB) format
    outstr = ''
    for item in coords:
        outstr = outstr + string.rjust(('%4.3f' % (item)),8)
    return outstr

def generatePDB(loadedPDB):
    """
    Output a loadedPDB in PDB format
    """
    outputList = []
    atomCounter = 1
    for res in loadedPDB:
        for atom in res.atomlist:
            recordType = 'ATOM'
            if res.groupType == 'het':
                recordType = 'HETATM'
            pdbline = string.ljust(recordType, 6) + \
                      string.rjust(str(atomCounter),5) + '  ' + \
                      string.ljust(atom, 4) + res.restype + ' ' + res.chain + \
                      string.rjust(str(res.resnum),4) + 4*' ' + \
                      outputCoords(res[atom]) + '\n'
            outputList.append(pdbline)
            atomCounter = atomCounter + 1
    return outputList

def generate_one_model(splitPDB):
    """
    Input:  a loadedPDB split by chains
    Output:  a list of PDB lines, with secondary structure at the beginning
    and TERs at the end of each chain

    Thanks to:
    DSSP: W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637
    (secondary structure determination)
    dssp2pdb (James Stroud, 2002) - converting DSSP secondary structure
    assignments to pdb format
    """
    outputPDB = []
    #generate pdb for each chain; add TER at the end of each chain
    for chain in splitPDB:
        outputPDB = outputPDB + generatePDB(chain)
        if chain[0].groupType == 'std':
            outputPDB = outputPDB + ['TER\n']

    #Generate secondary structure
    #write out temporary pdb
    writeToPDB(outputPDB, 'tmp.pdb')
    os.system('dssp tmp.pdb tmp.dssp')
    os.system('dssp2pdb tmp.dssp > tmp.pdbheader')
    pdbheader = open('tmp.pdbheader', 'r').readlines()
    os.system('rm -f tmp.pdb tmp.dssp tmp.pdbheader')

    return pdbheader + outputPDB

def writeToPDB(outputLines, filename):
    outputFile = open(filename, 'w')
    for line in outputLines:
        outputFile.write(line)

def writeSinglePDB(loadedPDB, filename):
    #Write a list of PDB lines to a file (from generatePDB)
    #split by chains
    parsedPDB = parseByChain(loadedPDB)
    #generate and write a pdb model, including TERs and secondary structure
    outputPDB = generate_one_model(parsedPDB)
    writeToPDB(outputPDB, filename)

def writeMultipleModels(loadedPDBlist, filename):
    """
    For writing multiple models of the same protein to a pdb file.

    1) Each model must be bracketed by 'MODEL' and 'ENDMDL' keywords so
    that viewer programs do not bond distinct models together
    2) Each chain must have a unique chainID so that they can be
    distinguished by viewer programs
    3) secondary structure and TERs (see generate_one_model)
    """
    outputPDB = []
    chainList = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
                 'P','Q','R','S','T','U','V','W','X','Y','Z']
    chainNum = 0
    for model in loadedPDBlist:
        #add model keyword
        outputPDB = outputPDB + ['MODEL\n']
        parsedModel = parseByChain(model)
        #re-chain the model (so each chain in the pdb will be unique)
        for chain in parsedModel:
            if chain[0].groupType == 'std':
                resetChain(chain, chainList[chainNum])
                chainNum = chainNum + 1
        outputPDB = outputPDB + generate_one_model(parsedModel)
        outputPDB = outputPDB + ['ENDMDL\n']
    writeToPDB(outputPDB, filename)


#*******************************************************************************
#Analysis functions for loadedPDBs and PDBres residues

def getCentroid(atomlist):
    coordlist=[]
    for item in atomlist:
        coordlist.append(item.coords)
    return coordlib.vavg(coordlist)

def addIndices(loadedPDB,start=0):
    """
    Adds numerical indices to residues in a loadedPDB to facilitate
    later data processing.
    """
    Nres = len(loadedPDB)
    for i in range(0, Nres):
        setattr(loadedPDB[i], 'index', i+start)

def popLigands(loadedPDB):
    """
    removes ligands from a loadedPDB object and returns them.
    Input:   a loadedPDB object (output of loadProt)
    Output:  ligands are removed from loadedPDB
             function returns the ligands

    Intended usage:  Some protein structure calculations use ligands (e.g.
    contacts), while others must ignore them (protein torsion angles).
    For the latter kind of calculation, this allows you to remove the ligands
    during the calculation but have them handy to add back in at the end to
    preserve the pdb.
    """
    #Store protein and ligand atoms
    protPDB = []
    ligandPDB = []
    for res in loadedPDB:
        if res.groupType != 'het':
            protPDB.append(res)
        else:
            ligandPDB.append(res)
            #sys.stderr.write('removing hetero group ' + res.printid + '\n')
    #reset loadedPDB (change in place) to protPDB
    loadedPDB[:] = protPDB
    #return the ligands
    return ligandPDB

def getLigands(loadedPDB):
    #return all ligand residues from loadedPDB
    ligands = []
    for res in loadedPDB:
        if res.groupType == 'het':
            ligands.append(res)
    if len(ligands) == 0:
        sys.stderr.write('warning: found no ligands in this pdb\n')
    return ligands

def extractChains(loadedPDB, chains):
    #return any chain in chains (more than one allowed) from loadedPDB
    chainPDB = []
    for res in loadedPDB:
        if res.chain in chains:
            chainPDB.append(res)
    if len(chainPDB) == 0:
        print 'error:  cannot find chain', chains, 'in this PDB'
        print 'Check your chains and try again'
        sys.exit(1)
    return chainPDB

def parseByChain(loadedPDB):
    """
    Parse a loadedPDB into chains.
    Input:  loadedPDB (output of loadProt())
    Output:  a nested list of chains corresponding to loadedPDB

    Just a standard parsing routine.
    """
    ligands = popLigands(loadedPDB)
    chainList = [[loadedPDB[0]]]
    Nres = len(loadedPDB)
    for i in range(1, Nres):
        if loadedPDB[i].chain == loadedPDB[i-1].chain:
            chainList[-1].append(loadedPDB[i])
        else:
            chainList.append([loadedPDB[i]])
    if ligands != []:
        chainList.append(ligands)
    return chainList

def resetChain(loadedPDB, newChain):
    #Resets the chain of all residues of a loadedPDB to a new chain
    #Don't do anything if the chainID is already correct
    if loadedPDB[0].chain == newChain:
        return
    for res in loadedPDB:
        res.chain = newChain

def combine_identical_chains(chainParsedPDB):
    """
    Input:  chainParsedPDB (output of parseByChain(loadedPDB))
    Method:  using Nres of each chain, detect identical chain.  For ex., if I
    had four chains with following Nres:

    [141,146,141,146]

    I would split the list as follows:

    [ [141,141], [146,146] ]

    This allows me to compare identical monomers
    """
    Nchains = len(chainParsedPDB)
    #make NresList
    NresList = []
    for chain in chainParsedPDB:
        NresList.append(len(chain))
    #Lists to keep track of matches
    matchedList = [0] * Nchains
    combinedList = []
    for i in range(Nchains):
        #skip chains already matched.
        if matchedList[i] == 1:
            continue
        imatchList = [chainParsedPDB[i]]
        for j in range(i+1,Nchains):
            if NresList[i] == NresList[j]:
                imatchList.append(chainParsedPDB[j])
                matchedList[j] = 1
        combinedList.append(imatchList)
        matchedList[i] = 1
    return combinedList

def reset_chain(loadedPDB, chain):
    for res in loadedPDB:
        res.chain = chain
        res.printid = res.restype + ' ' + res.chain + ' ' + res.resid

def extractResidues(loadedPDBchain, resnum1, resnum2):
    """
    inputs:
    loadedPDBchain - a single chain of a pdb
    resnum1, resnum2 - the range to extract

    No protection against multiple chains
    """
    residueList = []
    for res in loadedPDBchain:
        if res.resnum >= resnum1:
            residueList.append(res)
            if res.resnum == resnum2:
                return residueList
    print 'residues', str(resnum1) + '-' + str(resnum2), 'not found'
    sys.exit(1)

def getSequence(loadedPDB):
    seqStr= ''
    for res in loadedPDB:
        if res.groupType == 'het':
            continue
        seqStr = seqStr + res.resaa
    return seqStr

def addasa(loadedPDB, pdbfile):
    """
    Call naccess on a pdb and then add the information into each residue of
    the loadedPDB
    """
    pdbLigands = popLigands(loadedPDB)
    pdbid = string.split(pdbfile, '.')[0]
    os.system('naccess ' + pdbfile)
    rsafname = pdbid + '.rsa'
    rsafile = open(rsafname, 'r')
    i=0
    for line in rsafile:
        if line[0:3] == 'RES':
            res = loadedPDB[i]
            setattr(loadedPDB[i], 'asatot', float(string.strip(line[22:28])))
            setattr(loadedPDB[i], 'asasc', float(string.strip(line[35:41])))
            setattr(loadedPDB[i], 'asabb', float(string.strip(line[48:54])))
            i=i+1
    loadedPDB[:] = loadedPDB + pdbLigands
    os.system('rm -f ' + rsafname)
    os.system('rm -f ' + pdbid + '.asa')

def getCAlist(loadedPDB):
    CAlist = []
    for res in loadedPDB:
        try:
            CAlist.append(res.CA)
        except AttributeError:
            continue
    return CAlist

def calcCAcentroid(loadedPDB):
    #returns centroid of the CA atoms from a loadedPDB
    CAlist = getCAlist(loadedPDB)
    return coordlib.vavg(CAlist)
