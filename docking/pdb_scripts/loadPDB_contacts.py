##############################################
#Module loadPDB.py
#created 10/13/03, Mike Daily
#last modified 3/22/04, Mike Daily
##############################################
import string
import vector

def extractCoords(PDBline):
    return [float(PDBline[30:38]), float(PDBline[38:46]), float(PDBline[46:54])] 

class PDBatom:
    #reduce a pdb atom line (columns of text) to a set of objects
    #describing that atom.
    def __init__(self, PDBline):
        recordType=PDBline[0:4]
        if recordType == 'ATOM':
            self.groupType='std'
        elif recordType == 'HETA':
            self.groupType='het'
        self.atomtype = string.strip(PDBline[13:17])
        self.restype = PDBline[17:20]
        self.chain = PDBline[21]
        self.resid = string.strip(PDBline[22:27])
        #Please note resid is a STRING not an int (works better
        #if resid is something like '14A')
        self.coords = [float(PDBline[30:38]), float(PDBline[38:46]), float(PDBline[46:54])]  

class reducedAtom:
    #reduces a pdbatom to just that atom's type and coordinates,
    #discarding all other information
    def __init__(self, pdbatom):
        self.atomtype = pdbatom.atomtype
        self.coords = pdbatom.coords

def loadAtoms(pdb, keepHetero):
    """
    takes opened pdb file 'pdb' and reduces each ATOM entry to a
    PDBatom as defined
    in the PDBatom class above.
    the output is a list of PDBatoms
    """
    atomlist = []
    for item in pdb:
        if item[0:4] == 'ATOM' or (keepHetero and item[0:6] == 'HETATM'):
            newatom = PDBatom(item)
            if newatom.restype != 'HOH':
                atomlist.append(newatom)
    return atomlist

def getResidues(pdbatoms):
    """
    takes a list of pdbatoms (output of loadAtoms above) and splits
    it into sublists, one for each residue
    The output looks like the following:
    [atom1, atom2, ... , Natoms(res1)] , [atom1, atom2, ... ,
    Natoms(res2)] , ... , Nres]
    The format of the pdbatom entries does not change, only the list
    structure.
    Each of the residue atom lists can then be reduced to a residue
    'class'
    """
    outputlist = []
    resatomlist = [pdbatoms[0]]
    for i in range (1, len(pdbatoms)):
        if pdbatoms[i].atomtype[0] != 'H':   #filter out H atoms
            if pdbatoms[i].resid == pdbatoms[i-1].resid and pdbatoms[i].chain == pdbatoms[i-1].chain:
                resatomlist.append(pdbatoms[i])
            else:
                outputlist.append(resatomlist)
                resatomlist = [pdbatoms[i]]
    outputlist.append(resatomlist)
    return outputlist

def getCAcoords(atomlist):
    for item in atomlist:
        if item.atomtype=='CA':
            return item.coords

def getCentroid(atomlist):
    coordlist=[]
    for item in atomlist:
        coordlist.append(item.coords)
    return vector.avg(coordlist)

class PDBres:
    """
    take a resatomlist (output of getResidues) and reduces it to a
    residue class, consisting of a type, chain, resid, and a list of
    atoms.  Each atom has been reduced from a pdbatom (output of class
    PDBatom) to a reduced atom (output of reducedAtom).  The pdbatom
    information is no longer needed, since it is the same for every residue.
    """
    def __init__(self, resatomlist):
        self.groupType=resatomlist[0].groupType
        self.restype = resatomlist[0].restype
        self.chain = resatomlist[0].chain
        self.resid = resatomlist[0].resid
        self.atomlist = []
        for item in resatomlist:
            self.atomlist.append(reducedAtom(item))
        if self.groupType=='std':
            self.CAcoords = getCAcoords(self.atomlist)
        elif self.groupType=='het':
            self.CAcoords = getCentroid(self.atomlist)

def reduceReslist(reslist):
    """
    reslist is the output of getResidues.  Convert each residue in
    reslist from a list of pdbatoms to a PDBres class (see class
    PDBres above).  The final output will be a list of PDBres.
    """
    PDBreslist = []
    for item in reslist:
        resclass=PDBres(item)
        if resclass.CAcoords != None:
            PDBreslist.append(resclass)
        else:
            print 'rejecting residue with incomplete backbone'
    return PDBreslist

def loadProt(pdb, keepHetero):
    """
    combines above functions
    Input is an opened pdb file (the pdb file must be opened
    separately)
    Output is a the output of function reduceReslist above.
    Function proceeds in several steps:
    1) Extract ATOM lines and reduce them to pdbatom classes (function
    loadAtoms)
    2) Parse pdbatom list (result of (1)) into residues (function
    getResidues)
    3) Reduce each residue in the residue list (result of (2)) into a
    PDBres (output of reduceReslist)
    """
    pdbatomlist = loadAtoms(pdb, keepHetero)
    resatomlist = getResidues(pdbatomlist)
    return reduceReslist(resatomlist)

##############################################################
"""
Additional functionalities
These functionalities will return raw data from a selected portion
of the pdb file, rather than pdb data in class format.  They are
most useful when the goal is to input some raw data from the pdb,
make some minor changes, and re-output the pdb or part of it.
However, some of the class functionalities are used in conditions.
"""
##############################################################

def getResList(pdblist):
    #pdblist is an opened pdb, not a pdbfile
    reslist = []
    pdbResidues = loadProt(pdblist, 0)
    for item in pdbResidues:
        reslist.append([item.restype, item.chain, item.resid])
    return reslist

def extractRes(pdb, chain, resnum):
    """
    pdb is the read PDB file, not the raw PDB file.
    extracts residue 'resnum' of chain 'chain' from pdb file 'PDBfile'
    returns atomlist, the list of atoms corresponding to that residue
    """
    atomlist=[]
    for i in range(0, len(pdb)):
        if len(pdb[i]) > 4 and pdb[i][0:4] == 'ATOM':
            pdbatom=PDBatom(pdb[i])
            if pdbatom.chain == chain and pdbatom.resid == str(resnum) and pdbatom.atomtype != 'OXT':
                atomlist.append(pdb[i])
        elif atomlist != []:
            i = len(pdb) + 5
    return atomlist

def extractHetero(pdb, heterochain, heteroname):
    #pdb is the read PDB file, not the raw PDB file
    #extracts all atoms with hetero group name "heteroname" of chain "chain" from PDB file "PDBfile"
    #returns atomlist, the list of atoms corresponding to that hetero group(s).
    atomlist=[]
    for i in range(0, len(pdb)):
        if len(pdb[i]) > 4 and pdb[i][0:6] == 'HETATM':
            hetatm=PDBatom(pdb[i])
            if hetatm.chain == heterochain and hetatm.restype == string.upper(heteroname):
                atomlist.append(pdb[i])
        elif atomlist != []:
            i = len(pdb) + 5
    return atomlist

#New subroutine added 3/22/04, Mike Daily
def getChain(loadedPDB, chain):
    #extracts chain 'chain' from a loaded PDB
    pdbChain = []
    for item in loadedPDB:
        if item.chain == chain or (item.groupType=='het' and item.chain == ' '):
            item.chain = chain
            pdbChain.append(item)
    return pdbChain
