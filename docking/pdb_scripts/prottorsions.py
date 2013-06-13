###################################
#Module prottorsions
#created 6/18/03
#last modified 7/17/03
#calculates phi,psi,omega, and chi angles for every residue in a protein
#replaces PhiPsi.py
#Mike Daily
####################################
import vector
import loadPDB
import string
import sys

class torsionRes:
    def __init__(self, PDBres):
        self.restype = PDBres.restype
        self.chain = PDBres.chain
        self.resid = PDBres.resid
        #Some variables to keep track of important atom types
        #These will be used later.
        bbatomtypelist = ['N', 'CA', 'C', 'CB']
        bbcounter = 0
        bblist = []
        scatomtypelist = [Gresname(self.restype), Dresname(self.restype), Eresname(self.restype), Zresname(self.restype)]
        sccounter = 0
        sclist = []
        """
        Explanation:
        The loading algorithm shown below saves time by (indirectly) keeping track
        of what has already been read in.  Only looks for one BB atom and one SC atom
        in each pdb line.
        """
        for item in PDBres.atomlist:
            if bbcounter < 4 and item.atomtype == bbatomtypelist[bbcounter]:
                bblist.append(item)
                bbcounter = bbcounter + 1
            if sccounter < 4 and item.atomtype == scatomtypelist[sccounter]:
                sclist.append(item)
                sccounter = sccounter + 1
        """
        if not bb_complete(bblist):
            sys.stderr.write('rejecting residue with incomplete backbone\n')
            return None
        """
        self.N = bblist[0]
        self.CA = bblist[1]
        self.C = bblist[2]
        bblength = len(bblist)
        sclength = len(sclist)
        if bblength > 3:
            self.CB = bblist[3]
            if sclength > 0:
                self.Gres = sclist[0]
                if sclength > 1:
                    self.Dres = sclist[1]
                    if sclength > 2:
                        self.Eres = sclist[2]
                        if sclength > 3:
                            self.Zres = sclist[3]

class emptyRes:
    def __init__(self, restype, chain, resid):
        self.restype = restype
        self.chain = chain
        self.resid = resid

class angleRes:
    #Converts a list of torsion angles into a more usable class format.
    def __init__(self, restype, chain, resid, omega, phi, psi, chilist, rotstring):
        self.restype = restype
        self.chain = chain
        self.resid = resid
        self.omega = omega
        self.phi = phi
        self.psi = psi
        self.chilist = chilist
        self.rotstring = rotstring

def bb_complete(bblist):
    testlist=[]
    for item in bblist:
        testlist.append(item.atomtype)
    print testlist
    return ('N' in testlist)*('CA' in testlist)*('C' in testlist)

def num(resstring):
    #For dealing with noncanonical residue numbers (e.g. '134A')
    nums = ''
    for item in resstring:
        if item >= '0' and item <= '9':
            nums = nums + item
    return int(nums)

def Gresname(restype):
    #returns an EXACT atom name to look for for CG of a residue, given the residue type
    if restype == 'CYS':
        return 'SG'
    elif restype in ['ILE', 'VAL']:
        return 'CG1'
    elif restype == 'SER':
        return 'OG'
    elif restype == 'THR':
        return 'CG2'
    return 'CG'

def Dresname(restype):
    #returns an EXACT atom name to look for for CD of a residue, given the residue type
    if restype == 'MET':
        return 'SD'
    if restype == 'HIS':
        return 'CD2'
    elif restype in ['ASP','ASN']:
        return 'OD1'
    elif restype in ['PHE', 'ILE', 'LEU', 'TRP', 'TYR']:
        return 'CD1'
    return 'CD'

def Eresname(restype):
    #returns an EXACT atom name to look for for CE of a residue, given the residue type
    if restype in ['GLU', 'GLN']:
        return 'OE1'
    elif restype == 'ARG':
        return 'NE'
    return 'CE'

def Zresname(restype):
    if restype == 'LYS':
        return 'NZ'
    return 'CZ'

def evalChi1(N, CA, CB, CG):
    #vector.torsion evaluates between -180 and 180
    rawTorsion = vector.torsion(N, CA, CB, CG)
    #Dunbrack & Cohen (1997) define chi1 between -120 and 240, so correct angles between -180 and -120
    if rawTorsion >= -180 and rawTorsion <= -120:          
        chi1 = rawTorsion + 360                             
    else:
        chi1 = rawTorsion
    #Use Dunbrack and Cohen definition of rotamers gauche+ (+), trans(T), and gauche- (-)
    if chi1 >= 0 and chi1 <= 120:                           
        return [chi1, '+']
    elif chi1 >= 120 and chi1 <= 240:
        return [chi1, 'T']
    return [chi1, '-']

def evalChi2(restype, CA, CB, CG, CD):
    """
    This function takes the raw torsion from vector.torsion and alters it in one of four ways.
    The simplest is to not modify rawTorsion (trp case).  The next simplest case is to use evalChi1
    (for the 8 chi2 residues not listed here.  The asp/asn and phe/tyr/his cases use the same basic
    principle, because carboxylates (and essentially amides too) and phenyl rings, respectively are
    planar and symmetrical about CG.  Thus, chi2 and chi2+180 are identical orientations.  The
    asp/asn and phe/tyr cases differ in the bounds used by Dunbrack, which are -90 to +90 and -30 to
    +150, respectively. 
    """
    if restype in ['ASP', 'ASN']:
        rawTorsion = vector.torsion(CA, CB, CG, CD)
        if rawTorsion <= -90:
            chi2 = rawTorsion + 180
        elif rawTorsion >= 90:
            chi2 = rawTorsion - 180
        else:
            chi2 = rawTorsion
        #Use Dunbrack and Cohen definition of rotamers gauche+ (+), trans(T), and gauche- (-) for the asp/asn chi2 case
        if chi2 >= 30 and chi2 <= 90:                                  
            return [chi2, '+']
        elif chi2 >= -30 and chi2 <= 30:
            return [chi2, 'T']
        return [chi2, '-']
    if restype in ['PHE', 'TYR', 'HIS']:
        rawTorsion = vector.torsion(CA, CB, CG, CD)
        if rawTorsion <= -30:
            chi2 = rawTorsion + 180
        elif rawTorsion >= 150:
            chi2 = rawTorsion - 180
        else:
            chi2 = rawTorsion
        #Use Dunbrack and Cohen definition of rotamers gauche+ (+), trans(T), and gauche- (-) for the phe/tyr/his chi2 case
        if chi2 >= 30 and chi2 <= 150:                                
            return [chi2, 'G']
        return [chi2, 'T']
    #Use Dunbrack and Cohen definition of rotamers gauche+ (+), trans(T), and gauche- (-) for the trp chi2 case
    if restype == 'TRP':
        chi2 = vector.torsion(CA, CB, CG, CD)
        if chi2 >= -180 and chi2 <= -60:                              
            return [chi2, '+']
        elif chi2 >= -60 and chi2 <= 60:
            return [chi2, 'T']
        return [chi2, '-']
    return evalChi1(CA, CB, CG, CD)

def evalChi3(restype, CB, CG, CD, CE):
    """
    This function does not evaluate chi3 directly but instead uses evalChi2
    or evalChi1.  Glu and gln can use the asp case of evalChi2; all other residues use evalChi1.
    """
    if restype == 'GLU' or restype == 'GLN':
        return evalChi2('ASP', CB, CG, CD, CE)
    return evalChi1(CB, CG, CD, CE)



def makeTorsionPDB(loadedPDB, cap_chains):
    """
    Converts each residue in loadedPDB (output of loadPDB.loadProt) to a torsionRes (see class above)
    Also, puts N- and C-terminal caps at the end of each chain.  This makes it easier to calculate
    torsion angles for terminal residues, since some torsion angles require atoms from residue i-1
    or i+1.
    """
    #put an N-terminal cap on the N-terminus of the first chain
    torsionPDB = []
    if cap_chains:
        torsionPDB.append(emptyRes('NTM', loadedPDB[0].chain, str(int(loadedPDB[0].resid) - 1)))
    for i in range(len(loadedPDB)):
        torsion_Res=torsionRes(loadedPDB[i])
        if torsion_Res != 'None':
            torsionPDB.append(torsion_Res)
        if cap_chains and i < len(loadedPDB)-1 and loadedPDB[i].chain != loadedPDB[i+1].chain:
            #cap C-terminus of internal chains
            torsionPDB.append(emptyRes('CTM', loadedPDB[i].chain, str(int(loadedPDB[i].resid) + 1)))
            #cap N terminus of internal chains
            torsionPDB.append(emptyRes('NTM', loadedPDB[i+1].chain, str(int(loadedPDB[i+1].resid) - 1)))
    #cap C terminus of last chain
    if cap_chains:
        torsionPDB.append(emptyRes('CTM', loadedPDB[-1].chain, str(int(loadedPDB[-1].resid) + 1)))
    return torsionPDB

def processPDB(pdb, cap_chains):
    loadedPDB = loadPDB.loadProt(pdb, 0)
    torsionPDB = makeTorsionPDB(loadedPDB, cap_chains)
    return torsionPDB

def extractResidues(torsionPDB, chain, minres, maxres):
    """
    extract residues minres-1 to maxres+1 from torsionPDB
    (calculation of phi and psi depends on residues i-1 and i+1.
    This list will then be passed into the torsion angle calculator to generate
    a torsionProfile for the residues of interest.
    """
    torsionResidues = []
    for item in torsionPDB:
        resnum = int(item.resid)
        if item.chain == chain and resnum >= minres -1 and resnum <= maxres + 1:
            torsionResidues.append(item)
    return torsionResidues


def calctorsions(torsionReslist):
    """
    Input:  a list of torsionResidues (any contiguous segment of the protein)
    Output:  A list of angleResidues (see class angleRes above)
    """
    angleReslist = []
    for i in range(1, len(torsionReslist) - 1):
        previous = torsionReslist[i-1]
        current = torsionReslist[i]
        next = torsionReslist[i+1]
        prev_resnum = num(previous.resid)
        resnum = num(current.resid)
        next_resnum = num(next.resid)
        Ngap = resnum - prev_resnum > 1
        Cgap = next_resnum - resnum > 1
        termrestypes = ['NTM', 'CTM']
        beginres = previous.restype == 'NTM'
        endres = next.restype == 'CTM'
        dummyres = current.restype in termrestypes
        #print current.restype, current.chain, current.resid
        #print Ngap, Cgap, beginres, endres
        """
        Remember:
        omega is defined by Ca(i-1), C(i-1), N, CA
        phi is defined by C(i-1), N, CA, C
        psi is defined by N, CA, C, N(i+1)
        """
        if not(dummyres):
            if beginres or Ngap:
                omega = phi = 0
            else:
                omega = vector.torsion(previous.CA.coords, previous.C.coords, current.N.coords, current.CA.coords)
                phi = vector.torsion(previous.C.coords, current.N.coords, current.CA.coords, current.C.coords)
            if endres or Cgap:
                psi = 0
            else:
                psi = vector.torsion(current.N.coords, current.CA.coords, current.C.coords, next.N.coords)
            #Now, evaluate the chi angles
            #put chi values in a list
            chilist = []
            #rotstring-for calculating rotamers.  See the evalChi functions for more details.
            rotstring = ''
            #Based on what atoms the residue contains (see class torsionRes for more info),
            #can determine which chi angles to evaluate
            #Use the most rigorous criterion - does the residue contain the appropriate atoms?
            #(see chi1res, chi2res, etc. variables below)
            currentatoms = dir(current)
            chi1res = 'Gres' in currentatoms and current.restype != 'PRO'
            #Use nested 'if' for greater efficiency
            #see functions evalChi1, evalChi2, and evalChi3 for more on how chi angles are determined.
            """
            Remember:
            Chi1 is defined by N, CA, CB, CG
            Chi2 is defined by CA, CB, CG, CD
            and so on down the sidechainfor other chi angles
            """
            if chi1res:
                chi1 = evalChi1(current.N.coords, current.CA.coords, current.CB.coords, current.Gres.coords)
                chilist.append(chi1[0])
                rotstring = rotstring + chi1[1]
                chi2res = 'Dres' in currentatoms
                if chi2res:
                    chi2 = evalChi2(current.restype, current.CA.coords, current.CB.coords, current.Gres.coords, current.Dres.coords)
                    chilist.append(chi2[0])
                    rotstring = rotstring + chi2[1]
                    chi3res = 'Eres' in currentatoms
                    if chi3res:
                        chi3 = evalChi3(current.restype, current.CB.coords, current.Gres.coords, current.Dres.coords, current.Eres.coords)
                        chilist.append(chi3[0])
                        rotstring = rotstring + chi3[1]
                        chi4res = 'Zres' in currentatoms
                        if chi4res:
                            chi4 = evalChi1(current.Gres.coords, current.Dres.coords, current.Eres.coords, current.Zres.coords)
                            chilist.append(chi4[0])
                            rotstring = rotstring + chi4[1]
            angleReslist.append(angleRes(current.restype, current.chain, current.resid, omega, phi, psi, chilist, rotstring))
    return angleReslist

def getTorsions(fname, chain, minres, maxres):
    """
    Input:  filename and selected range for which to calculate torsions
    Output: A list of angleResidues (see class angleRes above)
    """
    pdb = open(fname, 'r').readlines()
    torsionPDB = processPDB(pdb, 1)
    if chain == '+':
        return calctorsions(torsionPDB)
    elif minres < 0 and maxres < 0:
        #get the whole chain
        ressection = extractResidues(torsionPDB, chain, -5, 10000)
        return calctorsions(ressection)
    else:
        ressection = extractResidues(torsionPDB, chain, minres, maxres)
        return calctorsions(ressection)
