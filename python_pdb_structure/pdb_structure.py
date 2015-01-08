# some foundation classes for representing the contents of a PDB
# from cgkit.cgtypes import vec3
from vector3d import vector3d
import copy

from amino_acids import *


def is_integer(s):
    try:
        int(s)
        return True
    except:
        return False


class Atom:
    def __init__(self):
        self.xyz = vector3d()
        self.name = ""
        self.pdb_index = 0
        self.occupancy = 1.0
        self.bfactor = 0
        self.het = False
        self.conformer = ""
        self.ishet = False
        self.isprot = False

class Residue:
    def __init__(self):
        self.atoms = []
        self.atmap = {}  #name to atom map
        self.stripped_atmap = {}
        self.resstring = ""  # residue index + insertion code -- a string
        self.resname = ""
        self.insertion_code = ""
        self.chain = None  #pointer to the containing chain

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.atmap[atom.name] = atom
        self.stripped_atmap[atom.name.strip()] = atom

    def atom(self, atname):
        if atname in self.atmap:
            return self.atmap[atname]
        elif atname.strip() in self.stripped_atmap:
            return self.stripped_atmap[atname.strip()]
        else:
            print "Error in looking up atom", atname, "in residue", self.resname
            sys.exit(1)

    def has_atom(self, atname):
        #print "has atom? ", self.resname, atname, atname in self.atmap, atname.strip() in self.stripped_atmap
        return (atname in self.atmap) or (atname.strip() in self.stripped_atmap)

    # insert this residue into a Chain

    def claim(self, containing_chain):
        self.chain = containing_chain

    #return the chain-resstring tuple for this residue that uniquely identifies it
    #Requires that the residue has already been inserted into a chain

    def resid(self):
        assert ( self.chain )
        return self.chain.chain_name + " " + self.resstring


class Chain:
    def __init__(self):
        self.chain_name = ""
        self.residues = []
        self.resmap = {}

    def add_residue(self, residue):
        self.residues.append(residue)
        self.resmap[residue.resstring] = residue
        residue.claim(self)

    def replace_residue(self, newresidue):
        # in place replacement; keep the original location of the residue in the self.residues array
        copyres = copy.copy(newresidue)
        copyres.claim(self)
        if newresidue.resstring not in self.resmap:
            print "Could not replace residue", newresidue.resstring
            print len(self.resmap)
        assert ( newresidue.resstring in self.resmap )
        for i in xrange(len(self.residues)):
            if self.residues[i].resstring == newresidue.resstring:
                self.residues[i] = copyres
        self.resmap[newresidue.resstring] = copyres

    def residue(self, resstring):
        return self.resmap[resstring]


class PDBStructure:
    def __init__(self):
        self.chains = []
        self.chainmap = {}
        self.name = ""

    def set_name(self,name):
        self.name = name

    def residue(self, chnm, resstring):
        return self.chainmap[chnm].residue(resstring)

    def add_chain(self, chain):
        self.chains.append(chain)
        self.chainmap[chain.chain_name] = chain

    def read_from_lines(self, lines):
        last_chain = Chain()
        last_residue = Residue()
        for line in lines:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                chnm = self.chain_name_from_pdbline(line)
                if last_chain.chain_name != "" and chnm != last_chain.chain_name:
                    last_chain.add_residue(last_residue)
                    last_residue = Residue()
                    if last_chain.chain_name not in self.chainmap: self.add_chain(last_chain)
                    if chnm not in self.chainmap:
                        last_chain = Chain()
                    else:
                        # this chain has already been added, but now we have more residues
                        last_chain = self.chainmap[chnm]
                if last_chain.chain_name == "":
                    last_chain.chain_name = chnm
                resstring = self.resstring_from_pdbline(line)
                if last_residue.resname != "" and last_residue.resstring != resstring:
                    last_chain.add_residue(last_residue)
                    last_residue = Residue()
                if last_residue.resname == "":
                    last_residue.resname = self.resname_from_pdbline(line)
                    last_residue.resstring = resstring
                    #print "Read residue", last_residue.resname, last_chain.chain_name, last_residue.resstring
                atom = Atom()
                atom.xyz = self.xyz_from_pdbline(line)
                atom.name = self.atname_from_pdbline(line)
                atom.pdb_index = self.atnum_from_pdbline(line)
                atom.occupancy = self.occupancy_from_pdbline(line)
                atom.bfactor = self.bfactor_from_pdbline(line)

                if line[0:6] == "HETATM":
                    atom.ishet =  True
                else:
                    atom.ishet = False

                if line[0:4] == "ATOM":
                    atom.isprot = True
                else:
                    atom.isprot = False
                atom.conformer = line[16]
#                atom.conformer = self.conformer_from_pdbline(line)

                last_residue.add_atom(atom)
        if last_residue.resname != "":
            last_chain.add_residue(last_residue)
        if last_chain.chain_name != "" and last_chain.chain_name not in self.chainmap:
            self.add_chain(last_chain)

    def pdb_atname_range(self):
        return ( 12, 16 )

    def pdb_atnum_range(self):
        return (  6, 11 )

    def pdb_conformer_range(self):
        return 16
    def pdb_resname_range(self):
        return ( 17, 20 )

    def pdb_chain_name_range(self):
        return ( 21, 22 )

    def pdb_resstring_range(self):
        return ( 22, 27 )

    def pdb_xcoord_range(self):
        return ( 30, 38 )

    def pdb_ycoord_range(self):
        return ( 38, 46 )

    def pdb_zcoord_range(self):
        return ( 46, 54 )

    def pdb_occupancy_range(self):
        return ( 56, 60 )

    def pdb_bfactor_range(self):
        return ( 61, 66 )

    def atname_from_pdbline(self, line):
        return line[self.pdb_atname_range()[0]:self.pdb_atname_range()[1]]

    def atnum_from_pdbline(self, line):
        return line[self.pdb_atnum_range()[0]:self.pdb_atnum_range()[1]]

    def resname_from_pdbline(self, line):
        return line[self.pdb_resname_range()[0]:self.pdb_resname_range()[1]]

    def chain_name_from_pdbline(self, line):
        return line[self.pdb_chain_name_range()[0]:self.pdb_chain_name_range()[1]]

    def resstring_from_pdbline(self, line):
        return line[self.pdb_resstring_range()[0]:self.pdb_resstring_range()[1]].strip()

    def xyz_from_pdbline(self, line):
        if len(line) < 50:
            return None
        xstr = line[self.pdb_xcoord_range()[0]:self.pdb_xcoord_range()[1]]
        ystr = line[self.pdb_ycoord_range()[0]:self.pdb_ycoord_range()[1]]
        zstr = line[self.pdb_zcoord_range()[0]:self.pdb_zcoord_range()[1]]
        return vector3d(float(xstr), float(ystr), float(zstr))

    def occupancy_from_pdbline(self, line):
        return float(line[self.pdb_occupancy_range()[0]:self.pdb_occupancy_range()[1]])

    def conformer_from_pdbline(self,line):
        return str( line[self.pdb_conformer_range() ])

    def bfactor_from_pdbline(self, line):
        return float(line[self.pdb_bfactor_range()[0]:self.pdb_bfactor_range()[1]])

    def pdb_lines(self):
        lines = []
        count_atoms = 0
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    line = " " * 80
                    count_atoms += 1
                    if atom.het:
                        line = line[:0] + "HETATM" + line[6:]
                    else:
                        line = line[:0] + "ATOM" + line[4:]
                    line = line[:self.pdb_atname_range()[0]] + atom.name + line[self.pdb_atname_range()[1]:]
                    line = line[:self.pdb_atnum_range()[0]] + ( "%5d" % count_atoms ) + line[self.pdb_atnum_range()[1]:]
                    line = line[:self.pdb_resname_range()[0]] + res.resname + line[self.pdb_resname_range()[1]:]
                    line = line[:self.pdb_chain_name_range()[0]] + chain.chain_name + line[
                                                                                      self.pdb_chain_name_range()[1]:]
                    if is_integer(res.resstring):
                        line = line[:self.pdb_resstring_range()[0]] + ("%5s" % (res.resstring + " ") ) + line[
                                                                                                         self.pdb_resstring_range()[
                                                                                                             1]:]
                    else:
                        line = line[:self.pdb_resstring_range()[0]] + ("%5s" % (res.resstring    ) ) + line[
                                                                                                       self.pdb_resstring_range()[
                                                                                                           1]:]
                    line = line[:self.pdb_xcoord_range()[0]] + ( " %7.3f" % atom.xyz.x() ) + line[
                                                                                             self.pdb_xcoord_range()[
                                                                                                 1]:]
                    line = line[:self.pdb_ycoord_range()[0]] + ( " %7.3f" % atom.xyz.y() ) + line[
                                                                                             self.pdb_ycoord_range()[
                                                                                                 1]:]
                    line = line[:self.pdb_zcoord_range()[0]] + ( " %7.3f" % atom.xyz.z() ) + line[
                                                                                             self.pdb_zcoord_range()[
                                                                                                 1]:]
                    line = line[:self.pdb_occupancy_range()[0]] + ( "%4.2f" % atom.occupancy ) + line[
                                                                                                 self.pdb_occupancy_range()[
                                                                                                     1]:]
                    line = line[:self.pdb_bfactor_range()[0]] + ( "%5.2f" % atom.bfactor ) + line[
                                                                                             self.pdb_bfactor_range()[
                                                                                                 1]:]

                    lines.append(line + "\n")
            lines.append("TER\n")
        return lines

    def take_a_selfie(self):
        for chain in self.chains:
            for res in chain.residues:
                print res.resstring

    def take_a_clean_selfie(self):
        for atom in self.cleanprotein:
            print atom.name
        for atom in self.cleanhetatm:
            print atom.name

    def remove_water(self): ## need to remove the residues, not make a new list
        numwaterremoved = 0
        self.remove_these_HOHS_res =[]
        for chain in self.chains:
            for res in chain.residues:
                    print res.resstring
                    if res.resname == "HOH":
                        print "is HOH"
                        self.remove_these_HOHS_res.append(res)
                        numwaterremoved +=1

        print "Removed a total of %s waters " %numwaterremoved


    def remove_hydrogens(self):
        from re import findall
        print "Removing Hydrogens"
        self.remove_these_hydrogens_atom = []
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    #print str(atom.name)
                    regex = findall('[^A-Z][H][A-Z0-9]?[ ]?[0-9]?',atom.name)
                    if regex != []:
                        print "Found a Hydrogen, removing %s " % str(atom.name)
                        self.remove_these_hydrogens_atom.append(atom)


    def split_hetatm_protein(self):

        print "Splitting structure into hetatm and protein "
        self.cleanhetatm=[]
        self.cleanprotein=[]

        for chain in self.chains:
            for res in chain.residues:
                if res not in self.remove_these_HOHS_res:
                     for atom in res.atoms:
                         if atom not in self.remove_these_hydrogens_atom:
                             if atom not in self.remove_these_confs_atom:

                                if (atom.ishet):
                                    self.cleanhetatm.append(atom)
                                else:
                                    self.cleanprotein.append(atom)

        print "Found %s hetatms in structure " % len(self.cleanhetatm)

    def remove_other_conformations(self):
        self.remove_these_confs_atom = []
        for chain in self.chains:
            for res in chain.residues:
                for atom in res.atoms:

                    if (atom.conformer !=" "):
                        if (atom.conformer !="A"):
                            print "Removing alternate conformations atom name %s conf %s" %(atom.name, atom.conformer)
                            self.remove_these_confs_atom.append(atom)


    def get_active_site_res(self,cutoff):
        from math import sqrt
        self.activesiteproteinatoms = []

        for atom in self.cleanprotein:
            for ligand_atom in self.cleanhetatm:
                dist = sqrt(atom.xyz.distance_squared(ligand_atom.xyz))
                print "Distance is %s" %dist
                if (dist < cutoff ):
                    print "This is less than the cutoff (%s) " %cutoff
                    self.activesiteproteinatoms.append(atom)

        print "Found %s atoms in the active site with cutoff = %s" %(len(self.activesiteproteinatoms),cutoff)

def nresidues(pdb):
    count = 0
    for chain in pdb.chains:
        count += len(chain.residues)
    return count


# create a "subset" of the pdb that contains all the residues in the structure
# this subset is in the standard subset format of a list of tuples, where tuple[0]
# is the chain (string) and tuple[1] is the residue (string).
def subset_from_pdb(pdb):
    subset = []
    for ch in pdb.chains:
        for res in ch.residues:
            subset.append((ch.chain_name, res.resstring ))
    return subset


def pdbstructure_from_file(fname):
    pdb = PDBStructure()
    pdb.read_from_lines(open(fname).readlines())
    return pdb


# crude -- looks at all positions
def sequence_identity_rate_for_two_pdbstructures(pdb1, pdb2):
    tot = 0.0
    rec = 0.0
    for ch in pdb1.chainmap:
        assert ( ch in pdb2.chainmap );
        for res in pdb1.chainmap[ch].resmap:
            assert ( res in pdb2.chainmap[ch].resmap )
            tot += 1.0
            if pdb1.residue(ch, res).resname == pdb2.residue(ch, res).resname:
                rec += 1.0
    return rec / tot


# compute an axis-aligned bounding box for the given pdb structure
def xyz_limits_for_pdb(pdb):
    lower_xyz = vector3d()
    upper_xyz = vector3d()
    first = True
    count = 0
    for chain in pdb.chains:
        for res in chain.residues:
            for atom in res.atoms:
                count += 1
                if first:
                    #print "first", count
                    first = False
                    lower_xyz.copy(atom.xyz)
                    upper_xyz.copy(atom.xyz)
                else:
                    lower_xyz.min(atom.xyz)
                    upper_xyz.max(atom.xyz)
    #print "xyz from", count, "atoms"
    return lower_xyz, upper_xyz

class ActiveSiteStructure(PDBStructure):
    def __init__(self,cutoff):
        PDBStructure.__init__(self)
        self.cutoff = cutoff

    def clean(self):

        self.remove_water()
        self.remove_hydrogens()
        self.remove_other_conformations()
        self.split_hetatm_protein()

        self.get_active_site_res(self.cutoff)

def activesitestructure_from_file(fname,cutoff,name):
    activesite = ActiveSiteStructure(cutoff)
    activesite.read_from_lines(open(fname).readlines())
    activesite.set_name(name)
    activesite.clean()
    return activesite

# how to slice when removing atoms from the list
# should replace all for chain in self.chians, for res in chian.res, for atom in res.atoms
#somelist[:] = [tup for tup in somelist if determine(tup)]