import re
from string import Template
class Bunch:
    """A simple class for bundling together arbitrary labeled data"""
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def convert_or_blank(string,convert_func):
    """sometimes fields are blank in sdf files.
    if a string cannot be converted, return an empty string"""
    try:
        return convert_func(string)
    except ValueError:
        return ""

def parse_counts(line):
    """Parse a counts line
    Returns a Bunch object with the following members:
    n_atoms - Number of atoms in the molecule
    n_bonds - Number of bonds in the molecule
    n_lists - Number of lists in the molecule
    unused_a - An unused field from characters 9:12
    chiral_flag - Chiral Flag
    n_stext - Number of STEXT entriess in the molecule
    unused_b - An unused field from characters 18:21
    unused_c - An unused field from characters 21:24
    unused_d - An unused field from characters 24:27
    unused_e - An unused field from characters 27:30
    n_prop - Number of properties in the molecule (obselete)
    version - Version flag"""
    return Bunch(
        n_atoms = convert_or_blank(line[0:3],int),
        n_bonds = convert_or_blank(line[3:6],int),
        n_lists = convert_or_blank(line[6:9],int),
        unused_a = convert_or_blank(line[9:12],int),
        chiral_flag = convert_or_blank(line[12:15],int),
        n_stext = convert_or_blank(line[15:18],int),
        unused_b = convert_or_blank(line[18:21],int),
        unused_c = convert_or_blank(line[21:24],int),
        unused_d = convert_or_blank(line[24:27],int),
        unused_e = convert_or_blank(line[27:30],int),
        n_prop = convert_or_blank(line[30:33],int),
        version = line[33:39]
    )

def counts_to_string(data):
    """turn a counts object back into a string"""
    n_atom_string = str(data.n_atoms).rjust(3," ")
    n_bond_string = str(data.n_bonds).rjust(3," ")
    n_lists_string = str(data.n_lists).rjust(3," ")
    unused_a_string = str(data.unused_a).rjust(3," ")
    chiral_flag_string = str(data.chiral_flag).rjust(3," ")
    n_stext_string = str(data.n_stext).rjust(3," ")
    unused_b_string = str(data.unused_b).rjust(3," ")
    unused_c_string = str(data.unused_c).rjust(3," ")
    unused_d_string = str(data.unused_d).rjust(3," ")
    unused_e_string = str(data.unused_e).rjust(3," ")
    n_prop_string = str(data.n_prop).rjust(3," ")
    version_string = str(data.version).rjust(6," ")
    return n_atom_string +\
          n_bond_string  +\
          n_lists_string +\
          unused_a_string +\
          chiral_flag_string +\
          n_stext_string +\
          unused_b_string +\
          unused_c_string +\
          unused_d_string +\
          unused_e_string +\
          n_prop_string +\
          version_string

def parse_atom(line,atom_index):
    """Parse an atom line, given an atom_index
    Returns a Bunch object with the following members:
    x_coord - x coordinate
    y_coord - y coordinate
    z_coord - z coordinate
    symbol - atom symbol
    mass_diff - mass difference
    charge - charge as represented in the sdf format (5 = -1).
    stereo_parity - Stereo Parity
    hydrogen_count - Hydrogen Count (Query)
    stereo_care - Stereo Care Box (Query)
    valence - Valence
    h0_designator - H0 Designator
    unused_a - An unused datafield for characters 54:57
    unused_b - An unused datafield for characters 57:60
    atom_atom_map - Atom-Atom Mapping (Reaction)
    invert_flag - Inversion/Retention flag (Reaction)
    exact_change - Exact Change Flag (Reaction)
    index - Atom Index
    """
    
    return Bunch(
        x_coord = convert_or_blank(line[0:10],float),
        y_coord = convert_or_blank(line[10:20],float),
        z_coord = convert_or_blank(line[20:30],float),
        symbol = line[31:34],
        mass_diff = convert_or_blank(line[34:36],int),
        charge = convert_or_blank(line[36:39],int),
        stereo_parity = convert_or_blank(line[39:42],int),
        hydrogen_count = convert_or_blank(line[42:45],int),
        stereo_care = convert_or_blank(line[45:48],int),
        valence = convert_or_blank(line[48:51],int),
        h0_designator = convert_or_blank(line[51:54],int),
        unused_a = convert_or_blank(line[54:57],int),
        unused_b = convert_or_blank(line[57:60],int),
        atom_atom_map = convert_or_blank(line[60:63],int),
        invert_flag = convert_or_blank(line[63:66],int),
        exact_change = convert_or_blank(line[66:69],int),
        index = convert_or_blank(atom_index,int)
    )

def atom_to_string(data):
    """turn an atom object back into a string"""
    x_coord_string = ('%(val)0.4f' % {"val" : data.x_coord } ).rjust(10," ")
    y_coord_string = ('%(val)0.4f' % {"val" : data.y_coord } ).rjust(10," ")
    z_coord_string = ('%(val)0.4f' % {"val" : data.z_coord} ).rjust(10," ")
    symbol_string = str(data.symbol).rjust(3," ")
    mass_diff_string = str(data.mass_diff).rjust(2," ")
    charge_string = str(data.charge).rjust(3," ")
    stereo_parity_string = str(data.stereo_parity).rjust(3," ")
    hydrogen_count_string = str(data.hydrogen_count).rjust(3," ")
    stereo_care_string = str(data.stereo_care).rjust(3," ")
    valence_string = str(data.valence).rjust(3," ")
    h0_designator_string = str(data.h0_designator).rjust(3," ")
    unused_a_string = str(data.unused_a).rjust(3," ")
    unused_b_string = str(data.unused_b).rjust(3," ")
    atom_atom_map_string = str(data.atom_atom_map).rjust(3," ")
    invert_flag_string = str(data.invert_flag).rjust(3," ")
    exact_change_string = str(data.exact_change).rjust(3," ")
    
    return x_coord_string+ \
        y_coord_string+ \
        z_coord_string+ \
        " "+ \
        symbol_string+ \
        mass_diff_string+ \
        charge_string+ \
        stereo_parity_string+ \
        hydrogen_count_string+ \
        stereo_care_string+ \
        valence_string+ \
        h0_designator_string+ \
        unused_a_string+ \
        unused_b_string+ \
        atom_atom_map_string+ \
        invert_flag_string+ \
        exact_change_string

def parse_bond(line):
    """Parse a bond line
    Returns a Bunch object with the following members:
    first_atom - First Atom index
    second_atom - Second Atom index
    bond_type - Bond Type
    unused_a - An unused datafield for characters 12:15
    Bond Topology - The bond topology
    Reaction Center - Reacting center topology
    """
    
    return Bunch(
        first_atom = convert_or_blank(line[0:3],int),
        second_atom = convert_or_blank(line[3:6],int),
        bond_type = convert_or_blank(line[6:9],int),
        bond_stereo = convert_or_blank(line[9:12],int),
        unused_a = convert_or_blank(line[12:15],int),
        bond_topology = convert_or_blank(line[15:18],int),
        reaction_center = convert_or_blank(line[18:21],int)
    )

def bond_to_string(data):
    first_atom_string = str(data.first_atom).rjust(3," ")
    second_atom_string = str(data.second_atom).rjust(3," ")
    bond_type_string = str(data.bond_type).rjust(3," ")
    bond_stereo_string = str(data.bond_stereo).rjust(3," ")
    unused_a_string = str(data.unused_a).rjust(3," ")
    bond_topology_string = str(data.bond_topology).rjust(3," ")
    reaction_center_string = str(data.reaction_center).rjust(3," ")
    
    return first_atom_string+ \
        second_atom_string+ \
        bond_type_string+ \
        bond_stereo_string+ \
        unused_a_string+ \
        bond_topology_string+ \
        reaction_center_string

def parse_header(line_list):
    """Parse a header block
    Returns a Bunch object with the following members:
    name - Molecule name
    program - Program data
    comments - Comment data
    """
    
    return Bunch(
        name = line_list[0],
        program = line_list[1],
        comments = line_list[2],
    )

def header_to_strings(data):
    """Parse a header block to a list of strings"""
    return [
        data.name,
        data.program,
        data.comments
    ]

def parse_data_block(line_list):
    """Parse a data block
    Returns a map with the Data element name as key and data lines as value"""
    
    data_map = {}
    inside_block = False
    data_list = []
    header_key = ""
    for line in line_list:
        if len(line) > 0 and line[0] == ">":
            header_match = re.match(">.*<(.*)>",line)
            
            #re.match returns None if nothing matches
            if header_match:
                header_key = header_match.group(1)
                inside_block = True
                continue
        
        #this is the end of a block, indidcated by a blank line
        #except, you can have an empty block indicated by two blank lines
        elif line == "" and len(data_list) != 0:
            inside_block = False
            assert(header_key != "")
            data_map[header_key] = data_list
            data_list = []
            continue
        else:
            data_list.append(line)
    return data_map

def data_block_to_strings(data):
    data_string_list = []
    for header_name in data:
        data_lines = data[header_name]
        header_string = ">  <"+header_name+">"
        data_string_list.append(header_string)
        for line in data_lines:
            data_string_list.append(line)
        data_string_list.append("")
    return data_string_list


class MolFile:
    """A class for parsing, storing and outputting molfiles"""
    def __init__(self, line_list):
        """Create a molfile given a list of lines"""
        
        #Parse Header lines
        header_lines = line_list[0:3]
        self.header_data = parse_header(header_lines)
        
        #Parse counts lines
        count_line = line_list[3]
        self.count_data = parse_counts(count_line)
        
        #identify start and end of atom and bond blocks
        atom_block_start = 4
        atom_block_end = atom_block_start + self.count_data.n_atoms
        
        bond_block_start = atom_block_end
        bond_block_end = bond_block_start + self.count_data.n_bonds
        
        #parse atom lines
        atom_lines = line_list[atom_block_start:atom_block_end]
        self.atom_data = [parse_atom(line,index) for index,line in enumerate(atom_lines)]
        
        #parse bond lines
        bond_lines = line_list[bond_block_start:bond_block_end]
        self.bond_data = [parse_bond(line) for line in bond_lines]
        
        #identify start and end of property block
        property_block_start = bond_block_end
        property_block_end = line_list.index("M  END")
        
        #for now, just store raw property lines
        self.property_lines = line_list[property_block_start:property_block_end]
        
        #parse data map
        data_item_start = property_block_end+1
        data_lines = line_list[data_item_start:]
        self.data_map = parse_data_block(data_lines)
        
    
    def data_to_string_list(self):
        """Output the stored molfile data as a list of strings"""
        string_list = []
        
        #output header
        string_list.extend(header_to_strings(self.header_data))
        #output counts
        string_list.append(counts_to_string(self.count_data))
        #output atoms
        for data in self.atom_data:
            string_list.append(atom_to_string(data))
        #output bonds
        for data in self.bond_data:
            string_list.append(bond_to_string(data))
        #output properties
        string_list.extend(self.property_lines)
        string_list.append("M  END")
        #output data
        string_list.extend(data_block_to_strings(self.data_map))
        
        return string_list
        
    def add_data_entry(self,key,value_list):
        """add a new data item to molfile"""
        self.data_map[key] = value_list
        
    def get_molecule_name(self):
        return self.header_data.name

    