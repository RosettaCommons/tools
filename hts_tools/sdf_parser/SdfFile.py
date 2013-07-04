import MolFile

class SdfFile:
    """Parse an SDF file into a set of MolFile objects and store in a dictionary by molecule name"""
    def __init__(self):
        self.molecule_map = {}
        
    def add_file(self,file,only_first_conformer=False):
        """Add an SDF file to the molecule map"""
        file_handle = open(file,'r')
        line_list = []
        for line in file_handle:
            line = line.rstrip()
            if line == "$$$$":
                mol_record = MolFile.MolFile(line_list)
                name = mol_record.get_molecule_name()
                self.molecule_map[name] = mol_record
                line_list = []
                if only_first_conformer:
                    break
                else:
                    continue
            
            line_list.append(line)
        file_handle.close()
            
    def get_record(self,name):
        """given a molecule name, return the associated record"""
        return self.molecule_map[name]
        
    def write_sdf_file(self,file):
        """Write all the records in the sdf file to the given """
        outfile = open(file,'w')
        for name in self.molecule_map:
            line_list = self.molecule_map[name].data_to_string_list()
            for line in line_list:
                outfile.write(line+"\n")
            outfile.write("$$$$\n")
                
        outfile.close()
        
def parse_sdf(filename):
    """A generator that reads in sdf files and emits a MolFile object for each entry"""
    file = open(filename,'r')
    line_list = []
    for line in file:
        line = line.rstrip()
        if line == "$$$$":
            mol_record = MolFile.MolFile(line_list)
            line_list = []
            yield mol_record
            continue
        
        line_list.append(line)
    file.close()