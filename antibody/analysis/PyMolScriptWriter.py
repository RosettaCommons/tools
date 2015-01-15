#!/usr/bin/python
#Author: Jared Adolf-Bryfogle

import os
import sys
import re
from collections import defaultdict

class PyMolScriptWriter:
    """
    Class to help build PyMol scripts using arbitrary lists of PDBs.
    """
    def __init__(self, outdir):
        self.base_dir = os.path.split(os.path.abspath(__file__))[0]
        self.set_outdir(outdir)

        self.reset_script()



        self.colors = []
        self.color_types = defaultdict()

        self.vis_options = ["cartoon", "spheres", "lines", "dots", "sticks", "surface", "mesh", "nonbonded"]
        self._read_colors(self.base_dir+"/"+"simple_pymol_colors.txt")

        self.pdbs = []
        self.final_names = []

    def _read_colors(self, path):
        """
        Reads PyMOl colors text file.  Loads colors.
        """
        INFILE = open(path, 'r')
        color_type = ""
        for line in INFILE:
            line = line.strip()
            if not line: continue
            if line.startswith("#"): continue

            lineSP = line.split()
            if lineSP[0] == "TYPE":
                color_type = lineSP[1].lower()
                self.color_types[color_type] = []
                continue

            self.colors.append(lineSP[0])
            self.color_types[color_type].append(lineSP[0])
        INFILE.close()
        print "Done reading PyMol color types"

    def set_outdir(self, outdir):
        self.output_dir = outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    def write_script(self, fname):
        OUTFILE = open(self.output_dir+"/"+fname, 'w')
        for line in self.script_lines:
            OUTFILE.write(line+"\n")
        OUTFILE.close()

    def save_script(self, fname):
        self.write_script(fname)

    def reset_script(self):
        self.script_lines = []

    def clear(self):
        self.reset_script()
        self.pdbs = []
        self.final_names = []


    ####################################################################################################
    ## Helpful Functions
    ###################################################################################################

    def get_color_types(self):
        return self.color_types.keys()

    def get_vis_types(self):
        return self.vis_options

    def get_colors_of_type(self, color_type):
        color_type = color_type.lower()
        return self.color_types[color_type]

    def get_final_names(self):
        """
        Get the final names PyMOL will use after loading PDBs.
        """
        return self.final_names

    def get_sele(self, chain, resid_array):
        """
        Get a a selection from an array of residue IDs and a particular chain.
        If the residue Id is a two-element tupple, then add a selection between the first and last element
        """
        def get_entry(resid):
            if type(resid) == tuple:
                if len(resid) == 2:
                    start = resid[0]
                    end = resid[1]
                    entry = "resi "+repr(start)+"-"+repr(end)
                    return entry
                else:
                    raise Exception("Tuple for PyMolScriptWriter must be length 2!")
            else:
                entry = "resi "+repr(resid)
                return entry

        if len(resid_array) == 1:
            sele = "chain "+chain+" & "+get_entry(resid_array[0])
        else:
            sele = "chain "+chain+" & ( "+get_entry(resid_array[0])
            for resi in resid_array[1:]:
                sele = sele +" | "+get_entry(resi)
            sele = sele+" )"

        return sele



    ####################################################################################################
    ## Build PyMol Script:
    ###################################################################################################

    def add_line(self, line):
        """
        Add an arbitrary line to the script
        """
        self.script_lines.append(line)

    def add_save_session(self, session_path):
        """
        Add a line to save the session to a path
        """
        self.script_lines.append("cmd.save('"+session_path+"')")

    def add_show(self, vis_type, sele=""):
        """
        Show a representation.  Optionally with a particular selection
        """
        if not vis_type in self.vis_options:
            print("Type "+vis_type+" not a known vis_option.  Options are: \n"+repr(self.vis_options))

        if not sele:
            self.script_lines.append("show "+vis_type)
        else:
            self.script_lines.append("show "+vis_type+", "+sele)

    def add_center(self, sele = None):
        if sele:
            self.add_line("center "+sele)
        else:
            self.add_line("center")

    def add_hide(self, vis_type, sele=""):
        """
        Hide a representation.  Optionally with a particular selection.
        """
        if not type in self.vis_options:
            print("Type "+vis_type+" not a known vis_option.  Options are: \n"+repr(self.vis_options))

        if not sele:
            self.script_lines.append("show "+vis_type)
        else:
            self.script_lines.append("show "+vis_type+", "+sele)

    def add_load_pdb(self, pdb_path, load_as):
        """
        Add line to load a PDB Path into PyMol
        Optionally load them as a particular name
        Will then set the final names PyMol uses to the object.
        """
        self.pdbs.append(pdb_path)
        name = os.path.basename(pdb_path)
        name = "".join(name.split(".")[0:-1])
        basenameSP = name.split('.')
        if re.search(".pdb.gz", name):
            basename = "".join(basenameSP[:len(basenameSP)-2])
        else:
            basename = "".join(basenameSP[:len(basenameSP)-1])


        if not load_as:
            self.final_names.append(name)
            self.script_lines.append("load "+pdb_path)
        else:
            self.final_names.append(load_as)
            self.script_lines.append("load "+pdb_path+", "+load_as)

    def add_load_pdbs(self, pdb_paths, load_as = ""):
        """
        Add lines to load the list of PDB paths into PyMol
        Optionally load them as a particular name
        Will then set the final names PyMol uses to the object.
        """
        i = 0
        for path in pdb_paths:
            print path
            if load_as:
                self.add_load_pdb(path, load_as[i])
            else:
                self.add_load_pdb(path)
            i+=1

    def add_align_all(self, sele1 = "", sele2="", limit_to_bb=True, pair_fit = False):
        """
        Align all to the first model
        """
        for name in self.final_names[1:]:
            self.add_align_to(name, self.final_names[0], sele1, sele2, limit_to_bb, pair_fit)


    def add_align_all_to(self, model, sele1 = "", sele2="", limit_to_bb=True, pair_fit = False):
        """
        Align all to a particular model
        """
        for name in self.final_names:
            if name !=model:
                self.add_align_to(name, model, sele1, sele2)


    def add_align_to(self, model1, model2, sele1="", sele2 = "", limit_to_bb = True, pair_fit = False):
        """
        Align one model to another, optionally specifying a selection.
        """

        align = "align "
        if pair_fit:
            align = "pair_fit "

        bb = ""
        if limit_to_bb:
            bb = " & name n+ca+c+o "

        if not sele1:
            self.script_lines.append(align+model1+bb+","+model2+bb)
        else:
            self.script_lines.append(align+model1+" & "+sele1+bb+", "+model2+" &"+sele2+bb)

    def add_group_objects(self, names, new_group_name):
        """
        Group a set of pre-loaded names to the new group.
        """
        names_str = " ".join(names)

        self.script_lines.append("group "+new_group_name+", "+names_str)

    def add_color(self, name, color):
        if not color in self.colors:

            sys.exit("Color not understood by PyMol: "+color+" See simple_pymol_colors for list of acceptable colors")

        self.script_lines.append("color "+color+", "+name)




########################################################################################################################
## Helper Functions
########################################################################################################################

def run_pymol_script(script_path, run_gui = False, delete_script = False):
    """
    Run the script of the given path.
    """

    if not os.path.exists(script_path):
        raise Exception(script_path +" does not exist...")

    if run_gui:
        os.system("pymol "+script_path)
    else:
        os.system("pymol -c "+script_path)

    if delete_script:
        os.remove(script_path)