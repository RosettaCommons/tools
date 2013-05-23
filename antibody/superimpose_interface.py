#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   superimpose_interface.py
## @brief  part of antibody.py with PyRosetta dependency, replaces ProFit
## @author JKLeman

import os, sys, re, json, commands, shutil
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )
from rosetta import *
rosetta.init()

LH_info = {'ref_tmpl': "template.light_heavy.pdb", 
            'ref':       {'L': "template.light_heavy_onlyL.pdb", 
                          'H': "template.light_heavy_onlyH.pdb"}, 
            'mob_tmpl':  {'L': "template.threaded.FRL.pdb", 
                          'H': "template.threaded.FRH.pdb"}, 
            'mob':       {'L': "template.threaded.FRL_onlyL.pdb", 
                          'H': "template.threaded.FRH_onlyH.pdb"},
            'out_indiv': {'L': "template.superimposed.FRL.pdb", 
                          'H': "template.superimposed.FRH.pdb"}, 
            'out': "FR.pdb", 
#            contains CDR loops, this fails when superimposition before grafting
#            'interface': {'L': [34, 36, 38, 43, 44, 46, 87, 89, 98, 100], 
#                          'H': [35, 37, 39, 44, 45, 47, 91, 93, 103, 105]}
            'interface': {'L': [36, 38, 43, 44, 46, 87, 98, 100], 
                          'H': [37, 39, 44, 45, 47, 91, 93, 103, 105]}
}

def main(args):
    ''' Script for superimposing templates within the antibody.py script.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('--prefix',
      action="store", default='output/',
      help="Prefix for output files. Should be dir name. Default is ./ string.",
    )

    #parse options
    (options, args) = parser.parse_args(args=args[1:])
    prefix = options.prefix

    for chain in ['L', 'H']:

      #get reference pose
      command = "grep ATOM " + prefix + "/" + LH_info['ref_tmpl'] + " | awk '{if($5 == \"" + chain + "\") print}' > " + prefix + LH_info['ref'][chain]
      status, output = commands.getstatusoutput(command)
      ref_pose = pose_from_pdb(prefix + LH_info['ref'][chain])

      #get mobile pose
      command = "grep ATOM " + prefix + "/" + LH_info['mob_tmpl'][chain] + " | awk '{if($5 == \"" + chain + "\") print}' > " + prefix + LH_info['mob'][chain]
      status, output = commands.getstatusoutput(command)
      mob_pose = pose_from_pdb(prefix + LH_info['mob'][chain])

      #define atommap
      atom_map = core.id.AtomID_Map_AtomID()
      print mob_pose.total_residue(), ref_pose.total_residue()

      #initialize atommap with longer pose
      if (mob_pose.total_residue() > ref_pose.total_residue()):
        core.pose.initialize_atomid_map(atom_map, mob_pose)
      else:
        core.pose.initialize_atomid_map(atom_map, ref_pose)

      #go through interface residues, get residue IDs
      for residue in LH_info['interface'][chain]:
        ref_resID = ref_pose.pdb_info().pdb2pose(chain, residue)
        mob_resID = mob_pose.pdb_info().pdb2pose(chain, residue)
        print chain, residue, ref_resID, mob_resID

        #go through atoms, get atomIDs
        for atom in ['N', 'CA', 'C', 'O']:
          mob_atomID = core.id.AtomID(mob_pose.residue(mob_resID).atom_index(atom), mob_resID)
          ref_atomID = core.id.AtomID(ref_pose.residue(ref_resID).atom_index(atom), ref_resID)

          #sets key, value pairs in atommap using atomIDs
          atom_map.set(mob_atomID, ref_atomID)

      #for checking, is atom_map populated?
      for i in range(1, len(mob_pose.chain_sequence(1))+1):
        for j in ['N', 'CA', 'C', 'O']:
          atomID = core.id.AtomID(mob_pose.residue(i).atom_index(j), i)
          print i, j, atom_map.get(atomID)

      #superimpose mobile pose onto ref pose using atommap
      core.scoring.superimpose_pose(mob_pose, ref_pose, atom_map)
      out = mob_pose.dump_pdb(prefix + '/' + LH_info['out_indiv'][chain])    

    #combine L and H chains into FR.pdb
    command = "cat " + prefix + "/" + LH_info['out_indiv']['L'] + " " + prefix +  "/" + LH_info['out_indiv']['H'] + " > " + prefix + "/" + LH_info['out']
    os.system(command)


if __name__ == "__main__": main(sys.argv)

