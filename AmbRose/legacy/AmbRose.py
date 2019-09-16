# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.

"""File:  AmbRose.py

Brief: This module enables pose to amber-trajetory interconversion and performs minimization of PDBs using the Amber FF or Amber minimization engine (sander).

Details: Requires pytraj, sander, and tLeap, which are available in AMBER (version >= 16): ambermd.org. 

Author:  Kristin Blacklock, Hai Nguyen
Blacklock <kristin.blacklock@rutgers.edu>

"""

import os
import sys
import time
import numpy
from glob import glob

## For MPI
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank, n_cores = comm.rank, comm.size

## For Amber energies and methods
import pytraj as pt
from pytraj.trajectory import Trajectory
import sander

## For per-residue energy decomposition
from MMPBSA_mods import API as MMPBSA_API

## For minimization using Amber FF (not sander engine).
from scipy.optimize import minimize

## Rosetta methods
from rosetta import *
init()

#######################################################################
##### Pose to Traj Methods ############################################

def initial_pose_to_traj( pose ):
    '''
    Generates a trajectory and coordinate map from a pose.
    '''

    ## Use tLeap to make topology and coordinates files. 
    parm_file, rst7_file = pose_to_parm( pose )

    ## Save the parm_file in the pose object.
    if pose.pdb_info():
        pose.pdb_info().modeltag( parm_file )
    else:
        new_pose_pdbinfo = rosetta.core.pose.PDBInfo()
        pose.pdb_info( new_pose_pdbinfo )
        pose.pdb_info().modeltag( parm_file )


    ## Parse the parm atom list and rst7 coordinate list by residue.
    coordinates_by_residue, atoms_by_residue = parse_by_residue( pose, 
                                                                 parm_file, 
                                                                 rst7_file
                                                                 )

    ## Use the parsed atom and coordinate lists to generate a 
    ## pose-to-parm map.
    pose_ambercoordinate_array, \
    pose_ambercoordinate_map = map_pose_to_parm( pose, 
                                                 coordinates_by_residue, 
                                                 atoms_by_residue 
                                                 )
    
    ## Use pytraj to load the coordinates and parm file into an
    ## Amber trajectory.
    pose_amber_trajectory = pt.Trajectory( xyz=pose_ambercoordinate_array, 
                                           top=parm_file 
                                           )

    return pose_amber_trajectory, pose_ambercoordinate_map


def subsequent_pose_to_traj( pose, coord_map ):
    ''' 
    Generates a trajectory using a previously generated coordinate map.
    '''

    ## Use the previously generated Amber-Pose coordinate map to get 
    ## the pose's coordinates in parm order.
    pose_ambercoodinate_arry = get_pose_coordinates_using_ambermap( pose, 
                                                                    coord_map 
                                                                    )
    ## Get the parm filename, stored in the pdb_info().modeltag().
    parm_file = pose.pdb_inf().modeltag()

    ## Use pytraj to load the coordinates and parm file into an
    ## Amber trajectory.
    pose_amber_trajectory = pt.Trajectory( xyz=pose_ambercoordinate_array, 
                                           top=parm_file )

    return pose_amber_trajectory


def pose_to_parm( pose ):
    ''' 
    Generates parm7 and rst7 files from a pose
    '''

    ## Name output files according to date and time (to keep from 
    ## overwriting).
    date_time = time.strftime("%y%m%d.%H%M%S")
    pose_name = "pose.{DT}".format(DT=date_time)

    ## Dump pose as PDB.
    pose.dump_pdb('{PNAME}.pdb'.format(PNAME=pose_name))

    ## Remove all hydrogens.
    os.system("sed -i '/     H  /d' {PNAME}.pdb".format(PNAME=pose_name))

    ## Write temporary tleap.in file, saving PDB coordinates as parm7 
    ## and rst7 files.
    parm_file = pose_name + ".parm7"
    rst7_file = pose_name + ".rst7"
    write_tleapfile( parm_file, rst7_file, pose_name )

    ## Run tLeap using the temporary tleap.in file.
    os.system('tleap -f {PNAME}.tleap.in'.format(PNAME=pose_name))

    ## Clean up unnecessary files.
    os.system('rm {PNAME}.tleap.in'.format(PNAME=pose_name))
    os.system('rm leap.log')

    return parm_file, rst7_file


def write_tleapfile( parm_file, rst7_file, pose_name ):
    '''
    Writes input file for tLeap program.
    The force field and radii used here are best for assessing the
        GB energy of proteins.
    '''
    
    ## Open a writable .tleap.in file
    with open('{PNAME}.tleap.in'.format( PNAME=pose_name ),'w') as tfile:

        ## Load the ff14SBonlysc force field.
        tfile.write("source leaprc.protein.ff14SBonlysc\n" )

        ## Load the dumped PDB.
        tfile.write("m = loadpdb {PNAME}.pdb\n".format( PNAME=pose_name ))

        ## Set the GB radii set to mbondi3, which is best for proteins.
        tfile.write("set default pbradii mbondi3\n" )

        ## Save the molecule as a parm file and coordinates file.
        tfile.write("saveamberparm m {parm} {rst}\n".format( parm=parm_file, rst=rst7_file ))
        tfile.write("quit" )

def parse_by_residue( pose, parm_file, rst7_file ):
    ''' 
    Separates the atoms in the parm file and the coordinates in 
    the rst7 file into residues.
    '''
    
    parm_atom_list = parse_parm( parm_file )
    coordinate_list = parse_rst7( rst7_file )

    ## TODO: assert coordinate_list and parm_atom_list have same number
    ##of items.

    ## Instantiate empty lists to hold residues lists.
    coords_by_residue = []
    atoms_by_residue = []

    ## Instantiate empty residue list.
    residue_atoms = []
    residue_coords = []

    ## For each atom in the atom list (from the parm_file):
    for atom_index in range( len( parm_atom_list ) ):

        ## If the first atom of a residue is reached:
        if parm_atom_list[ atom_index ] == 'N':

            ## Append this residue list to coordinate list and the 
            ## atom list.
            coords_by_residue.append( residue_coords )
            atoms_by_residue.append( residue_atoms )
            
            ## Instantiate new lists for this residue.
            residue_coords = []
            residue_atoms = []
        
        ## Add this atom's coordinates to the residue list.
        residue_coords.append( coordinate_list[ atom_index ] )

        ## Add this atom's identifier to the atoms list.
        residue_atoms.append( parm_atom_list[ atom_index ] )

    ## Once all the atoms have been gone through, append this residue's
    ## coordinate list and this residue's atom list to the bigger lists.
    coords_by_residue.append( residue_coords )
    atoms_by_residue.append( residue_atoms )
    
    ## Remove any instances of empty lists from the bigger lists.
    coords_by_residue[:] = [x for x in coords_by_residue if x != []]
    atoms_by_residue[:] = [x for x in atoms_by_residue if x != []]

    ## coords_by_residue is of the form: [ res1_numpy_array, 
    ##                                     res2_numpy_array,
    ##                                     ...,
    ##                                     resN_numpy_array ]

    ## atoms_by_residue is of the formt: [ [res1_atomlist], 
    ##                                     [res2_atomlist],
    ##                                     ...,
    ##                                     [resN_atomlist] ] 

    return coords_by_residue, atoms_by_residue


def parse_rst7( rst7_filename ):
    '''
    Parses the rst7 file into  a numpy array of coordinates.
    '''

    ## Read contents of rst7_file
    with open( rst7_filename ) as rst7_file:
        rst7_lines = rst7_file.readlines()

    ## Instantiate a new numpy array of 1 row and 3 columns.
    coordinates = numpy.zeros((1,3))

    ## For each line in the file:
    ## (except for the first two, which are "Name" and "No. of atoms" lines.
    for rst7_line in rst7_lines[2:]:

        ## Split the line into individual coordinates.
        coords_array = rst7_line.split()

        ## The first 3 coordinates correspond to one atom. Add these
        ## the numpy array of coordinates.
        coordinates = numpy.vstack([ coordinates, 
                                     numpy.array([ float(coords_array[0]), 
                                                   float(coords_array[1]), 
                                                   float(coords_array[2]) 
                                                   ]) 
                                     ])
        
        ## If there is a second set of coordinates, add them to the
        ## coordinates array as well.
        if len(coords_array) > 3:

            coordinates = numpy.vstack([ coordinates, 
                                         numpy.array([ float(coords_array[3]), 
                                                       float(coords_array[4]), 
                                                       float(coords_array[5]) 
                                                       ])
                                         ])

    ## Delete the first row in the coordinates array (it's zeroes).
    coordinates = numpy.delete( coordinates, (0), axis=0 )

    ## coordinates is an Nx3 numpy array, where N = total residues.
    return coordinates


def parse_parm( parm_file_name ):
    '''
    Parsed the parm file into a list of atom names.
    '''

    ## Read the contents of the parm_file.
    with open( parm_file_name ) as parm_file:
        parm_lines = parm_file.readlines()
    
    ## Instantiate an empty list to hold the list of atoms.
    atoms = []

    #### Parse only the lines between the "%FLAG ATOM_NAME" line and the 
    #### "%FLAG CHARGE" line, without parsing the "%FORMAT" line after 
    #### "%FLAG ATOM_NAME" line.
 
    ## Start with parsing off.
    parsing = False

    ## For each line in the parm_file:
    for line in parm_lines:

        ## If the end of the atoms block is reached,
        if "%FLAG CHARGE" in line:

            ## Turn parsing off.
            parsing = False
        
        ## If we are still parsing,
        if parsing:

            if "%FORMAT" not in line:
                
                ## Split the line into 4-character blocks (atom) and 
                ## append each to the atom list.
                for atom in [line[i:i+4] for i in range(0, len(line), 4)]:
                    atoms.append(atom)

        ## If the start of the block is reached,
        if "%FLAG ATOM_NAME" in line:
            
            ## Turn parsing on
            parsing = True
    
    ## Once all lines have been gone through, remove any  newline 
    ## characters from the list of atoms.
    atoms[:] = [x for x in atoms if '\n' not in x]

    ## Remove extra whitespaces, too.
    atoms[:] = [x.strip(' ') for x in atoms]

    ## atoms is a list of strings.
    return atoms


def map_pose_to_parm( pose, coords_by_residue, atoms_by_residue ):
    '''
    Maps Amber atoms to Pose atoms via atom-distances.
    '''

    ## Instantiate a numpy array of 1 row and 3 columns.
    pose_coords_in_parm_order = numpy.zeros((1,3))

    ## Instantiate a list that will hold residue-dictionaries that map
    ## Amber atoms to pose atoms.
    amber_to_pose_residue_dicts = []

    ## For each residue in the atoms_by_residue list,
    for amber_residue_index in range( len( atoms_by_residue ) ):

        ## Instantiate a list that will hold pose atoms that have
        ## already been mapped to an Amber atom.
        pose_atoms_used = []
        
        ## Instantiate a dictionary for this residue that will hold 
        ## the (Amber_atom : pose_atom) map.
        residue_dict = {}

        ## For each atom in this residue,
        for amber_atom_index in range( len( atoms_by_residue[ 
                                                amber_residue_index ] ) ):

            ## Make ridiculously large number as the starting minimum 
            ## for the Amber_atom-pose_atom distances.
            minimum = 100

            ## For each pose atom in this pose residue,
            pose_residue = pose.residue( amber_residue_index + 1)
            for pose_atom_index in range(1, pose_residue.natoms() + 1):

                ## If we haven't already mapped this atom to an Amber atom,
                if pose_atom_index not in pose_atoms_used:

                    ## If this atom isn't a virtual atom 
                    ## (Amber doesn't use those),
                    if not pose_residue.atom_type(pose_atom_index).is_virtual():


                        ## Gather the pose's atom coordinates in a numpy array.
                        pose_atom_coords = numpy.array([ pose_residue.atom( pose_atom_index ).xyz().x,
                                                         pose_residue.atom( pose_atom_index ).xyz().y, 
                                                         pose_residue.atom( pose_atom_index ).xyz().z ])
                        
                        ## Calculate the distance between this atom and
                        ## the Amber atom.
                        distance = numpy.linalg.norm( coords_by_residue[ amber_residue_index ][ amber_atom_index ]
                                                      - pose_atom_coords )
                        if distance < minimum:
                            minimum = distance
                            closest_pose_atom = pose_atom_index
                            closest_pose_atom_coords = pose_atom_coords
                
            ## Append the final closest atom to the "Used" list.
            pose_atoms_used.append( closest_pose_atom )

            ## Update the list of coordinates with the closest atom's 
            ## coordinates.
            pose_coords_in_parm_order = numpy.vstack([ pose_coords_in_parm_order, 
                                                       closest_pose_atom_coords ])

            ## Add this (Amber-atom : Pose-atom) mapping to the 
            ## residue dictionary.
            residue_dict[ amber_atom_index ] = closest_pose_atom

        ## After all of the Amber-residue's atoms have been mapped,
        ## append this residue dictionary to the list.
        amber_to_pose_residue_dicts.append( residue_dict )

    ## Delete the first row in the coordinates list(it's zeroes).
    pose_coords_in_parm_order = numpy.delete( pose_coords_in_parm_order, (0), axis=0 )

    ## pytraj wants a (No. Frames) x (No. Atoms) x 3 array, so make
    ## set of coordinates into a pytraj set of 1 frame.
    pose_coords_in_parm_order = numpy.array( [pose_coords_in_parm_order] )

    ## pose_coords_in_parm_order is of the form:
    ##                      frame  [
    ##                        atoms  [ 
    ##                          coords  [pose_coords_amberatom1],
    ##                                  [pose_coords_amberatom2],
    ##                                  ...,
    ##                                  [pose_coords_amberatomZ]
    ##                               ]
    ##                             ]

    ## amber_to_pose_residue_dicts is of the form: 
    ##     [ {res1_dict}, {res2_dict}, ..., {resN_dict} ]

    return pose_coords_in_parm_order, amber_to_pose_residue_dicts


def get_pose_coordinates_using_ambermap( pose, coord_map ):
    '''
    Gathers pose coordinates in parm file order using a previously
    generated coordinate map.
    '''

    ## Instantiate a numpy array of 1 row and 3 columns.
    pose_coords_in_parm_order = numpy.zeros((1,3))

    ## For each residue in the coordinate map,
    for amber_residue_index in range(len( coord_map ) ):
        pose_residue = pose.residue( amber_residue_index+1 )

        ## Get and order the dictionary keys (Amber atoms) as a list.
        amber_atom_keys = coord_map[ amber_residue_index ].keys()
        amber_atom_keys.sort()

        ## For each atom in the list of atoms,
        for amber_atom_index in amber_atom_keys:

            ## Use the dictionary to find the pose-atom index.
            pose_atom_index = coord_map[ amber_residue_index ][ amber_atom_index ]

            ## Gather the pose-atom's coordinates in a numpy array.
            pose_atom_coords = numpy.array([ pose_residue.atom( pose_atom_index ).xyz().x, 
                                             pose_residue.atom( pose_atom_index ).xyz().y, 
                                             pose_residue.atom( pose_atom_index ).xyz().z 
                                             ])

            ## Add this atom's coordinates to the list of coordinates.
            pose_coords_in_parm_order = numpy.vstack([ pose_coords_in_parm_order, 
                                                       pose_atom_coords ])

    ## Delete the first row of the coordinate list (it's zeroes).
    pose_coords_in_parm_order = numpy.delete( pose_coords_in_parm_order, (0), axis=0 )

    ## Make a 1-frame numpy array from the coordinate list.
    pose_coords_in_parm_order = numpy.array([ pose_coords_in_parm_order ])

    return pose_coords_in_parm_order


#######################################################################
##### PDB to Parm/RST7 Methods #########################################
def convert_pdbs_to_rst7_parm7_files( pdblist, same_parm ):
    
    tleap_template = '''logFile log.{pdbfile_root}
source leaprc.ff14SBonlysc
m = loadpdb NoH_{pdbfile_root}.pdb
set default pbradii mbondi3
saveamberparm m {parmfile_name}.parm7 NoH_{pdbfile_root}.rst7
quit
'''

    rst7_list = []
    parmfile_list = []

    for pdb in pdblist:
        if '/' in pdb:
            pdbfile_dir = '/'.join(pdb.split('/')[:-1])
            pdbfile = pdb.split('/')[-1]
        else:
            pdbfile_dir = os.getcwd()
            pdbfile = pdb
        
        pdbfile_root = pdbfile.replace('.pdb','')
        pdbfile_root = pdbfile_root.replace('\n','')

        os.system("sed '/     H  /d' {PDB_file} > NoH_{PDB}".format( PDB_file=pdbfile_dir + '/' + pdbfile_root+'.pdb', PDB=pdbfile_root+'.pdb' ))
	
        if same_parm == 1:
            parmfile = 'parm'
        else:
	    parmfile = pdbfile_root

        tleap_command = tleap_template.format(pdbfile_root=pdbfile_root, parmfile_name=parmfile)
        
        leap_in = 'tleap.{}.in'.format(pdbfile_root)
        with open(leap_in, 'w') as leap_in_file:
            leap_in_file.write(tleap_command)
    
        os.system('tleap -f {}'.format(leap_in))
        os.remove(leap_in)
        os.remove('log.{}'.format(pdbfile_root))

        rst7_list.append('NoH_{}.rst7'.format(pdbfile_root))
        if same_parm and len(parmfile_list) > 0:
            continue;
        else:
            parmfile_list.append( '{}.parm7'.format(parmfile) )
    
    return rst7_list, parmfile_list	

#######################################################################
##### Traj to Pose Methods ############################################

def traj_to_pose_version1( pose, amber_traj, coord_map, frame=0 ):
    '''
    Converts a pose to Amber coordinates using a previously generated
    coordinate map.
    '''
    
    ## Instantiate an empty vector1 of AtomIDs to hold list of pose atom
    ## ids in parm order.
    atom_ids = utility.vector1_AtomID()

    ## Instantiate an empty vector1 of xyzVectors of doubles to hold the
    ## amber trajectory atom coordinates.
    coordinates = utility.vector1_xyzVector_double()

    ## Starting atom index, 
    amber_atom_index = -1

    ## For each residue in the coordinate map,
    for amber_residue_index in range( len( coord_map ) ):

        ## For each Amber atom in this residue's list of atoms,
        for amber_residue_atom in range( len( coord_map[ amber_residue_index ] ) ):

            ## Increment the amber atom index by 1
            amber_atom_index += 1

            ## Get the pose's atom ID for this Amber atom according 
            ## to this residue's coordinate map.
            atom_id = core.id.AtomID( coord_map[ amber_residue_index ][ amber_residue_atom ], 
                                      amber_residue_index+1 )
            
            ## Gather the Amber atom's coordinates in an xyzVector.
            point = numeric.xyzVector_double( amber_traj[ frame ][ amber_atom_index ][0], 
                                              amber_traj[ frame ][ amber_atom_index ][1], 
                                              amber_traj[ frame ][ amber_atom_index ][2] )
            
            ## Append the AtomID to the list of AtomIDs.
            atom_ids.append( atom_id )

            ## Append the coordinates to the list of coordinates.
            coordinates.append( point )

    ## Set the pose's atom coordinates to the Amber atom coordinates.
    pose.batch_set_xyz( atom_ids, coordinates )


def get_pose_coords( pose ):
    '''
    Gathers the pose coordinates into a 1-frame numpy array.
    '''

    ## Instantiate a numpy array of 1 row and 3 columns.
    coordinates = numpy.zeros((1,3))
    
    ## For each residue in the pose,
    for residue_index in range(1, pose.total_residue()+1 ):
        pose_residue = pose.residue( residue_index )

        ## For each atom in the residue,
        for atom_index in range(1, pose.residue( residue_index ).natoms()):

            ## Make a 1x3 numpy array of the atom coordinates
            pose_atom_coords = numpy.array([ pose_residue.atom( atom_index ).xyz().x, 
                                             pose_residue.atom( atom_index ).xyz().y, 
                                             pose_residue.atom( atom_index ).xyz().z 
                                             ])

            coordinates = numpy.vstack([ coordinates, 
                                          pose_atom_coords ])

    ## Delete the first row of the coordinates array (it's zeroes).
    coordinates = numpy.delete( coordinates, (0), axis=0 )
    
    ## Make a 1-frame numpy array from the coordinate list.
    coordinates = numpy.array( [coordinates] )

    return coordinates


#######################################################################
##### Energy-Getters ##################################################

def get_energies_MPI( trajectory, structure_names, igb ):
    '''
    Calculates Amber energy from trajectory using MPI.
        igb = 8 for implicit solvent
    Note: Does not work when trajectory generated with pt.iterload
    '''

    energy_data = pt.pmap_mpi(pt.energy_decomposition, trajectory, igb )

    ## Add the structure name to the energy dictionary to make
    ## printing the energies easier.
    energy_data['description'] = structure_names

    return energy_data;

def get_energies( trajectory, structure_names ):
    '''
    Calculates Amber energy from trajectory.
    '''

    energy_data = pt.energy_decomposition( trajectory, igb=8 )

    ## Add the structure name to the energy dictionary to make
    ## printing the energies easier.
    energy_data['description'] = structure_names

    return energy_data;


def get_energy_term( traj, term ):
    '''
    Calculates Amber energy from trajectory and returns a specific term
    only.
    '''

    energy_data = pt.energy_decomposition( trajectory, igb=8 )
    
    return energy_data[term];


def get_energies_per_residue( coordinates, parm_file ):
    '''
    Calculates per-residue energies.
    '''

    os.system('MMPBSA.py -i mmgbsa.in -cp {parm} -y {rst7}'.format( parm=parm_file, 
                                                                    rst7=coordinates 
                                                                    ))
    ## Parse the data
    data = MMPBSA_API.load_mmpbsa_info('_MMPBSA_info')

    decomp_data = data['decomp']['gb']['complex']['TDC']

    return decomp_data


def print_energies( energy_data ):
    '''
    Prints energies to logfile or screen in a nice way.
    '''

    energy_keys = energy_data.keys()
    energy_keys.sort()

    for frame in range(len(energy_data[ energy_keys[0] ])):
        for key in energy_keys:
            print("{KEY}{SPACER}{VALUE}".format( KEY=key, SPACER='.'*(15-len(key)), VALUE=energy_data[ key ][frame] ) )
        print('\n')


def write_energies( energy_data, outfile ):
    '''
    Writes the energies to the outfile.
    '''

    energy_keys = energy_data.keys()
    energy_keys.sort()

    with open( outfile, 'w' ) as scorefile:

        header = ''
        for e_key in energy_keys:
            header += e_key + '\t'
        scorefile.write(header+'\n')
        
        for frame in range(len(energy_data[ energy_keys[0] ])):
            scoreline = ''
            for e_key in energy_keys:
                scoreline += str(energy_data[ e_key ][frame]) + '\t'
            scorefile.write(scoreline+'\n')


def get_ca_rmsd( traj, reference=0 ):
    
    ##TODO: Pick reference structure.
    ca_rmsd = pt.rmsd(traj, ref=reference, mask=['@CA'])
    return ca_rmsd;


def update_energies_with_rmsds( energies, rmsd ):
    energies['rmsd'] = rmsd
    return energies


#######################################################################
##### For SciPy Minimization ##########################################

def target( coordinates ):
    sander.set_positions( coordinates )
    e, f = sander.energy_forces()
    return e.tot, -numpy.array(f)


# This isn't finished, so it shouldn't be part of the API. -- mszegedy
# def minimize_with_amberff( trajectory, parm_file ):
#     '''
#     Minimized using SciPy Minimizer and Amber FF
#     '''
#     inp = sander.gas_input(8)
#     for index, frame in enumerate( trajectory ):
#         with sander.setup( parm_file, frame.xyz, frame.box, inp ):
#             res = minimize( target, frame.xyz.flatten(), method='L-BFGS-B', jac=True, tol=1e-8, options=dict(maxiter=200, disp=True ))
#         trajectory[ index ] = res.x.reshape( trajectory.n_atoms, 3 )
#         #new_frame = res.x.reshape( trajectory.n_atoms, 3 )
#         #trajectory.__setitem__( index, new_frame )

#######################################################################
##### Minimize PDBs with Sander #######################################
def run_each_core(cmlist):
    '''run a chunk of total_commands in each core

    Parameters
    ----------
    cmlist : a list of commands
    '''
    for cm in cmlist:
        os.system(cm)

def get_commands_for_my_rank(total_commands, rank, n_cores):
    import numpy as np
    arr = np.array(total_commands)

    sub_commands = np.array_split(arr, n_cores)[rank]
    if len(sub_commands) > 0:
        return sub_commands
    else:
        return ['echo nothing',]

def get_total_commands( parm7, rst7_files, overwrite, COMMAND_TEMPLATE):
    commands = []
    for rst7_file in rst7_files:
        # make sure rst7 is relative path
        abspath_rst7 = os.path.abspath(rst7_file)
        rst7 = rst7_file.split('/')[-1]
        restart_ext = '.' + rst7.split('.')[-1]
        command = COMMAND_TEMPLATE.format(
            sander = 'sander',
            overwrite = overwrite,
            minin= "min.in",
            prmtop = parm7,
            rst7 = rst7,
            abspath_rst7 = abspath_rst7,
            rst7_no_ext = rst7.strip(restart_ext))

        commands.append(command)
    return commands

def write_min_script(): 

    minimization_script = '''
&cntrl
    imin = 1, maxcyc = 1000,
    ntx = 1,
    ntxo = 2,
    ntwr = 100, ntpr = 100,
    cut = 999.0,
    ntb = 0, igb = 8,
    ntmin = 3, drms = 0.01,
/
'''
    with open("min.in", "w") as minfile:
        minfile.write( minimization_script )

def batch_minimize_with_sander( parm7_file, rst7_list, num_cores, overwrite ):
    '''
    Minimize using Amber's sander program.
    
    '''
    write_min_script()
    COMMAND_TEMPLATE = '{sander} {overwrite} -i {minin} -p {prmtop} -c {abspath_rst7} -r min_{rst7} -o out/min_{rst7_no_ext}.out -ref {abspath_rst7}'
    
    if overwrite == 1:
        overwrite_option = '-O'
    else: 
        overwrite_option = ''
    commands = get_total_commands( parm7_file, rst7_list, overwrite_option, COMMAND_TEMPLATE )

    if rank == 0:
        pwd = '/'.join(os.getcwd().split('/')[-4:])
        print('number of runs = {} in {}'.format(len(rst7_list), pwd))

    try:
        os.mkdir('out')
    except OSError:
        pass
    
    myrank_cmlist = get_commands_for_my_rank(commands, rank, num_cores)
    
    n_structures = len(myrank_cmlist)
    x = comm.gather(n_structures, root=0)
    if rank == 0:
        print('max structures per node = {}'.format(max(x)))
    run_each_core(myrank_cmlist)



#######################################################################
##### Miscellaneous Methods ###########################################

def chunks(l,n):
    '''
    Separates a list (l) into a lists of lists, where each sublist is
    of length n (or less)
    '''
    n = max(1,n)
    return [l[i:i+n] for i in range(0, len(l), n)]


#######################################################################
##### For MD Simulations ##############################################

def write_explicit_tleapfile( parm_file, rst7_file, pose_name ):
    '''
    Writes input file for tLeap program.
    The force field used here is best for running explicit
        solvent MD simulations.
    '''
    
    ## Open a writable .tleap.in file
    with open('{PNAME}.tleap.in'.format( PNAME=pose_name ),'w') as tfile:

        ## Load the ff14SBonlysc force field.
        tfile.write("source leaprc.ff14SBonlysc\n" )

        ## Library for adding Na+ and Cl- ions.
        tfile.write("loadAmberParams frcmod.ionsjc_tip3p\n")

        ## Load the dumped PDB.
        tfile.write("m = loadpdb {PNAME}.pdb\n".format( PNAME=pose_name ))

        ## Add Ions.
        tfile.write("addions m Cl- 0\n")
        tfile.write("addions m Na+ 0\n")

        ## Solvate Box.
        tfile.write("solvatebox m TIP3PBOX 10.0\n")

        ## Save the molecule as a parm file and coordinates file.
        tfile.write("saveamberparm m {parm} {rst}\n".format( parm=parm_file, rst=rst7_file ))
        tfile.write("quit" )

def change_hmass( pose_name, rst7_file ):
    '''
    Changes the hydrogen mass for 4-femtosecond time steps.
    '''
    with open("hmass.in","w") as hmass_file:
        hmass_file.write("HMassRepartition\n")
        hmass_file.write("outparm {pose}.newHmass.rst7\n".format( pose=pose_name ))
    
    os.system("$AMBERHOME/bin/parmed.py {rst7} hmass.in".format( rst7=rst7_file ))
