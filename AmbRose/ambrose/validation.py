'''Functions for validating objects and files produced by AMBRose.'''

import os

from . import errors
from . import traj_to_poses

def check_file_existence(path, file_type='AMBER file'):
    '''Checks whether a file exists. If it doesn't exist, this raises an
    `AMBERFileError`_ with the path of the file and a description of it.

    Parameters
    ----------
    path : str
        Path to file to check existence of.
    file_type : str
        Noun phrase describing file, to be inserted into error message.
    '''

    if not os.path.isfile(path):
        raise errors.AMBERFileError(
            f'Failed to produce {file_type} at {path}; cannot continue.')

def check_pose_convertibility(pose, crd_path, top_path):
    '''Checks whether a Pose changes its topology upon being converted to AMBER
    coordinates and back, given the pose and the paths to the coordinates and
    topology files it's been converted to.

    Parameters
    ----------
    pose : rosetta.core.Pose
        Pose to check.
    crd_path : str
        Path to coordinates file Pose has been converted to.
    top_path : str
        Path to topology file Pose has been converted to.
    '''

    ambered_pose = traj_to_poses.pose_from_amber_params(crd_path, top_path)
    # pylint: disable=no-member
    if ambered_pose.size() > pose.size():
        bad_residue = None
        for i in range(1, ambered_pose.size()+1):
            if ambered_pose.residue(i).type() != pose.residue(i).type():
                bad_residue = f'{i} {ambered_pose.residue(i).name()}'
                break
        raise errors.TopologySizeError(
            'Pose changed size when piped through AMBER! It somehow got '
            'bigger, which should be impossible. Please open an issue on the '
            'Github repo and describe your error there in detail. Here\'s the '
            '(first) residue that sprouted out of thin air:\n  ' + bad_residue)
    if ambered_pose.size() < pose.size():
        bad_residue = None
        for i in range(1, pose.size()+1):
            if ambered_pose.residue(i).type() != pose.residue(i).type():
                bad_residue = f'{i} {pose.residue(i).name()}'
                break
        raise errors.TopologySizeError(
            'Pose changed size when piped through AMBER! Here\'s the (first) '
            'residue that went missing:\n  ' + bad_residue)
    for i in range(1, pose.size()+1):
        bad_residue = None
        other_bad_residue = None
        if ambered_pose.residue(i).type() != pose.residue(i).type():
            bad_residue       = f'{i} {pose.residue(i).name()}'
            other_bad_residue = f'{i} {ambered_pose.residue(i).name()}'
            break
    if bad_residue is not None:
        raise errors.TopologyTypeError(
            'Pose changed at least one residue\'s type when piped through '
            'AMBER! Here\'s the offending residue:\n  ' + bad_residue + '\n'
            'It changed to:\n  ' + other_bad_residue)
    # pylint: enable=no-member
