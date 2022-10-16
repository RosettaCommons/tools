'''Functions for selecting Poses from (indexable) sequences of Poses, to be used
for the ``pose_selector`` attribute of `AMBERMover`_s that have one.'''

import pyrosetta as pr # pylint: disable=import-error

def last(l):
    '''Returns the last Pose in a sequence.

    Parameters
    ----------
    l : list or TrajToPoses
        A sequence (indexable object) containing Poses that you want to select a
        Pose from.

    Returns
    -------
    rosetta.core.Pose
        The selected Pose from the given list.
    '''
    return l[-1]

def lowest_energy_from_end(l):
    '''Returns the lowest energy Pose in roughly the rightmost half of a
    sequence. The region that the Pose is selected from will contain at minimum
    the rightmost 1/e of the list, and is extended leftwards until a Pose is
    reached that is more than 2 standard deviations from the average of the
    rightmost 1/e (this Pose is excluded).

    Parameters
    ----------
    l : list or TrajToPoses
        A sequence (indexable object) containing Poses that you want to select a
        Pose from.

    Returns
    -------
    rosetta.core.Pose
        The selected Pose from the given list.
    '''

    if not l:
        raise ValueError('Input list must have at least one element.')
    if len(l) == 1:
        return l[0]

    scorefxn = pr.get_fa_scorefxn()
    if len(l) == 2:
        energies = tuple(scorefxn.score(pose) for pose in l)
        return l[energies.index(min(energies))]
    start_index = int(len(l)*0.63212)
    energies = []
    # let's not assume l is sliceable, although both lists and TrajToPoses
    # objects are. also we're building the list backwards because we're gonna
    # append Poses with earlier indices to it later
    for i in range(len(l)-1, start_index-1, -1):
        energies.append(scorefxn.score(l[i]))
    # this would be easier with numpy but I don't want to introduce a
    # dependency just to do a mean and a standard deviation
    mean = sum(energies)/len(energies)
    std = (sum((e-mean)**2 for e in energies) / len(energies))**0.5
    start_index -= 1 # guaranteed >= 0 since len(l) > 1
    while True:
        next_energy = scorefxn.score(l[start_index])
        if abs(next_energy-mean)/std > 2:
            break
        energies.append(next_energy)
        if start_index == 0:
            break
        start_index -= 1
    # the lowest energy is basically guaranteed to be unique, so this is okay:
    min_index = energies.index(min(energies))
    return l[start_index + min_index]
