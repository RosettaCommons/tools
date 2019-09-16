'''Enums for specifying options (currently solvent options) for AMBERMovers.'''

import enum

class Solvent(enum.Enum):
    '''A solvent choice for anything that requires you to set or provide a
    solvent, like `pose_to_amber_params`_ or `AMBERSimulateMover`_. Can be:

    SPCE_WATER
        SPC/E (extended simple point charge) waters.
    TIP3P_WATER
        TIP3P (transferable intramolecular potential with 3 points) waters.
    METHANOL
        Methanol.
    '''

    SPCE_WATER = 'SPCE_WATER'
    TIP3P_WATER = 'TIP3P_WATER'
    METHANOL = 'METHANOL'

class SolventShape(enum.Enum):
    '''The choice of shape for the periodic boundaries you're using for an
    explicit-solvent simulation, for anything that requires you to set or
    provide it, like `pose_to_amber_params`_ or `AMBERSimulateMover`_. Can be:

    BOX
        Rectangular prism.
    OCT
        Truncated octahedron.
    '''

    BOX = 'BOX'
    OCT = 'OCT'
