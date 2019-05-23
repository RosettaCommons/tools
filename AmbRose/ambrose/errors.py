'''AMBRose-specific exceptions, to do with errors during conversion between
Rosetta and AMBER.'''

class TopologyError(RuntimeError):
    '''An error encountered when there's something wrong with a topology.'''

class TopologySizeError(TopologyError):
    '''An error encountered when a topology has the wrong number of residues or
    atoms.'''

class TopologyTypeError(TopologyError):
    '''An error encountered when a topology has the wrong kind of residues.'''
