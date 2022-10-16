'''Utility functions applicable to any submodule of AMBRose.'''

import sys
from . import consts

def debug_print(value, *more_values, sep=' ', end='\n', file=sys.stdout,
                flush=True):
    if consts.DEBUG:
        print(value, *more_values, sep=sep, end=end, file=file, flush=flush)
