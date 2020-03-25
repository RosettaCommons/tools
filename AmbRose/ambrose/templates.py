'''Functions that produce template ``mdin_dict``s for various common use cases,
with a nicer API than AMBER's default. To be used with either of the top-level
functions run_md_ or `dict_to_namelist_file`_, if you are looking to directly
run a simulation or to just create a NAMELIST file, respectively.'''

def minimize_wet(*, cutoff=15.):
    '''Template of dict of mdin parameters for an explicit solvent minimization.

    Parameters
    ----------
    cutoff : int or float, optional
        Van der Waals cutoff in angstroms.

    Returns
    -------
    dict
    '''

    return {
        'imin':   1, # minimization on
        'maxcyc': 500,
        'ncyc':   250,
        'cut':    cutoff,
    }

def minimize_dry(*, cutoff=15.):
    '''Template of dict of mdin parameters for an implicit solvent minimization.

    Parameters
    ----------
    cutoff : int or float, optional
        Van der Waals cutoff in angstroms.

    Returns
    -------
    dict
    '''

    return {
        'imin':   1, # minimization on
        'maxcyc': 500,
        'ncyc':   250,
        'ntb':    0, # periodicity off
        'igb':    1, # Born solvent on
        'cut':    cutoff,
    }

def simulate_wet(*, duration=None, temperature=None, starting_temperature=0,
                 cutoff=15., output_interval=1.):
    '''Template of dict of mdin parameters for an explicit solvent simulation.

    Parameters
    ----------
    duration : int or float
        Duration of simulation in picoseconds.
    temperature : int or float
        Temperature of simulation in kelvins. 0 by default.
    starting_temperature : int or float, optional
        Starting temperature of simulation in kelvins.
    cutoff : int or float
        Van der Waals cutoff in angstroms.
    output_interval : int or float
        Amount of time to wait in between writing frames and log data, in
        picoseconds.

    Returns
    -------
    dict
    '''

    assert duration is not None, 'Duration must be specified.'
    assert temperature is not None, 'Temperature must be specified.'
    return {
        'nstlim': int(duration*1000),
        'tempi':  starting_temperature,
        'temp0':  temperature,
        'cut' :   cutoff,
        'ntt':    1, # weak coupling thermostat
        'ntp':    1, # constant pressure dynamics
        'pres0':  1.013, # 1 atm pressure
        'ntwx':   int(output_interval*1000),
        'ntpr':   int(output_interval*1000),
    }

def simulate_dry(*, duration=None, temperature=None, starting_temperature=0,
                 cutoff=15., output_interval=1., seed=-1):
    '''Template of dict of mdin parameters for an implicit solvent simulation.

    Parameters
    ----------
    duration : int or float
        Duration of simulation in picoseconds.
    temperature : int or float
        Temperature of simulation in kelvins.
    starting_temperature : int or float, optional
        Starting temperature of simulation in kelvins. 0 by default.
    cutoff : int or float, optional
        Van der Waals cutoff in angstroms. 15 by default.
    output_interval : int or float, optional
        Amount of time to wait in between writing frames and log data, in
        picoseconds. 1 by default.
    seed : int, optional
        Seed to use for random operations (e.g. Langevin thermostat).
        Nonnegative integers are valid seeds, while negative integers cause
        the seed to be determined by the hardware.

    Returns
    -------
    dict
    '''

    assert duration is not None, 'Duration must be specified.'
    assert temperature is not None, 'Temperature must be specified.'
    return {
        'nstlim':   int(duration*1000),
        'tempi':    starting_temperature,
        'temp0':    temperature,
        'cut':      cutoff,
        'igb':      1, # Born solvent
        'ntt':      3, # Langevin thermostat
        'gamma_ln': 1.0, # Langevin collision frequency
        'ntwx':     int(output_interval*1000),
        'ntpr':     int(output_interval*1000),
        'ig':       seed
    }
