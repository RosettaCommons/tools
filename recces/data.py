from abc import ABCMeta, abstractmethod, abstractproperty
from collections import defaultdict
import glob
import os
import re

import numpy as np
from scipy.misc import logsumexp

import util

N_SCORE_TERMS = 10
# KT_IN_KCAL = 0.593
KT_IN_KCAL = 0.61633135471  # 37 Celcius

# WHAM parameters
# TODO: This assumes the energy mostly falls in -100 to 800 RU. Should be
# able to use a general automatic way to figure out the bin ranges.
BIN_SIZE = 0.1
MIN_ENERGY = -100.05
MAX_ENERGY = 800.05
WHAM_BINS = np.arange(MIN_ENERGY, MAX_ENERGY + BIN_SIZE, BIN_SIZE)
N_BINS = len(WHAM_BINS)
WHAM_ENERGIES = (WHAM_BINS[1:] + WHAM_BINS[:-1]) * 0.5
WHAM_CONVERGE_LIMIT = 0.000001
WHAM_NORMALIZE_FILE_NAME = "norm_factor.npy"


class BaseMinFunc(object):
    """Base class for a function that can be evaluated to get its value and
    derivative, and be reweighted with a new set of parameters.
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def value(self):
        """Current value of the function."""
        raise NotImplementedError

    @abstractproperty
    def deriv(self):
        """Current derivative of the function."""
        raise NotImplementedError

    @abstractmethod
    def reweight(self, new_weight):
        """Reweight the parameter with the new_weight."""
        raise NotImplementedError


class SingleSimulation(BaseMinFunc):
    """API for data in a single Turner-rule style simulation
    (one-strand or duplex).
    """

    def __init__(self, data_folder, orig_weight=None,
                 use_existing_norm_factor=True, down_sampling_ratio=None,
                 bootstrap=False, legacy_output=False):
        """Create a SingleSimulation object.

        Parameters
        ----------
        data_folder : Folder name where the data is stored. The name of the
                      folder should correspond to the sequence being simulated.
        orig_weight : The original score term weights used in the simulation.
        use_existing_norm_factor : Whether to use existing normalization factor
                                   computed by WHAM stored in folder.
        down_sampling_ratio : If not None, the raw data will be down sampled
                              with the given ratio.
        bootstrap : If True, resample the data using bootstrapping.
        legacy_output : For analyzing deprecated old output format.
        """
        # TODO: Folder name should correspond to the sequence here.
        # Can be made more general.
        self.name = os.path.basename(data_folder)
        self.curr_weight = (np.array(orig_weight) if orig_weight is not None
                            else np.ones(N_SCORE_TERMS))

        # Load data
        kt_to_fn = _get_kt_to_file_names(data_folder, ".bin.gz")
        raw_data = []
        kt_list = []
        for kt, file_list in kt_to_fn.iteritems():
            raw_data.append(
                np.vstack(util.load_2d_bin_gz(fn) for fn in file_list))
            kt_list.append(kt)

        # Bootstraping
        if bootstrap:
            raw_data_bootstrap = []
            for d in raw_data:
                resample_index = np.random.randint(len(d), size=len(d))
                raw_data_bootstrap.append(d[resample_index])
            raw_data = raw_data_bootstrap

        # Compute normalization with WHAM
        normalization_file = os.path.join(data_folder,
                                          WHAM_NORMALIZE_FILE_NAME)
        if use_existing_norm_factor and os.path.isfile(normalization_file):
            normalization = np.load(normalization_file)
        else:
            histograms = [self._energy2hist(data[:, 1], data[:, 0])
                          for data in raw_data]
            normalization, _ = _wham(histograms, kt_list)
            np.save(normalization_file, normalization)

        # Down sampling
        if down_sampling_ratio is not None:
            resampled_data = []
            while raw_data:
                data = raw_data.pop()
                new_data = _rand_choice(data, down_sampling_ratio)
                resampled_data.append(new_data)
            raw_data = list(reversed(resampled_data))
            del resampled_data
            del data
        histograms = [self._energy2hist(data[:, 1], data[:, 0])
                      for data in raw_data]

        dos_raw = np.sum(histograms, axis=0) / normalization
        seq = self.name.split('_')
        torsion_volume = util.torsion_volume(*seq)
        normalization *= np.sum(dos_raw) / torsion_volume

        full_data = np.vstack(raw_data)
        energy = full_data[:, 1]
        normalization_idx = self._energy2binidx(energy)
        data_weight = full_data[:, 0] / normalization[normalization_idx]
        self._energy0 = energy
        self._data_weight = data_weight
        if legacy_output:
            self._score_terms = full_data[:, 2:] / self.curr_weight
        else:
            self._score_terms = full_data[:, 2:]

        self._update_free_energy_and_deriv()

    def reweight(self, new_weight):
        if np.any(new_weight != self.curr_weight):
            # Only reweight when weights have changed
            self.curr_weight = new_weight
            self._update_free_energy_and_deriv()

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv

    def get_ST_weights(self, kt_list):
        """Get optimized Simulated Tempering weights for a list of kT.

        Parameters
        ----------
        kt_list : List of kT (float).

        Returns
        -------
        List of optimized simulated tempering weights.
        """
        new_energy = self._get_curr_energy()
        dos = self._energy2hist(new_energy, self._data_weight)
        weights = [0]
        for kt1, kt2 in zip(kt_list[:-1], kt_list[1:]):
            hist1 = np.column_stack((
                WHAM_ENERGIES, dos * np.exp(-WHAM_ENERGIES / kt1)))
            hist2 = np.column_stack((
                WHAM_ENERGIES, dos * np.exp(-WHAM_ENERGIES / kt2)))
            delta_weight = util.get_ST_delta(hist1, hist2, kt1, kt2)
            weights.append(weights[-1] + delta_weight)
        return weights

    def _update_free_energy_and_deriv(self):
        """Update the free energy and the derivatives with current parameters.
        """
        new_energy = self._get_curr_energy()
        # Use logsumexp here to reduce numerical errors
        free_energy = -logsumexp(-new_energy, b=self._data_weight)

        # Derivatives for each score term
        deriv = np.sum(self._data_weight[:, None] * self._score_terms *
                       np.exp(-new_energy + free_energy)[:, None], axis=0)
        self._value = free_energy
        self._deriv = deriv

    def _get_curr_energy(self):
        # score_weight_diff = self.curr_weight - self._orig_weight
        # energy = self._energy0 + self._score_terms.dot(score_weight_diff)
        return self._score_terms.dot(self.curr_weight)

    @staticmethod
    def _energy2hist(energy, weights):
        binidx = SingleSimulation._energy2binidx(energy)
        hist = np.bincount(binidx, weights, len(WHAM_BINS) - 1)
        return hist

    @staticmethod
    def _energy2binidx(energies):
        binidx = (energies - MIN_ENERGY) / BIN_SIZE
        binidx[binidx < 0] = 0
        binidx[binidx > len(WHAM_BINS) - 2] = len(WHAM_BINS) - 2
        return binidx.astype(np.uint32)


class SingleHistSimulation(BaseMinFunc):
    """Histogram based simulation data. For computing free energies.
    Cannot be minimized.
    """
    def __init__(self, data_folder):
        """Create a SingleHistSimulation object.

        Parameters
        ----------
        data_folder : Folder name where the data is stored. The name of the
                      folder should correspond to the sequence being simulated.
        """
        self.name = os.path.basename(data_folder)

        # Scores for Hist bins.
        # Assumes all the hist_scores.gz files are the same
        hist_scores = util.load_1d_bin_gz(
            _get_file_names(data_folder, '_hist_scores.gz')[0],
            dtype=np.float64)

        kt_to_fn = _get_kt_to_file_names(data_folder, '.hist.gz')
        hist_list = []
        kt_list = []
        for kt, fn_list in kt_to_fn.iteritems():
            all_sub_hists = [util.load_1d_bin_gz(fn, dtype=np.uint64)
                             for fn in fn_list]
            hist_counts = np.sum(all_sub_hists, axis=0)
            hist_list.append(hist_counts)
            kt_list.append(kt)

        _, dos_raw = _wham(hist_list, kt_list, hist_scores)
        seq = self.name.split('_')
        torsion_volume = util.torsion_volume(*seq)
        normalization = sum(dos_raw) / torsion_volume
        self._dos = dos_raw / normalization
        self._dos_scores = hist_scores

    def get_free_energy(self, kt=1):
        """Free energy of the system.

        Parameters
        ----------
        kt : current temperature in Rosetta Unit. Default to 1.

        Returns
        -------
        Computed free energy.
        """
        return -kt * logsumexp(-self._dos_scores / kt, b=self._dos)

    @property
    def value(self):
        return self.get_free_energy()

    @property
    def deriv(self):
        return np.zeros(N_SCORE_TERMS)

    def reweight(self, new_weight):
        pass


class SimulationCombination(BaseMinFunc):
    """Linear combination of multiple simulations.
    """
    def __init__(self, simualtions, sim_weights='mean',
                 start_weight=None):
        """Create a SimulationCombination object.

        Parameters
        ----------
        simulations : List of simulations to be combined.
        sim_weights : Weights for each simulation in the linear combination.
                      It should be a list of float of the same size as the
                      simulations. Two special cases are 'mean' and 'sum'.
        start_weight : Initial weight of the simulation. If None it will use
                       the weight of the first one in positive_simulation.
        """
        self.simulations = simualtions
        if sim_weights == 'mean':
            sim_weights = [1. / len(simualtions)] * len(simualtions)
        elif sim_weights == 'sum':
            sim_weights = [1.] * len(simualtions)
        if len(simualtions) != len(sim_weights):
            raise ValueError("Length of simulations and "
                             " sim_weights are unequal")
        self.sim_weights = sim_weights

        start_weight = (simualtions[0].curr_weight
                        if start_weight is None else start_weight)
        self.curr_weight = start_weight
        self._value = self._deriv = None
        self.reweight(start_weight)

    def reweight(self, new_weight):
        self.curr_weight = new_weight
        self._value = 0
        self._deriv = np.zeros(N_SCORE_TERMS)
        for s, w in zip(self.simulations, self.sim_weights):
            s.reweight(new_weight)
            self._value += s.value * w
            self._deriv += s.deriv * w

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv


class TurnerRuleMinimizeFunc(BaseMinFunc):
    """This object combines all simulations needed for Turner Rule computation
    and returns the square error of the simulations. Useful for minimizing the
    parameters to achieve maximum fit to the data.
    """
    def __init__(self, single_simulations, dangling_weight=0.1,
                 start_weight=None):
        """Create a TurnerRuleMinimizeFunc object.

        Parameters
        ----------
        all_single_simulations : Dict of single simulations needed for Turner
                rule calculation. Use the name of the simulation as the key.
        dangling_weight : Weight for dangling end data in the minimization
                (duplex data is weighted as 1).
        start_weight : Initial weight of the simulation. If None it will use
                the weight of the first one in all_single_simulations.
        """
        self._simulation_dict = {}
        self._target_values = {}
        self.dangling_weight = dangling_weight

        def update_simulation_and_target_value(seq, is_duplex=True):
            self._simulation_dict[seq] = get_simulation_from_seq(
                seq, single_simulations, symmetrize=False)

            if seq in util.DG_ALL:
                self._target_values[seq] = util.DG_ALL[seq]
            else:
                self._target_values[seq] = util.DG_ALL[(seq[1], seq[0])]
            if is_duplex:
                s1 = util.seq_parse(seq[0])
                s2 = util.seq_parse(seq[1])
                terminal1 = tuple(sorted((s1[0], s2[1])))
                terminal2 = tuple(sorted((s1[1], s2[0])))
                if terminal1 == ('c', 'g') and terminal2 in util.DG_TERMINAL:
                    self._target_values[seq] += util.DG_TERMINAL[terminal2]
                elif terminal2 == ('c', 'g') and terminal1 in util.DG_TERMINAL:
                    self._target_values[seq] -= util.DG_TERMINAL[terminal1]
            self._target_values[seq] /= KT_IN_KCAL

        for seq in util.DG_CANONICAL:
            update_simulation_and_target_value(seq)
            if seq[0] != seq[1]:
                update_simulation_and_target_value((seq[1], seq[0]))

        if dangling_weight != 0:
            for seq in util.DG_DANGLING:
                update_simulation_and_target_value(seq, False)
        start_weight = (
            single_simulations.values()[0].curr_weight
            if start_weight is None else start_weight)
        self._value = self._deriv = None
        self.reweight(start_weight)

    def reweight(self, new_weight):
        self._value = 0
        self._deriv = np.zeros(N_SCORE_TERMS)
        for seq, sim in self._simulation_dict.iteritems():
            sim.reweight(new_weight)
            expt_free_energy = self._target_values[seq]
            weight = self.dangling_weight if seq in util.DG_DANGLING else 1
            self._value += (sim.value - expt_free_energy) ** 2 * weight
            self._deriv += (
                2 * (sim.value - expt_free_energy) * sim.deriv * weight)

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv


def load_raw_simulations(folder, *args, **kwargs):
    """Load all the simulations in folder.

    Parameters
    ----------
    folder : Name to the folder containing data for all sequences.
    optional arguments : see SingleSimulation docstring

    Returns
    -------
    simulations : dict of seq_name: SingleSimulation.
    """
    simulations = {}
    for fn in glob.glob(folder + "/*"):
        sim = SingleSimulation(fn, *args, **kwargs)
        simulations[sim.name] = sim
    return simulations


def get_simulation_from_seq(seq, raw_simulations, symmetrize=True):
    '''Generate composite simulation from sequence.

    Parameters
    ----------
    seq : sequence string, e.g. auu_aau for a duplex
    raw_simulations : simulation dict from load_raw_simulations
    symmetrize : symmetrize the sequence

    Returns
    -------
    sim : composite simulation for the imput seq
    '''
    def combine_sim(positive, negative):
        sim_weight = [1] * len(positive) + [-1] * len(negative)
        return SimulationCombination(positive + negative, sim_weight)

    def _get_duplex(seq1, seq2):
        positive, negative = _get_duplex_seqs(seq1, seq2)
        postive_simu = [raw_simulations[name] for name in positive]
        negative_simu = [raw_simulations[name] for name in negative]
        return combine_sim(postive_simu, negative_simu)

    def _get_dangling(seq1, seq2):
        positive, negative = _get_dangling_seqs(seq1, seq2)
        postive_simu = [raw_simulations[name] for name in positive]
        negative_simu = [raw_simulations[name] for name in negative]
        return combine_sim(postive_simu, negative_simu)

    seq1, seq2 = util.seq_parse(seq[0]), util.seq_parse(seq[1])
    if len(seq1) == len(seq2):
        sim = _get_duplex(seq1, seq2)
        if symmetrize and seq1 != seq2:
            sim_sym = _get_duplex(seq2, seq1)
            sim = SimulationCombination([sim, sim_sym], sim_weights='mean')
        return sim
    else:
        sim = _get_dangling(seq1, seq2)
        return sim


def get_g_from_seq(seq, raw_g, symmetrize=True):
    '''Get the delta G from free energies of each construct.

    Parameters
    ----------
    seq : sequence string, e.g. auu_aau for a duplex
    raw_g : dict of sequence: G
    symmetrize : symmetrize the sequence

    Returns
    -------
    Free energy (delta-G) for the sequence.
    '''
    def _get_duplex(seq1, seq2):
        positive, negative = _get_duplex_seqs(seq1, seq2)
        return (np.sum([raw_g[name] for name in positive], axis=0) -
                np.sum([raw_g[name] for name in negative], axis=0))

    def _get_dangling(seq1, seq2):
        positive, negative = _get_dangling_seqs(seq1, seq2)
        return (np.sum([raw_g[name] for name in positive], axis=0) -
                np.sum([raw_g[name] for name in negative], axis=0))

    seq1, seq2 = util.seq_parse(seq[0]), util.seq_parse(seq[1])
    if len(seq1) == len(seq2):
        g = _get_duplex(seq1, seq2)
        if symmetrize and seq1 != seq2:
            g_sym = _get_duplex(seq2, seq1)
            g = 0.5 * (g + g_sym)
        return g
    else:
        return _get_dangling(seq1, seq2)


def get_terminal_g(base1, base2, raw_g):
    '''Get the terminal constribution from
    free energies of each construct.

    Parameters
    ----------
    base1 : base string for the first base
    base1 : base string for the second base
    raw_g : dict of sequence: G

    Returns
    -------
    Free energy (delta-G) for the terminal contribution.
    '''
    combinations = (
        ('g' + base1, base2 + 'c'),
        ('g' + base2, base1 + 'c'),
        ('c' + base1, base2 + 'g'),
        ('c' + base2, base1 + 'g'),
    )

    terminal_g = []
    for seq1, seq2 in combinations:
        g1 = get_g_from_seq((seq1, seq2), raw_g, symmetrize=False)
        g2 = get_g_from_seq((seq2, seq1), raw_g, symmetrize=False)
        terminal_g.append((g1 - g2) / 2)
    return terminal_g


def _wham(histograms, kt_list, energies=WHAM_ENERGIES):
    """Weighted Histogram Analysis Method (WHAM).
    Returns the the normalization factor for each energy bin,
    and the full density of states (DoS).
    """
    kt_list = np.asarray(kt_list)

    hist_sum = np.sum(histograms, axis=0)
    n_sample_per_hist = np.sum(histograms, axis=1)

    exp_f_list = np.ones_like(kt_list)
    boltzmann_factor = np.exp(
        -np.tile(energies, (len(kt_list), 1)).T / kt_list)
    bin_size = energies[1] - energies[0]

    while True:
        exp_f_list_old = exp_f_list
        weight_per_bin = np.sum(n_sample_per_hist * bin_size * exp_f_list *
                                boltzmann_factor, axis=1)

        dos_raw = hist_sum / weight_per_bin  # Density of states
        exp_f_list = 1.0 / np.sum(dos_raw * bin_size * boltzmann_factor.T,
                                  axis=1)
        diff = (np.abs(exp_f_list - exp_f_list_old) /
                (exp_f_list + exp_f_list_old))
        if np.all(diff <= WHAM_CONVERGE_LIMIT):
            return weight_per_bin, dos_raw


def _rand_choice(data, ratio):
    weight = data[:, 0]
    cumsum_weight = np.cumsum(weight)
    sum_weight = np.sum(weight)
    n_sample_new = sum_weight * ratio
    select_idx = np.searchsorted(
        cumsum_weight, np.random.rand(n_sample_new) * sum_weight)
    new_data = data[select_idx]
    new_data[:, 0] = 1
    return new_data


def _get_duplex_seqs(seq1, seq2):
    s0 = [util.END_BP[0]] + list(seq1)
    s1 = list(seq2) + [util.END_BP[1]]
    s2 = s0[:-1]
    s3 = s1[1:]
    positive = ''.join(s0) + '_' + ''.join(s1), ''.join(s2), ''.join(s3)
    negative = ''.join(s2) + '_' + ''.join(s3), ''.join(s0), ''.join(s1)
    return positive, negative


def _get_dangling_seqs(seq1, seq2):
    s0 = [util.END_BP[0]] + list(seq1)
    s1 = list(seq2) + [util.END_BP[1]]
    dup0 = ''.join(s0) + '_' + ''.join(s1)
    if len(s0) > len(s1):
        dup1 = ''.join(s0[:-1]) + '_' + ''.join(s1)
        strand0 = ''.join(s0)
        strand1 = ''.join(s0[:-1])
    else:
        dup1 = ''.join(s0) + '_' + ''.join(s1[1:])
        strand0 = ''.join(s1)
        strand1 = ''.join(s1[1:])
    positive = dup0, strand1
    negative = dup1, strand0
    return positive, negative


def _get_kt_to_file_names(folder_name, fn_suffix):
    kt_to_fn = defaultdict(list)
    for fn in _get_file_names(folder_name, fn_suffix):
        basename = os.path.basename(fn)
        # The names should look like "ST_1_3.00FN_SUFFIX"
        kt = float(basename.split('_')[-1].replace(fn_suffix, ''))
        if kt < 0:
            kt = np.inf
        kt_to_fn[kt].append(fn)
    return kt_to_fn


def _get_file_names(folder_name, fn_suffix):
    folder_pattern = folder_name[:]
    folder_pattern = re.sub(r'\[', '[[]', folder_pattern)
    folder_pattern = re.sub(r'(?<!\[)\]', '[]]', folder_pattern)
    full_pattern = os.path.join(folder_pattern, "*%s" % fn_suffix)
    return glob.glob(full_pattern)
