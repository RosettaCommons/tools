from abc import ABCMeta, abstractmethod, abstractproperty
from collections import defaultdict
import glob
import os
import re

import numpy as np
from scipy.misc import logsumexp

import turner_util1 as util


N_SCORE_TERMS = 10
KT_IN_KCAL = 0.593

# WHAM parameters
# TODO: This assumes the energy mostly falls in -100 to 500 RU. Should be
# able to use a general automatic way to figure out the bin ranges.
BIN_SIZE = 0.1
MIN_ENERGY = -100.25
MAX_ENERGY = 500.25
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

    def __init__(self, data_folder, orig_weight=np.ones(N_SCORE_TERMS),
                 use_existing_norm_factor=True, down_sampling_ratio=None):
        """Create a SingleSimulation object.

        Parameters
        ----------
        data_folder : Folder name where the data is stored. The name of the
                      folder should correspond to the sequence being simulated.
        orig_weight : The original score term weights used in the simulation.
        ru_in_kt : The initial guess of Rosetta Unit in the unit of kT
                   (at room temperature).
        use_existing_norm_factor : Whether to use existing normalization factor
                                   computed by WHAM stored in folder.
        down_sampling_ratio : If not None, the raw data will be down sampled
                              with the given ratio.
        """
        self._folder = data_folder
        # TODO: Folder name should correspond to the sequence here.
        # Can be made more general.
        self.name = os.path.basename(data_folder)
        self.curr_weight = self._orig_weight = np.array(orig_weight)

        # Load data
        files = defaultdict(list)
        pattern = data_folder[:]
        pattern = re.sub(r'\[', '[[]', pattern)
        pattern = re.sub(r'(?<!\[)\]', '[]]', pattern)
        file_names = os.path.join(pattern, "*.bin.gz")

        for fn in glob.glob(file_names):
            basename = os.path.basename(fn)
            # The names should look like "ST_1_3.00.bin.gz"
            kt = float(basename.split('_')[-1].replace('.bin.gz', ''))
            if kt < 0:
                kt = np.inf
            files[kt].append(fn)
        raw_data = []
        kt_list = []
        for kt, file_list in files.iteritems():
            raw_data.append(
                np.vstack(util.load_2d_bin_gz(fn) for fn in file_list))
            kt_list.append(kt)

        # Compute normalization with WHAM
        normalization_file = os.path.join(data_folder, WHAM_NORMALIZE_FILE_NAME)
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
        normalization *= sum(dos_raw) / torsion_volume

        full_data = np.vstack(raw_data)
        energy = full_data[:, 1]
        normalization_idx = self._energy2binidx(energy)
        data_weight = full_data[:, 0] / normalization[normalization_idx]
        self._energy0 = energy
        self._data_weight = data_weight
        self._score_terms = full_data[:, 2:] / self._orig_weight
        self._update_free_energy_and_deriv()

    def reweight(self, new_weight):
        if np.any(new_weight != self.curr_weight):
            self.curr_weight = new_weight
            self._update_free_energy_and_deriv()

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv

    def _update_free_energy_and_deriv(self):
        """Update the free energy and the derivatives with current parameters.
        """
        # score_weight_diff = self.curr_weight - self._orig_weight
        # new_energy = self._energy0 + self._score_terms.dot(score_weight_diff)
        new_energy = self._score_terms.dot(self.curr_weight)
        # Use logsumexp here to reduce numerical errors
        free_energy = -logsumexp(-new_energy, b=self._data_weight)

        # Derivatives for each score term
        deriv = np.sum(self._data_weight[:, None] * self._score_terms *
                       np.exp(-new_energy + free_energy)[:, None], axis=0)
        self._value = free_energy
        self._deriv = deriv

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


class SimulationCombination(BaseMinFunc):
    """Combine multiple simulations by adding and/or substracting them from
     each other.
    """
    def __init__(self, positive_simualtions, negative_simulations,
                 start_weight=None):
        """Create a SimulationCombination object.

        Parameters
        ----------
        positive_simulations : Simulations that contribute positively to the
                               combination.
        negative_simulations : Simulations that contribute negatively.
        start_weight : Initial weight of the simulation. If None it will use the
                       weight of the first one in positive_simulation.
        """
        self.positive_simulations = positive_simualtions
        self.negative_simulations = negative_simulations
        start_weight = positive_simualtions[0].curr_weight \
            if start_weight is None else start_weight
        self.curr_weight = start_weight
        self._value = self._deriv = None
        self.reweight(start_weight)

    def reweight(self, new_weight):
        self.curr_weight = new_weight
        self._value = 0
        self._deriv = np.zeros(N_SCORE_TERMS)
        for simulation in self.positive_simulations:
            simulation.reweight(new_weight)
            self._value += simulation.value
            self._deriv += simulation.deriv

        for simulation in self.negative_simulations:
            simulation.reweight(new_weight)
            self._value -= simulation.value
            self._deriv -= simulation.deriv

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv


class MeanSimulation(BaseMinFunc):
    """Combine multiple simulations by taking the mean of them.
    """
    def __init__(self, simulations, start_weight=None):
        """Create a MeanSimulation object.

        Parameters
        ----------
        simulations : Simulations to be averaged.
        start_weight : Initial weight of the simulation. If None it will use the
                       weight of the first one in simulations.
        """
        self._simulations = simulations
        start_weight = simulations[0].curr_weight if start_weight is None \
            else start_weight
        self.curr_weight = start_weight
        self._value = self._deriv = None
        self.reweight(start_weight)

    def reweight(self, new_weight):
        for s in self._simulations:
            s.reweight(new_weight)
        self._value = np.mean([s.value for s in self._simulations])
        self._deriv = np.mean([s.deriv for s in self._simulations], axis=0)

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
    def __init__(self, single_simulations, dangling_weight=0,
                 start_weight=None):
        """Create a TurnerRuleMinimizeFunc object.

        Parameters
        ----------
        all_single_simulations : Dict of single simulations needed for Turner rule calculation. Use the
                                 name of the simulation as the key.
        dangling_weight : Weight for dangling end data in the minimization (duplex data is weighted as 1).
        start_weight : Initial weight of the simulation. If None it will use the
                       weight of the first one in all_single_simulations.
        """
        self._simulation_dict = {}
        self.dangling_weight = dangling_weight
        for seq in util.DG_CANONICAL:
            self._simulation_dict[seq] = get_simulation_from_seq(seq, single_simulations)

        if dangling_weight != 0:
            for seq in util.DG_DANGLING:
                self._simulation_dict[seq] = get_simulation_from_seq(seq, single_simulations)
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
            expt_free_energy = util.DG_ALL[seq] / KT_IN_KCAL
            weight = self.dangling_weight if seq in util.DG_DANGLING else 1
            self._value += (sim.value - expt_free_energy) ** 2 * weight
            self._deriv += 2 * (sim.value - expt_free_energy) * sim.deriv * weight

    @property
    def value(self):
        return self._value

    @property
    def deriv(self):
        return self._deriv


def load_raw_simulations(folder, *args, **kwargs):
    simulations = {}
    for fn in glob.glob(folder + "/*"):
        sim = SingleSimulation(fn, *args, **kwargs)
        simulations[sim.name] = sim
    return simulations


def get_simulation_from_seq(seq, raw_simulations):
    def _get_duplex(seq1, seq2):
        s0 = [util.END_BP[0]] + list(seq1)
        s1 = list(seq2) + [util.END_BP[1]]
        s2 = s0[:-1]
        s3 = s1[1:]
        positive = ''.join(s0) + '_' + ''.join(s1), ''.join(s2), ''.join(s3)
        negative = ''.join(s2) + '_' + ''.join(s3), ''.join(s0), ''.join(s1)
        postive_simu = [raw_simulations[name] for name in positive]
        negative_simu = [raw_simulations[name] for name in negative]
        return SimulationCombination(postive_simu, negative_simu)

    def _get_dangling(seq1, seq2):
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
        postive_simu = [raw_simulations[name] for name in positive]
        negative_simu = [raw_simulations[name] for name in negative]
        return SimulationCombination(postive_simu, negative_simu)

    seq1, seq2 = util.seq_parse(seq[0]), util.seq_parse(seq[1])
    if len(seq1) == len(seq2):
        sim = _get_duplex(seq1, seq2)
        sym_seq = seq2, seq1
        if sym_seq != (seq1, seq2):
            sim_sym = _get_duplex(sym_seq[0], sym_seq[1])
            sim = MeanSimulation([sim, sim_sym])
        return sim
    else:
        sim = _get_dangling(seq1, seq2)
        return sim


def _wham(histograms, kt_list):
    """Use Weighted Histogram Analysis Method (WHAM) to compute the full
    density of states (DoS) and the normalization factor for each
    energy bin.
    """
    kt_list = np.asarray(kt_list)

    hist_sum = np.sum(histograms, axis=0)
    n_sample_per_hist = np.sum(histograms, axis=1)

    exp_f_list = np.ones_like(kt_list)
    boltzmann_factor = np.exp(
        -np.tile(WHAM_ENERGIES, (len(kt_list), 1)).T / kt_list)

    while True:
        exp_f_list_old = exp_f_list
        weight_per_bin = np.sum(n_sample_per_hist * BIN_SIZE * exp_f_list *
                                boltzmann_factor, axis=1)
        dos_raw = hist_sum / weight_per_bin  # Density of states
        exp_f_list = 1.0 / np.sum(dos_raw * BIN_SIZE * boltzmann_factor.T,
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

