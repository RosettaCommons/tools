import util
import data
import numpy as np

DU_SEQ = [
    ('cu', 'Z[1AP]g'),
    ('gZ[1AP]', 'uc'),
    ('gu', 'Z[1AP]c'),
    ('cZ[1AP]', 'ug'),
]

DU_EXPT = [
    -5.33 * 0.5,
    -5.085 * 0.5,
    -5.77 * 0.5,
    -5.25 * 0.5,
    -2.24,
    -2.66,
]

DU_EXPT_SEQ = [
    (DU_SEQ[0], DU_SEQ[1]),
    (DU_SEQ[0], DU_SEQ[3]),
    (DU_SEQ[1], DU_SEQ[2]),
    (DU_SEQ[2], DU_SEQ[3]),
    (DU_SEQ[0],),
    (DU_SEQ[2],),
]


class MinFunc(data.BaseMinFunc):
    """Minimization function for all experimental data including the DU ones.
    """
    def __init__(self, single_simulations, dangling_weight=0.1,
                 start_weight=None):
        self._simulation_dict = {}
        self._target_values = {}
        self.dangling_weight = dangling_weight

        def update_simulation_and_target_value(seq, is_duplex=True):
            self._simulation_dict[seq] = data.get_simulation_from_seq(
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
            self._target_values[seq] /= data.KT_IN_KCAL

        for seq in util.DG_ALL:
            if seq in util.DG_DANGLING:
                if dangling_weight != 0:
                    update_simulation_and_target_value(seq, False)
            else:
                update_simulation_and_target_value(seq)
                if seq[0] != seq[1]:
                    update_simulation_and_target_value((seq[1], seq[0]))

        for i, seq in enumerate(DU_EXPT_SEQ):
            if i < 4:
                sim1 = data.get_simulation_from_seq(seq[0],
                                                    single_simulations)
                sim2 = data.get_simulation_from_seq(seq[1],
                                                    single_simulations)
                self._simulation_dict[seq] = data.SimulationCombination(
                    [sim1, sim2], sim_weights='mean')
            else:
                self._simulation_dict[seq] = data.get_simulation_from_seq(
                    seq[0], single_simulations, symmetrize=False)
            self._target_values[seq] = DU_EXPT[i]

        start_weight = (
            single_simulations.values()[0].curr_weight
            if start_weight is None else start_weight)
        self._value = self._deriv = None
        self.reweight(start_weight)

    def reweight(self, new_weight):
        self._value = 0
        self._deriv = np.zeros(data.N_SCORE_TERMS)
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
