from glob import glob
import numpy as np
from data import SingleSimulation
import os

folder = "raw_data1"
orig_wt = np.array([0.5679, 0.3123, 2.963, 7.407, 0.06173, 0.3086, 0.3951, 1.136])
epsilon = 1e-5

for fn in glob(folder + "/*"):
    sim = SingleSimulation(fn, orig_weight=orig_wt, down_sampling_ratio=0.001)

    new_wt = np.random.rand(9)
    sim.reweight(new_wt)
    curr_val = sim.value
    deriv_analytical = sim.deriv
    for i in xrange(9):
        test_wt = new_wt.copy()
        test_wt[i] += epsilon
        sim.reweight(test_wt)
        deriv_numerical = (sim.value - curr_val) / epsilon
        err = (deriv_numerical - deriv_analytical[i]) / deriv_analytical[i]
        print deriv_numerical, deriv_analytical[i], err
