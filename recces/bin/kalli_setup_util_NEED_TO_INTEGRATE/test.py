from glob import glob
import numpy as np
from data import SingleSimulation
import os

folder = "raw_data"
orig_wt = np.array([0.5679, 0.3123, 2.963, 7.407, 0.06173, 0.3086, 0.3951, 1.136])
epsilon = 1e-5

fn = folder + "/gcc_ggc"
# sim = SingleSimulation(fn, orig_weight=orig_wt, use_existing_norm_factor=False)
sim = SingleSimulation(fn, orig_weight=orig_wt, down_sampling_ratio=0.001, use_existing_norm_factor=False)
print sim.value

