from glob import glob

from data import SingleSimulation
import numpy as np

folder = "raw_data"
orig_wt = np.array([0.3651, 0.2008, 0.003492, 0.00001, 9.5238, 3.8095, 0.5079, 0.3968, 1.746, 1.5873])

for fn in glob(folder + "/*"):
    print fn
    sim = SingleSimulation(fn, use_existing_norm_factor=False, orig_weight=orig_wt)
    del sim
