from glob import glob
import numpy as np
import sys
from data import SingleSimulation, TurnerRuleMinimizeFunc, load_raw_simulations
from scipy.optimize import minimize

folder = "raw_data"
orig_wt = np.array([0.3651, 0.2008, 0.003492, 0.00001, 9.5238, 3.8095, 0.5079, 0.3968, 1.746, 1.5873])

simulations = load_raw_simulations(folder, orig_weight=orig_wt, down_sampling_ratio=0.01)

turner_min = TurnerRuleMinimizeFunc(simulations, start_weight=orig_wt, dangling_weight=0.1)

start_weight = np.array([1, 1, 0.01, 1, 5, 5, 1, 1, 1, 1])
n_cycles = 100

bounds = [(0, 20) for _ in start_weight]
bounds[0] = bounds[1] = (0.1, 20)
bounds[2] = (0, 0.01)
bounds[4] = bounds[5] = (0.5, 20)

def min_func(wt):
    turner_min.reweight(wt)
    return turner_min.value, turner_min.deriv

out = open(sys.argv[1], 'w')

import time
for i in xrange(n_cycles):
    time0 = time.time()
    # Randomly start with 1/10 to 10x of the start_weight
    init_guess = start_weight * 10 ** np.random.uniform(-1, 1, len(start_weight))
    result = minimize(min_func, x0=init_guess, method="TNC", jac=True, bounds=bounds)
    min_weight = result.x
    min_val, _ = min_func(min_weight)
    out_str = repr(min_val) + "; "
    for i in min_weight:
        out_str += repr(i) + ", "
    out_str = out_str[:-2]
    out.write(out_str + '\n')
    out.flush()
    print time.time() - time0
out.close()
