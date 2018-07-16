import sys
import data
import numpy as np

folder = sys.argv[1]

outfile =folder + "/g_list.out"
orig_wt = np.array([0.3651, 0.2008, 0.003492, 0.00001, 9.5238, 3.8095, 0.5079, 0.3968, 1.746, 1.5873])

sim = data.SingleSimulation(folder, orig_weight=orig_wt)

scores = []
weights = []
for line in open("min_result.out"):
    a, b = line.strip().split('; ')
    scores.append(float(a))
    weights.append([float(i) for i in b.split(', ')])

g_list = []
for wt in weights:
    sim.reweight(wt)
    g_list.append(sim.value)

np.savetxt(outfile, g_list)
