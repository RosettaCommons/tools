import numpy as np

# got this from running analyze_04-26.py
#results = np.genfromtxt('results_stripped_no_elec.txt')
results = np.genfromtxt('results_stripped.txt')
expt = np.array([-0.93, -1.10, -2.24, -2.08, -1.33, -2.35, -2.11, -3.26, -2.36, -3.42])
rmse = np.sqrt(np.sum((results - expt) ** 2)/expt.size)
print rmse
