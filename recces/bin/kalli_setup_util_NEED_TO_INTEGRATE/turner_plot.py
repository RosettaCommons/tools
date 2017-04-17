import numpy as np
import matplotlib.pyplot as plt
import turner_util


def turner_eval(seqs, filename='freeE_orig.npy'):
    freeE_list = []
    for seq in seqs:
        kT, freeE = turner_util.get_delta_G(*seq, filename=filename)
        freeE_list.append(freeE)
    freeE_list = np.column_stack(freeE_list)
    return freeE_list[70]

seqs = [
    ('aa', 'uu'), ('au', 'au'), ('ac', 'gu'), ('ag', 'cu'), ('ua', 'ua'),
    ('uc', 'ga'), ('ug', 'ca'), ('cc', 'gg'), ('cg', 'cg'), ('gc', 'gc')]
G = turner_eval(seqs)
expt = np.array([
    -0.93, -1.10, -2.24, -2.08, -1.33, -2.35, -2.11, -3.26, -2.36, -3.42])
print np.sqrt(np.sum((G - expt) ** 2) / G.shape[0])
plt.plot(G, expt, 'ro')


seqs1 = []
seqs2 = []
for i in 'aucg':
    for j in 'aucg':
        seqs1.append((i + j, turner_util.convert_wcna(i)))
        seqs2.append((turner_util.convert_wcna(j), i + j))
seqs = seqs1 + seqs2
expt = np.array([
    -0.8, -0.5, -0.8, -0.6, -1.7, -0.8, -1.7, -1.2, -1.1, -0.4, -1.3,
    -0.6, -0.7, -0.1, -0.7, -0.1, -0.3, -0.3, -0.4, -0.2, -0.5, -0.3,
    -0.2, -0.1, -0.2, -0.3, 0, 0, -0.3, -0.1, -0.2, -0.2])
G = turner_eval(seqs)
print np.sqrt(np.sum((G - expt) ** 2) / G.shape[0])
plt.plot(G, expt, 'mx')


seqs = [
    ('ag', 'uu'), ('au', 'gu'), ('cg', 'ug'), ('cu', 'gg'), ('gg', 'uc'),
    ('gu', 'gc'), ('ga', 'uu'), ('ug', 'ua')]
expt = np.array([
    -0.55, -1.36, -1.41, -2.11, -1.53, -2.51, -1.27, -1.00])
G = turner_eval(seqs)
print np.sqrt(np.sum((G - expt) ** 2) / G.shape[0])
plt.plot(G, expt, 'bs')

plt.plot([-3.8, 0.3], [-3.8, 0.3], 'k-')
plt.ylabel('EXPT')
plt.xlabel('PRED')
plt.show()
