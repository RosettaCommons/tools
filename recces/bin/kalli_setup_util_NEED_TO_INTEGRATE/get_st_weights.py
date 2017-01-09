import turner_util_new as util
from glob import glob
import os.path as op

out = open('st_weights.txt', 'w')

for folder in glob('raw_data/*'):
    seq = op.basename(folder)
    prerun_folder = op.join(folder, 'prerun')
    scores_file = op.join(prerun_folder, 'prerun_hist_scores.gz')
    temp, weights = util.weight_evaluate(prerun_folder, scores_file)
    out_str = [seq]
    for w in weights:
        out_str.append(repr(w))
    out.write(' '.join(out_str) + '\n')
out.close()
