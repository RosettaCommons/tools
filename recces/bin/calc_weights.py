#!/usr/bin/env python

from recces.util import weight_evaluate
from os.path import exists

hist_scores_file = 'prerun_hist_scores.gz'
if not exists( hist_scores_file): hist_scores_file = 'ST_hist_scores.gz'
assert( exists( hist_scores_file ) )

(temps, wts) = weight_evaluate('./', hist_scores_file )

for wt in wts: print wt,
print
