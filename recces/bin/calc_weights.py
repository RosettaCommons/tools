#!/usr/bin/env python

from recces.util import weight_evaluate
(temps, wts) = weight_evaluate('./', 'prerun_hist_scores.gz')
for wt in wts: print wt,
print
