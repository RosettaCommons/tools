#!/bin/sh

cat subdirs/*/scorefile | sort | uniq > all_scores
awk '{print $59,$53}' all_scores | sort -n -k2 > allscores_nres
grep " 0.0" allscores_nres > zero_nres
wc -l all_scores zero_nres
echo "note all_scores is +3 from formatting"
