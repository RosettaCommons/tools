import glob
out = open("reweight.job_list", 'w')

for f in glob.glob('raw_data/*'):
    out.write('./\tpython weight_g_eval.py %s\n' % f)

