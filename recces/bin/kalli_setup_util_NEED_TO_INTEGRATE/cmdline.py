out = open("turner_wt_min.job_list", 'w')

for i in xrange(100):
    cmdline = "./\tpython min_test1.py turner_min%s.out\n" % i
    out.write(cmdline)
out.close()
