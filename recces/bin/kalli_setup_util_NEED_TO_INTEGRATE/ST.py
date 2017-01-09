import os
from os.path import abspath
import turner_util_new as util
import pickle

n_repeat = 1
job_file = open("ST.job_list", 'w')

n_cycle_inf = 300000
n_cycle_ST = 9000000
cmdline_common = util.get_base_cmdline()
one_strand, duplex = pickle.load(open("prerun.seq"))
kTs = [0.8, 1.0, 1.4, 1.8, 3.0, 7.0, 30.0]
seq_weights = {}
for line in open('st_weights.txt'):
    seq, weights = line.strip().split(' ', 1)
    seq_weights[seq] = weights

root_path = abspath('./raw_data/')
data_folder = abspath('./')
os.chdir(root_path)

for seq in one_strand:
    os.chdir(seq)

    cmdline1 = cmdline_common
    cmdline1 += '-n_cycle %d ' % n_cycle_inf
    cmdline1 += "-seq1 %s " % seq
    cmdline1 += "-temps -1 "
    cmdline1 += "-out_prefix " + abspath("./kT_inf")
    job_file.write(cmdline1 + '\n')
    #job_file.write("./\t" + cmdline1 + '\n')

    cmdline2 = cmdline_common
    cmdline2 += '-n_cycle %d ' % n_cycle_ST
    cmdline2 += "-seq1 %s " % seq

    weights = seq_weights[os.path.basename(seq)]
    cmdline2 += "-temps "
    for kT in kTs:
        cmdline2 += '%s ' % kT
    cmdline2 += '-st_weights '
    cmdline2 += '%s ' % weights

    cmdline2 += "-out_prefix "
    for i in xrange(n_repeat):
        cmdline = cmdline2 + abspath("./ST_%d" % i)
        job_file.write(cmdline + '\n')
        #job_file.write("./\t" + cmdline + '\n')

    os.chdir(root_path)

for seq1, seq2 in duplex:
    seq = '%s_%s' % (seq1, seq2)
    os.chdir(seq)

    cmdline0 = cmdline_common
    cmdline0 += "-seq1 %s " % seq1
    cmdline0 += "-seq2 %s " % seq2

    cmdline1 = cmdline0 + '-n_cycle %d ' % n_cycle_inf
    cmdline1 += "-temps -1 "
    cmdline1 += "-out_prefix " + abspath("./kT_inf")
    job_file.write(cmdline1 + '\n')
    #job_file.write("./\t" + cmdline1 + '\n')

    cmdline2 = cmdline0 + '-n_cycle %d ' % n_cycle_ST

    weights = seq_weights[seq]
    cmdline2 += "-temps "
    for kT in kTs:
        cmdline2 += '%s ' % kT
    cmdline2 += '-st_weights '
    cmdline2 += '%s ' % weights

    cmdline2 += "-out_prefix "
    for i in xrange(n_repeat):
        cmdline = cmdline2 + abspath("./ST_%d" % i)
        job_file.write(cmdline + '\n')
        #job_file.write("./\t" + cmdline + '\n')

    os.chdir(root_path)
