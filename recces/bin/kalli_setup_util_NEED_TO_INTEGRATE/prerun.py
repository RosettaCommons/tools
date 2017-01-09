import turner_util_new as util
import os
from os.path import abspath, isdir, isfile, join
import pickle


job_file = open("prerun.job_list", 'w')

temps = util.temp_list()
print temps

one_strand, duplex = util.design_seq_from_file('seq_list')
used_strand = []
used_duplex = []

n_cycle = 300000
#cmdline_common = '/home/fcchou/Rosetta/main/source/bin/'
cmdline_common = util.get_base_cmdline()
cmdline_common += '-n_cycle %d ' % n_cycle
print cmdline_common

root_path = abspath('./raw_data/')
if not isdir(root_path):
    os.mkdir(root_path)
data_folder = abspath('./')
os.chdir(root_path)

for seq in one_strand:
    if isfile(join(seq, 'kT_inf_-1.00.bin.gz')):
        continue
    used_strand.append(seq)
    os.mkdir(seq)
    os.chdir(seq)
    os.mkdir('prerun')
    os.chdir('prerun')

    for T in temps:
        cmdline = cmdline_common
        cmdline += "-seq1 %s " % seq
        cmdline += "-temps %s " % T
        cmdline += "-out_prefix " + abspath("./prerun")
        job_file.write(cmdline + '\n')
        #job_file.write('./\t' + cmdline + '\n')
    os.chdir(root_path)

for seq1, seq2 in duplex:
    folder_name = '%s_%s' % (seq1, seq2)
    if isfile(join(folder_name, 'kT_inf_-1.00.bin.gz')):
        continue
    used_duplex.append((seq1, seq2))
    os.mkdir(folder_name)
    os.chdir(folder_name)
    os.mkdir('prerun')
    os.chdir('prerun')

    for T in temps:
        cmdline = cmdline_common
        cmdline += "-seq1 %s " % seq1
        cmdline += "-seq2 %s " % seq2
        cmdline += "-temps %s " % T
        cmdline += "-out_prefix " + abspath("./prerun")
        job_file.write(cmdline + '\n')
        #job_file.write('./\t' + cmdline + '\n')
    os.chdir(root_path)
job_file.close()

pickle.dump((used_strand, used_duplex), open("../prerun.seq", 'w'))
