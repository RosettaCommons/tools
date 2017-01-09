import turner_util
from glob import glob


seq_list = glob('raw_data/*')
for seq in seq_list:
    print seq
    turner_util.folder_analyze(seq)
