#!/usr/bin/python

import os
import pfpd_const as protocol

EXTRACT_MODEL = os.path.join(protocol.ROSETTA_BIN, 'extract_pdbs.linuxgccrelease') + \
                ' -database ' + protocol.ROSETTA_DB + ' @extract_flags > log'


def extracting_flags(lowest_sc_struct):
    with open('extract_flags', 'w') as extract:
        extract.write('-in:file:silent decoys.silent\n'
                      '-in:file:tags {}\n'.format(lowest_sc_struct))


def extract_pdb():
    with open('score.sc', 'r') as score_file:
        scores = score_file.readlines()
        header = scores[1].split()
        scores = scores[2:]  # SEQUENCE line + header
        reweighted_column = header.index('reweighted_sc')
        description_column = header.index('description')

    sorted_sc = sorted(scores, key=lambda sc_line: float(sc_line.split()[reweighted_column]))
    lowest_sc_struct = sorted_sc[0].split()[description_column]
    extracting_flags(lowest_sc_struct)
    os.system(EXTRACT_MODEL)

if __name__ == "__main__":
    extract_pdb()
