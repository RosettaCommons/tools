#!/usr/bin/python

import sys
import os
import pfpd_const as protocol


def rescoring_flags(lowest_sc_struct):
    with open('rescore_flags', 'w') as rescore:
        rescore.write('-jd2:mpi_file_buf_job_distributor\n'
                      '-in:file:silent decoys.silent\n'
                      '-scorefile rescore.sc\n'
                      '-out:file:silent_struct_type binary\n'
                      '-out:file:silent decoys_rescored.silent\n'
                      '-flexPepDocking:flexpep_score_only\n'
                      '-use_input_sc\n'
                      '-unboundrot {receptor}\n'
                      '-native {nat}\n'
                      '-mute all\n-unmute protocols.flexPepDocking'.format(nat=lowest_sc_struct, receptor=receptor))


def rescoring():
    with open('score.sc', 'r') as score_file:
        scores = score_file.readlines()
        header = scores[1].split()  # scores[0] = SEQUENCE:
        scores = scores[2:]
        reweighted_column = header.index('reweighted_sc')
        description_column = header.index('description')

    sorted_sc = sorted(scores, key=lambda sc_line: float(sc_line.split()[reweighted_column]))
    lowest_sc_struct = sorted_sc[0].split()[description_column] + '.pdb'

    rescoring_flags(lowest_sc_struct)

    if sc_func == 'ref2015':
        os.system(protocol.FPD.format(flags='rescore_flags'))
    elif sc_func == 'talaris14':
        os.system(protocol.FPD_TALARIS.format(flags='rescore_flags'))


if __name__ == "__main__":
    sc_func = sys.argv[1]
    receptor = sys.argv[2]
    rescoring()
