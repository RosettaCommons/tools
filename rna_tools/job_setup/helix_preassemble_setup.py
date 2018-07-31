#!/usr/bin/env python
import subprocess
import argparse


parser = argparse.ArgumentParser(
    description='Setup scripts to preassemble helices for FARFAR.')
parser.add_argument(
    '-secstruct', required=True, help='Secondary structure file.')
parser.add_argument('-fasta', required=True, help='Fasta file.')
parser.add_argument(
    '-offset', type=int, default=0, help='Residue number offset')
parser.add_argument(
    '-out_prefix', default='helix', help='Prefix of the output scripts.')
parser.add_argument(
    '-out_cmdlines', default='CMDLINES',
    help='Output file for the cmdlines for running the jobs.')
args = parser.parse_args()


# Get the helix stems
def get_helix_stems(secstruct_file):
    def read_sectruct(secstruct_file):
        with open(secstruct_file) as f:
            secstruct = f.readline().strip()
        return secstruct

    def get_bp_list(secstruct):
        left_base_in_bp = []
        bp = []
        for i, char in enumerate(secstruct):
            if char == '(':
                left_base_in_bp.append(i + 1)
            elif char == ')':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        for i, char in enumerate(secstruct):
            if char == '[':
                left_base_in_bp.append(i + 1)
            elif char == ']':
                if not left_base_in_bp:
                    raise ValueError("Invalid secstruct!")
                bp.append((left_base_in_bp.pop(), i + 1))
        if left_base_in_bp:
            raise ValueError("Invalid secstruct!")
        return bp

    def get_stems(bp):
        stems = []
        stem = []
        for i, j in bp:
            if (not stem) or (i - 1 == stem[-1][0] and j + 1 == stem[-1][1]):
                stem.append((i, j))
            else:
                stems.append(stem)
                stem = [(i, j)]
        if stem:
            stems.append(stem)
        return stems

    secstruct = read_sectruct(secstruct_file)
    print secstruct
    bp = sorted(get_bp_list(secstruct))
    return get_stems(bp)

helix_stems = get_helix_stems(args.secstruct)
offset = args.offset

n_struct = 100
n_cycle = 1000
cmdline_base = bytearray(
    "rna_denovo -score:weights stepwise/rna/rna_res_level_energy4.wts -minimize_rna true -use_legacy_job_distributor true -restore_talaris_behavior true "
)
cmdline_base += "-nstruct %d " % n_struct
cmdline_base += "-cycles %d " % n_cycle
cmdline_base += "-fasta %s " % args.fasta
cmdline_base += "-secstruct_file %s " % args.secstruct
cmdline_base += "-offset %d " % offset

out = open(args.out_cmdlines, 'w')
print >> out, "# cmdlines:"
stem_idx = 0
res_list = bytearray()
out_silent_files = bytearray()
min_helix_len = 2
for stem in helix_stems:
    if len(stem) < min_helix_len:
        continue
    res_num = (
        '%d-%d ' % (stem[0][0] + offset, stem[-1][0] + offset) +
        '%d-%d' % (stem[-1][1] + offset, stem[0][1] + offset)
    )

    tag = "%s%d" % (args.out_prefix, stem_idx)
    out_script = "%s.RUN" % tag
    res_list += res_num + ' '
    out_silent_files += '%s.out ' % tag
    stem_idx += 1

    cmdline = cmdline_base[:]
    cmdline += "-working_res %s " % res_num
    #cmdline += "-out_script %s " % out_script
    cmdline += "-tag %s " % tag
    subprocess.check_call([str(elem) for elem in cmdline.split()])
    print >> out, "source %s" % out_script
print >> out, (
    "# cmdlines args for FARFAR:\n"
    "# -silent %s -input_silent_res %s"
    % (out_silent_files, res_list)
)
