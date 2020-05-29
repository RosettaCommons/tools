import os
import struct
import gzip
import glob
import math
import itertools
import warnings

import numpy as np
import scipy.optimize


# Canonical AUCG RNA bases and their base-pair-partners
CANONICAL_RNA_BP = {'a': 'u', 'g': 'c', 'u': 'a', 'c': 'g'}

# Ending base-pair used in the computation
END_BP = 'g', 'c'

# Experiment measured free energies
DG_CANONICAL = {('aa', 'uu'): -0.93,
                ('au', 'au'): -1.1,
                ('ac', 'gu'): -2.24,
                ('ag', 'cu'): -2.08,
                ('ua', 'ua'): -1.33,
                ('uc', 'ga'): -2.35,
                ('ug', 'ca'): -2.11,
                ('cc', 'gg'): -3.26,
                ('cg', 'cg'): -2.36,
                ('gc', 'gc'): -3.42}

DG_DANGLING = {('au', 'u'): -0.5,
               ('ug', 'a'): -1.2,
               ('uu', 'a'): -0.8,
               ('cg', 'g'): -0.6,
               ('u', 'aa'): -0.3,
               ('ag', 'u'): -0.6,
               ('ga', 'c'): -0.7,
               ('cc', 'g'): -1.3,
               ('a', 'uu'): -0.3,
               ('c', 'gg'): -0.2,
               ('c', 'ag'): -0.2,
               ('u', 'ga'): -0.3,
               ('aa', 'u'): -0.8,
               ('gg', 'c'): -0.1,
               ('u', 'ca'): -0.2,
               ('u', 'ua'): -0.5,
               ('ua', 'a'): -1.7,
               ('a', 'cu'): -0.3,
               ('g', 'cc'): 0,
               ('gc', 'c'): -0.7,
               ('g', 'gc'): -0.2,
               ('gu', 'c'): -0.1,
               ('g', 'uc'): -0.2,
               ('a', 'gu'): -0.1,
               ('c', 'ug'): -0.1,
               ('uc', 'a'): -1.7,
               ('a', 'au'): -0.3,
               ('ca', 'g'): -1.1,
               ('ac', 'u'): -0.8,
               ('g', 'ac'): -0.4,
               ('c', 'cg'): 0,
               ('cu', 'g'): -0.4}

DG_GU = {('ag', 'uu'): -0.35,
         ('au', 'gu'): -0.90,
         ('cg', 'ug'): -1.25,
         ('cu', 'gg'): -1.77,
         ('gg', 'uc'): -1.80,
         ('gu', 'gc'): -2.15,
         ('ga', 'uu'): -0.51,
         ('ug', 'ua'): -0.39}

DG_IGIC = {('cZ[IGU]', 'Z[ICY]g'): -2.46,
           ('Z[ICY]Z[IGU]', 'Z[ICY]Z[IGU]'): -2.45,
           ('cZ[ICY]', 'Z[IGU]g'): -3.46,
           ('gZ[IGU]', 'Z[ICY]c'): -3.07,
           ('Z[IGU]Z[IGU]', 'Z[ICY]Z[ICY]'): -3.30,
           ('gZ[ICY]', 'Z[IGU]c'): -4.00,
           ('Z[IGU]Z[ICY]', 'Z[IGU]Z[ICY]'): -4.61}

DG_ALL = dict(DG_CANONICAL.items() + DG_DANGLING.items() +
              DG_GU.items() + DG_IGIC.items())

DG_TERMINAL = {('a', 'u'): 0.45,
               ('g', 'u'): 0,
               ('Z[ICY]', 'Z[IGU]'): -0.19}


def get_base_cmdline(save_terms=True, save_scores=True,
                     weight_file='turner.wts'):
    """Return the base Rosetta command line for Turner rule computation."""
    cmdline_common = "recces_turner "
    #cmdline_common = "turner_new "
    cmdline_common += "-rna::corrected_geo "
    cmdline_common += "-score:rna_torsion_potential RNA11_based_new "
    cmdline_common += "-chemical::enlarge_H_lj "
    cmdline_common += "-analytic_etable_evaluation false "
    if save_terms or save_scores:
         cmdline_common += "-save_score_terms "
    #if save_terms:
    #    cmdline_common += "-save_terms "
    #if save_scores:
    #    cmdline_common += "-save_scores "
    cmdline_common += "-score:weights %s " % os.path.abspath(weight_file)
    return cmdline_common


def design_seq(base1, base2=None):
    """Design sequences to be simulated for a base/base-pair.

    Parameters
    ----------
    base1 : Name of the first base.
    base2 : Name of the second base paired with base1. If not given base1 is
            assumed to be and dangling base.

    Returns
    -------
    one_strand, duplex : Set of single- and double-strand sequences needed
                         to be simulated.
    """
    one_strand = set()
    duplex = set()

    def add_seq_to_set(seq_set, seq1, seq2=None):
        if seq2 is None:
            seq_set.add(''.join(seq1))
        else:
            seq_set.add((''.join(seq1), ''.join(seq2)))

    def add_seq_duplex(seq1, seq2):
        add_seq_to_set(one_strand, seq1[:-1])
        add_seq_to_set(one_strand, seq2[1:])
        add_seq_to_set(duplex, seq1[:-1], seq2[1:])
        add_seq_to_set(one_strand, seq1)
        add_seq_to_set(one_strand, seq2)
        add_seq_to_set(duplex, seq1, seq2)

    if base2 is None:  # Dangling end
        for na1, na2 in CANONICAL_RNA_BP.iteritems():
            seq1 = END_BP[0], na1, base1
            seq2 = na2, END_BP[1]
            add_seq_to_set(one_strand, seq1)
            add_seq_to_set(one_strand, seq2)
            add_seq_to_set(duplex, seq1, seq2)
            add_seq_to_set(one_strand, seq1[:-1])
            add_seq_to_set(duplex, seq1[:-1], seq2)

            seq1 = END_BP[0], na1
            seq2 = base1, na2, END_BP[1]
            add_seq_to_set(one_strand, seq1)
            add_seq_to_set(one_strand, seq2)
            add_seq_to_set(duplex, seq1, seq2)
    else:  # Base-pair
        seq1 = END_BP[0], base1, base1
        seq2 = base2, base2, END_BP[1]
        add_seq_duplex(seq1, seq2)

        seq1 = END_BP[0], base2, base2
        seq2 = base1, base1, END_BP[1]
        add_seq_duplex(seq1, seq2)

        seq1 = END_BP[0], base1, base2
        seq2 = base1, base2, END_BP[1]
        add_seq_duplex(seq1, seq2)

        seq1 = END_BP[0], base2, base1
        seq2 = base2, base1, END_BP[1]
        add_seq_duplex(seq1, seq2)
        for na1, na2 in CANONICAL_RNA_BP.iteritems():
            seq1 = END_BP[0], na1, base1
            seq2 = base2, na2, END_BP[1]
            add_seq_duplex(seq1, seq2)

            seq1 = END_BP[0], na1, base2
            seq2 = base1, na2, END_BP[1]
            add_seq_duplex(seq1, seq2)

            seq1 = END_BP[0], base1, na1
            seq2 = na2, base2, END_BP[1]
            add_seq_duplex(seq1, seq2)

            seq1 = END_BP[0], base2, na1
            seq2 = na2, base1, END_BP[1]
            add_seq_duplex(seq1, seq2)
    return one_strand, duplex


def design_seq_from_list(seq_list):
    """Design sequences for simulation from list of sequences.

    Parameters
    ----------
    seq_list : List of sequences. Elements must be list/tuple of bases,
               or single strings for dangling end, e.g. [('A', 'U'),
               'Z[IGU]', ('Z[INO]',)].

    Returns
    -------
    one_strand, duplex : Set of single- and double-strand sequences needed
                         to be simulated.
    """
    one_strand = set()
    duplex = set()
    for seq in seq_list:
        if isinstance(seq, str):
            base1, base2 = seq, None
        elif len(seq) == 1:
            base1, base2 = seq[0], None
        else:
            base1, base2 = seq
        one_strand_new, duplex_new = design_seq(base1, base2)
        one_strand |= one_strand_new
        duplex |= duplex_new
    return one_strand, duplex


def design_seq_from_file(filename):
    """Design sequences for simulation from an input file.

    Parameters
    ----------
    filename : Name of the input file. Each line contains the name of one or
               two bases, separated by space.

    Returns
    -------
    one_strand, duplex : Set of single- and double-strand sequences needed
                         to be simulated.
    """
    seq_list = [line.split() for line in open(filename)]
    one_strand, duplex = design_seq_from_list(seq_list)
    return one_strand, duplex


def load_2d_bin_gz(filename, dtype=np.float32):
    """Reads *.bin.gz file output by Rosetta into a 2d numpy array.

    Parameters
    ----------
    filename : Name of input file.
    dtype : Data type of the array. Default to np.float32, corresponds to
            float in C.

    Returns
    -------
    data : 2D numpy array.
    """
    all_data = gzip.open(filename, 'rb').read()
    dim = struct.unpack('QQ', all_data[0:16])
    data = np.fromstring(all_data[16:], dtype=dtype).reshape(dim)
    return data


def load_1d_bin_gz(filename, dtype=np.float64):
    """Reads *.bin.gz file output by Rosetta into a 1d numpy array.

    Parameters
    ----------
    filename : Name of input file.
    dtype : Data type of the array. Default to np.float32, corresponds to
            float in C.

    Returns
    -------
    data : 1d numpy array.
    """
    all_data = gzip.open(filename, 'rb').read()
    data = np.fromstring(all_data, dtype=dtype)
    return data


def weight_evaluate(folder, hist_score):
    """Evaluate the weight of simulated tempering from prerun data.

    Parameters
    ----------
    folder : Folder name containing the prerun histogram data.
    hist_score : File name storing the score for each histogram bin.

    Returns
    -------
    kT_list : List of kT in Rosetta unit.
    weight_list : List of weights correspond to each kT
    """
    scores = load_1d_bin_gz(hist_score, dtype=np.float64)

    def get_hist(filename):
        hist = load_1d_bin_gz(filename, dtype=np.uint64)
        return np.column_stack((scores, hist))

    working_dir = os.getcwd()
    os.chdir(folder)
    file_list = glob.glob('*.hist.gz')
    data_list = [
        [name, float(name.split('_')[-1].replace('.hist.gz', ''))]
        for name in file_list]
    data_list = sorted(data_list, key=lambda x: x[1])
    data_list = zip(*data_list)
    kT_list = data_list[1]
    name_list = data_list[0]
    hist_list = map(get_hist, name_list)
    weight_list = [0]

    for hist1, kT1, hist2, kT2 in itertools.izip(
            hist_list[:-1], kT_list[:-1], hist_list[1:], kT_list[1:]):
        delta_weight = get_ST_delta(hist1, hist2, kT1, kT2)
        weight_list.append(weight_list[-1] + delta_weight)

    os.chdir(working_dir)
    return kT_list, weight_list


def get_ST_delta(hist1, hist2, kt1, kt2):
    """Compute the optimal weight difference for simulated tempering.

    Parameters
    ----------
    hist1 : Histogram for temperature 1
    hist2 : Histogram for temperature 2
    kt1 : kT for temperature 1
    kt2 : kT for temperature 2

    Returns
    -------
    delta_final : Optimal simulated tempering weight delta.
    """
    beta1 = 1.0 / kt1
    beta2 = 1.0 / kt2
    sum1 = np.sum(hist1[:, 1])
    sum2 = np.sum(hist2[:, 1])

    log_prob_base1 = -(beta2 - beta1) * hist1[:, 0]
    log_prob_base2 = -(beta1 - beta2) * hist2[:, 0]

    def get_acpt_rate(delta, log_prob_base, hist, sum_hist):
        log_prob = log_prob_base + delta
        log_prob[log_prob > 0] = 0
        return np.sum(np.exp(log_prob) * hist[:, 1]) / sum_hist

    def target_func(delta):
        acpt_rate1 = get_acpt_rate(delta, log_prob_base1, hist1, sum1)
        acpt_rate2 = get_acpt_rate(-delta, log_prob_base2, hist2, sum2)
        return acpt_rate1 - acpt_rate2

    # Simple heuristic for the init weight based on avg. energy
    E1 = np.sum(hist1[:, 0] * hist1[:, 1]) / sum1
    E2 = np.sum(hist2[:, 0] * hist2[:, 1]) / sum2
    delta_start = (beta2 - beta1) * (E1 + E2) / 2
    search_width = 5
    delta_range = (delta_start - search_width), (delta_start + search_width)

    delta_final = scipy.optimize.brentq(
        target_func, delta_range[0], delta_range[1])
    acpt_rate = get_acpt_rate(delta_final, log_prob_base1, hist1, sum1)
    if acpt_rate < 0.1:
        warnings.warn(" Acceptance rate (%s) lower than 0.1" % acpt_rate)
    return delta_final


def seq_parse(seq):
    """Parse Rosetta style sequence into list

    Parameters
    ----------
    seq : Rosetta-style sequence, e.g. auZ[IGU]g

    Returns
    -------
    parsed_seq : parsed sequence in Python tuple, e.g. ('a', 'u', 'Z[IGU]',
    'g')
    """
    in_brac = False
    parsed_seq = []
    if not seq:
        return parsed_seq
    i = 0
    for j, c in enumerate(seq):
        if j == 0:
            continue
        if in_brac:
            if c == ']':
                in_brac = False
        elif c == '[':
            in_brac = True
        else:
            parsed_seq.append(seq[i:j])
            i = j
    parsed_seq.append(seq[i:])
    return tuple(parsed_seq)


def torsion_volume(seq1, seq2='', aform_torsion_range=2 * math.pi / 3):
    """Compute phase space volume of a sequence.

    Parameters
    ----------
    seq1 : Sequence 1 for duplex or the full sequence for single strand.
    seq2 : Sequence 2 for duplex. Skip it for single strand.
    aform_torsion_range : Range for A form torsion angles, for duplex only.
                          Default value is 120 degrees (+/- 60).

    Returns
    -------
    Phase space volume of the sequence.
    """
    from math import pi

    def seq_len(seq):
        return len(seq_parse(seq))

    len1 = seq_len(seq1)
    len2 = seq_len(seq2)
    min_len = min(len1, len2)
    diff_len = abs(len1 - len2)
    if min_len == 0:  # One-strand
        #print "applying phase space volume for single strand with length ", diff_len
        return (2 * pi) ** (6 * diff_len - 5) * (2 ** diff_len)
    else:
        #print "applying phase space volume for two strand with helix length ", min_len, " and dangle length ", diff_len
        volume = aform_torsion_range ** (12 * min_len - 10)
        volume *= (2 * pi) ** (6 * diff_len) * (2 ** diff_len)
        return volume

def temp_list():
    return [0.8, 1, 1.4, 1.8, 3, 7, 30]
