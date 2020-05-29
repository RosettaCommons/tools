import os
import struct
import numpy as np
import gzip
import glob
import math
import sys
import itertools


def convert_wcna(base):
    if base == 'a':
        return 'u'
    if base == 'u':
        return 'a'
    if base == 'g':
        return 'c'
    if base == 'c':
        return 'g'


def get_cmdline():
    cmdline_common = "turner_new "
    cmdline_common += "-rna::corrected_geo "
    cmdline_common += "-save_terms "
    cmdline_common += "-save_scores "
    cmdline_common += "-score:weights %s " % os.path.abspath('turner.wts')
    cmdline_common += "-score:rna_torsion_potential RNA11_based_new "
    cmdline_common += "-chemical::enlarge_H_lj "
    cmdline_common += "-analytic_etable_evaluation false "
    return cmdline_common


def design_seq(base1, base2=None):
    one_strand = set()
    duplex = set()
    na_list = 'aucg'

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

    # Always use GC as first base-pair
    if base2 is None:  # Dangling end
        for na in na_list:
            seq1 = ['g', na, base1]
            seq2 = [convert_wcna(na), 'c']
            add_seq_to_set(one_strand, seq1)
            add_seq_to_set(one_strand, seq2)
            add_seq_to_set(duplex, seq1, seq2)
            add_seq_to_set(one_strand, seq1[:-1])
            add_seq_to_set(duplex, seq1[:-1], seq2)

            seq1 = ['g', na]
            seq2 = [base1, convert_wcna(na), 'c']
            add_seq_to_set(one_strand, seq1)
            add_seq_to_set(one_strand, seq2)
            add_seq_to_set(duplex, seq1, seq2)
    else:  # Base-pair
        seq1 = ['g', base1, base1]
        seq2 = [base2, base2, 'c']
        add_seq_duplex(seq1, seq2)

        seq1 = ['g', base2, base2]
        seq2 = [base1, base1, 'c']
        add_seq_duplex(seq1, seq2)

        seq1 = ['g', base1, base2]
        seq2 = [base1, base2, 'c']
        add_seq_duplex(seq1, seq2)

        seq1 = ['g', base2, base1]
        seq2 = [base2, base1, 'c']
        add_seq_duplex(seq1, seq2)
        for na in na_list:
            seq1 = ['g', na, base1]
            seq2 = [base2, convert_wcna(na), 'c']
            add_seq_duplex(seq1, seq2)

            seq1 = ['g', na, base2]
            seq2 = [base1, convert_wcna(na), 'c']
            add_seq_duplex(seq1, seq2)

            seq1 = ['g', base1, na]
            seq2 = [convert_wcna(na), base2, 'c']
            add_seq_duplex(seq1, seq2)

            seq1 = ['g', base2, na]
            seq2 = [convert_wcna(na), base1, 'c']
            add_seq_duplex(seq1, seq2)
    return one_strand, duplex


def design_seq_from_file(filename):
    one_strand = set()
    duplex = set()
    for line in open(filename):
        seqs = line.split()
        a, b = design_seq(*seqs)
        one_strand |= a
        duplex |= b
    return one_strand, duplex


def temp_list():
    return [0.8, 1, 1.4, 1.8, 3, 7, 30]


def load_bin_gz(filename, dtype=np.float32):
    all_data = gzip.open(filename, 'rb').read()
    dim = struct.unpack('QQ', all_data[0:16])
    data = np.fromstring(all_data[16:], dtype=dtype).reshape(dim)
    return data


def load_1d_bin_gz(filename, dtype=np.float32):
    all_data = gzip.open(filename, 'rb').read()
    data = np.fromstring(all_data, dtype=dtype)
    return data


def weight_evaluate(folder):
    def weight_optimize(hist1, hist2, beta1, beta2):
        sum1 = np.sum(hist1[:, 1])
        sum2 = np.sum(hist2[:, 1])
        E1 = np.sum(hist1[:, 0] * hist1[:, 1]) / sum1
        E2 = np.sum(hist2[:, 0] * hist2[:, 1]) / sum2

        # Simple heuristic for the init weight based on avg. energy
        delta_weight_start = (beta2 - beta1) * (E1 + E2) / 2

        log_prob_base1 = -(beta2 - beta1) * hist1[:, 0]
        log_prob_base2 = -(beta1 - beta2) * hist2[:, 0]

        def get_acpt_rate(delta_weight, log_prob_base, hist, sum_hist):
            log_prob = log_prob_base + delta_weight
            log_prob[log_prob > 0] = 0
            return np.sum(np.exp(log_prob) * hist) / sum_hist

        def get_score(r1, r2):
            return np.fmin(r1, r2)

        def best_weight(delta_weights):
            acpt_rate1 = np.array([
                get_acpt_rate(weight, log_prob_base1, hist1[:, 1], sum1)
                for weight in delta_weights])
            acpt_rate2 = np.array([
                get_acpt_rate(-weight, log_prob_base2, hist2[:, 1], sum2)
                for weight in delta_weights])
            score = get_score(acpt_rate1, acpt_rate2)
            idx = np.argmax(score)
            lowest_weight = delta_weights[idx]
            acpt_rate = math.sqrt(acpt_rate1[idx] * acpt_rate2[idx])
            return lowest_weight, acpt_rate

        # Do a coarse grid search
        delta_weights = np.arange(
            delta_weight_start - 5, delta_weight_start + 5, 0.1)
        lowest_weight, acpt_rate = best_weight(delta_weights)

        # Now do a finer search using the best coarse result
        delta_weights = np.arange(
            lowest_weight - 0.5, lowest_weight + 0.5, 0.01)
        lowest_weight, acpt_rate = best_weight(delta_weights)

        if acpt_rate < 0.1:
            sys.stderr.write("Warning: Avg accept rate < 0.1!!!\n")
            sys.stderr.write(
                "T1: %s, T2: %s, rate: %s\n" %
                (1 / beta1, 1 / beta2, acpt_rate))
        return lowest_weight
    ####################################
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
        delta_weight = weight_optimize(hist1, hist2, 1.0 / kT1, 1.0 / kT2)
        weight_list.append(weight_list[-1] + delta_weight)

    os.chdir(working_dir)
    return kT_list, weight_list


def get_hist(filename):
    scores = load_1d_bin_gz('/work/01937/fcchou/projects/turner/run5/hist_scores.gz', dtype=np.float64)
    hist = load_1d_bin_gz(filename, dtype=np.uint64)
    print hist, scores
    return np.column_stack((scores, hist))


def get_hist_from_folder(folder):
    cwd = os.getcwd()
    os.chdir(folder)
    file_list = glob.glob('*.hist.gz')
    file_list = [
        [name, float(name.split('_')[-1].replace('.hist.gz', ''))]
        for name in file_list]
    kT_list = []
    for i in file_list:
        if i[1] not in kT_list:
            kT_list.append(i[1])

    kT_list = np.sort(kT_list)
    score_list = load_1d_bin_gz('/home/fcchou/projects/turner_rule/run5/hist_scores.gz', dtype=np.float64)
    n_bin = len(score_list)
    n_kT = kT_list.shape[0]
    hist_list = np.zeros((n_kT, n_bin))

    for filename, kT in file_list:
        hist = load_1d_bin_gz(filename, dtype=np.uint64)
        hist_list[np.searchsorted(kT_list, kT)] += hist
    N_list = np.sum(hist_list, axis=1)

    os.chdir(cwd)

    kT_list[kT_list < 0] = np.inf
    return kT_list, N_list, score_list, hist_list


def wham(kT_list, N_list, scores, hist_list):
    hist_sum = np.sum(hist_list, axis=0)

    nonzero_idx = hist_sum.nonzero()[0]
    first_non_zero = nonzero_idx[0]
    last_non_zero = nonzero_idx[-1] + 1

    hist_sum = hist_sum[first_non_zero:last_non_zero]
    scores = scores[first_non_zero:last_non_zero]

    DoS = np.zeros_like(hist_sum)
    norm_factor = np.zeros_like(DoS)
    exp_f_list_pre = np.ones_like(N_list)
    exp_f_list = np.ones_like(N_list)
    bin_size = scores[1] - scores[0]
    n_kT = kT_list.shape[0]
    boltzmann_factor = np.exp(-np.tile(scores, (n_kT, 1)).T / kT_list)

    is_converge = False
    converge_limit = 0.00000001

    while not is_converge:
        norm_factor = np.sum(
            N_list * bin_size * exp_f_list * boltzmann_factor, axis=1)
        DoS = hist_sum / norm_factor
        exp_f_list = 1.0 / np.sum(DoS * bin_size * boltzmann_factor.T, axis=1)
        diff = np.abs(
            (exp_f_list - exp_f_list_pre) /
            (exp_f_list + exp_f_list_pre) * 0.5)
        if np.all(diff <= converge_limit):
            is_converge = True
        else:
            exp_f_list_pre = exp_f_list
    return DoS, norm_factor, scores


def compute_freeE(DoS, scores, total_volume):
    kT_in_kcal = 0.593
    kT_list = np.arange(0.3, 2, 0.01)
    n_kT = kT_list.shape[0]
    boltzmann_factor = np.exp(-np.tile(scores, (n_kT, 1)).T / kT_list).T
    DoS_volume = np.sum(DoS)

    Z_list = total_volume * np.sum(DoS / DoS_volume * boltzmann_factor, axis=1)
    A_list = -kT_in_kcal * np.log(Z_list)
    E_list = (
        np.sum(DoS * scores * boltzmann_factor, axis=1) /
        np.sum(DoS * boltzmann_factor, axis=1) / kT_list * kT_in_kcal)
    minus_TS_list = A_list - E_list
    return kT_list, A_list, E_list, minus_TS_list


def compute_folder_freeE(DoS_filename, folder):
    data = np.load(os.path.join(folder, DoS_filename))
    kT_list, A_list, E_list, minus_TS_list = compute_freeE(
        data[:, 1], data[:, 0],
        torsion_volume(*os.path.basename(folder).split('_')))
    np.save(
        '%s/%s' % (folder, DoS_filename.replace('DoS', 'freeE')),
        np.column_stack((kT_list, A_list, E_list, minus_TS_list)))


def folder_analyze(folder):
    kT_list, N_list, scores, hist_list = get_hist_from_folder(folder)
    for kT, hist in itertools.izip(kT_list, hist_list):
        nonzero_idx = hist.nonzero()[0]
        first = nonzero_idx[0]
        last = nonzero_idx[-1] + 1
        hist_final = np.column_stack((
            scores[first:last], hist[first:last]))
        np.save('%s/hist_raw_kT_%.3f.npy' % (folder, kT), hist_final)
    DoS, norm_factor, scores = wham(kT_list, N_list, scores, hist_list)
    np.save('%s/DoS_orig.npy' % folder, np.column_stack((scores, DoS)))
    np.save(
        '%s/WHAM_norm_factors.npy' % folder,
        np.column_stack((scores, norm_factor)))
    kT_list, A_list, E_list, minus_TS_list = compute_freeE(
        DoS, scores, torsion_volume(*os.path.basename(folder).split('_')))
    np.save(
        '%s/freeE_orig.npy' % folder,
        np.column_stack((kT_list, A_list, E_list, minus_TS_list)))


def seq_parse(seq):
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
    return parsed_seq


def torsion_volume(seq1, seq2=''):
    from math import pi

    def seq_len(seq):
        return len(seq_parse(seq))
    len1 = seq_len(seq1)
    len2 = seq_len(seq2)
    min_len = min(len1, len2)
    diff_len = abs(len1 - len2)
    # print seq1, seq2, min_len, diff_len
    if min_len == 0:  # One-strand
        return (2 * pi) ** (6 * diff_len - 5) * (2 ** diff_len)
    else:
        volume = (2 * pi / 3) ** (12 * min_len - 10)
        volume *= (2 * pi) ** (6 * diff_len) * (2 ** diff_len)
        return volume


def get_delta_G(seq1, seq2, filename='freeE_orig.npy'):
    seq1_parsed = seq_parse(seq1)
    seq2_parsed = seq_parse(seq2)
    len1 = len(seq1_parsed)
    len2 = len(seq2_parsed)

    def freeE_from_seq(seq):
        data = np.load('raw_data/%s/%s' % (seq, filename))
        return data[:, 1]
    kT_list = np.load('raw_data/g%s/%s' % (seq1, filename))[:, 0]
    assert(0 < len1 <= 2 and 0 < len2 <= 2)
    if len1 != len2:  # Dangling end case
        dup1 = 'g%s_%sc' % (seq1, seq2)
        if len1 > len2:
            strand0 = ''.join(['g'] + seq1_parsed[:-1])
            strand1 = 'g' + seq1
            dup0 = '%s_%sc' % (strand0, seq2)
        else:
            strand0 = ''.join(seq2_parsed[1:] + ['c'])
            strand1 = seq2 + 'c'
            dup0 = 'g%s_%s' % (seq1, strand0)
        freeE = (
            freeE_from_seq(dup1) - freeE_from_seq(strand1) -
            freeE_from_seq(dup0) + freeE_from_seq(strand0))
    else:  # Duplex case
        def freeE_duplex():
            strand0 = ''.join(['g'] + seq1_parsed[:-1])
            strand1 = ''.join(seq2_parsed[1:] + ['c'])
            strand2 = 'g' + seq1
            strand3 = seq2 + 'c'
            dup0 = '%s_%s' % (strand0, strand1)
            dup1 = 'g%s_%sc' % (seq1, seq2)
            freeE = (
                freeE_from_seq(dup1) - freeE_from_seq(strand2) -
                freeE_from_seq(strand3) - freeE_from_seq(dup0) +
                freeE_from_seq(strand0) + freeE_from_seq(strand1))
            return freeE
        freeE = freeE_duplex()
        seq1_parsed, seq2_parsed = seq2_parsed, seq1_parsed
        seq1 = ''.join(seq1_parsed)
        seq2 = ''.join(seq2_parsed)
        freeE += freeE_duplex()
        freeE *= 0.5
    return kT_list, freeE


def turner_duplex_eval(filename='freeE_orig.npy'):
    seqs = [
        ('aa', 'uu'), ('au', 'au'), ('ac', 'gu'), ('ag', 'cu'), ('ua', 'ua'),
        ('uc', 'ga'), ('ug', 'ca'), ('cc', 'gg'), ('cg', 'cg'), ('gc', 'gc')]
    expt = np.array([
        -0.93, -1.10, -2.24, -2.08, -1.33, -2.35, -2.11, -3.26, -2.36, -3.42])
    freeE_list = []
    for seq in seqs:
        kT, freeE = get_delta_G(*seq, filename=filename)
        freeE_list.append(freeE)
    freeE_list = np.column_stack(freeE_list)
    rmse = np.sqrt(np.sum((freeE_list - expt) ** 2, axis=1) / expt.size)
    min_id = np.argmin(rmse)
    return kT[min_id], rmse[min_id], kT, rmse


def turner_dangle_eval(filename='freeE_orig.npy'):
    seqs1 = []
    seqs2 = []
    for i in 'aucg':
        for j in 'aucg':
            seqs1.append((i + j, convert_wcna(i)))
            seqs2.append((convert_wcna(j), i + j))
    seqs = seqs1 + seqs2
    expt = np.array([
        -0.8, -0.5, -0.8, -0.6, -1.7, -0.8, -1.7, -1.2, -1.1, -0.4, -1.3,
        -0.6, -0.7, -0.1, -0.7, -0.1, -0.3, -0.3, -0.4, -0.2, -0.5, -0.3,
        -0.2, -0.1, -0.2, -0.3, 0, 0, -0.3, -0.1, -0.2, -0.2])
    freeE_list = []
    for seq in seqs:
        kT, freeE = get_delta_G(*seq, filename=filename)
        freeE_list.append(freeE)
    freeE_list = np.column_stack(freeE_list)
    rmse = np.sqrt(np.sum((freeE_list - expt) ** 2, axis=1) / expt.size)
    min_id = np.argmin(rmse)
    return kT[min_id], rmse[min_id], kT, rmse
