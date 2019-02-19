from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap
from scipy.stats import gaussian_kde
from functools import partial

random.seed(42)

abundance_range, total_partners = load(open('data_stats_dump.dmp', 'r'))
ploidy_vs_size = load(open('ploidy_vs_size.dmp'))

abundance_range = abundance_range[abundance_range > 1]
total_partners = np.array(total_partners)
total_partners = total_partners[total_partners > 1]
total_partners = total_partners[total_partners < 45]
total_partners = total_partners.tolist()
total_partners_old = total_partners

total_partners = np.random.gamma(1.538, 2.016, len(total_partners))
total_partners = np.ceil(total_partners).astype(np.int16)
total_partners = total_partners[total_partners > 1]
total_partners_new = total_partners

alpha_factor = 1.0

deviation_from_ideality = 1.0

# plt.title('Complex Size distribution')
# data = np.array(total_partners)
# density = gaussian_kde(data.flatten())
# xs = np.linspace(data.min(), data.max(), 50)
# plt.plot(xs, density(xs), 'k', label='gamma distro fitted to mean/std')
#
# data = np.array(total_partners_old)
# density = gaussian_kde(data.flatten())
# xs = np.linspace(data.min(), data.max(), 50)
# plt.plot(xs, density(xs), 'r', label='interaction partners proxy')
#
# plt.xlabel('complex size')
# plt.ylabel('distribution density')
# plt.legend()
# plt.show()

total_partners = total_partners_old

# print total_partners
print np.any(np.array(total_partners) == 0)

abundance_range = np.random.permutation(abundance_range)


def calculate_van_Hoeff(aneuploidy_factor, complex_size):
    # We assume we have 1000 complexes that absorb 6000 genes total
    genes = range(0, 6000)
    # complexes are consecutive genes runs (simplifying assumption)
    shuffle(genes)

    genes = genes[:int(len(genes) * aneuploidy_factor)]
    genes = sorted(genes)
    genes = np.array(genes)

    complex_counter = 0
    free_protein_counter = len(genes)
    for a, b in zip(range(0, 6000 - complex_size, complex_size), range(complex_size, 6000, complex_size)):
        if np.sum(np.logical_and(genes > a, genes < b)) + 1 == complex_size:
            complex_counter += 1
            free_protein_counter -= complex_size

    total_van_Hoeff = complex_counter + free_protein_counter

    return total_van_Hoeff


def generate_complex_ids(total_partners):
    complex_contents = []
    i = 0
    local_partners = (np.array(total_partners) + 1).tolist()

    while i < len(abundance_range) - 45:  # -45 is here for manual abundant complex injection
        # while i < len(abundance_range):
        current_complex_size = random.choice(local_partners)
        complex = []
        for j in range(i, i + current_complex_size):
            complex.append(j)
        complex_contents.append(complex)
        # print complex
        i += current_complex_size

    last_complex = complex_contents[-1]
    last_complex = np.array(last_complex[last_complex < len(abundance_range)]).tolist()
    if type(last_complex) is int:
        del complex_contents[-1]
    else:
        complex_contents[-1] = last_complex

    # Manual large, abundant complex injection:
    complex_contents.append(range(complex_contents[-1][-1], len(abundance_range)))

    return complex_contents


def align_complex_abundances(complex_contents, abundance_correlation=0.7):
    aligned_abundances = np.copy(abundance_range)
    sorted_abundances = np.sort(abundance_range)

    for complex in complex_contents[:-1]:
        # print complex
        # print aligned_abundances[complex]

        # # Manual injection boosting of small complexes abundances:
        # if len(complex) < 45 and len(complex) > 40:
        #     aligned_abundances[complex] *= 100
        #     loc_ab_corr = 0.99
        #     average_abundance = np.mean(aligned_abundances[complex])
        #     aligned_abundances[complex] = aligned_abundances[complex]*(1 - loc_ab_corr) +\
        #                                 average_abundance*loc_ab_corr

        # else:
        #     average_abundance = np.mean(aligned_abundances[complex])
        #     aligned_abundances[complex] = aligned_abundances[complex]*(1 - abundance_correlation) +\
        #                                 average_abundance*abundance_correlation

        average_abundance = np.mean(aligned_abundances[complex])
        aligned_abundances[complex] = aligned_abundances[complex] * (1 - abundance_correlation) + \
                                      average_abundance * abundance_correlation
        # print aligned_abundances[complex]

    # Manual large, abundant complex injection:
    sorted_abundances = np.sort(abundance_range)
    print sorted_abundances
    print len(complex_contents[-1])
    average_abundance = np.mean(sorted_abundances[complex_contents[-1]])
    min_abundance = np.min(sorted_abundances[complex_contents[-1]])
    print average_abundance, min_abundance
    aligned_abundances[complex_contents[-1]] = sorted_abundances[complex_contents[-1]] * (1 - abundance_correlation) + \
                                               average_abundance * abundance_correlation

    return aligned_abundances


def calculate_free_mol_entities(aneuploidy_factor, complex_contents, abundances, buckets={}):
    multiplication_factor = np.random.binomial(1, aneuploidy_factor, len(abundances)) + 1
    aneup_abundance = abundances * multiplication_factor

    total_proteins = sum(aneup_abundance)
    total_complexes = 0
    bound_proteins = 0

    for complex in complex_contents:
        # print complex
        # print aneup_abundance[complex]
        complex_abundances = aneup_abundance[complex]
        loc_min = np.min(complex_abundances)
        total_complexes += loc_min
        bound_proteins += loc_min * len(complex)

        complex_size = len(complex)
        # yup, we are doing it with pointers and pointers are uncool and unpythonic.
        if buckets and complex_size in buckets.keys():
            buckets[complex_size][0] += loc_min
            buckets[complex_size][1] += np.sum(aneup_abundance[complex])
            buckets[complex_size][2] += np.sum(aneup_abundance[complex] - loc_min)

    total_molecules = total_proteins - bound_proteins + total_complexes

    return total_molecules


def core_sim_loop(base, total_partners, abundance_correlation=0.7, repeats=5, buckets=[]):
    # type: (object, object, object, object, object) -> object
    complex_contents = generate_complex_ids(total_partners)
    aligned_abundances = align_complex_abundances(complex_contents, abundance_correlation)
    re_runs = []

    bucket_dict = {}
    for key in buckets:
        bucket_dict[key] = [[] for _ in range(repeats)]

    buckets = bucket_dict

    for i in range(0, repeats):
        read_out = []

        for aneuploidy_factor in base:
            sub_buckets = {}
            if buckets:
                for key in buckets.keys():
                    sub_buckets[key] = [0, 0, 0]
            read_out.append(calculate_free_mol_entities(aneuploidy_factor, complex_contents,
                                                        aligned_abundances, sub_buckets))

            for key in buckets.keys():
                buckets[key][i].append(sub_buckets[key])

        read_out = np.array(read_out)
        re_runs.append(read_out)

    re_runs = np.array(re_runs)
    means = np.mean(re_runs, 0)
    stds = np.std(re_runs, 0)
    means, stds = (means / means[0], stds / means[0])

    for key, value in buckets.iteritems():
        val = np.array(value)
        buckets[key] = np.hstack((np.mean(val, axis=0), np.std(val, axis=0)))

    return means, stds, buckets


base = np.linspace(0.0, 1.0, 20).tolist()

arr_base = np.array(base)
arr_base += 1

means, stds, pre_buckets = core_sim_loop(base, total_partners,
                                         0.7,
                                          20,
                                         [3, 15])
np.savetxt('base.csv', arr_base, delimiter='\t', fmt='%f')
np.savetxt('means.csv', means, delimiter='\t', fmt='%f')
np.savetxt('stds.csv', stds, delimiter='\t', fmt='%f')


