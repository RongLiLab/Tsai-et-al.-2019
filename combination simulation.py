from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import gaussian_kde

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

    genes = genes[:int(len(genes)*aneuploidy_factor)]
    genes = sorted(genes)
    genes = np.array(genes)

    complex_counter = 0
    free_protein_counter = len(genes)
    for a, b in zip(range(0, 6000-complex_size, complex_size), range(complex_size, 6000, complex_size)):
        if np.sum(np.logical_and(genes > a, genes < b))+1 == complex_size:
            complex_counter += 1
            free_protein_counter -= complex_size

    total_van_Hoeff = complex_counter + free_protein_counter

    return total_van_Hoeff


def generate_complex_ids():
    complex_contents = []
    i = 0
    local_partners = (np.array(total_partners)+1).tolist()

    # while i < len(abundance_range)-45:  # -45 is here for manual abundant complex injection
    while i < len(abundance_range):
        current_complex_size = random.choice(local_partners)
        complex = []
        for j in range(i, i+current_complex_size):
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

    # # Manual large, abundant complex injection:
    # complex_contents.append(range(complex_contents[-1][-1], len(abundance_range)))

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
        aligned_abundances[complex] = aligned_abundances[complex]*(1 - abundance_correlation) +\
                                        average_abundance*abundance_correlation
        # print aligned_abundances[complex]

    # # Manual large, abundant complex injection:
    # sorted_abundances = np.sort(abundance_range)
    # print sorted_abundances
    # print len(complex_contents[-1])
    # average_abundance = np.mean(sorted_abundances[complex_contents[-1]])
    # print average_abundance
    # aligned_abundances[complex_contents[-1]] = sorted_abundances[complex_contents[-1]]*(1-abundance_correlation) +\
    #                                             average_abundance*abundance_correlation


    return aligned_abundances


def calculate_free_mol_entities(aneuploidy_factor, complex_contents, abundances, buckets={}):
    multiplication_factor = np.random.binomial(1, aneuploidy_factor, len(abundances))+1
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


def core_sim_loop(base, abundance_correlation=0.7, repeats=5, buckets=[]):
    complex_contents = generate_complex_ids()
    aligned_abundances = align_complex_abundances(complex_contents, abundance_correlation)
    re_runs = []

    bucket_dict = {}
    for key in buckets:
        bucket_dict[key] = [[]for _ in range(repeats)]

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


if __name__ == "__main__":

    corrfactor = ploidy_vs_size[0, 1]
    ploidy_vs_size[:, 1] /= corrfactor
    ploidy_vs_size[:, 2] /= corrfactor

    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base)
    arr_base += 1
    cmap = get_cmap('gnuplot')

    # plt.figure(num=None, figsize=(5, 6), dpi=300, facecolor='w', edgecolor='k')

    plt.title('Molecules abundance vs ploidy vs abundance correlation')
    # plt.plot(arr_base, arr_base, 'ok')
    plt.plot(arr_base, corrfactor*np.cbrt(arr_base), '--k', label='0', lw=2)

    ax = plt.gca()
    ax.set_axisbelow(True)

    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

    plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    plt.yticks(np.linspace(8, 15, 8), [8, '', 10, '', 12, '', 14, ''])

    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")

    # filtered = lowess(ploidy_vs_size[:, 1], ploidy_vs_size[:, 0], is_sorted=True, frac=0.05, it=0)
    # plt.errorbar(ploidy_vs_size[:, 0], ploidy_vs_size[:, 1], yerr=ploidy_vs_size[:, 2], fmt='or', label='exp')
    # plt.plot(filtered[:, 1], filtered[:, 0], 'ob', label='exp')

    corr_min = 0.75
    corr_max = 0.95

    for abundance_correlation in np.linspace(corr_min, corr_max, 5):
        abundance_correlation = abundance_correlation*10
        c = cmap((abundance_correlation*0.1-0.7)/(corr_max-corr_min)*0.7)
        means, stds, buckets = core_sim_loop(base, abundance_correlation*0.1, 20, [3, 7, 15])

        print buckets

        # plt.errorbar(arr_base, means, yerr=stds, fmt='o', color=c, label=abundance_correlation*0.1)
        plt.plot(arr_base,
                 corrfactor*np.cbrt(means),
                 '--',
                 color=c,
                 label=abundance_correlation*10, lw=2)
        plt.fill_between(arr_base,
                         corrfactor*np.cbrt(means-stds),
                         corrfactor*np.cbrt(means+stds),
                         color=c,
                         alpha=.3)


    plt.legend(loc='best', ncol=3, title='complex abundance corr. (%)')
    # plt.axis([0.95, 2.05, 0.9, 3.2])
    plt.axis([0.95, 2.05, 8, 15])
    plt.show()

    for key, value_matrix in buckets.iteritems():

        plt.figure()
        plt.title('Molecules abundance vs ploidy % s for complex size %s' % (corr_max, key))

        print value_matrix

        norm = value_matrix[0, 0] + value_matrix[0, 2]
        value_matrix /= norm

        plt.plot(arr_base,
                 value_matrix[:, 0],
                 '-k',
                 label='complexes',
                 lw=2
                 )

        plt.plot(arr_base,
                 value_matrix[:, 1],
                 '-g',
                 label='total proteins',
                 lw=2
                 )

        plt.plot(arr_base,
                 value_matrix[:, 2],
                 '-b',
                 label='free proteins',
                 lw=2
                 )

        plt.plot(arr_base,
                 value_matrix[:, 0] + value_matrix[:, 2],
                 '-r',
                 label='total molecules',
                 lw=2
                 )

        plt.fill_between(arr_base,
                         value_matrix[:, 0]+value_matrix[:, 3],
                         value_matrix[:, 0]-value_matrix[:, 3],
                         color='k',
                         alpha=.3)

        plt.fill_between(arr_base,
                         value_matrix[:, 1]+value_matrix[:, 4],
                         value_matrix[:, 1]-value_matrix[:, 4],
                         color='g',
                         alpha=.3)

        plt.fill_between(arr_base,
                         value_matrix[:, 2]+value_matrix[:, 5],
                         value_matrix[:, 2]-value_matrix[:, 5],
                         color='b',
                         alpha=.3)

        plt.fill_between(arr_base,
                         value_matrix[:, 0] + value_matrix[:, 2] + np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
                         value_matrix[:, 0] + value_matrix[:, 2] - np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
                         color='r',
                         alpha=.3)

        plt.legend(loc='best', ncol=1, title='')
        plt.axis([0.95, 2.05, 0, 20])

        plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
        plt.yticks(np.linspace(0, 20, 11))

        ax = plt.gca()
        ax.set_axisbelow(True)
        ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
        ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

        plt.xlabel("Ploidy")
        plt.ylabel("Abundance relative to total free molecule abundance in haploid")

        plt.show()
