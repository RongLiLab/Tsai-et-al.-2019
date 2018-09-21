from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import gaussian_kde
from functools import partial
from scipy.interpolate import interp1d

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

large_complex_boost = 150
large_complex_correlation = 0.95
base_abundance_correlation = 0.85

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


def generate_complex_ids(total_partners):
    complex_contents = []
    i = 0
    local_partners = (np.array(total_partners)+1).tolist()

    while i < len(abundance_range) - 45:  # -45 is here for manual abundant complex injection
    # while i < len(abundance_range):
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

    # Manual large, abundant complex injection:
    complex_contents.append(range(complex_contents[-1][-1], len(abundance_range)))

    return complex_contents


def align_complex_abundances(complex_contents, abundance_correlation=0.7):
    aligned_abundances = np.copy(abundance_range)
    sorted_abundances = np.sort(abundance_range)

    for complex in complex_contents[:-1]:
        # print complex
        # print aligned_abundances[complex]

        # Manual injection boosting of small complexes abundances:
        if len(complex) < 45 and len(complex) > 40:
            aligned_abundances[complex] *= large_complex_boost
            loc_ab_corr = large_complex_correlation
            average_abundance = np.mean(aligned_abundances[complex])
            aligned_abundances[complex] = aligned_abundances[complex]*(1 - loc_ab_corr) +\
                                        average_abundance*loc_ab_corr

        else:
            average_abundance = np.mean(aligned_abundances[complex])
            aligned_abundances[complex] = aligned_abundances[complex]*(1 - abundance_correlation) +\
                                        average_abundance*abundance_correlation

        average_abundance = np.mean(aligned_abundances[complex])
        aligned_abundances[complex] = aligned_abundances[complex]*(1 - abundance_correlation) +\
                                        average_abundance*abundance_correlation
        # print aligned_abundances[complex]

    # Manual large, abundant complex injection:
    sorted_abundances = np.sort(abundance_range)
    print sorted_abundances
    print len(complex_contents[-1])
    average_abundance = np.mean(sorted_abundances[complex_contents[-1]])
    min_abundance = np.min(sorted_abundances[complex_contents[-1]])
    print average_abundance, min_abundance
    aligned_abundances[complex_contents[-1]] = sorted_abundances[complex_contents[-1]]*(1-abundance_correlation) +\
                                                average_abundance*abundance_correlation


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


def core_sim_loop(base, total_partners, abundance_correlation=0.7, repeats=5, buckets=[], alpha=1, ideality_correction=1):

    complex_contents = generate_complex_ids(total_partners)
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
    means, stds = (means / means[0]*alpha*ideality_correction, stds / means[0]*alpha*ideality_correction)

    for key, value in buckets.iteritems():
        val = np.array(value)
        buckets[key] = np.hstack((np.mean(val, axis=0), np.std(val, axis=0)))

    return means, stds, buckets


# reference (black line)
# variated thing
# variated thing values
#
def diameter_plot_array(defautlt_params, base_array,
                        names_array=None, base_val=None, base_name=None,
                        array_vals_name='', y_min=8, y_max=15, yticks=8):
    """
    Generates a plot of states based on a variation of a single parameter.

    :param base_array: 1D numpy array
    :param names_array: 1D numpy string arrays

    :return:
    """
    plt.title('Cell diameter vs ploidy vs %s' % array_vals_name)

    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base)
    arr_base += 1

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

    plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    base_yticks_array = np.linspace(y_min, y_max, yticks)
    base_yticks_array_names = [int(val) if i % 2 == 0 else '' for i, val in enumerate(base_yticks_array)]
    plt.yticks(base_yticks_array, base_yticks_array_names)

    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")

    cmap = get_cmap('brg')

    if names_array is None:
        names_array = base_array

    if base_val is not None and base_name is None:
        base_name = base_val

    # base array is assumed to be a 1D numpy array

    color_remap = np.linspace(0.5, 0.9, len(base_array))[np.argsort(base_array)[np.argsort(base_array)]]

    defautlt_params['base'] = base
    defautlt_params['total_partners'] = total_partners

    # for key, value in defautlt_params.iteritems():
    #     print key, value

    running_partial = core_sim_loop

    free_kwarg = ''
    for key, value in defautlt_params.iteritems():
        if value is not None:
            running_partial = partial(running_partial, **{key: value})
        else:
            free_kwarg = key

    # partial_function = partial(core_sim_loop, **defautlt_params)
    partial_function = running_partial

    # todo: form a partial function
    # => kwargs

    # reference plot
    if base_val is not None:
        # means, stds, pre_buckets = core_sim_loop(base,
        #                                          total_partners,
        #                                          base_val,
        #                                          20,
        #                                          [3, 15])

        means, stds, pre_buckets = partial_function(**{free_kwarg: base_val})

        plt.plot(arr_base,
                 corrfactor * np.cbrt(means),
                 '--k',
                 label=base_name,
                 lw=2)

        plt.fill_between(arr_base,
                         corrfactor * np.cbrt(means - stds),
                         corrfactor * np.cbrt(means + stds),
                         color='k',
                         alpha=.3)



    # the variation loop
    for i, (var_value, color, name) in enumerate(zip(base_array, color_remap, names_array)):

        # color generation from the 0.5-0.9 space
        color = cmap((color - 0.5) / (0.9 - 0.5) * 0.8 + 0.1)

        # generation of the curve proper
        # means, stds, pre_buckets = core_sim_loop(base,
        #                                          total_partners,
        #                                          var_value,
        #                                          20,
        #                                          [3, 15])

        means, stds, pre_buckets = partial_function(**{free_kwarg: var_value})

        plt.plot(arr_base,
                 corrfactor * np.cbrt(means),
                 '--',
                 color=color,
                 label=name,
                 lw=2)

        plt.fill_between(arr_base,
                         corrfactor * np.cbrt(means - stds),
                         corrfactor * np.cbrt(means + stds),
                         color=color,
                         alpha=.3)

    # here what is plotted is varied as well.
    plt.legend(loc='best', ncol=3, title=array_vals_name)

    # and y-axis here needs to be adjusted dynamically as well
    plt.axis([0.95, 2.05, y_min, y_max])
    plt.show()


def get_osmotic_pressure(y_min=8, y_max=18, yticks=8):

    plt.title('Cell diameter vs ploidy')

    base = np.linspace(0.0, 1.0, 20).tolist()
    arr_base = np.array(base)
    arr_base += 1

    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)

    plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    base_yticks_array = np.linspace(y_min, y_max, yticks)
    base_yticks_array_names = [int(val) if i % 2 == 0 else '' for i, val in enumerate(base_yticks_array)]
    plt.yticks(base_yticks_array, base_yticks_array_names)

    plt.xlabel("Ploidy")
    plt.ylabel("Equivalent Diameter")

    means, stds, pre_buckets = core_sim_loop(base, total_partners,
                                             abundance_correlation=base_abundance_correlation)

    interp_prediction = interp1d(arr_base, corrfactor*np.cbrt(means), kind='cubic')

    plt.plot(arr_base,
             corrfactor * np.cbrt(means),
             '--k',
             label='',
             lw=2)

    plt.fill_between(arr_base,
                     corrfactor * np.cbrt(means - stds),
                     corrfactor * np.cbrt(means + stds),
                     color='k',
                     alpha=.3)

    sorting_index = np.argsort(ploidy_vs_size[:, 0])

    plt.errorbar(ploidy_vs_size[sorting_index, 0],
                 ploidy_vs_size[sorting_index, 1]*corrfactor,
                 yerr=ploidy_vs_size[sorting_index, 2]*corrfactor,
                 fmt='ro--')

    undervolume = (ploidy_vs_size[sorting_index, 1]*corrfactor)**3/(interp_prediction(
        ploidy_vs_size[sorting_index, 0]))**3

    alpha_0 = np.exp(-0.05e6*(18.03e-3)/8.314/293)

    press = []
    for value in undervolume:
        alpha = 1 - (1 - alpha_0)/value
        pressure = -8.314*293./18.03e-3*np.log(alpha)
        print value, alpha, pressure, pressure/0.05e6
        if not np.isnan(pressure):
            press.append(pressure)

    press = np.array(press)

    aneuploid_means = press[1:-1]
    euploid_means = press[[0, -1],]

    # print alpha_0, alpha_0*undervolume, -8.314*293./18.03e-3*np.log(undervolume)

    plt.legend(loc='best', ncol=3)

    # and y-axis here needs to be adjusted dynamically as well
    plt.axis([0.95, 2.05, y_min, y_max])
    plt.show()

    true_euploids_op, true_aneuploid_op = load(open('osmotic_pressure.dmp', 'r'))

    normalizer_1 = np.median(euploid_means)
    normalizer_2 = np.median(true_euploids_op)

    true_euploids_op = np.array(true_euploids_op)/normalizer_2
    true_aneuploid_op = np.array(true_aneuploid_op)/normalizer_2

    euploid_means = euploid_means/normalizer_1
    aneuploid_means = aneuploid_means/normalizer_1

    plt.suptitle("large complex boost: x%.2f; large complex corr: %.2f, overall corr: %.2f\n"
                 "simulations medians ratio: 1:%.2f |"
                 " true medians ratio: 1:%.2f" % (large_complex_boost,
                                               large_complex_correlation,
                                               base_abundance_correlation,
                 np.median(aneuploid_means)/np.median(euploid_means),
                 np.median(true_aneuploid_op)/np.median(true_euploids_op)))
    # print "simulations medians ratio: 1:%f" % (np.median(aneuploid_means)/np.median(euploid_means))
    # print "true medians ratio: 1:%f" % (np.median(true_aneuploid_op)/np.median(true_euploids_op))

    vibration = 0.05*(2*np.random.rand(len(euploid_means))-1)
    plt.plot(1+vibration, euploid_means, 'ko')
    vibration = 0.05*(2*np.random.rand(len(aneuploid_means))-1)
    plt.plot(2+vibration, aneuploid_means, 'ko')

    vibration = 0.05 * (2 * np.random.rand(len(true_euploids_op)) - 1)
    plt.plot(3 + vibration, true_euploids_op, 'ko')
    vibration = 0.05 * (2 * np.random.rand(len(true_aneuploid_op)) - 1)
    plt.plot(4 + vibration, true_aneuploid_op, 'ko')

    bplot1 = plt.boxplot([euploid_means, aneuploid_means, true_euploids_op, true_aneuploid_op],
                         patch_artist=True
                         )

    colors = ['lightblue', 'pink', 'lightblue', 'pink']
    for bplot in [bplot1]:
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)


    plt.ylabel('Turgor pressure\nrelative to the euploid median')

    # add x-tick labels
    plt.xticks([y + 1 for y in range(4)], ['simulated\neuploid\nturgor',
                                           'simulated\naneuploid\nturgor',
                                           'measured\neuploid\nturgor',
                                           'measured\naneuploid\nturgor'])

    plt.show()


if __name__ == "__main__":

    corrfactor = ploidy_vs_size[0, 1]
    ploidy_vs_size[:, 1] /= corrfactor
    ploidy_vs_size[:, 2] /= corrfactor

    # plt.figure(num=None, figsize=(5, 6), dpi=300, facecolor='w', edgecolor='k')

    get_osmotic_pressure()

    ################################################################################################
    ## Swipe for the abundance correlation swipe
    ################################################################################################

    # corr_vars = np.linspace(0.5, 0.9, 5)
    # kwargs = {'abundance_correlation': None, 'repeats': 20, 'buckets': [5, 13], 'alpha': 1, 'ideality_correction': 1}
    # diameter_plot_array(kwargs, corr_vars, corr_vars*100, 0, 0, 'complex abundance corr. (%)', 8, 15, 8)

    # corr_vars = np.linspace(0.15, 0.95, 5)
    # kwargs = {'abundance_correlation': 0.7, 'repeats': 20, 'buckets': [5, 13], 'alpha': None, 'ideality_correction': 1}
    # diameter_plot_array(kwargs, corr_vars, corr_vars*100, 1.0, 100, 'alpha value. (%)', 4, 14, 11)
    #
    # corr_vars = np.linspace(0.5, 1.5, 6)
    # kwargs = {'abundance_correlation': 0.7, 'repeats': 20, 'buckets': [5, 13], 'alpha': 1, 'ideality_correction': None}
    # diameter_plot_array(kwargs, corr_vars, corr_vars*100, 1.0, 100, 'ideality correction factor. (%)', 8, 15, 8)

    # ################################################################################################
    # ## Rendering routines for the abundance
    # ################################################################################################
    #
    # for key, value_matrix in buckets.iteritems():
    #
    #     plt.figure()
    #     plt.title('Molecules abundance vs ploidy % s for complex size %s' % (curr_corr, key))
    #
    #     print value_matrix
    #
    #     norm = value_matrix[0, 0] + value_matrix[0, 2]
    #     value_matrix /= norm
    #
    #     plt.plot(arr_base,
    #              value_matrix[:, 0],
    #              '-k',
    #              label='complexes',
    #              lw=2
    #              )
    #
    #     # plt.plot(arr_base,
    #     #          value_matrix[:, 1]/value_matrix[0, 1],
    #     #          '-g',
    #     #          label='total proteins',
    #     #          lw=2
    #     #          )
    #
    #     plt.plot(arr_base,
    #              value_matrix[:, 2],
    #              '-b',
    #              label='free proteins',
    #              lw=2
    #              )
    #
    #     plt.plot(arr_base,
    #              value_matrix[:, 0] + value_matrix[:, 2],
    #              '-r',
    #              label='total molecules',
    #              lw=2
    #              )
    #
    #     plt.fill_between(arr_base,
    #                      value_matrix[:, 0]+value_matrix[:, 3],
    #                      value_matrix[:, 0]-value_matrix[:, 3],
    #                      color='k',
    #                      alpha=.3)
    #
    #     # plt.fill_between(arr_base,
    #     #                  (value_matrix[:, 1]+value_matrix[:, 4])/value_matrix[0, 1],
    #     #                  (value_matrix[:, 1]-value_matrix[:, 4])/value_matrix[0, 1],
    #     #                  color='g',
    #     #                  alpha=.3)
    #
    #     plt.fill_between(arr_base,
    #                      value_matrix[:, 2]+value_matrix[:, 5],
    #                      value_matrix[:, 2]-value_matrix[:, 5],
    #                      color='b',
    #                      alpha=.3)
    #
    #     plt.fill_between(arr_base,
    #                      value_matrix[:, 0] + value_matrix[:, 2] + np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
    #                      value_matrix[:, 0] + value_matrix[:, 2] - np.sqrt(np.power(value_matrix[:, 3], 2) + np.power(value_matrix[:, 5], 2)),
    #                      color='r',
    #                      alpha=.3)
    #
    #     plt.legend(loc='best', ncol=1, title='')
    #     plt.axis([0.95, 2.05, 0, 4])
    #
    #     plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    #     plt.yticks(np.linspace(0, 4, 9))
    #
    #     ax = plt.gca()
    #     ax.set_axisbelow(True)
    #     ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    #     ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    #
    #     plt.xlabel("Ploidy")
    #     plt.ylabel("Abundance relative to total free molecule abundance in haploid")
    #
    #     plt.show()

    # ################################################################################################
    # ## Rendering routines for the complex size swipe
    # ################################################################################################
    #
    # plt.title('Cell diameter vs ploidy vs complex size')
    # plt.plot(arr_base, corrfactor*np.cbrt(arr_base), '--k', label='1', lw=2)
    #
    # ax = plt.gca()
    # ax.set_axisbelow(True)
    #
    # ax.yaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    # ax.xaxis.grid(color='gray', linestyle='solid', alpha=0.5)
    #
    # plt.xticks(np.linspace(1, 2, 9), ['1.00', '', '1.25', '', '1.50', '', '1.75', '', '2.00'])
    # plt.yticks(np.linspace(8, 15, 8), [8, '', 10, '', 12, '', 14, ''])
    #
    # plt.xlabel("Ploidy")
    # plt.ylabel("Equivalent Diameter")
    #
    # complex_sizes = [2, 5, 10, 20, 40]
    #
    # for c_size in complex_sizes:
    #     c = cmap(c_size/40.*0.8+0.1)
    #
    #     means, stds, buckets = core_sim_loop(base,[c_size], 0.75, 20, [])
    #
    #     plt.plot(arr_base,
    #              corrfactor*np.cbrt(means),
    #              '--',
    #              color=c,
    #              label=c_size, lw=2)
    #
    #     plt.fill_between(arr_base,
    #                      corrfactor*np.cbrt(means-stds),
    #                      corrfactor*np.cbrt(means+stds),
    #                      color=c,
    #                      alpha=.3)
    #
    #
    # plt.legend(loc='best', ncol=3, title='complex size')
    # plt.axis([0.95, 2.05, 8, 15])
    # plt.show()
