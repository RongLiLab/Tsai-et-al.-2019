from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random
from matplotlib.cm import get_cmap

random.seed(42)

abundance_range, total_partners = load(open('data_stats_dump.dmp', 'r'))
abundance_range = abundance_range[abundance_range > 1]
total_partners = np.array(total_partners)
total_partners = total_partners[total_partners > 1]
total_partners = total_partners[total_partners < 45]
total_partners = total_partners.tolist()



plt.hist(total_partners)
plt.show()

total_partners = np.random.gamma(1.538, 2.016, len(total_partners))
total_partners = np.ceil(total_partners).astype(np.int16)
total_partners = total_partners[total_partners > 1]

plt.hist(total_partners)
plt.show()

print total_partners
print np.any(np.array(total_partners) == 0)

abundance_range = np.random.permutation(abundance_range)

#TODO: further improvement idea: use actual experimental data for protein perturbation by aneuploidy from Rancati


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

    return complex_contents


def align_complex_abundances(complex_contents, abundance_correlation=0.7):
    aligned_abundances = np.copy(abundance_range)

    for complex in complex_contents:
        # print complex
        # print aligned_abundances[complex]
        average_abundance = np.mean(aligned_abundances[complex])
        aligned_abundances[complex] = aligned_abundances[complex]*(1 - abundance_correlation) +\
                                      average_abundance*abundance_correlation
        # print aligned_abundances[complex]

    return aligned_abundances


def calculate_free_mol_entities(aneuploidy_factor, complex_contents, abundances):
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

        # we have formed stable complexes

    total_molecules = total_proteins - bound_proteins + total_complexes

    return total_molecules


def core_simulation_loop():
    pass


def core_sim_loop(base, abundance_correlation=0.7, repeats=5):
    complex_contents = generate_complex_ids()
    aligned_abundances = align_complex_abundances(complex_contents, abundance_correlation)
    re_runs = []
    for _ in range(0, repeats):
        read_out = []

        for aneuploidy_factor in base:
            # read_out.append(calculate_van_Hoeff(aneuploidy_factor, 6))
            read_out.append(calculate_free_mol_entities(aneuploidy_factor, complex_contents,
                                                        aligned_abundances))

        read_out = np.array(read_out)
        re_runs.append(read_out)

    re_runs = np.array(re_runs)
    means = np.mean(re_runs, 0)
    stds = np.std(re_runs, 0)
    means, stds = (means / means[0], stds / means[0])

    return means, stds


if __name__ == "__main__":

    base = np.arange(0.00, 1.05, 0.05).tolist()
    arr_base = np.array(base)
    arr_base += 1
    cmap = get_cmap('Reds')

    plt.title('Molecules abundance vs ploidy vs abundance correlation')
    # plt.plot(arr_base, arr_base, 'ok')
    plt.plot(arr_base, np.cbrt(arr_base), '-k', label='0')

    for abundance_correlation in range(5, 10):
        c = cmap(abundance_correlation*0.1)
        means, stds = core_sim_loop(base, abundance_correlation*0.1, 10)
        # plt.errorbar(arr_base, means, yerr=stds, fmt='o', color=c, label=abundance_correlation*0.1)
        plt.plot(arr_base, np.cbrt(means), color=c, label=abundance_correlation*0.1)
        plt.fill_between(arr_base, np.cbrt(means-stds), np.cbrt(means+stds), color=c, alpha=.3)

    plt.xlabel("ploidy")
    # plt.ylabel("pressure")
    plt.ylabel("size")
    plt.legend(loc=8, ncol=3, title='complex member abundance correlation')
    # plt.axis([0.95, 2.05, 0.9, 3.2])
    plt.axis([0.95, 2.05, 0.9, 1.5])
    plt.show()

    # plt.title('Cell size vs ploidy at complex memeber correlation of %s' % abundance_correlation)
    # plt.plot(base, np.cbrt(base), 'ok', label='euploidy linear interpolation')
    # plt.errorbar(base, np.cbrt(means), yerr=np.cbrt(means+stds)-np.cbrt(means-stds), fmt='or', label='simulation results')
    # plt.xlabel("ploidy")
    # plt.ylabel("size")
    # plt.axis([0.95, 2.05, 0.9, 1.5])
    # plt.show()
