from random import shuffle
import numpy as np
from matplotlib import pyplot as plt
from pickle import load
import random


abundance_range, total_partners = load(open('data_stats_dump.dmp', 'r'))
abundance_range = abundance_range[abundance_range > 1]

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

    while i < len(abundance_range):
        current_complex_size = random.choice(total_partners+1)
        complex = []
        for j in range(i, i+current_complex_size):
            complex.append(j)
        print complex
        complex_contents.append(complex)
        i += current_complex_size

    return complex_contents


def calculate_free_mol_entities(aneuploidy_factor, complex_contents):
    multiplication_factor = np.random.binomial(1, 1+aneuploidy_factor, len(abundance_range))




    # we have formed stable complexes


if __name__ == "__main__":

    base = np.arange(0.05, 1., 0.05).tolist()
    read_out = []

    for aneuploidy_factor in base:
        read_out.append(calculate_van_Hoeff(aneuploidy_factor, 6))

    base.insert(0, 0)
    base.append(1)
    read_out.insert(0, 0)
    read_out.append(1000)

    base = np.array(base)
    read_out = np.array(read_out)

    print base
    print read_out

    base = base + 1
    read_out = read_out + 1000

    plt.plot(base, read_out)
    plt.xlabel("ploidy")
    plt.ylabel("pressure")
    plt.axis([1, 2, 0, 7000])
    plt.show()

    plt.title('Cell size vs ploidy')
    plt.plot(base, np.cbrt(read_out), label='simulation results')
    plt.plot(base, np.cbrt(base*1000), label='linear volume incerease')
    plt.xlabel("ploidy")
    plt.ylabel("size")
    plt.axis([1, 2, 5, 20])
    plt.show()