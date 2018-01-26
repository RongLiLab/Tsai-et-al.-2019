from random import shuffle
import numpy as np
from matplotlib import pyplot as plt


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

    plt.plot(base, np.sqrt(read_out))
    plt.xlabel("ploidy")
    plt.ylabel("size")
    plt.axis([1, 2, 0, 100])
    plt.show()