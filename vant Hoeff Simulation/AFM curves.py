from csv import reader as csv_reader
import numpy as np
from matplotlib import pyplot as plt


aneuploids_file = "aneuploids.csv"
haploids_file = "haploid.csv"


def read_data(source_fle):
    data_accumulator = []
    with open(source_fle, 'rb') as source:
        reader = csv_reader(source)
        header = reader.next()
        for line in reader:
            data_accumulator.append(line)

        return np.array(data_accumulator).astype(np.float)


def show_curves(data_array, color='r'):
    for value in np.unique(data_array[:, 0]):
        xy = data_array[data_array[:, 0] == value, 1:]
        x = xy[:, 0]
        y = xy[:, 1]
        if x[-1] < 45:
            x /= x[-1]
            plt.plot(x, y, color, alpha=0.5)


if __name__ == "__main__":
    data = read_data(aneuploids_file)
    show_curves(data, 'r')

    data = read_data(haploids_file)
    show_curves(data, 'k')

    plt.show()