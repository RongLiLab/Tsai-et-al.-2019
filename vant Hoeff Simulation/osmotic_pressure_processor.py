from csv import reader as csv_reader
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde, t, sem
from pickle import dump

type_index = []
master_table = []


def smooth_histogram(data, color='k'):
    fltr = np.logical_not(np.isnan(data))
    density = gaussian_kde(data[fltr].flatten())
    xs = np.linspace(data[fltr].min(), data[fltr].max(), 200)
    plt.plot(xs, density(xs), color)


def clear_outlier(quadruplet, FDR=0.05, no_clear=True):

    potential_outliers = []

    if no_clear:
        return quadruplet

    for i in range(0, 3):
        potential_outlier = quadruplet[i]
        others_index = np.ones((4,), bool)
        others_index[i] = False
        all_others = quadruplet[others_index]
        lower, upper = t.interval(1-FDR, 2, np.mean(all_others), np.std(all_others))
        if potential_outlier > upper or potential_outlier < lower:
            potential_outliers.append(i)

    if len(potential_outliers) == 0:
        print 'no one cleared in %s' % quadruplet
        return quadruplet

    if len(potential_outliers) == 1:
        others_index = np.ones((4,), bool)
        others_index[potential_outliers[0]] = False
        print 'cleared %s in %s' % (quadruplet[potential_outliers[0]], quadruplet)
        return quadruplet[others_index]

    raise Exception('something went wrong with outlier calculation for the set %s; %s' %
                    (quadruplet, potential_outliers))


with open('20171107_yeast.csv', 'rb') as source_file:
    reader = csv_reader(source_file)
    header = reader.next()
    for line in reader:
        strain_type = line[1][:1]
        strain_index = line[1][1:]
        type_index.append(strain_type)
        master_table.append((float(strain_index), float(line[2])))

type_index = np.array(type_index)
master_table = np.array(master_table)

aneuploid_selector = type_index == 'a'
euploid_selector = type_index == 'e'

euploids_raw = master_table[euploid_selector][:, 1]
aneuploids_raw = master_table[aneuploid_selector][:, 1]

aneuploid_clusters = np.unique(master_table[aneuploid_selector, 0])
euploid_clusters = np.unique(master_table[euploid_selector, 0])

aneuploid_cleared_clusters = []
euploid_cleared_clusters = []

for cluster in aneuploid_clusters:
    aneuploid_cleared_clusters.append(clear_outlier(master_table[np.logical_and(master_table[:,
                                                                         0] == cluster,
                                               aneuploid_selector), 1]))

for cluster in euploid_clusters:
    euploid_cleared_clusters.append(clear_outlier(master_table[np.logical_and(master_table[:,
                                                                        0] == cluster,
                                             euploid_selector), 1]))

aneuploid_means = []
euploid_means = []

for i, cluster in enumerate(aneuploid_cleared_clusters):
    aneuploid_means.append(np.mean(cluster))
    plt.errorbar(i, np.mean(cluster), yerr=sem(cluster), fmt='rs--')

for i, cluster in enumerate(euploid_cleared_clusters):
    euploid_means.append(np.mean(cluster))
    plt.errorbar(i, np.mean(cluster), yerr=sem(cluster), fmt='ks--')

plt.show()

print "means ratio: 1:%f" % (np.mean(aneuploid_means)/np.mean(euploid_means))

vibration = 0.05*(2*np.random.rand(len(euploid_means))-1)
plt.plot(1+vibration, euploid_means, 'ko')
vibration = 0.05*(2*np.random.rand(len(aneuploid_means))-1)
plt.plot(2+vibration, aneuploid_means, 'ko')
plt.boxplot([euploid_means, aneuploid_means])

plt.show()

dump((euploid_means, aneuploid_means), open('osmotic_pressure.dmp', 'w'))

# print euploids_raw

# smooth_histogram(euploids_raw)
# smooth_histogram(aneuploids_raw, 'r')
#
# plt.show()


