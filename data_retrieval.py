from csv import reader as csv_reader
import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
from pickle import dump

from scipy import stats

abundance_range = []
abundance_table = defaultdict(int)
interaction_table = []


with open('4932-GPM_201408.txt', 'r') as source:
    reader = csv_reader(source, delimiter='\t')
    for line in reader:
        abundance = float(line[2])
        abundance_range.append(abundance)
        name = line[1].split('.')[1]
        abundance_table[name] = abundance


with open('CervBinaryHQ.txt', 'r') as source:
    reader = csv_reader(source, delimiter='\t')
    header = reader.next()
    for line in reader:
        prot1 = line[0]
        prot2 = line[1]
        interaction_table.append([prot1, prot2])


interaction_table = np.array(interaction_table)
total_partners = []
partner_abundances = []

for unique_gene in np.unique(interaction_table):
    c1 = np.count_nonzero(interaction_table[:, 0] == unique_gene)
    c2 = np.count_nonzero(interaction_table[:, 1] == unique_gene)
    interacting_partners = set(interaction_table[interaction_table[:, 0] == unique_gene, 1]) | set(interaction_table[interaction_table[:, 1] == unique_gene, 0])
    interacting_partners = list(interacting_partners)
    total_partners.append(len(interacting_partners))
    for partner in interacting_partners:
        partner_abundances.append([abundance_table[unique_gene], abundance_table[partner]])


abundance_range = np.array(abundance_range)
log_range = np.log(abundance_range[abundance_range > 0])

partner_abundances = np.array(partner_abundances)

plt.hist(log_range, bins=50)
plt.show()

plt.hist(np.log(total_partners), bins=50)
plt.show()

plt.loglog(partner_abundances[:, 0], partner_abundances[:, 1], '.k')
plt.show()

slope, intercept, r_value, p_value, std_err = stats.linregress(partner_abundances[:, 0], partner_abundances[:, 1])
print slope, intercept, r_value, p_value, std_err

sel = np.logical_and(partner_abundances[:, 0] > 0, partner_abundances[:, 1] > 0)
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(partner_abundances[sel, 0]), np.log(partner_abundances[sel, 1]))
print slope, intercept, r_value, p_value, std_err


dump((abundance_range, total_partners), open('data_stats_dump.dmp', 'w'))
