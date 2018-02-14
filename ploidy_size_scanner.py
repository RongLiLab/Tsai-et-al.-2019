from csv import reader as csv_reader
from os import path, walk
import numpy as np
import pickle

directory = "C:\Users\Andrei\Dropbox\workspaces\JHU\Hung-Ji paper\size data"
ploidy_dict = {}
master_table = []


def line_reduce(line):
    return [elt for elt in line if elt]


def per_file_loop(fle):
    collector = []
    with open(path.join(directory, fle), 'r') as source:
        reader = csv_reader(source, delimiter='\t')
        reader.next()
        for line in reader:
            # print line
            line = line_reduce(line)
            if line:
                collector.append(float(line[5]))

    collector = np.array(collector)

    return np.mean(collector), np.std(collector)


with open(path.join(directory, 'ploidy_size.tsv'), 'r') as source:
    reader = csv_reader(source, delimiter='\t')
    reader.next()
    for line in reader:
        ploidy_dict[line[0].split('.')[0]] = float(line[1])

# print ploidy_dict

for loc, sub_folders, files in walk(directory):
    for fle in files:
        if fle[-4:] == '.txt':
            # print fle
            ploidy = ploidy_dict[fle[:-4]]
            mean, std = per_file_loop(fle)
            master_table.append([ploidy, mean, std])

master_table = np.array(master_table)
pickle.dump(master_table, open('ploidy_vs_size.dmp', 'w'))
print master_table