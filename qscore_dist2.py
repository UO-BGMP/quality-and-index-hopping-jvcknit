#!/usr/bin/env python3.6

import numpy as np
import matplotlib.pyplot as plt
import argparse
import gzip
import os
from collections import Counter


def get_arguments():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(prog='qscore_dist',
        description= 'This program plots a distribution of quality scores over the base position of sequencing reads in a fastq file.')

    parser.add_argument('-f', '--infile', help='specify input file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    return parser.parse_args()


def convert_phred(letter):
    """Converts a single ASCII character into a phred score,
    based on phred+33"""
    score = ord(letter) - 33    # illumina 1.8 phred+33
    return score

args = get_arguments()
file = args.infile.name

all_qscores = np.zeros(101)

mean_readscores = Counter()



with open(file, 'r') as fh:
    # index line number to select only phred lines
    count = 0
    NR = 1
    for line in fh:
        line = str(line).strip('\n')

        count = 2
        # select only phred string
        if NR % 4 == 0:
            # index base position in the string
            base_pos = 0
            readscore = 0
            for base in line:
                # an array of all the scores
                all_qscores[base_pos] = ((int(convert_phred(base)) + all_qscores[base_pos]) / (count))

                readscore += int(convert_phred(base))
                read_avg = int(readscore//len(line))
                # increment base position
                base_pos += 1
            count += 1

            mean_readscores[read_avg] += 1
        # increment line number
        NR += 1

# make new directory for plot files
try:
    os.mkdir('Plots')
    os.chdir('Plots')
except:
    os.chdir('Plots')


# determine the quality score distribution across the read length
xdata = np.arange(0, len(all_qscores[all_qscores>0]), 1)
plt.figure()
plt.bar(xdata, all_qscores[all_qscores>0], width = 0.5)
print(all_qscores[all_qscores>0])

# label the axes
plt.title(str(file)+'\n'+'Mean quality score for each Base Position')
plt.xlabel('Base Position')
plt.ylabel('Mean Quality Score')
plt.savefig(file+'_dist1.png')



# determine the quality score distribution across reads
x = []
y = []
for key in mean_readscores:
    x.append(key)
    y.append(mean_readscores[key])


plt.figure()
plt.bar(x, y, width = 0.5)
plt.title(str(file)+'\n'+'Frequency of average quality scores over all reads')
plt.xlabel('Avg QS/read')
plt.ylabel('Frequency')
plt.savefig(file+'_dist2.png')

# go back to directory script was run in
os.chdir('../')