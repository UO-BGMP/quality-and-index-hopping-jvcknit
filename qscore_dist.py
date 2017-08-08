#!/usr/bin/env python3.6

import numpy as np
import matplotlib.pyplot as plt
import argparse
import gzip


def get_arguments():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(prog='qscore_dist',
        description= 'This program plots a distribution of quality scores over the base position of sequencing reads in a fastq file.')

    parser.add_argument('-f', '--infile', help='specify input file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-p', '--prefix', help='specify plot file output name',
                            required=True,
                            type=str)

    return parser.parse_args()


def convert_phred(letter):
    """Converts a single ASCII character into a phred score,
    based on phred+33"""
    score = ord(letter) - 33    # illumina 1.8 phred+33
    return score


args = get_arguments()
file = args.infile.name
pre = args.prefix

# determine the array dimentions
with gzip.open(file, 'r') as dim:
    count_lines = 0
    for line in dim:
        line = str(line, 'utf-8').strip('\n')

        if count_lines % 4 == 1:
            read_length = (len(line))

        count_lines += 1

sum_scores = np.zeros(101)
mean_scores = np.zeros(101)
all_qscores = np.zeros((int(count_lines/4), read_length), dtype=float)


with gzip.open(file, 'r') as fh:
    # index line number to select only phred lines
    NR = 1
    for line in fh:
        line = str(line, 'utf-8').strip('\n')

        # select only phred string
        if NR % 4 == 2:
            # index base position in the string
            base_pos = 0
            for base in line:

                # an array of all the scores
                all_qscores[NR//4, base_pos] = convert_phred(base)

                # increment base position
                base_pos += 1

        # increment line number
        NR += 1

print(all_qscores)

# get the variance for each base position
varr_array = np.var(all_qscores, axis=0)
# print(varr_array)

# get the stdev for each base position
stdev_array = np.std(all_qscores, axis=0)
# print(stdev_array)

# get the median for each position
med_array = np.median(all_qscores, axis=0)
# print(med_array)

# ######################## Plot These Results ##########################

xdata = np.arange(0, read_length, 1)
plt.figure()
# plot the variance as a green line
plt.plot(xdata, varr_array, 'g-', label='Variance')

# plot the standard deviation as a blue line
plt.plot(xdata, stdev_array, 'b-', label='Standard Deviation')

# plot median as a red line
plt.plot(xdata, med_array, 'r-', label='Median')

# legend the plot according to the plot labels
plt.legend()

# label the axes
plt.xlabel('Base Position')
plt.ylabel('Statistics')
plt.savefig(pre+'_dist.png')
