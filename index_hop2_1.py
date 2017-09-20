#!/usr/bin/env python3.6

'''
Parses four fastq files from paired end sequencing output, Fwd,
IndexFwd, IndexRev, Rev. by a quality score cutoff -c, and determines
the level of index swapping observed by counting paired indeces and
swapped indices. Possible indexes are determined from index file -i
'''

import argparse
import os
import itertools
from collections import Counter


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
        return complement(s[::-1])


def get_arguments():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(prog='Index_hop2',
                                     description='Demultiplexing of \
                           paired end sequences by index pariing')

    parser.add_argument('-1', '--R1',
                        help='Read1 file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-2', '--R2',
                        help='Read2 file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-3', '--R3',
                        help='Read3 file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-4', '--R4',
                        help='Read4 file (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-i', '--Ind',
                        help='Indeces File (path)',
                        required=True,
                        type=argparse.FileType('rt', encoding='UTF-8 '))

    parser.add_argument('-c', '--qcutoff',
                        help='specify the qscore cutoff',
                        required=True,
                        type=int)

    return parser.parse_args()


args = get_arguments()

# set min quality score
min_qscore = args.qcutoff

# open some files
R1 = args.R1.name
R2 = args.R2.name
R3 = args.R3.name
R4 = args.R4.name
Index = args.Ind.name

# initialize index dictonary
Index_dict = {}

# populate index dictonary with known indices
with open(Index, 'r') as ind:
    for line in ind:
        line = line.strip('\n').split('\t')
        Index_dict[line[3]] = line[4]

# all possible indices
all_indices = []

# build list of all indices
for key in Index_dict:
    all_indices.append(Index_dict[key])

# all combinations of dual indices
all_combinations_list = list(
    itertools.combinations_with_replacement(all_indices, 2)
    )

# dict of all unique combinations of indices
all_combinations_dict = {}

# initialize the count to zero for each combination
for comb in all_combinations_list:
    all_combinations_dict[comb] = 0

# initialize empty dictonary for undetermined reads (contain N's)
undetermined_dict = Counter()

with open(R1, 'r') as r1,\
     open(R2, 'r') as r2,\
     open(R3, 'r') as r3,\
     open(R4, 'r') as r4:

    # Initialize line number
    NL = 0

    # loop over the lines in all files
    for line in zip(r1, r2, r3, r4):

        # strip each 'line' in the tuple of lines
        line = [l.strip() for l in line]

        # name the subsequent lines for later access
        if NL % 4 == 0:
            head = line

        if NL % 4 == 1:
            seq = line

        if NL % 4 == 2:
            plus = line

        if NL % 4 == 3:
            qline = line

            for j in qline[1:2]:
                qlist = []
                qtotal = 0
                for qual in j:
                    qtotal += ord(qual) - 33
                    qavg = qtotal/len(j)

                    if qavg >= min_qscore:

                        qlist.append(qual)
                        # print(qlist)

                        # identify that index pair
                        pair = seq[1], reverse_complement(seq[2])

                        # check if the index pair is a possible combination
                        if pair in all_combinations_dict:

                            # increment count of that combination
                            all_combinations_dict[pair] += 1

                        else:
                            undetermined_dict[pair] += 1

        # increment line number
        NL += 1

# make new directory for plot files
try:
    os.mkdir('Results'+'_'+str(min_qscore))
    os.chdir('Results'+'_'+str(min_qscore))
except:
    os.chdir('Results'+'_'+str(min_qscore))

# open a file to hold the matched index pair counts
# and write a first line
match = open('match_out.tsv', 'w')
match.write('Matched Index Pair\tCounts\n')

# open a file to hold the swapped index pair counts
# and write a first line
swap = open('swapped_out.tsv', 'w')
swap.write('Swapped Index Pair\tCounts\n')

# open file to hold the undetermined reads
# and write a first line
und = open('undetermined_out.tsv', 'w')
und.write('Undetermined Index Pair\tCounts\n')


# wite matched and swapped indecies to their respective files
swapped = 0
for comb in all_combinations_dict:
    if comb[0] == comb[1]:
        match.write(str(comb[0]) + '_' + str(comb[1]) + '\t' +
                    str(all_combinations_dict[comb]) + '\n')
    else:
        swapped += 1
        swap.write(str(comb[0]) + '_' +
                   str(comb[1]) + '\t' +
                   str(all_combinations_dict[comb]) + '\n')

for com in undetermined_dict:
    und.write((str(com[0]) + '_' +
               str(com[1]) + '\t' + str(undetermined_dict[com]) + '\n'))

# how many combinations were recorded
match = 0
for comb in all_combinations_dict:
    match += all_combinations_dict[comb]

undetermined = 0
for com in undetermined_dict:
    undetermined += undetermined_dict[com]

stats = open('stats.tsv', 'w')
stats.write('Filtered' + '\t' + 'Matched' + '\t' +
            'Swapped' + '\t' + 'Undetermined' + '\n')
stats.write(str(match+swapped+undetermined) + '\t' +
            str(match) + '\t' +
            str(swapped) + '\t' + str(undetermined)+'\n')

# Go back into parent
os.chdir('../')
