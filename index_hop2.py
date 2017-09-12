#!/usr/bin/env python3.6

'''
Parses four fastq files, Read1, Read2, Read3, Read4 by a quality score
cutoff -c, and determines the level of index swapping observed by counting paired indeces and swapped indices, determined from an index file -i.
'''

import argparse
import gzip
import os
import itertools

def get_arguments():
    '''Define and return command line options.'''
    parser = argparse.ArgumentParser(prog='Index_hop2',
        description=
         'Demultiplexing of paired end sequences by index pariing')

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

    parser.add_argument('-c', '--qcutoff', help='specify the qscore cutoff',
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
	itertools.combinations_with_replacement(all_indices,2)
	)

# dict of all unique combinations of indices
all_combinations_dict = {}

# initialize the count to zero for each combination
for comb in all_combinations_list:
	all_combinations_dict[comb] = 0

with open(R1 ,'r') as r1,\
	 open(R2, 'r') as r2,\
	 open(R3, 'r') as r3,\
	 open(R4, 'r') as r4:

	# Initialize line number
	NL = 0

	# loop over the lines in all files
	for line in zip(r1,r2,r3,r4):

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
				for qual in j:
					if (ord(qual) - 33)  >= min_qscore:
						qlist.append(qual)
						#print(qlist)

				# check that all quality scores for that index are
				# above quality limit
				if len(qlist) == len(j):

					# identify that index pair
					pair = seq[1],seq[2]

					# check if the index pair is a possible combination
					if pair in all_combinations_dict:

						# increment count of that combination
						all_combinations_dict[pair] += 1

		# increment line number
		NL += 1

# open a file to hold the matched index pair counts
# and write a first line
match = open('ind_out.tsv','w')
match.write('Matched Index Pair\tCounts\n')

# open a file to hold the swapped index pair counts
# and write a first line
swap = open('swapped_out.tsv', 'w')
swap.write('Swapped Index Pair\tCounts\n')


# wite matched and swapped indecies to their respective files
for comb in all_combinations_dict:
	if comb[0]==comb[1]:
		match.write(str(comb)+'\t'+
			str(all_combinations_dict[comb])+'\n')
	else:
		swap.write(str(comb)+'\t'+str(all_combinations_dict[comb])+'\n')

# how many combinations were recorded
total =  0
for comb in all_combinations_dict:
	total += all_combinations_dict[comb]

print('Total combinations: ', str(total))