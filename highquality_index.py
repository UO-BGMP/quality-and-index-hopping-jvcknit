#!/usr/bin/env python3.6

import argparse
import gzip
import itertools

# def get_arguments():
#     '''Define and return command line options.'''
#     parser = argparse.ArgumentParser(prog='qscore_dist',
#         description= 'This program plots a distribution of quality scores over the base position of sequencing reads in a fastq file.')

#     parser.add_argument('-f', '--infile', help='specify input file (path)',
#                         required=True,
#                         type=argparse.FileType('rt', encoding='UTF-8 '))

#     parser.add_argument('-c', '--qcutoff', help='specify the qscore cutoff',
#                             required=True,
#                             type=str)

#     return parser.parse_args()

# define working files

min_qscore = 30

R1 = 'R1.fastq'
R2 = 'R2.fastq'
R3 = 'R3.fastq'
R4 = 'R4.fastq'
Index = 'index.tsv'

# initialize index dictonary
Index_dict = {}

# populate index dictonary with known indices
with open(Index, 'r') as ind:
	for line in ind:
		line = line.strip('\n').split('\t')
		Index_dict[line[0]] = line[1]

# all possible indices
all_indices = []

# build list of all indices
for key in Index_dict:
	all_indices.append(Index_dict[key])

# all combinations of dual indices
all_combinations_list = list(itertools.combinations_with_replacement(all_indices,2))

# dict of all unique combination of indices
all_combinations_dict = {}
for comb in all_combinations_list:
	all_combinations_dict[comb] = 0

with open(R1 ,'r') as r1, open(R2, 'r') as r2, open(R3, 'r') as r3, open(R4, 'r') as r4:
	NR = 0
	for line in zip(r1,r2,r3,r4):
		line = [l.strip() for l in line]
		if NR % 4 == 0:
			head = line
			# print(head)
		if NR % 4 == 1:
			seq = line
		if NR % 4 == 2:
			plus = line
		if NR % 4 == 3:
			qline = line

			for j in qline[1:2]:
				qlist = []
				for qual in j:
					if (ord(qual) - 33)  >= min_qscore:
						qlist.append(qual)
						#print(qlist)
				if len(qlist) == len(j):
					pair = seq[1],seq[2]
					if pair in all_combinations_dict:
						all_combinations_dict[pair] += 1
						head = [f.split(' ')[0] for f in head]
						h1 = head[0]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h2 = head[1]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h3 = head[2]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h4 = head[3]+ '_'+ str(pair[0]+ '_'+ pair[1])

						heads = [h1,h2,h3,h4]
						#print(heads)


		NR += 1

print(all_combinations_dict)

total =  0
for comb in all_combinations_dict:
	total += all_combinations_dict[comb]

print(total)




