#!/usr/bin/env python3.6

'''
Parses four fastq files, Read1, Read2, Read3, Read4 by a quality score
cutoff. And determines the level of index swapping observed by sorting
the records into new fastq files, if they were part of a possible
index combination.

'''

import argparse
import gzip
import os
import itertools

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    bases = list(seq)
    #bases.remove('\n')
    bases = [complement[base] for base in bases]
    return ''.join(bases)


def reverse_complement(s):
        return complement(s[::-1])

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

# set min quality score
min_qscore = 10

# open some files
R1 = 'R1.fastq'
R2 = 'R2.fastq'
R3 = 'R3.fastq'
R4 = 'R4.fastq'
Index = 'indexes.txt'

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

					pair = seq[1],reverse_complement(seq[2])


					# make new directory for output files
					try:
						os.removedirs('Output')
						os.mkdir('Output')
						os.chdir('Output')
					except:
						os.chdir('Output')

					# check if the index pair is a possible combination
					if pair in all_combinations_dict:

						# increment count of that combination
						all_combinations_dict[pair] += 1

						# rename head lines with their index pair
						head = [f.split(' ')[0] for f in head]
						h0 = head[0]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h1 = head[1]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h2 = head[2]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h3 = head[3]+ '_'+ str(pair[0]+ '_'+ pair[1])

						heads = [h0,h1,h2,h3]
						Rs = ['R1','R2','R3','R4']
						for e,r in enumerate(Rs):
							r = open(str(pair[0]+ '_'+ pair[1]+'_'+r+'.fastq'), 'w')
							r.write(heads[e]+'\n'+seq[e]+'\n'+plus[e]+'\n'+qline[e]+'\n')
							r.close()

					# if pair isn't a possible combination
					# output to an undetermined records file
					else:
						ud1 = open(
							'Undetermined_index_pair_R1.fastq', 'w'
							)
						ud2 = open(
							'Undetermined_index_pair_R2.fastq', 'w'
							)

						head = [f.split(' ')[0] for f in head]
						h0 = head[0]+ '_'+ str(pair[0]+ '_'+ pair[1])
						h3 = head[3]+ '_'+ str(pair[0]+ '_'+ pair[1])

						ud1.write(
							h0+'\n'+seq[0]+'\n'+plus[0]+'\n'+qline[0]+'\n'
							)
						ud2.write(
							h3+'\n'+seq[3]+'\n'+plus[3]+'\n'+qline[3]+'\n'
							)
						# close undetermined files
						ud1.close()
						ud2.close()


					# go back to parent directory
					os.chdir('../')

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
		match.write(str(comb[0])+'_'+str(comb[1])+'\t'+
			str(all_combinations_dict[comb])+'\n')
	else:
		swap.write(str(comb[0])+'_'+str(comb[1])+'\t'+str(all_combinations_dict[comb])+'\n')

# how many combinations were recorded
total =  0
for comb in all_combinations_dict:
	total += all_combinations_dict[comb]

print(total)




