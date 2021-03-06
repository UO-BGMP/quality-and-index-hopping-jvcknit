# Quality-and-index-hopping-jvcknit

To asses the quality of the reads from the following fastq files of foward and reverse reads, and their associated indices:

```
1294_S1_L008_R1_001.fastq
1294_S1_L008_R2_001.fastq
1294_S1_L008_R3_001.fastq
1294_S1_L008_R4_001.fastq
```

Quality score distributions were generated for each file by running `qscore_dist2.py` on each of the files. The program was run on the talapas high performance cluster using the following SLURM script for each file.

```
#!/bin/bash
#SBATCH --partition=fat     ### partition
#SBATCH --job-name=QD1      ### Job Name
#SBATCH --output=QD1.out    ### Job output
#SBATCH --error=QD1.err     ### Job error
#SBATCH --time=0-24:00:00   ### Time limit in Days-HH:MM:SS
#SBATCH --nodes=1              ### Node count required for the job
#SBATCH --ntasks-per-node=28   ### Nuber of tasks/Node
#SBATCH --mail-user=jvancamp@uoregon.edu     ### Email
#SBATCH --mail-type=ALL                      ### When start/stop/error

# load the latest version of zlib
ml zlib/1.2.11

# load the latest python
ml python3/3.6.1

# run qscore_dist2.py on read1
./qscore_dist2.py -f 1294_S1_L008_R1_001.fastq
```

# R1 plot
![R1](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R1_001.fastq_dist1.png)

# R1 plot
![R1](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R1_001.fastq_dist2.png)

# R2 plot
![R2](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R2_001.fastq_dist1.png)

# R2 plot
![R2](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R2_001.fastq_dist2.png)

# R3 plot
![R3](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R3_001.fastq_dist1.png)

# R3 plot
![R3](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R3_001.fastq_dist2.png)

# R4 plot
![R4](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R4_001.fastq_dist1.png)

# R4 plot
![R4](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R4_001.fastq_dist2.png)


The above plots show the largest proportion of quality scores above an average quality score of 30 for the whole read. This seems like a reasonable starting quality score threshold for quality filtering of the reads. A more strict threshold of 35 could be applied and many of the reads would still be retained. Here, reads were filtered using both a quality score threshold of 30, and 35.


## Quality Score Filtering, and index hopping

The program `index_hop2.py` was run on the fastq files. The program was run on the talapas high performance computer using the following SLURM script with a quality cutoff of 30 (flag `-c 30` ).

```
#!/bin/bash
#SBATCH --partition=fat        ### partition
#SBATCH --job-name=ihop        ### Job Name
#SBATCH --output=ihop.out      ### Job output
#SBATCH --error=ihop.err       ### Job error
#SBATCH --time=0-24:00:00      ### time limit, Days-HH:MM:SS
#SBATCH --nodes=1              ### Node count
#SBATCH --ntasks-per-node=56   ### Nuber of tasks/node

cd ../

ml Python/3.6.1
ml matplotlib/2.0.1-Python-3.6.1

# hop the indexes
./index_hop2_1.py -1 1294_S1_L008_R1_001.fastq \
-2 1294_S1_L008_R2_001.fastq \
-3 1294_S1_L008_R3_001.fastq \
-4 1294_S1_L008_R4_001.fastq \
-i indexes.txt -c 30
```

The reults of the output at are shown below:

```
$cat match_out.tsv
Matched Index Pair	Counts
GTAGCGTA_GTAGCGTA	15170760
CGATCGAT_CGATCGAT	10424338
GATCAAGG_GATCAAGG	12168683
AACAGCGA_AACAGCGA	16933442
TAGCCATG_TAGCCATG	19843352
CGGTAATC_CGGTAATC	9246665
CTCTGGAT_CTCTGGAT	65953935
TACCGGAT_TACCGGAT	146346803
CTAGCTCA_CTAGCTCA	31980971
CACTTCAC_CACTTCAC	7737998
GCTACTCT_GCTACTCT	13378747
ACGATCAG_ACGATCAG	14671490
TATGGCAC_TATGGCAC	20887249
TGTTCCGT_TGTTCCGT	29163416
GTCCTAAG_GTCCTAAG	16322301
TCGACAAG_TCGACAAG	7194025
TCTTCGAC_TCTTCGAC	79631086
ATCATGCG_ATCATGCG	18464234
ATCGTGGT_ATCGTGGT	12680245
TCGAGAGT_TCGAGAGT	21649368
TCGGATTC_TCGGATTC	8023188
GATCTTGC_GATCTTGC	6790575
AGAGTCCA_AGAGTCCA	20254740
AGGATAGC_AGGATAGC	16128732
```


```
$ cat swapped_out.tsv
Swapped Index Pair	Counts
GTAGCGTA_CGATCGAT	307
GTAGCGTA_GATCAAGG	385
GTAGCGTA_AACAGCGA	513
GTAGCGTA_TAGCCATG	458
GTAGCGTA_CGGTAATC	227
GTAGCGTA_CTCTGGAT	1344
GTAGCGTA_TACCGGAT	2869
GTAGCGTA_CTAGCTCA	1293
GTAGCGTA_CACTTCAC	173
GTAGCGTA_GCTACTCT	2789
GTAGCGTA_ACGATCAG	1380
GTAGCGTA_TATGGCAC	537
GTAGCGTA_TGTTCCGT	777
GTAGCGTA_GTCCTAAG	621
GTAGCGTA_TCGACAAG	175
GTAGCGTA_TCTTCGAC	1910
GTAGCGTA_ATCATGCG	495
GTAGCGTA_ATCGTGGT	407
GTAGCGTA_TCGAGAGT	1214
GTAGCGTA_TCGGATTC	236
GTAGCGTA_GATCTTGC	224
GTAGCGTA_AGAGTCCA	778
GTAGCGTA_AGGATAGC	435
CGATCGAT_GATCAAGG	281
CGATCGAT_AACAGCGA	256
CGATCGAT_TAGCCATG	366
CGATCGAT_CGGTAATC	349
CGATCGAT_CTCTGGAT	2098
CGATCGAT_TACCGGAT	3369
CGATCGAT_CTAGCTCA	643
CGATCGAT_CACTTCAC	215
CGATCGAT_GCTACTCT	327
CGATCGAT_ACGATCAG	465
CGATCGAT_TATGGCAC	317
CGATCGAT_TGTTCCGT	785
CGATCGAT_GTCCTAAG	355
CGATCGAT_TCGACAAG	181
CGATCGAT_TCTTCGAC	1877
CGATCGAT_ATCATGCG	317
CGATCGAT_ATCGTGGT	269
CGATCGAT_TCGAGAGT	471
CGATCGAT_TCGGATTC	182
CGATCGAT_GATCTTGC	204
CGATCGAT_AGAGTCCA	436
CGATCGAT_AGGATAGC	516
GATCAAGG_AACAGCGA	306
GATCAAGG_TAGCCATG	566
GATCAAGG_CGGTAATC	188
GATCAAGG_CTCTGGAT	1098
GATCAAGG_TACCGGAT	2385
GATCAAGG_CTAGCTCA	469
GATCAAGG_CACTTCAC	153
GATCAAGG_GCTACTCT	373
GATCAAGG_ACGATCAG	933
GATCAAGG_TATGGCAC	385
GATCAAGG_TGTTCCGT	543
GATCAAGG_GTCCTAAG	657
GATCAAGG_TCGACAAG	229
GATCAAGG_TCTTCGAC	33103
GATCAAGG_ATCATGCG	443
GATCAAGG_ATCGTGGT	692
GATCAAGG_TCGAGAGT	422
GATCAAGG_TCGGATTC	171
GATCAAGG_GATCTTGC	340
GATCAAGG_AGAGTCCA	425
GATCAAGG_AGGATAGC	379
AACAGCGA_TAGCCATG	557
AACAGCGA_CGGTAATC	274
AACAGCGA_CTCTGGAT	1961
AACAGCGA_TACCGGAT	3836
AACAGCGA_CTAGCTCA	987
AACAGCGA_CACTTCAC	290
AACAGCGA_GCTACTCT	377
AACAGCGA_ACGATCAG	852
AACAGCGA_TATGGCAC	746
AACAGCGA_TGTTCCGT	886
AACAGCGA_GTCCTAAG	555
AACAGCGA_TCGACAAG	213
AACAGCGA_TCTTCGAC	3542
AACAGCGA_ATCATGCG	1022
AACAGCGA_ATCGTGGT	650
AACAGCGA_TCGAGAGT	673
AACAGCGA_TCGGATTC	639
AACAGCGA_GATCTTGC	1078
AACAGCGA_AGAGTCCA	872
AACAGCGA_AGGATAGC	826
TAGCCATG_CGGTAATC	278
TAGCCATG_CTCTGGAT	1575
TAGCCATG_TACCGGAT	3988
TAGCCATG_CTAGCTCA	1012
TAGCCATG_CACTTCAC	429
TAGCCATG_GCTACTCT	332
TAGCCATG_ACGATCAG	639
TAGCCATG_TATGGCAC	699
TAGCCATG_TGTTCCGT	865
TAGCCATG_GTCCTAAG	672
TAGCCATG_TCGACAAG	333
TAGCCATG_TCTTCGAC	2212
TAGCCATG_ATCATGCG	651
TAGCCATG_ATCGTGGT	311
TAGCCATG_TCGAGAGT	606
TAGCCATG_TCGGATTC	308
TAGCCATG_GATCTTGC	248
TAGCCATG_AGAGTCCA	594
TAGCCATG_AGGATAGC	490
CGGTAATC_CTCTGGAT	963
CGGTAATC_TACCGGAT	12917
CGGTAATC_CTAGCTCA	478
CGGTAATC_CACTTCAC	333
CGGTAATC_GCTACTCT	197
CGGTAATC_ACGATCAG	263
CGGTAATC_TATGGCAC	388
CGGTAATC_TGTTCCGT	367
CGGTAATC_GTCCTAAG	193
CGGTAATC_TCGACAAG	113
CGGTAATC_TCTTCGAC	1844
CGGTAATC_ATCATGCG	248
CGGTAATC_ATCGTGGT	594
CGGTAATC_TCGAGAGT	284
CGGTAATC_TCGGATTC	234
CGGTAATC_GATCTTGC	188
CGGTAATC_AGAGTCCA	272
CGGTAATC_AGGATAGC	461
CTCTGGAT_TACCGGAT	23014
CTCTGGAT_CTAGCTCA	4845
CTCTGGAT_CACTTCAC	1135
CTCTGGAT_GCTACTCT	4358
CTCTGGAT_ACGATCAG	2076
CTCTGGAT_TATGGCAC	2674
CTCTGGAT_TGTTCCGT	3961
CTCTGGAT_GTCCTAAG	1984
CTCTGGAT_TCGACAAG	912
CTCTGGAT_TCTTCGAC	9650
CTCTGGAT_ATCATGCG	1911
CTCTGGAT_ATCGTGGT	1803
CTCTGGAT_TCGAGAGT	3690
CTCTGGAT_TCGGATTC	1231
CTCTGGAT_GATCTTGC	951
CTCTGGAT_AGAGTCCA	3217
CTCTGGAT_AGGATAGC	1720
TACCGGAT_CTAGCTCA	7085
TACCGGAT_CACTTCAC	2317
TACCGGAT_GCTACTCT	2789
TACCGGAT_ACGATCAG	4159
TACCGGAT_TATGGCAC	6747
TACCGGAT_TGTTCCGT	8312
TACCGGAT_GTCCTAAG	4538
TACCGGAT_TCGACAAG	1807
TACCGGAT_TCTTCGAC	20721
TACCGGAT_ATCATGCG	3765
TACCGGAT_ATCGTGGT	4817
TACCGGAT_TCGAGAGT	7115
TACCGGAT_TCGGATTC	2504
TACCGGAT_GATCTTGC	1912
TACCGGAT_AGAGTCCA	4103
TACCGGAT_AGGATAGC	3767
CTAGCTCA_CACTTCAC	527
CTAGCTCA_GCTACTCT	3075
CTAGCTCA_ACGATCAG	1247
CTAGCTCA_TATGGCAC	1346
CTAGCTCA_TGTTCCGT	1496
CTAGCTCA_GTCCTAAG	1021
CTAGCTCA_TCGACAAG	27104
CTAGCTCA_TCTTCGAC	3977
CTAGCTCA_ATCATGCG	3085
CTAGCTCA_ATCGTGGT	713
CTAGCTCA_TCGAGAGT	1713
CTAGCTCA_TCGGATTC	639
CTAGCTCA_GATCTTGC	544
CTAGCTCA_AGAGTCCA	1720
CTAGCTCA_AGGATAGC	962
CACTTCAC_GCTACTCT	160
CACTTCAC_ACGATCAG	633
CACTTCAC_TATGGCAC	530
CACTTCAC_TGTTCCGT	405
CACTTCAC_GTCCTAAG	186
CACTTCAC_TCGACAAG	115
CACTTCAC_TCTTCGAC	1478
CACTTCAC_ATCATGCG	215
CACTTCAC_ATCGTGGT	156
CACTTCAC_TCGAGAGT	238
CACTTCAC_TCGGATTC	221
CACTTCAC_GATCTTGC	154
CACTTCAC_AGAGTCCA	260
CACTTCAC_AGGATAGC	341
GCTACTCT_ACGATCAG	569
GCTACTCT_TATGGCAC	466
GCTACTCT_TGTTCCGT	905
GCTACTCT_GTCCTAAG	485
GCTACTCT_TCGACAAG	150
GCTACTCT_TCTTCGAC	1607
GCTACTCT_ATCATGCG	463
GCTACTCT_ATCGTGGT	327
GCTACTCT_TCGAGAGT	624
GCTACTCT_TCGGATTC	201
GCTACTCT_GATCTTGC	246
GCTACTCT_AGAGTCCA	696
GCTACTCT_AGGATAGC	396
ACGATCAG_TATGGCAC	670
ACGATCAG_TGTTCCGT	859
ACGATCAG_GTCCTAAG	783
ACGATCAG_TCGACAAG	542
ACGATCAG_TCTTCGAC	2013
ACGATCAG_ATCATGCG	665
ACGATCAG_ATCGTGGT	356
ACGATCAG_TCGAGAGT	776
ACGATCAG_TCGGATTC	267
ACGATCAG_GATCTTGC	279
ACGATCAG_AGAGTCCA	760
ACGATCAG_AGGATAGC	639
TATGGCAC_TGTTCCGT	164471
TATGGCAC_GTCCTAAG	1282
TATGGCAC_TCGACAAG	450
TATGGCAC_TCTTCGAC	7081
TATGGCAC_ATCATGCG	738
TATGGCAC_ATCGTGGT	492
TATGGCAC_TCGAGAGT	1488
TATGGCAC_TCGGATTC	1682
TATGGCAC_GATCTTGC	598
TATGGCAC_AGAGTCCA	769
TATGGCAC_AGGATAGC	1071
TGTTCCGT_GTCCTAAG	706
TGTTCCGT_TCGACAAG	373
TGTTCCGT_TCTTCGAC	4175
TGTTCCGT_ATCATGCG	1119
TGTTCCGT_ATCGTGGT	755
TGTTCCGT_TCGAGAGT	1513
TGTTCCGT_TCGGATTC	1558
TGTTCCGT_GATCTTGC	535
TGTTCCGT_AGAGTCCA	1270
TGTTCCGT_AGGATAGC	897
GTCCTAAG_TCGACAAG	290
GTCCTAAG_TCTTCGAC	1876
GTCCTAAG_ATCATGCG	601
GTCCTAAG_ATCGTGGT	337
GTCCTAAG_TCGAGAGT	525
GTCCTAAG_TCGGATTC	432
GTCCTAAG_GATCTTGC	318
GTCCTAAG_AGAGTCCA	530
GTCCTAAG_AGGATAGC	438
TCGACAAG_TCTTCGAC	1107
TCGACAAG_ATCATGCG	14458
TCGACAAG_ATCGTGGT	150
TCGACAAG_TCGAGAGT	550
TCGACAAG_TCGGATTC	179
TCGACAAG_GATCTTGC	97
TCGACAAG_AGAGTCCA	218
TCGACAAG_AGGATAGC	266
TCTTCGAC_ATCATGCG	2326
TCTTCGAC_ATCGTGGT	8023
TCTTCGAC_TCGAGAGT	4028
TCTTCGAC_TCGGATTC	2516
TCTTCGAC_GATCTTGC	1612
TCTTCGAC_AGAGTCCA	2408
TCTTCGAC_AGGATAGC	3358
ATCATGCG_ATCGTGGT	799
ATCATGCG_TCGAGAGT	692
ATCATGCG_TCGGATTC	309
ATCATGCG_GATCTTGC	408
ATCATGCG_AGAGTCCA	908
ATCATGCG_AGGATAGC	889
ATCGTGGT_TCGAGAGT	806
ATCGTGGT_TCGGATTC	235
ATCGTGGT_GATCTTGC	210
ATCGTGGT_AGAGTCCA	506
ATCGTGGT_AGGATAGC	417
TCGAGAGT_TCGGATTC	603
TCGAGAGT_GATCTTGC	333
TCGAGAGT_AGAGTCCA	750
TCGAGAGT_AGGATAGC	823
TCGGATTC_GATCTTGC	441
TCGGATTC_AGAGTCCA	291
TCGGATTC_AGGATAGC	439
GATCTTGC_AGAGTCCA	299
GATCTTGC_AGGATAGC	598
AGAGTCCA_AGGATAGC	778
```

```
$ cat undetermined_out.tsv | wc -l
209507
```

There were 209507 undetermined index pairs due to low quality base positions read from index pairs resulting in "N"s present in the index pairs. To determine the number of "N" present in sequence lines of the raw files, I used the following bash command:


```
$ for file in *.fastq; do cat $file | awk 'NR%4==2 {print $0}' | grep 'N' | wc -l ; done
2602560
3976613
3328051
3591851
```

This represents a range from


To asses how qualiity score filtering affected the level of index hoping `index_hop2.py` was run on the files again wtih a q-score threshold of 35 using the flag `-c 35` producing the following results:

```
$ cat match_out.tsv
Matched Index Pair	Counts
GTAGCGTA_GTAGCGTA	7208537
CGATCGAT_CGATCGAT	4933975
GATCAAGG_GATCAAGG	5796436
AACAGCGA_AACAGCGA	8077276
TAGCCATG_TAGCCATG	9532164
CGGTAATC_CGGTAATC	4392347
CTCTGGAT_CTCTGGAT	31425190
TACCGGAT_TACCGGAT	70030710
CTAGCTCA_CTAGCTCA	15215680
CACTTCAC_CACTTCAC	3684205
GCTACTCT_GCTACTCT	6237707
ACGATCAG_ACGATCAG	6972608
TATGGCAC_TATGGCAC	9981360
TGTTCCGT_TGTTCCGT	13628726
GTCCTAAG_GTCCTAAG	7798019
TCGACAAG_TCGACAAG	3438640
TCTTCGAC_TCTTCGAC	38401642
ATCATGCG_ATCATGCG	8700512
ATCGTGGT_ATCGTGGT	5823446
TCGAGAGT_TCGAGAGT	9533559
TCGGATTC_TCGGATTC	3711662
GATCTTGC_GATCTTGC	3236404
AGAGTCCA_AGAGTCCA	9418815
AGGATAGC_AGGATAGC	7670823
```

```
$ cat swapped_out.tsv
Swapped Index Pair	Counts
GTAGCGTA_CGATCGAT	136
GTAGCGTA_GATCAAGG	174
GTAGCGTA_AACAGCGA	230
GTAGCGTA_TAGCCATG	202
GTAGCGTA_CGGTAATC	103
GTAGCGTA_CTCTGGAT	588
GTAGCGTA_TACCGGAT	1286
GTAGCGTA_CTAGCTCA	543
GTAGCGTA_CACTTCAC	81
GTAGCGTA_GCTACTCT	1322
GTAGCGTA_ACGATCAG	635
GTAGCGTA_TATGGCAC	234
GTAGCGTA_TGTTCCGT	345
GTAGCGTA_GTCCTAAG	269
GTAGCGTA_TCGACAAG	81
GTAGCGTA_TCTTCGAC	834
GTAGCGTA_ATCATGCG	225
GTAGCGTA_ATCGTGGT	183
GTAGCGTA_TCGAGAGT	554
GTAGCGTA_TCGGATTC	108
GTAGCGTA_GATCTTGC	99
GTAGCGTA_AGAGTCCA	356
GTAGCGTA_AGGATAGC	188
CGATCGAT_GATCAAGG	121
CGATCGAT_AACAGCGA	112
CGATCGAT_TAGCCATG	163
CGATCGAT_CGGTAATC	161
CGATCGAT_CTCTGGAT	928
CGATCGAT_TACCGGAT	1519
CGATCGAT_CTAGCTCA	277
CGATCGAT_CACTTCAC	99
CGATCGAT_GCTACTCT	142
CGATCGAT_ACGATCAG	208
CGATCGAT_TATGGCAC	143
CGATCGAT_TGTTCCGT	343
CGATCGAT_GTCCTAAG	154
CGATCGAT_TCGACAAG	82
CGATCGAT_TCTTCGAC	825
CGATCGAT_ATCATGCG	139
CGATCGAT_ATCGTGGT	123
CGATCGAT_TCGAGAGT	208
CGATCGAT_TCGGATTC	82
CGATCGAT_GATCTTGC	93
CGATCGAT_AGAGTCCA	189
CGATCGAT_AGGATAGC	237
GATCAAGG_AACAGCGA	136
GATCAAGG_TAGCCATG	254
GATCAAGG_CGGTAATC	84
GATCAAGG_CTCTGGAT	494
GATCAAGG_TACCGGAT	1086
GATCAAGG_CTAGCTCA	204
GATCAAGG_CACTTCAC	69
GATCAAGG_GCTACTCT	166
GATCAAGG_ACGATCAG	430
GATCAAGG_TATGGCAC	172
GATCAAGG_TGTTCCGT	240
GATCAAGG_GTCCTAAG	293
GATCAAGG_TCGACAAG	104
GATCAAGG_TCTTCGAC	15743
GATCAAGG_ATCATGCG	202
GATCAAGG_ATCGTGGT	318
GATCAAGG_TCGAGAGT	187
GATCAAGG_TCGGATTC	75
GATCAAGG_GATCTTGC	145
GATCAAGG_AGAGTCCA	194
GATCAAGG_AGGATAGC	168
AACAGCGA_TAGCCATG	252
AACAGCGA_CGGTAATC	119
AACAGCGA_CTCTGGAT	871
AACAGCGA_TACCGGAT	1682
AACAGCGA_CTAGCTCA	442
AACAGCGA_CACTTCAC	128
AACAGCGA_GCTACTCT	168
AACAGCGA_ACGATCAG	381
AACAGCGA_TATGGCAC	322
AACAGCGA_TGTTCCGT	379
AACAGCGA_GTCCTAAG	243
AACAGCGA_TCGACAAG	94
AACAGCGA_TCTTCGAC	1597
AACAGCGA_ATCATGCG	451
AACAGCGA_ATCGTGGT	277
AACAGCGA_TCGAGAGT	287
AACAGCGA_TCGGATTC	289
AACAGCGA_GATCTTGC	499
AACAGCGA_AGAGTCCA	389
AACAGCGA_AGGATAGC	365
TAGCCATG_CGGTAATC	119
TAGCCATG_CTCTGGAT	711
TAGCCATG_TACCGGAT	1797
TAGCCATG_CTAGCTCA	444
TAGCCATG_CACTTCAC	199
TAGCCATG_GCTACTCT	151
TAGCCATG_ACGATCAG	290
TAGCCATG_TATGGCAC	318
TAGCCATG_TGTTCCGT	386
TAGCCATG_GTCCTAAG	307
TAGCCATG_TCGACAAG	151
TAGCCATG_TCTTCGAC	991
TAGCCATG_ATCATGCG	296
TAGCCATG_ATCGTGGT	137
TAGCCATG_TCGAGAGT	254
TAGCCATG_TCGGATTC	133
TAGCCATG_GATCTTGC	112
TAGCCATG_AGAGTCCA	264
TAGCCATG_AGGATAGC	223
CGGTAATC_CTCTGGAT	434
CGGTAATC_TACCGGAT	6139
CGGTAATC_CTAGCTCA	210
CGGTAATC_CACTTCAC	156
CGGTAATC_GCTACTCT	89
CGGTAATC_ACGATCAG	119
CGGTAATC_TATGGCAC	174
CGGTAATC_TGTTCCGT	156
CGGTAATC_GTCCTAAG	87
CGGTAATC_TCGACAAG	51
CGGTAATC_TCTTCGAC	841
CGGTAATC_ATCATGCG	108
CGGTAATC_ATCGTGGT	274
CGGTAATC_TCGAGAGT	125
CGGTAATC_TCGGATTC	103
CGGTAATC_GATCTTGC	85
CGGTAATC_AGAGTCCA	120
CGGTAATC_AGGATAGC	203
CTCTGGAT_TACCGGAT	10346
CTCTGGAT_CTAGCTCA	2151
CTCTGGAT_CACTTCAC	507
CTCTGGAT_GCTACTCT	2041
CTCTGGAT_ACGATCAG	921
CTCTGGAT_TATGGCAC	1180
CTCTGGAT_TGTTCCGT	1750
CTCTGGAT_GTCCTAAG	861
CTCTGGAT_TCGACAAG	413
CTCTGGAT_TCTTCGAC	4194
CTCTGGAT_ATCATGCG	839
CTCTGGAT_ATCGTGGT	776
CTCTGGAT_TCGAGAGT	1651
CTCTGGAT_TCGGATTC	560
CTCTGGAT_GATCTTGC	416
CTCTGGAT_AGAGTCCA	1451
CTCTGGAT_AGGATAGC	772
TACCGGAT_CTAGCTCA	3063
TACCGGAT_CACTTCAC	1012
TACCGGAT_GCTACTCT	1244
TACCGGAT_ACGATCAG	1790
TACCGGAT_TATGGCAC	2806
TACCGGAT_TGTTCCGT	3631
TACCGGAT_GTCCTAAG	1835
TACCGGAT_TCGACAAG	774
TACCGGAT_TCTTCGAC	8865
TACCGGAT_ATCATGCG	1639
TACCGGAT_ATCGTGGT	1967
TACCGGAT_TCGAGAGT	2815
TACCGGAT_TCGGATTC	1071
TACCGGAT_GATCTTGC	825
TACCGGAT_AGAGTCCA	1798
TACCGGAT_AGGATAGC	1668
CTAGCTCA_CACTTCAC	234
CTAGCTCA_GCTACTCT	1429
CTAGCTCA_ACGATCAG	555
CTAGCTCA_TATGGCAC	593
CTAGCTCA_TGTTCCGT	661
CTAGCTCA_GTCCTAAG	465
CTAGCTCA_TCGACAAG	12820
CTAGCTCA_TCTTCGAC	1764
CTAGCTCA_ATCATGCG	1444
CTAGCTCA_ATCGTGGT	319
CTAGCTCA_TCGAGAGT	780
CTAGCTCA_TCGGATTC	287
CTAGCTCA_GATCTTGC	239
CTAGCTCA_AGAGTCCA	770
CTAGCTCA_AGGATAGC	434
CACTTCAC_GCTACTCT	73
CACTTCAC_ACGATCAG	290
CACTTCAC_TATGGCAC	237
CACTTCAC_TGTTCCGT	177
CACTTCAC_GTCCTAAG	79
CACTTCAC_TCGACAAG	55
CACTTCAC_TCTTCGAC	659
CACTTCAC_ATCATGCG	96
CACTTCAC_ATCGTGGT	68
CACTTCAC_TCGAGAGT	102
CACTTCAC_TCGGATTC	103
CACTTCAC_GATCTTGC	69
CACTTCAC_AGAGTCCA	113
CACTTCAC_AGGATAGC	155
GCTACTCT_ACGATCAG	255
GCTACTCT_TATGGCAC	201
GCTACTCT_TGTTCCGT	394
GCTACTCT_GTCCTAAG	210
GCTACTCT_TCGACAAG	63
GCTACTCT_TCTTCGAC	696
GCTACTCT_ATCATGCG	198
GCTACTCT_ATCGTGGT	145
GCTACTCT_TCGAGAGT	276
GCTACTCT_TCGGATTC	88
GCTACTCT_GATCTTGC	114
GCTACTCT_AGAGTCCA	311
GCTACTCT_AGGATAGC	174
ACGATCAG_TATGGCAC	303
ACGATCAG_TGTTCCGT	390
ACGATCAG_GTCCTAAG	355
ACGATCAG_TCGACAAG	250
ACGATCAG_TCTTCGAC	910
ACGATCAG_ATCATGCG	289
ACGATCAG_ATCGTGGT	162
ACGATCAG_TCGAGAGT	346
ACGATCAG_TCGGATTC	117
ACGATCAG_GATCTTGC	121
ACGATCAG_AGAGTCCA	335
ACGATCAG_AGGATAGC	290
TATGGCAC_TGTTCCGT	78436
TATGGCAC_GTCCTAAG	594
TATGGCAC_TCGACAAG	203
TATGGCAC_TCTTCGAC	3253
TATGGCAC_ATCATGCG	333
TATGGCAC_ATCGTGGT	216
TATGGCAC_TCGAGAGT	672
TATGGCAC_TCGGATTC	782
TATGGCAC_GATCTTGC	268
TATGGCAC_AGAGTCCA	340
TATGGCAC_AGGATAGC	491
TGTTCCGT_GTCCTAAG	303
TGTTCCGT_TCGACAAG	167
TGTTCCGT_TCTTCGAC	1798
TGTTCCGT_ATCATGCG	493
TGTTCCGT_ATCGTGGT	335
TGTTCCGT_TCGAGAGT	649
TGTTCCGT_TCGGATTC	717
TGTTCCGT_GATCTTGC	248
TGTTCCGT_AGAGTCCA	564
TGTTCCGT_AGGATAGC	382
GTCCTAAG_TCGACAAG	132
GTCCTAAG_TCTTCGAC	836
GTCCTAAG_ATCATGCG	269
GTCCTAAG_ATCGTGGT	147
GTCCTAAG_TCGAGAGT	235
GTCCTAAG_TCGGATTC	200
GTCCTAAG_GATCTTGC	143
GTCCTAAG_AGAGTCCA	234
GTCCTAAG_AGGATAGC	195
TCGACAAG_TCTTCGAC	488
TCGACAAG_ATCATGCG	6898
TCGACAAG_ATCGTGGT	70
TCGACAAG_TCGAGAGT	240
TCGACAAG_TCGGATTC	82
TCGACAAG_GATCTTGC	43
TCGACAAG_AGAGTCCA	98
TCGACAAG_AGGATAGC	117
TCTTCGAC_ATCATGCG	1066
TCTTCGAC_ATCGTGGT	3814
TCTTCGAC_TCGAGAGT	1796
TCTTCGAC_TCGGATTC	1144
TCTTCGAC_GATCTTGC	746
TCTTCGAC_AGAGTCCA	1095
TCTTCGAC_AGGATAGC	1546
ATCATGCG_ATCGTGGT	347
ATCATGCG_TCGAGAGT	308
ATCATGCG_TCGGATTC	137
ATCATGCG_GATCTTGC	180
ATCATGCG_AGAGTCCA	412
ATCATGCG_AGGATAGC	387
ATCGTGGT_TCGAGAGT	346
ATCGTGGT_TCGGATTC	105
ATCGTGGT_GATCTTGC	95
ATCGTGGT_AGAGTCCA	218
ATCGTGGT_AGGATAGC	181
TCGAGAGT_TCGGATTC	253
TCGAGAGT_GATCTTGC	132
TCGAGAGT_AGAGTCCA	299
TCGAGAGT_AGGATAGC	323
TCGGATTC_GATCTTGC	196
TCGGATTC_AGAGTCCA	125
TCGGATTC_AGGATAGC	188
GATCTTGC_AGAGTCCA	133
GATCTTGC_AGGATAGC	282
AGAGTCCA_AGGATAGC	339
```

```
$ cat undetermined_out.tsv | wc -l
141636
```


It is evident that more reads were thrown out because of the more strict cutoff by looking at the number of undetermined reads from 209507 undetermined reads at a cutoff of 30 to 141636 at a cutoff of 35.


A deeper look into the swapping of reads can be found in swapping.pdf and swapping.Rmd



