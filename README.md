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

![R1 Quality Score Distributuion](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R1_001.fastq_dist1.png)

# R1 plot
![R1 Quality Score distribution plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R1_001.fastq_dist2.png)

# R2 plot
![R2 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R2_001.fastq_dist1.png)

# R2 plot
![R2 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R2_001.fastq_dist2.png)

# R3 plot
![R3 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R3_001.fastq_dist1.png)

# R3 plot
![R3 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R3_001.fastq_dist2.png)

# R4 plot
![R4 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R4_001.fastq_dist1.png)

# R4 plot
![R4 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/Plots/1294_S1_L008_R4_001.fastq_dist2.png)


The above plots show the largest proportion of quality scores above an average quality score of 30 for the whole read making this a reasonable quality score threshold for quality filtering of the reads. A more strict threshold of 35 could be applied and many of the reads would still be retained.

The program `index_hop2.py` was run on the fastq files. The program was run on the talapas high performance computer with the following SLURM script.

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

To asses how qualiity score filtering affected the level of index hoping `index_hop2.py` was run on the files again wtih a q-score threshold of 35.

The reults from the run produced the following results.


To determine the number of reads containing "N's"

$ for file in *.fastq; do cat $file | awk 'NR%4==2 {print $0}' | grep 'N' | wc -l ; done
2602560
3976613
3328051
3591851




