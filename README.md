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
#SBATCH --partition=fat            ### partition
#SBATCH --job-name=QD1     ### Job Name
#SBATCH --output=QD1.out    ### File in which to store job output
#SBATCH --error=QD1.err     ### File in which to store job error messages
#SBATCH --time=0-24:00:00           ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --ntasks-per-node=28        ### Nuber of tasks to be launched per Node
#SBATCH --mail-user=jvancamp@uoregon.edu        ### Notifyme
#SBATCH --mail-type=ALL                         ### All of it

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



