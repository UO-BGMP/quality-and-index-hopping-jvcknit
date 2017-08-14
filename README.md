# Quality-and-index-hopping-jvcknit

The following plots were generated by running the `qscore_dist.py` on the test files in the following slurm script:

```
#!/bin/bash
#SBATCH --partition=long            ### partition
#SBATCH --job-name=qscore_dist      ### Job Name
#SBATCH --output=qscore_dist.out    ### File in which to store job output
#SBATCH --error=qscore_dist.err     ### File in which to store job error messages
#SBATCH --time=0-02:00:00           ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   ### Node count required for the job
#SBATCH --ntasks-per-node=28        ### Nuber of tasks to be launched per Node
#SBATCH --mail-user=jvancamp@uoregon.edu        ### Notifyme
#SBATCH --mail-type=ALL                         ### All of it

# load the latest version of zlib
ml zlib/1.2.11

# cd into the plots directory
cd ~/index_hop/plots

# get quality score distributions for each file
~/index_hop/qscore_dist.py -f /projects/bgmp/test_data/Undetermined_S0_R1_001_trunc.fastq.gz -p R1

~/index_hop/qscore_dist.py -f /projects/bgmp/test_data/Undetermined_S0_R2_001_trunc.fastq.gz -p R2

~/index_hop/qscore_dist.py -f /projects/bgmp/test_data/Undetermined_S0_R3_001_trunc.fastq.gz -p R3

~/index_hop/qscore_dist.py -f /projects/bgmp/test_data/Undetermined_S0_R4_001_trunc.fastq.gz -p R4

```


# R1 plot

![R1 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/R1_dist.png)

# R2 plot
![R2 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/R2_dist.png)

# R3 plot
![R3 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/R3_dist.png)

# R4 plot
![R4 plot](https://github.com/UO-BGMP/quality-and-index-hopping-jvcknit/blob/index_hop.working/R4_dist.png)



TO DO:

* figure out why some files have indices below qscore
* add eac


