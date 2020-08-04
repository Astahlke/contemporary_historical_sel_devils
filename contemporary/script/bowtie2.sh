#!/bin/bash

#SBATCH --mem=50G
#SBATCH --tasks-per-node=20

cd $SLURM_SUBMIT_DIR

echo "running on "
hostname

source /usr/modules/init/bash

module load python
module load bowtie2/2.1.0
module load samtools

python /mnt/lfs2/soraia/TD/4_combine/Rapture_2018/PHAV5_HT30/bowtie2_lane.py

echo "finished"