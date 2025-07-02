#!/bin/bash
#SBATCH --partition=defq,hmemq,shortq,devq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --mem=0
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=jennifer.feron@nottingham.ac.uk
#SBATCH --output=log/%x.%j.o
#SBATCH --error=log/%x.%j.e

module load gcc-uoneasy
module load foss/2022a
module load GSL

cd /gpfs01/home/ppxjf3/peculiar_vel/

./Simfast21/simfast21 ./data


