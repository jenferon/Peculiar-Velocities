#! /usr/bin/bash
#SBATCH -p shortq,defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=0:30:00
#SBATCH --output=log/%x.%j.o
#SBATCH --error=log/%x.%j.e
#SBATCH --array=0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

source ~/.bash_profile

cd /gpfs01/home/ppxjf3/peculiar_vel/
activate jen

#python make_fg.py
#python noise_norm_SKA.py
#python run_GPR.py
python build_lc.py


n=$(echo "201" | bc)
echo $n

for ((  i = 0 ;  i <= n;  i++  ))
do
    freq=$(echo "75+0.1*$i" | bc)
    echo $freq
    rm -rf "/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/Noise_N512_FOV1.000_${freq}MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits"
    echo "/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/Noise_N512_FOV1.000_${freq}MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits"
done
