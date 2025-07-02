#!/bin/bash
#SBATCH --partition=defq,hmemq,shortq,devq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --mail-user=jennifer.feron@nottingham.ac.uk
#SBATCH --output=log/%x.%j.o
#SBATCH --error=log/%x.%j.e

module load gcc-uoneasy
module load foss/2022a
module load GSL

cd /gpfs01/home/ppxjf3/peculiar_vel/

echo
echo Lightcone
echo Generating the simulation in directory: $1
echo

DIR="${0%/*}"
nu1=109.0 #larger frequency
nu2=108.0 #smaller frequency
del_nu_lc=0.1
FoV_deg=1.0
zmax=13.0  #(zmax-zmin)/del_z MUST be an integer
zmin=11.0
del_z=0.1
Dim=64

date

#echo ------------------------------------------ Velocities ------------------------------------------
#n=$(echo "scale=0; (($zmax-$zmin)/$del_z)" | bc)
#echo $n

#for ((  i = 0 ;  i <= n;  i++  ))
#do
#    z=$(echo "scale=3; $zmax-$del_z*$i" | bc)
#    echo $z
#    mpirun -np ${SLURM_NTASKS} ./Simfast21/vel ./data $z 1
#    mpirun -np ${SLURM_NTASKS} ./Simfast21/vel ./data $z 2
#    mpirun -np ${SLURM_NTASKS} ./Simfast21/vel ./data $z 3
#done

echo ------------------------------------------ Lightcone ------------------------------------------
mpirun -np ${SLURM_NTASKS} ./Simfast21/RSD_LC ./data/ $nu1 $nu2 $del_nu_lc $FoV_deg $Dim $zmin $zmax $del_z

date

# . lightcone_script /gpfs01/home/ppxjf3/peculiar_vel/Simfast21/

#python dat_to_fits.py
