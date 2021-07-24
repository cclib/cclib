#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-10:10:00 # Runtime in D-HH:MM
#SBATCH --job-name=a
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --output=cfour.out
#SBATCH --cluster=invest
#SBATCH --partition=kjordan
#SBATCH --constraint=smp

export CFOUR_NUM_CORES=12
export OMP_NUM_THREADS=1
module purge
module load mkl/2017.1.132
module load gcc/9.2.0
module use /ihome/kjordan/shu8/software/privatemodules/smp/
module load cfour/cfour2_gcc_mkl2017
cd $SLURM_SUBMIT_DIR
for i in *_ZMAT
do
    cd $SLURM_SUBMIT_DIR
    filename="${i%_ZMAT}"
    echo $filename
    cp $i ZMAT
    cp $SLURM_SUBMIT_DIR/ZMAT $SLURM_SCRATCH
    cp $SLURM_SUBMIT_DIR/GENBAS $SLURM_SCRATCH
    cd $SLURM_SCRATCH
    xcfour &> output
    cp output $SLURM_SUBMIT_DIR
    cd $SLURM_SUBMIT_DIR
    mv output ${filename}_ZMAT.out
done
