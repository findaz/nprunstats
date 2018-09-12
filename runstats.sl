#!/bin/bash -l

#SBATCH -p regular
#SBATCH -t 01:30:00
#SBATCH -J stats.dir
#SBATCH -o runstats.dir.out
#SBATCH -e runstats.dir.err
#SBATCH -N 1 
#SBATCH -C haswell
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user
#SBATCH --account=desi

##SBATCH -p shared
##SBATCH -t 03:00:00
##SBATCH -J basschute
##SBATCH -n 1 
##SBATCH --mem=10GB
##SBATCH -C haswell

##SBATCH -p debug
##SBATCH -t 00:05:00
##SBATCH -N 1
##SBATCH -C haswell
##SBATCH -o runstats.out
##SBATCH -e runstats.err

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

if [ "$NERSC_HOST" == "edison" ]
then
	if [[ -z "$NPROC" ]]; then
	   export NPROC=24
	   fi
fi

if [ "$NERSC_HOST" == "cori" ]
then
	if [[ -z "$NPROC" ]]; then
	    export NPROC=32
	   fi
fi

srun -n 1 -c $NPROC make runstats

