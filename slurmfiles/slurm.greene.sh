#!/bin/bash
#SBATCH --job-name="multilayer"               # Job name
#SBATCH --mail-user=jiarong.wu@nyu.edu        # Email address    
#SBATCH --nodes=2                             # Number of nodes requested, 32 each
#SBATCH --time=1:00:00                        # Time limit request
#SBATCH --mail-type=END 
##SBATCH --ntasks-per-node=32
##SBATCH --cpus-per-task=1

# The executable name
EXE=field_PM
# Parameter value
NLAYER=15
LEVEL=10
TEND=200
nu=0.000025
RAND=2
L0=200
Htheta=0.503
# Initial field info
FIELD=P0.02

mkdir ./surface
mkdir ./field
# If the initial field files are in the ../spectra directory
cp ${SCRATCH}/multilayer/spectra/F_kxky_${FIELD}_200m ./F_kxky
cp ${SCRATCH}/multilayer/spectra/kx_200m ./kx
cp ${SCRATCH}/multilayer/spectra/ky_200m ./ky
echo srun ./$EXE NLAYER=$NLAYER LEVEL=$LEVEL TEND=$TEND nu=$nu RAND=${RAND} L0=$L0 Htheta=${Htheta}

module purge
srun --mpi=pmi2 \
    /scratch/work/public/singularity/run-openmpi-4.1.2-ubuntu-22.04.1.bash \
    ./$EXE $NLAYER $LEVEL $TEND $nu $RAND $L0 $Htheta > message 2>&1 

