#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=stokes_ml
#Parameter value
NLAYER=30
N=256
RE=40000
ak=0.35
thetaH=0.503 # A numerical parameter that controls the dissipation of non-barotropic component

mkdir surface
mkdir field
srun ./$EXE $RE $NLAYER $N $ak $thetaH

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
