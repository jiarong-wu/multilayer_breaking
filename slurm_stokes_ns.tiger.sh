#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=stokes_ns
LEVEL=9
ak=0.35
BO=200
RE=40000

mkdir eta
mkdir field
srun ./$EXE $LEVEL $ak $BO $RE

#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
