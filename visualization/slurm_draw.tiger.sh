#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jiarongw@princeton.edu

#The executable name
EXE=stokes_ns_draw
LOC=0 #location where to take the slice 

for i in `seq 3 5`
do 
    let t=0+i*1
    echo $t
    srun ./$EXE $t $LOC 
    convert slicey_ux.ppm slicey_ux.png
    convert slicex_uy.ppm slicex_uy.png
    convert slicey_uy.ppm slicey_uy.png
    mv slicey_ux.png slicey_ux_t$t.png
    mv slicex_uy.png slicex_uy_t$t.png
    mv slicey_uy.png slicey_uy_t$t.png
    for j in `seq 1 9`
    do 
	echo $t.${j}
	srun ./$EXE $t.${j} $LOC
	convert slicey_ux.ppm slicey_ux.png
	convert slicex_uy.ppm slicex_uy.png
	convert slicey_uy.ppm slicey_uy.png
	mv slicey_ux.png slicey_ux_t$t.${j}.png
	mv slicex_uy.png slicex_uy_t$t.${j}.png
	mv slicey_uy.png slicey_uy_t$t.${j}.png
   done
done


#To move the whole directory to /tigress
#cp -r $ScratchDir /tigress/jiarongw/
#rm -rf $ScratchDir
