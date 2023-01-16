The source code and minimal working environment for the multilayer breaking wave project.

Install basilisk using the tarball provided here to ensure the correct environment.
Install additional library ppr.

To compile, run 
```make -f Makefile.parallel OUTPUT=field_PM```
which will generate an executable under the current fold.

`slurm_example.sh` is an example slurm script.
Folder ./test/ contains files that are still under testing

output_mpi.h is necessary for outputting 2D slices.

## Single breaking wave simulations
* stokes_ns.c and slurm_stokes_ns.sh for the two-phase Navier-Stokes simulation.
* stokes_ml.c and slurm_stokes_ml.sh for the multilayer simulation.
Both are compiled using the Makefile parallel file.
Time in the files are normalized while time in the energy files are not.