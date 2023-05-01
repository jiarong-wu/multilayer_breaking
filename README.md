The source code and minimal working environment for the multilayer breaking wave project.

## Installing basilisk, compiling and running on tiger
Install basilisk using the tarball provided here to ensure the correct environment.
Install additional library ppr.

To compile, run 
```make -f Makefile.parallel OUTPUT=field_PM```
which will generate an executable under the current fold.

`slurm_example.sh` is an example slurm script.

## Header file 
output_mpi.h is necessary for outputting 2D slices.
`spectrum.h` enables the spectrum read-in process. Currently it is fixed number of wavenumbers (32\*33) with linear spacing. MPI compatible. The spectrum is in binary and generated with separate python code.
`spectrum_swell.h` temporarily expands the spectrum to 65\*65 shape with still linear spacing. 

## Field scale cases
field_PM.c is the code used for the JFM paper, which uses synthesized wind-sea spectrum. field_swell.c temporarily expands this to spectrum with swells (by including spectrum_swell.h). Need to better factorize the code later.

## Single breaking wave simulations
* stokes_ns.c and slurm_stokes_ns.sh for the two-phase Navier-Stokes simulation.
* stokes_ml.c and slurm_stokes_ml.sh for the multilayer simulation.
Both are compiled using the Makefile parallel file.
Time in the files are normalized while time in the energy files are not.

## Others 
Folder ./test/ contains files that are still under testing

