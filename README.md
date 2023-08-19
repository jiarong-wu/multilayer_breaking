The source code and minimal working environment for the multilayer breaking wave project.

## Installing basilisk, compiling and running on tiger
Install basilisk using the tarball provided here to ensure the correct environment.
Install additional library [ppr](http://basilisk.fr/src/ppr/).

To compile, run e.g.
```make -f Makefile.parallel OUTPUT=field_PM```
which will generate an executable under the current fold.

The `slurmfiles/` folder contains examples of slurm scripts (I was using the Tiger cluster).

## Preparation of the spectrum
The `specgen/` folder contains the python code that generates the spectra in a format that `spectrum.h` can read. The `specgen/spectra_examples/` folder contains some examples of such spectra but they are written in binary so can't be directly viewed.

## Working examples
The `examples/` folder contains the working examples.

### Header file 
* `output_mpi.h` is necessary for outputting 2D slices.
* `spectrum.h` enables the spectrum read-in process. Currently, it uses fixed number of wavenumbers (32\*33) with linear spacing. MPI compatible. The spectrum is in binary and generated with separate python code.
* `spectrum_swell.h` temporarily expands the spectrum to 65\*65 shape with still linear spacing. 
* `vorticity.h` computes and outputs vorticity fields and gradient fields.

### Field scale cases
* `field_PM.c` is the code used for the JFM paper, which uses synthesized wind-sea spectra. 
* `field_swell.c` temporarily expands this to spectrum with swells (by including spectrum_swell.h). Need to better factorize the code later.

### Single breaking wave simulations
* `stokes/stokes_ns.c` and `stokes/slurm_stokes_ns.sh` for the two-phase Navier-Stokes simulation.
* `stokes/stokes_ml.c` and `stokes/slurm_stokes_ml.sh` for the multilayer simulation.
Both are compiled using the Makefile parallel file.
(Time in the files is normalized while time in the energy files is not.)

## Tests
Folder `test/` contains files that are still under testing
* `test/capillary/` contains the extension to include surface tension with the [headerfiles from Clement Robert](http://basilisk.fr/sandbox/crobert/). It is working now and has been tested with a broadband spectrum as well.
* `test/smax/` contains tests on the effect of gradient-limiter threshold. At some point, I also tried to print out the order of event execution and at what point the gradient limiter is evoked. So there is a `hydro_test.h` and `nh_test.h` header file there. These checks can be activated by adding `-DEVENT=1` and `-DSMAX=1` respectively.
* `test/horidiff/` and `test/vertdiff/` contain tests on the effect of horizontal and vertical diffusion coefficients that is still ongoing.
* `test/vort/` contains vorticity computation files. The `vorticity.c` code has been re-written into `vorticity.h` so it's obsolete. `vorticity_ns.c` attempts to compute the vorticity field for the NS stokes case but was not working.
* `test/disp/` contains some tests on the basic dispersion relation that I did not have time to pursue.
* `test/current/` contains the first attempt to put current in the simulation.

## To-dos
* Replace the output_mpi.h since the field is already structured and does not require further interpolation. We should also use less precision (instead of floating number precision now) to save disk space.
* Find a better way to read in the spectrum and initialize instead of using the fixed length wavenumber array now.

