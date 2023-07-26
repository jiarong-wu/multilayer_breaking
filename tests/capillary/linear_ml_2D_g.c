/**
   # 3D breaking Stokes wave (multilayer solver)

   A steep, 3D, third-order Stokes wave is unstable and breaks. This is
   the 3D equivalent of this [test case](/src/test/stokes.c). The
   bathymetry is given by
   $$
   z_b(y) = - 0.5 + \sin(\pi y)/4
   $$
   i.e. it is shallower toward the back of the domain which causes the
   wave to break earlier there.

   ![Animation of the free-surface. The surface is coloured according to
   the $x$-component of the surface velocity.](breaking/movie.mp4)(
   width=100% )

   The solution is obtained using the layered model and demonstrates its
   robustness and a degree of realism even for this complex case. Note
   also the interesting longitudinal "scars". */

#include "grid/multigrid1D.h"
#include "hydro-tension-g.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "output_mpi.h" // Antoon's function for MPI compatible matrix output

/**
   The initial conditions are given by the wave steepness $ak$ and the
   Reynolds number $Re=c\lambda/\nu$ with $c$ the phase speed of the
   gravity wave and $\lambda$ its wavelength. */

double ak = 0.1;
double a0k0 = 0.1; // The long wave slope
double RE = 40000.;
double gpe_base = 0;
#define g_   9.8
#define k_  (10.*pi) // 5 Waves in the 1m box
#define h_   0.5
#define k0_ (0.25*pi) // The long wave is 8 times larger than the box
#define freq0_ (sqrt(g_*k0_)) // The long wave freq
// #define T0  (k_/sqrt(g_*k_)) // Time is not normalized anymore

/** Function for writing fields at time t.
*/
int writefields_1D (double t, const char *suffix) {
  char filename1[50], filename2[50], filename3[50], filename4[50];
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename2, "field/uz_%s_t%g_l%d", suffix, t, j);
    sprintf (filename3, "field/h_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename4, "field/phi_%s_t%g_l%d", suffix, t, j);
    FILE * fux = fopen (filename1, "w");
    FILE * fuz = fopen (filename2, "w");
    FILE * fh = fopen (filename3, "w");
    FILE * fphi = fopen (filename4, "w");
    foreach() {
      fprintf (fux, "%g %g \n", x, u.x[0,j]);
      fprintf (fuz, "%g %g \n", x, w[0,j]);
      fprintf (fh, "%g %g \n", x, h[0,j]);
      fprintf (fphi, "%g %g \n", x, phi[0,j]);
    }
    fflush (fux);
    fclose (fux);
    fflush (fuz);
    fclose (fuz);
    fflush (fh);
    fclose (fh);
    fflush (fphi);
    fclose (fphi);
  }
  return 0;
}

/**
   The domain is periodic in $x$ and resolved using 256$^2$
   points and 60 layers. */

int main(int argc, char * argv[])
{
  if (argc > 1)
    RE = atof (argv[1]);
  if (argc > 2)
    nl = atoi (argv[2]);
  if (argc > 3)
    N = atoi (argv[3]);
  if (argc > 4)
    ak = atof (argv[4]);
  if (argc > 5)
    theta_H = atof(argv[5]);
  /* max_slope = 0.4; */ // Change the slope limit
  origin (-L0/2.);
  periodic (right);
  TOLERANCE = 1e-4;
  G = g_;
  gpe_base = -0.5*sq(h_)*sq(L0)*g_;
  nu = 0.000001; // Overwrite the Reynolds number and use realistic value
  CFL = 0.1;
  CFL_H = 10;
  max_slope = 0.577; // Change the slope limiting value 
  run();
}

/**
   The initial conditions for the free-surface and velocity are given by
   a linear sinusoidal wave. */

double wave (double x, double y)
{
  double a_ = ak/k_;
  double eta = a_*cos(k_*x);
  return eta - y;
}

double u_x (double x, double y, double sigma) 
{ // k_ is a global variable
  return ak*sqrt(g_/k_+sigma*k_)*cos(k_*x)*exp(k_*y); 
}

double u_y (double x, double y, double sigma)
{
  return ak*sqrt(g_/k_+sigma*k_)*sin(k_*x)*exp(k_*y);
}

event init (i = 0)
{
  foreach() 
    sigma[] = 0.0000728;
  if (!restore("restart")) {
    geometric_beta (1./3., true);
    foreach() {
      /* zb[] = - 0.5 + sin(pi*y)/4.; */
      zb[] = -h_;
      eta[] = wave(x,y);
      double H = wave(x, 0) - zb[];
      foreach_layer() {
	h[] = H/nl;
      }
    }
    vertical_remapping (h,tracers);
    foreach() {
      double z = zb[];
      foreach_layer() {
	z += h[]/2.;
	u.x[] = u_x(x, z, sigma[]);
	w[] = u_y(x, z, sigma[]);
	z += h[]/2.;
      }
      foreach_layer() 
	G_var[] = G*(1. + a0k0*cos(k0_*x));
    }
  }
  else {
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method
    dtmax = 0.01;
    dt = dtnext (dtmax);
    char *suffix = "matrix";
    writefields_1D (t, suffix);
  }
}

// Modulation of long waves by varying gravity??
event var_g (i++) {
  foreach(){
    foreach_layer() 
      G_var[] = G*(1. + a0k0*cos(k0_*x-freq0_*t));
  }
}

// Horizontal diffusion was not added before!
event viscous_term (i++)
  horizontal_diffusion ((scalar *){u}, nu, dt);

/* event limiter (i = 0) { */
/*   gradient = minmod2; // test if this is able to remove the oscillation */
/*   theta = 1; // defined in util.h (theta=1 is minmod, the most dissipative, and theta=2 superbee, the least dissipative  */
/* } */


/**
   We log the evolution of the kinetic and potential energies.

   ~~~gnuplot Evolution of the kinetic, potential and total energy
   set xlabel 't/T0'
   plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:($3+0.0007) w l t 'potential', \
   '' u 1:(($2+$3+0.0007)/2.) w l t 'total/2'
   ~~~
*/

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    double zc = zb[];
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  static FILE * fp = fopen ("budget.dat", "w");
  fprintf (fp, "%g %g %g\n", t, ke/2., g_*gpe - gpe_base);
  fflush (fp);
}

/**
   And generate the movie of the free surface (this is quite expensive). */

event movieoutput (t += 0.1)
{
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uz_matrix_%g", t);
  FILE * feta = fopen (filename1, "w");
  FILE * fux = fopen (filename2, "w");
  FILE * fuz = fopen (filename3, "w");
  // In 1D directly print to file (need to write a mpi compatible version)
  foreach() {
    fprintf (feta, "%g %g \n", x, eta[]);
    fprintf (fux, "%g %g \n", x, u.x[0,nl-1]);
    fprintf (fuz, "%g %g \n", x, w[0,nl-1]);
  }
  fflush (feta);
  fclose (feta);
  fflush (fux);
  fclose (fux);
  fflush (fuz);
  fclose (fuz);
}

event field_log (t += 1) {
  char *suffix = "matrix";
  writefields_1D (t, suffix);
}

event snapshot (t += 5) {
  char dname[100];
  sprintf (dname, "dump%g", t);
  dump (dname);
}

event end (t = 100.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}


/**
   ## Parallel run

   The simulation was run in parallel on the [Occigen
   machine](https://www.cines.fr/calcul/materiels/occigen/) on 64 cores,
   using this script

   ~~~bash
   local% qcc -source -D_MPI=1 breaking.c
   local% scp _breaking.c user@occigen.cines.fr:
   ~~~

   ~~~bash
   #!/bin/bash
   #SBATCH -J breaking
   #SBATCH --constraint=HSW24
   #SBATCH --ntasks=64
   #SBATCH --threads-per-core=1
   #SBATCH --time=1:00:00
   #SBATCH --exclusive
   #SBATCH --mail-type=BEGIN,END,FAIL
   #SBATCH --mail-user=popinet@basilisk.fr

   module purge
   module load openmpi
   module load intel gcc

   NAME=breaking
   # mpicc -Wall -std=c99 -O2 _$NAME.c -o $NAME \
   #    -I/home/popinet/local -L/home/popinet/local/gl -L/home/popinet/local/lib \
   #    -lglutils -lfb_osmesa -lOSMesa -lGLU -lppr -lgfortran -lm

   export LD_LIBRARY_PATH=/home/popinet/local/lib:$LD_LIBRARY_PATH
   export PATH=$PATH:/home/popinet/local/bin
   rm -f *.ppm
   srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS $NAME 2> log > out
   ~~~

   The number of timesteps was 4450 and the runtime was 37 minutes with
   movie generation and 17 minutes without, corresponding to a
   computational speed of 277 000 point.timestep/sec/core (on 64 cores). */
