/**
# Breaking wave

We solve the two-phase Navier--Stokes equations with surface tension
and using a momentum-conserving transport of each phase. Gravity is
taken into account using the "reduced gravity approach" and the
results are visualised using Basilisk view. */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "tag.h"


/**
   We log some profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

/**
The primary parameters are the wave steepness $ak$, the Bond and
Reynolds numbers. */

double ak = 0.35;
double BO = 200.;
double RE = 40000.;

/**
   The default maximum level of refinement depends on the dimension. */

//int LEVEL = dimension == 2 ? 9 : 6;
int LEVEL = 9;
/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
   The density and viscosity ratios are those of air and water. */

#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
   Define if we want to use a Dirac viscous layer initialization. */
int DIRAC = 0;

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. */

#define k_  (2.*pi)
#define h_   0.5
#define g_   1.
#define T0  (k_/sqrt(g_*k_))

/**
The program takes optional arguments which are the level of
refinement, steepness, Bond and Reynolds numbers. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);    
  if (argc > 5)
    DIRAC = atof(argv[5]);

  /**
  The domain is a cubic box centered on the origin and of length
  $L0=1$, periodic in the x- and z-directions. */
   
  origin (-L0/2, -L0/2, -L0/2);
  periodic (right);
#if dimension > 2
  periodic (front);
#endif

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 1.;
  rho2 = RATIO;
  mu1 = 1.0/RE; //using wavelength as length scale
  mu2 = 1.0/RE*MURATIO;
  f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
#if TREE  
  N = 32;
#else
  N = 1 << LEVEL;
#endif
  run();
}

/**
## Initial conditions

These functions return the shape of a third-order Stokes wave with the
wavenumber and steepness given by the parameters above ($ak$ and
$_k_$). */

double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}
/**

We also calculate an approximation to a Dirac distribution on the wave surface.
This allows us to calculate a vortex sheet on the surface to provide a boundary
layer in the air above the water surface. */

double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}


/**
We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */
event init (i = 0)
{
  if (!restore ("restart")) {
    do {
      fraction (f, wave(x,y));

      /**
	 To initialise the velocity field, we first define the potential. */
      
      scalar Phi[];
      foreach() {
	double alpa = 1./tanh(k_*h_);
	double a_ = ak/k_;
	double sgma = sqrt(g_*k_*tanh(k_*h_)*
			   (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					      (sq(alpa) - 1.) + sq(alpa))));
	double A_ = a_*g_/sgma;
	double phi1 = A_*cosh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x);
	double phi2 = 3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
	  cosh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_);
	double phi3 = 1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
	  (9.*sq(alpa) - 13.)*
	  cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
	Phi[] = phi1 + ak*phi2 + ak*ak*phi3;
      } 
      boundary ({Phi});
      if (DIRAC){
	/** 
	      We calculate the vorticity in the Dirac layer. We need a separate
	      foreach here because we need the derivative of the potential phi.*/
	scalar vort2[];
	scalar psi[];
	foreach() {
	  vort2[] = -2.0*gaus(y,wave(x,y)+y,f[])*(Phi[1,0]-Phi[-1,0])/(2.*Delta);
	  psi[] = 0.0;
	}
	boundary({vort2,psi});
	psi[top] = dirichlet(0.);
	psi[bottom] = dirichlet(0.);
	/**
	   Solve the Poisson problem for the streamfunction psi given the vorticity field.*/
	poisson(psi, vort2);
	/**
	      And then define the velocity field using centered-differencing
	      of the streamfunction. */
      
	foreach()
	  {
	    u.x[] = (psi[0,1] - psi[0,-1])/(2.*Delta);
	    u.y[] = -(psi[1] - psi[-1])/(2.*Delta);
	  }
      }
      else{
	/**
	      If we choose not to use the Dirac layer, instead initialize
	      in the water only according to the potential already calculated.*/
	foreach(){
	  foreach_dimension()
	    u.x[] = (Phi[1] - Phi[-1])/(2.0*Delta) * f[];
	}
      }
      boundary ((scalar *){u});
      
    }

    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f,u},
			  (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5).nf);
#else
    while (0);
#endif
  }
}

/**
## Outputs

We are interested in the viscous dissipation rate. */

/**
## Outputs

We are interested in the viscous dissipation rate in both water and air. */

int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

/**
   We also want to count the drops and bubbles in the flow. */

/* event countDropsBubble(i++) */
/* { */
/*   scalar m1[]; //droplets */
/*   scalar m2[]; //bubbles */
/*   foreach(){ */
/*     m1[] = f[] > 0.5; //i.e. set m true if f[] is close to unity (droplets) */
/*     m2[] = f[] < 0.5; //m true if f[] close to zero (bubbles) */
/*   } */
/*   int n1 = tag(m1); */
/*   int n2 = tag(m2); */
/*   /\** */
/*      Having counted the bubbles, now we find their size. This example */
/*      is similar to the jet atomization problem. We are interested in */
/*      the volumes and positions of each droplet/bubble.*\/ */
/*   double v1[n1]; //droplet */
/*   coord b1[n1];  //droplet */
/*   double v2[n2]; //bubble */
/*   coord b2[n2];  //bubble */
/*   /\** */
/*      We initialize: *\/ */
/*   for (int j=0; j<n1; j++) */
/*     { */
/*       v1[j] = b1[j].x = b1[j].y = b1[j].z = 0.0; */
/*     } */
/*   for (int j=0; j<n2; j++) */
/*     { */
/*       v2[j] = b2[j].x = b2[j].y = b2[j].z = 0.0; */
/*     } */
/*   /\** */
/*      We proceed with calculation. *\/ */
/*   foreach_leaf() //droplets */
/*     { */
/*       if (m1[] > 0) { */
/* 	int j = m1[] - 1; */
/* 	v1[j] += dv()*f[]; //increment the volume of the droplet */
/* 	coord p = {x,y,z}; */
/* 	foreach_dimension() */
/* 	  b1[j].x += dv()*f[]*p.x; */
/*       } */
/*     } */
/*   foreach_leaf() //bubbles */
/*     { */
/*       if (m2[] > 0) { */
/* 	int j = m2[] - 1; */
/* 	v2[j] += dv()*(1.0-f[]); */
/* 	coord p = {x,y,z}; */
/* 	foreach_dimension() */
/* 	  b2[j].x += dv()*(1.0-f[])*p.x; */
/*       } */
/*     } */
/*   /\** */
/*      Reduce for MPI. *\/ */
/* #if _MPI */
/*   MPI_Allreduce (MPI_IN_PLACE, v1, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/*   MPI_Allreduce (MPI_IN_PLACE, b1, 3*n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/*   MPI_Allreduce (MPI_IN_PLACE, v2, n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/*   MPI_Allreduce (MPI_IN_PLACE, b2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); */
/* #endif */
/*   /\** */
/*      Output the volume and position of each droplet to file. *\/ */
/*   static FILE * fdrop = fopen("droplets.dat","w"); */
/*   static FILE * fbubb = fopen("bubbles.dat","w"); */
/*   for (int j=0; j<n1; j++) */
/*     { */
/*       fprintf (fdrop, "%d %g %d %g %g %g\n", i, t, */
/* 	       j, v1[j], b1[j].x/v1[j], b1[j].y/v1[j]); */
/*     } */
/*   for (int j=0; j<n2; j++) */
/*     { */
/*       fprintf (fbubb, "%d %g %d %g %g %g\n", i, t, */
/* 	       j, v2[j], b2[j].x/v2[j], b2[j].y/v2[j]); */
/*     } */
/* } */


/**
We log the evolution of the kinetic and potential energies and
dissipation rate as functions of the non-dimensional time. */

event graphs (i++) {
  static FILE * fpwater = fopen("budgetWater.dat", "a");
  static FILE * fpair = fopen("budgetAir.dat", "a");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
  }
  fprintf (fpwater, "%g %g %g %g\n",
	   t, ke/2., gpe + 0.125, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t, keAir/2., gpeAir + 0.125, dissAir);
  fprintf (ferr, "%g %g %g %g\n",
	   t, ke/2., gpe + 0.125, dissWater);
}

/** 
Output velocity field */
/** 
    Outputting slices on the fly. */
void sliceXY(char * fname,scalar s,double zp, int maxlevel){
  FILE *fpver = fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}


/**
   Output eta on the fly. */
void output_twophase_locate (double snapshot_time) {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,nx,ny,nz\n");
  foreach(){
    if (interfacial (point, f)){
      if (point.level == LEVEL) {
	coord n = mycs (point, f);
	double eta = pos[];
	if (point.level > 0) {
	  POINT_VARIABLES;
	  fprintf (feta, "%g,%g,%g,%g,%g,%g\n", x, z, eta, n.x, n.y, n.z);
	}
      }
    }
  }
  fflush (feta);
  fclose (feta);
}


event eta_output (t += 0.1*T0) {
    output_twophase_locate (t/T0);
}


/* event output_slice (t += 0.05)  */
/* { */
/*   char filename[100]; */
/*   double zslice = 0.; */
/*   sprintf (filename, "./field/ux_t%g_center", t); */
/*   sliceXY (filename,u.x,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/uy_t%g_center", t); */
/*   sliceXY (filename,u.y,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/uz_t%g_center", t); */
/*   sliceXY (filename,u.z,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/f_t%g_center", t); */
/*   sliceXY (filename,f,zslice,MAXLEVEL-1); */
/* } */

scalar pair[];
event turbulence_stat (t += 0.1*T0) {
  char filename[100];
  int Nslice = 256;
  double L0 = 2*pi;
  double zslice = -L0/2+L0/2./Nslice;
  for (int i=0; i<Nslice; i++) {
    zslice += L0/Nslice;
    sprintf (filename, "./field/ux_t%g_slice%d", t/T0, i);
    sliceXY (filename,u.x,zslice,9);
    sprintf (filename, "./field/uy_t%g_slice%d", t/T0, i);
    sliceXY (filename,u.y,zslice,9);
    sprintf (filename, "./field/uz_t%g_slice%d", t/T0, i);
    sliceXY (filename,u.z,zslice,9);
    sprintf (filename, "./field/f_t%g_slice%d", t/T0, i);
    sliceXY (filename,f,zslice,9);
  }
}

/**
~~gnuplot Evolution of kinetic, potential and total energies and dissipation rate.
set key auto col
plot 'log' u 1:2 w l, '' u 1:3 w l, '' u 1:(($2+$3)/2.) t 'total/2' w l, \
'' u 1:4 w l
~~

## Visualisation

We use Basilisk view (and output_ppm()) to display animations of the
results.

On some parallel systems, pipes tend to cause problems, so we switch
to simple uncompressed PPM outputs when running with MPI. Otherwise,
we use MPEG-4 file compression (which requires working pipes and
ffmpeg). */
#if _MPI
#  define POPEN(name, mode) fopen (name ".ppm", mode)
#else
#  define POPEN(name, mode) popen ("ppm2mp4 > " name ".mp4", mode)
#endif

event movies (t += 0.01*T0) {

  /**
  We first do simple movies of the volume fraction, level of
  refinement fields. In 3D, these are in a $z=0$ cross-section. */

  {
    static FILE * fp = POPEN ("f", "w");
    output_ppm (f, fp, min = 0, max = 1, n = 512);
  }

#if TREE
  {
    scalar l[];
    foreach()
      l[] = level;
    static FILE * fp = POPEN ("level", "w");
    output_ppm (l, fp, min = 5, max = LEVEL, n = 512);
  }
#endif

  /**
  <p><center>
  <video width="512" height="512" controls>
  <source src="wave/level.mp4" type="video/mp4">
  Your browser does not support the video tag.
  </video><br>
  Wave breaking. Animation of the level of refinement.
  </center></p>

  We use Basilisk view differently in 2D and 3D. */
  
  scalar omega[];
  vorticity (u, omega);
#if dimension == 2
  view (width = 800, height = 600, fov = 18.8);
  clear();

  /**
  We repeat the drawing periodically in the x-direction. */
  
  for (double x = -L0; x <= L0; x += L0)
    translate (x) {
      draw_vof ("f");
      squares ("omega", linear = true);
    }

  /**
  This gives the following movie.
  <p><center>
  <video width="800" height="600" controls>
  <source src="wave/movie.mp4" type="video/mp4">
  Your browser does not support the video tag.
  </video></center></p>
  */
#else // dimension == 3
  /**
  In 3D, we generate a first movie seen from below. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = -pi/6, fov = 20);
  clear();
  char s[80];
  sprintf (s, "t = %.2f T0", t/T0);
  draw_string (s, size = 80);
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
/*       for (double z = -3*L0; z <= L0; z += L0) */
/* translate (z = z) */
  draw_vof ("f");
    }
  {
    static FILE * fp = POPEN ("below", "w");
    save (fp = fp);
  }

  /**
  And a second movie, seen from above. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = pi/6, fov = 20);
  clear();

  /**
  In 3D, we are doubly-periodic (along x and z). */
  
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
      draw_vof ("f");
    }
#endif // dimension == 3
  {
    static FILE * fp = POPEN ("movie", "w");
    save (fp = fp);
  }
}

/**
## Dump/restore

To be able to restart, we dump the entire simulation at regular
intervals. */

event snapshot_dt (t += 0.1*T0) {
  char dumpname[100];
  sprintf(dumpname, "dump%g", t/T0);
  dump ( dumpname );
}

/**
## End 

The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
(alternatively 4) periods. */

event end (t = 6.*T0) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}
#endif

/**
## Running in parallel

This file will work in 2D or 3D, either with parallel multigrid
(without adaptivity), using for example:

~~bash
qcc -source -D_MPI=1 -grid=multigrid3D wave.c
scp _wave.c occigen.cines.fr:
~~

and then following a recipe similar to that of the
[atomisation](/src/examples/atomisation.c#on-occigen) example.

To use adaptivity, just do something like:

~~bash
qcc -source -D_MPI=1 -grid=octree wave.c
scp _wave.c occigen.cines.fr:
~~
*/
