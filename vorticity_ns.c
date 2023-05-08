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

/**
The primary parameters are the wave steepness $ak$, the Bond and
Reynolds numbers. */

double ak = 0.35;
double BO = 200.;
double RE = 40000.;

/**
   The default maximum level of refinement depends on the dimension. */

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

double TRESTORE; // restoring time

/**
   Re-write the vorticity function to compute all 3D components. */

vector omega[];
void vorticity3D (const vector u, vector omega)
{
  foreach() { 
      /* omega[] = ((fm.x[1] - fm.x[])*u.y[] + */
      /* 	       fm.x[1]*u.y[1] - fm.x[]*u.y[-1] - */
      /* 	       (fm.y[0,1] - fm.y[])*u.x[] + */
      /* 	       fm.y[]*u.x[0,-1] - fm.y[0,1]*u.x[0,1])/(2.*cm[]*Delta + SEPS); */
    omega.x[] = (0.5*(u.y[]+u.y[0,0,-1]) - 0.5*(u.y[]+u.y[0,0,1]) - 0.5*(u.z[]+u.z[0,-1]) + 0.5*(u.z[]+u.z[0,1]))/Delta;
    omega.y[] = (-0.5*(u.x[]+u.x[0,0,-1]) + 0.5*(u.x[]+u.x[0,0,1]) + 0.5*(u.z[]+u.z[-1]) - 0.5*(u.z[]+u.z[1]))/Delta;
    omega.z[] = (-0.5*(u.y[]+u.y[-1]) + 0.5*(u.y[]+u.y[1]) + 0.5*(u.x[]+u.x[0,-1]) - 0.5*(u.x[]+u.x[0,1]))/Delta;
    fprintf (stderr, "omega computed!\n");
  }
}

/**
We are interested in the viscous dissipation rate in both water and air. 
Here we modify the original function to output a field not just the integrated value. */

scalar rateWater[];
int dissipation_rate (scalar rateWater)
{
  foreach () {
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
    rateWater[] = mu1/rho[]*f[]*sqterm; //water
  }  
  fprintf (stderr, "Dissipation computed!\n");
  return 0;
}

/**
The program takes optional arguments which are the level of
refinement, steepness, Bond and Reynolds numbers. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    RE = atof(argv[3]);    
  if (argc > 4)
    TRESTORE = atof(argv[4]); // Restoring time for analysis, already normalized by T0

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
  mu1 = 1.0/RE; // using wavelength as length scale
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
We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */

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

event init (i = 0)
{
  char dumpname[100];
  sprintf (dumpname, "dump_t%g", TRESTORE);
  if (!restore ("restart")) {
    fprintf (stderr, "%s not found!\n", dumpname);
  }
  else {
    fprintf (stderr, "Computing vorticity! \n");
    /* vorticity3D (u, omega); */
    /* fprintf (stderr, "Computing dissipation! \n"); */
    /* dissipation_rate (rateWater); */
    /* char filename[100]; */
    /* int Nslice = 256; */
    /* double zslice = -L0/2+L0/2./Nslice; */
    /* for (int i=0; i<Nslice; i++) { */
    /*   zslice += L0/Nslice; */
    /*   sprintf (filename, "./field/omegax_t%g_slice%d", TRESTORE, i); */
    /*   sliceXY (filename,omega.x,zslice,9); */
    /*   sprintf (filename, "./field/omegay_t%g_slice%d", TRESTORE, i); */
    /*   sliceXY (filename,omega.y,zslice,9); */
    /*   sprintf (filename, "./field/omegaz_t%g_slice%d", TRESTORE, i); */
    /*   sliceXY (filename,omega.z,zslice,9); */
    /*   sprintf (filename, "./field/epsilon_t%g_slice%d", TRESTORE, i); */
    /*   sliceXY (filename,rateWater,zslice,9); */
   /* } */
 }
}



