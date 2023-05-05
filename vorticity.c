/**
# Field scale wave breaking (multilayer solver)
*/

#include "grid/multigrid.h"
#include "view.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "input.h" // Used for older input method (matrices of quantities)
#include "output_mpi.h" // Antoon's function for MPI compatible matrix output

/**
Definition of some controlling parameters. */

double TRESTORE = 50.; // t to restore
int NLAYER = 10; // number of layers
int LEVEL_data = 7; // horizontal resolution

int main (int argc, char * argv[])
{
  if (argc > 1)
    NLAYER = atoi(argv[1]); // # of layers
  if (argc > 2)
    LEVEL_data = atoi(argv[2]); // Horizontal resolution
  if (argc > 3)
    TRESTORE = atof(argv[3]); // Restoring time for analysis
  if (argc > 4)
    L0 = atof(argv[4]); // Box size 
  else
    L0 = 50.;
  origin (-L0/2., -L0/2.);
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data; 
  nl = NLAYER;
  run();
}

/**
## Computation of vorticity 
The following function computes vorticity vector omega. */

face vector hu, hf, ha;
vector omega;
scalar omegaz;
vector dzdx;

void vort ()
{
  hu = new face vector[nl];
  hf = new face vector[nl];
  ha = new face vector[nl];
  foreach_face (reduction (min:dtmax)) {
    double ax = a_baro (eta, 0);
    double H = 0.;
    foreach_layer() {

      /**
      The face velocity is computed as the height-weighted average of
      the cell velocities. */
      
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      hu.x[] = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;

      /**
      Different from the original face_field event in hydro.h, not going to bother with 
      advection computation. */
      
      double hff;
      hff = h[];
      // Should we do the following instead?
      /* double hff, un = pdt*(hu.x[] + pdt*ax)/Delta, a = sign(un); */
      /* int i = - (a + 1.)/2.; */
      /* double g = h.gradient ? h.gradient (h[i-1], h[i], h[i+1])/Delta : */
      /* 	(h[i+1] - h[i-1])/(2.*Delta); */
      /* hff = h[i] + a*(1. - a*un)*g*Delta/2.; */
      hf.x[] = fm.x[]*hff; // Why does hf need to be multiplied by fm?

      /**
      The flux and height-weighted accelerations are computed. */
            
      hu.x[] *= hf.x[];
      ha.x[] = hf.x[]*ax;
 
      H += hff;
    }
  }

  omega = new vector[nl];
  reset ({omega}, 0.);
  omegaz = new scalar[nl];
  reset ({omegaz}, 0.);
  foreach () { 
    foreach_layer () {
      omegaz[] = (0.5*(u.y[] + u.y[1])*fm.x[1] - 0.5*(u.y[] + u.y[-1])*fm.x[] \
		  + 0.5*(u.x[] + u.x[0,-1])*fm.y[] - 0.5*(u.x[] + u.x[0,1])*fm.y[0,1])/Delta;
      foreach_dimension () {
	if (point.l > 0) {
	  double area, circ;
	  area = Delta*(hf.y[] + hf.y[0,0,-1] + hf.y[0,1,0] + hf.y[0,1,-1])/4.; // Add fm later
	  circ = (-u.y[] + u.y[0,0,-1])*Delta - 0.5*(w[0,0,-1] + w[0,-1,-1])*0.5*(hf.y[] + hf.y[0,0,-1]) + \
	    0.5*(w[0,0,-1] + w[0,1,-1])*0.5*(hf.y[0,1,0] + hf.y[0,1,-1]);
	  omega.x[] = circ/area;      
	}
	else
	  omega.x[] = 0.; // not well defined for the bottom layer
      }
    }
  }
  
  // Analyze the slope
  dzdx = new vector[nl];
  reset ({dzdx}, 0.);
  foreach () {
    coord dz;
    foreach_dimension ()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer () {
      foreach_dimension () {
	dzdx.x[] = dz.x + hf.x[1] - hf.x[];
	dz.x += hf.x[1] - hf.x[];
      }
    }
  }
}

vector dzdxc;
void slope () {
  // Analyze the slope, test if it's different when taking a centered value
  dzdxc = new vector[nl];
  reset ({dzdxc}, 0.);
  foreach () {
    coord dz;
    foreach_dimension ()
      dz.x = (fm.x[1]*(zb[1] + zb[]) - fm.x[]*(zb[-1] + zb[]))/2.;
    foreach_layer () {
      foreach_dimension () {
	dzdxc.x[] = dz.x + h[1] - h[];
	dz.x += h[1] - h[];
      }
    }
  }
}
/**
## Write to files
A new writefields function with omega added. */

int writefields (double t, const char *suffix) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50], filename5[50], filename6[50], filename7[50], filename8[50], filename9[50], filename10[50], filename11[50];
  vector u_temp, omega_temp, dzdx_temp, dzdxc_temp;
  scalar w_temp, h_temp, omegaz_temp;
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%g_l%d", suffix, t, j);
    sprintf (filename2, "field/uy_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename3, "field/uz_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename4, "field/h_%s_t%g_l%d", suffix, t, j);
    sprintf (filename5, "field/omegax_%s_t%g_l%d", suffix, t, j);
    sprintf (filename6, "field/omegay_%s_t%g_l%d", suffix, t, j);
    sprintf (filename7, "field/omegaz_%s_t%g_l%d", suffix, t, j);
    sprintf (filename8, "field/dzdx_%s_t%g_l%d", suffix, t, j);
    sprintf (filename9, "field/dzdy_%s_t%g_l%d", suffix, t, j);
    sprintf (filename10, "field/dzdxc_%s_t%g_l%d", suffix, t, j);
    sprintf (filename11, "field/dzdyc_%s_t%g_l%d", suffix, t, j);
    if (j==0) {
      // The first layer is named u instead of u0
      sprintf (s, "u");
      u_temp = lookup_vector (s);
      sprintf (s, "w");
      w_temp = lookup_field (s);
      sprintf (s, "h");
      h_temp = lookup_field (s);
      sprintf (s, "omega");
      omega_temp = lookup_vector (s);
      sprintf (s, "omegaz");
      omegaz_temp = lookup_field (s);
      sprintf (s, "dzdx");
      dzdx_temp = lookup_vector (s);
      sprintf (s, "dzdxc");
      dzdxc_temp = lookup_vector (s);
    }
    else {
      sprintf (s, "u%d", j);
      u_temp = lookup_vector (s);
      sprintf (s, "w%d", j);
      w_temp = lookup_field (s);
      sprintf (s, "h%d", j);
      h_temp = lookup_field (s);
      sprintf (s, "omega%d", j);
      omega_temp = lookup_vector (s);
      sprintf (s, "omegaz%d", j);
      omegaz_temp = lookup_field (s);
      sprintf (s, "dzdx%d", j);
      dzdx_temp = lookup_vector (s);
      sprintf (s, "dzdxc%d", j);
      dzdxc_temp = lookup_vector (s);
    }
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi (u_temp.x, fux, N, linear = true);
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi (u_temp.y, fuy, N, linear = true);
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    output_matrix_mpi (w_temp, fuz, N, linear = true);
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    output_matrix_mpi (h_temp, fh, N, linear = true);
    fclose (fh);    
    FILE * fomegax = fopen (filename5, "w");
    output_matrix_mpi (omega_temp.x, fomegax, N, linear = true);
    fclose (fomegax);
    FILE * fomegay = fopen (filename6, "w");
    output_matrix_mpi (omega_temp.y, fomegay, N, linear = true);
    fclose (fomegay);
    FILE * fomegaz = fopen (filename7, "w");
    output_matrix_mpi (omegaz_temp, fomegaz, N, linear = true);
    fclose (fomegaz);
    FILE * fdzdx = fopen (filename8, "w");
    output_matrix_mpi (dzdx_temp.x, fdzdx, N, linear = true);
    fclose (fdzdx);
    FILE * fdzdy = fopen (filename9, "w");
    output_matrix_mpi (dzdx_temp.y, fdzdy, N, linear = true);
    fclose (fdzdy);
    FILE * fdzdxc = fopen (filename10, "w");
    output_matrix_mpi (dzdxc_temp.x, fdzdxc, N, linear = true);
    fclose (fdzdxc);
    FILE * fdzdyc = fopen (filename11, "w");
    output_matrix_mpi (dzdxc_temp.y, fdzdyc, N, linear = true);
    fclose (fdzdyc);
  }
  return 0;
}

/**
Read the dump file and compute vorticity and output. */
event init (i = 0)
{
  char dumpname[100];
  sprintf (dumpname, "dump_t%g", TRESTORE);
  if (!restore (dumpname)) {
    fprintf (stderr, "%s not found!\n", dumpname);
  }
  else {
    // We limit the first time step after the restart
    geometric_beta (1./3., true); // when restarting, remember to specify the grid mapping method, and this needs to match the original grid
    /* dtmax = 0.01; */
    /* dt = dtnext (dtmax); */
    /* char *suffix = "matrix"; */
    /* writefields (t, suffix); */
    vort ();
    slope ();
    char *suffix = "matrix";
    writefields (TRESTORE, suffix); // if I put t instead of TRESTORE here it doesn't work?? t=0. When does t get updated?
  }
}

event deletefield (t=end, last) {
  delete ((scalar *){hu,ha,hf,omega});
}
