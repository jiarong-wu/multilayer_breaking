/**
# Header file for vorticity computation.
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
## Computation of vorticity 
The following function computes vorticity vector omega. */

face vector hu, hf, ha;
vector Omega; // Use capital Omega to avoid conflict with the wave frequency omega
scalar Omegaz;
vector dzdx; // dzdx from face field. declared and then assigned space in the function

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

  Omega = new vector[nl];
  reset ({Omega}, 0.);
  Omegaz = new scalar[nl];
  reset ({Omegaz}, 0.);
  foreach () { 
    foreach_layer () {
      Omegaz[] = (0.5*(u.y[] + u.y[1])*fm.x[1] - 0.5*(u.y[] + u.y[-1])*fm.x[] \
		  + 0.5*(u.x[] + u.x[0,-1])*fm.y[] - 0.5*(u.x[] + u.x[0,1])*fm.y[0,1])/Delta;
      foreach_dimension () {
	if (point.l > 0) {
	  double area, circ;
	  area = Delta*(hf.y[] + hf.y[0,0,-1] + hf.y[0,1,0] + hf.y[0,1,-1])/4.; // Add fm later
	  circ = (-u.y[] + u.y[0,0,-1])*Delta - 0.5*(w[0,0,-1] + w[0,-1,-1])*0.5*(hf.y[] + hf.y[0,0,-1]) + \
	    0.5*(w[0,0,-1] + w[0,1,-1])*0.5*(hf.y[0,1,0] + hf.y[0,1,-1]);
	  Omega.x[] = circ/area;      
	}
	else
	  Omega.x[] = 0.; // not well defined for the bottom layer
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
	dzdx.x[] = dz.x + hf.x[1] - hf.x[]; // It's the same if we had used h[] here
	dz.x += hf.x[1] - hf.x[];
      }
    }
  }  
  delete ((scalar *){hu,ha,hf});
}

vector dzdxc; // centered dzdx
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
	dzdxc.x[] = dz.x + hf.x[1] - hf.x[];
	dz.x += hf.x[1] - hf.x[];
      }
    }
  }
}

/**
## Write to files
A new writefields function with Omega added. */

int writefields (double t, const char *suffix) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50], filename5[50], filename6[50], filename7[50], filename8[50], filename9[50], filename10[50], filename11[50];
  vector u_temp, Omega_temp, dzdx_temp, dzdxc_temp;
  scalar w_temp, h_temp, Omegaz_temp;
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%g_l%d", suffix, t, j);
    sprintf (filename2, "field/uy_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename3, "field/uz_%s_t%g_l%d", suffix, t, j);  
    sprintf (filename4, "field/h_%s_t%g_l%d", suffix, t, j);
    sprintf (filename5, "field/Omegax_%s_t%g_l%d", suffix, t, j);
    sprintf (filename6, "field/Omegay_%s_t%g_l%d", suffix, t, j);
    sprintf (filename7, "field/Omegaz_%s_t%g_l%d", suffix, t, j);
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
      sprintf (s, "Omega");
      Omega_temp = lookup_vector (s);
      sprintf (s, "Omegaz");
      Omegaz_temp = lookup_field (s);
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
      sprintf (s, "Omega%d", j);
      Omega_temp = lookup_vector (s);
      sprintf (s, "Omegaz%d", j);
      Omegaz_temp = lookup_field (s);
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
    FILE * fOmegax = fopen (filename5, "w");
    output_matrix_mpi (Omega_temp.x, fOmegax, N, linear = true);
    fclose (fOmegax);
    FILE * fOmegay = fopen (filename6, "w");
    output_matrix_mpi (Omega_temp.y, fOmegay, N, linear = true);
    fclose (fOmegay);
    FILE * fOmegaz = fopen (filename7, "w");
    output_matrix_mpi (Omegaz_temp, fOmegaz, N, linear = true);
    fclose (fOmegaz);
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
