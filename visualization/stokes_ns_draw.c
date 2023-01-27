#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  //reduced gravity
#include "view.h"
#include "lambda2.h"
#include "iso3D.h" // From Antoon, compute the 2D vof slice
//#include "save_data.h" // From Palas, for paraview


double snapshot_time = 0.; // Snapshot time
double xlocation = 0.; // Where we take the slice

#if _MPI
#  define POPEN(name, mode) fopen (name ".ppm", mode)
#else
#  define POPEN(name, mode) popen ("ppm2mp4 > " name ".mp4", mode)
#endif

void draw_wave1 () {
  scalar uwater[];
  foreach () {
    uwater[] = u.x[]*f[];
  }
  clear ();
  view (fov = 19.5, width = 600, height = 600, bg = {1., 1., 1.}, camera = "front");
  squares ("uwater", linear = true, n = {0,0,1}, alpha = xlocation, min = -0.624443, max = 0.624443, map = blue_white_red);
  cross_section ("f", alpha = xlocation, np = {0,0,1}, lc = {0.05,0.05,0.05}, lw = 10);
  save ("slicey_ux.ppm");
}

void draw_wave2 () {
  scalar vwater[];
  foreach () {
    vwater[] = u.z[]*f[];
  }
  clear ();
  view (fov = 19.5, width = 600, height = 600, bg = {1., 1., 1.}, camera = "right");
  squares ("vwater", linear = true, n = {1,0,0}, alpha = xlocation, min = -0.124888, max = 0.124888, map = blue_white_red);
  cross_section ("f", alpha = xlocation, np = {1,0,0}, lc = {0.05,0.05,0.05}, lw = 10);
  save ("slicex_uy.ppm");
}

void draw_wave3 () {
  scalar vwater[];
  foreach () {
    vwater[] = u.z[]*f[];
  }
  clear ();
  view (fov = 19.5, width = 600, height = 600, bg = {1., 1., 1.}, camera = "front");
  squares ("vwater", linear = true, n = {0,0,1}, alpha = xlocation, min = -0.124888, max = 0.124888, map = blue_white_red);
  cross_section ("f", alpha = xlocation, np = {0,0,1}, lc = {0.05,0.05,0.05}, lw = 10);
  save ("slicey_uy.ppm");
}

int main(int argc, char *argv[]) {
  if (argc > 1)
    snapshot_time = atof(argv[1]);
  if (argc > 2)
    xlocation = atof(argv[2]);
  fprintf(ferr, "Time at %gT, slice at z = %g\n", snapshot_time, xlocation);
  run();
}

event init(i=0)
{
  char targetname[100];
  sprintf (targetname, "dump%g", snapshot_time);
  if (!restore (targetname)) {
    fprintf(ferr, "Not restored!\n");
    return 1;
  }
  draw_wave1 ();
  draw_wave2 ();
  draw_wave3 ();
  /* char name[80]; */
  /* save_data ((scalar*){f,p}, (vector*){u}, (int)snapshot_time, t, "pvta"); */
  return 0;
}
