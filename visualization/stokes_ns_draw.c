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

// bool save_PS (struct _save p); //prototype

/** Draw the interface while taking a slice.  */

void draw_wave1 () {
  scalar uwater[];
  foreach () {
    uwater[] = u.x[]*f[];
  }
  clear ();
  view (fov = 19.5, width = 1200, height = 1200, bg = {1., 1., 1.}, camera = "front");
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
  view (fov = 19.5, width = 1200, height = 1200, bg = {1., 1., 1.}, camera = "right");
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
  view (fov = 19.5, width = 1200, height = 1200, bg = {1., 1., 1.}, camera = "front");
  squares ("vwater", linear = true, n = {0,0,1}, alpha = xlocation, min = -0.124888, max = 0.124888, map = blue_white_red);
  cross_section ("f", alpha = xlocation, np = {0,0,1}, lc = {0.05,0.05,0.05}, lw = 10);
  save ("slicey_uy.ppm");
}

/* void draw_command (void) { */
/*   scalar vwater[]; */
/*   foreach () { */
/*     vwater[] = u.z[]*f[]; */
/*   } */
/*   clear (); */
/*   squares ("vwater", linear = true, n = {0,0,1}, alpha = xlocation, min = -0.124888, max = 0.124888, map = blue_white_red); */
/*   cross_section ("f", alpha = xlocation, np = {0,0,1}, lc = {0.05,0.05,0.05}, lw = 10); */
/* } */

/* void draw_vector () { */
/*   view (fov = 19.5, width = 600, height = 600, bg = {1., 1., 1.}, camera = "front"); */
/*   save_PS ("slicey_uy.pdf"); */
/* }  */

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
  // draw_vector ();
  /* char name[80]; */
  /* save_data ((scalar*){f,p}, (vector*){u}, (int)snapshot_time, t, "pvta"); */
  return 0;
}


/** Saving the vector plot (from Antoon) */
/* bool save_PS (struct _save p) { */
/*   char ppm[] = "ppm"; */
/*   if (!p.format) { */
/*     p.format = ppm; */
/*     if (p.file) { */
/*       char * s = strchr (p.file, '.'), * dot = s; */
/*       while (s) { */
/* 	dot = s; */
/* 	s = strchr (s + 1, '.'); */
/*       } */
/*       if (dot) */
/* 	p.format = dot + 1; */
/*     } */
/*   } */
/*   bview * view = p.view ? p.view : get_view(); */
/*   if (p.file && (p.fp = fopen (p.file, "w")) == NULL) { */
/*     perror (p.file); */
/*     return false; */
/*   } */
/*   if (!p.fp) */
/*     p.fp = stdout; */
/*   if (!strcmp (p.format, "ps") || */
/*       !strcmp (p.format, "eps") || */
/*       !strcmp (p.format, "tex") || */
/*       !strcmp (p.format, "pdf") || */
/*       !strcmp (p.format, "svg") || */
/*       !strcmp (p.format, "pgf")) { */
/*     GLint format = (!strcmp (p.format, "ps") ? GL2PS_PS : */
/* 		    !strcmp (p.format, "eps") ? GL2PS_EPS : */
/* 		    !strcmp (p.format, "tex") ? GL2PS_TEX : */
/* 		    !strcmp (p.format, "pdf") ? GL2PS_PDF : */
/* 		    !strcmp (p.format, "svg") ? GL2PS_SVG : */
/* 		    !strcmp (p.format, "pgf") ? GL2PS_PGF : */
/* 		    -1); */
/*     GLint state = GL2PS_OVERFLOW; */
/*     GLint sort = p.sort ? p.sort : GL2PS_SIMPLE_SORT; */
/*     GLint options = p.options ? p.options : (GL2PS_SIMPLE_LINE_OFFSET | */
/* 					          GL2PS_SILENT | */
/* 					          GL2PS_BEST_ROOT | */
/* 					          GL2PS_OCCLUSION_CULL | */
/* 					          GL2PS_USE_CURRENT_VIEWPORT | */
/* 					     GL2PS_TIGHT_BOUNDING_BOX); */
/*     unsigned buffsize = 1 << 24; */
/*     while (state == GL2PS_OVERFLOW && buffsize <= MAXBUFFSIZE) { */
/*       gl2psBeginPage ("", "bview", */
/* 		      NULL, */
/* 		      format, sort, options,  */
/* 		      GL_RGBA, 0, NULL,  */
/* 		      0, 0, 0, */
/* 		      buffsize, p.fp, ""); */
      
/*       float res = view->res; */
/*       view->res = 0.; */
/*       view->vector=true; */
    
/*       draw_commands(); */
       
/*       glFinish (); */
/*       enable_fpe (FE_DIVBYZERO|FE_INVALID); */
/*       view->active = false; */
/*       view->vector = false; */
/*       view->res = res; */
/*       //draw(); */
/*       disable_fpe (FE_DIVBYZERO|FE_INVALID); */
/*       state = gl2psEndPage(); */
/*       enable_fpe (FE_DIVBYZERO|FE_INVALID); */
/*       buffsize *= 2; */
/*     } */
/*     if (state == GL2PS_OVERFLOW) */
/*       fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n"); */
/*   } */

/*   else { */
/*     fprintf (ferr, "save(): unknown format '%s'\n", p.format); */
/*     if (p.file) { */
/*       fclose (p.fp); */
/*       remove (p.file); */
/*     } */
/*     return false; */
/*   } */
  
/*   fflush (p.fp); */
/*   if (p.file) */
/*     fclose (p.fp); */
/*   return true; */
/* } */
