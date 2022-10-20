#line 1 "field_new_PM_restore-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "field_new_PM_restore-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif




#line 1 "/home/jiarongw/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#ifndef assert
# include <assert.h>
#endif
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif 1

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) do { type __tmp = a; a = b; b = __tmp; } while(0)
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout
#define sysfprintf fprintf

#if 1
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 84

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



static inline void qassert (const char * file, int line, const char * cond) {
  fprintf (ferr, "%s:%d: Assertion `%s' failed.\n", file, line, cond);
  abort();
}
#line 108 "/home/jiarongw/basilisk/src/common.h"
#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  if (!(d != NULL)) qassert ("/home/jiarongw/basilisk/src/common.h", 183, "d != NULL");
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free ()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  if (!(t->stack.len > 0)) qassert ("/home/jiarongw/basilisk/src/common.h", 455, "t->stack.len > 0");
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off ()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if 1
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  if (!(Trace.stack.len >= 2*sizeof(double))) qassert ("/home/jiarongw/basilisk/src/common.h", 559, "Trace.stack.len >= 2*sizeof(double)");
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if 1
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if 1
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if 1
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off ()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#elif 1

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  if (!(!in_prof)) qassert ("/home/jiarongw/basilisk/src/common.h", 669, "!in_prof"); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 671

#define prof_stop()\
  if (!(in_prof)) qassert ("/home/jiarongw/basilisk/src/common.h", 673, "in_prof"); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 676


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 685);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 686);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/jiarongw/basilisk/src/common.h", 687); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 695

#define mpi_all_reduce_array(v,type,op,elem) {\
  prof_start ("mpi_all_reduce");\
  type global[elem], tmp[elem];\
  for (int i = 0; i < elem; i++)\
    tmp[i] = (v)[i];\
  MPI_Datatype datatype;\
  if (!strcmp(#type, "double")) datatype = MPI_DOUBLE;\
  else if (!strcmp(#type, "int")) datatype = MPI_INT;\
  else if (!strcmp(#type, "long")) datatype = MPI_LONG;\
  else if (!strcmp(#type, "bool")) datatype = MPI_C_BOOL;\
  else {\
    fprintf (ferr, "unknown reduction type '%s'\n", #type);\
    fflush (ferr);\
    abort();\
  }\
  mpi_all_reduce0 (tmp, global, elem, datatype, op, MPI_COMM_WORLD);\
  for (int i = 0; i < elem; i++)\
    (v)[i] = global[i];\
  prof_stop();\
}\

#line 716


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init ()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_array(v,type,op,elem)

#endif

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 828 "/home/jiarongw/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  if (!(sizeof (int64_t) == sizeof (double))) qassert ("/home/jiarongw/basilisk/src/common.h", 838, "sizeof (int64_t) == sizeof (double)");
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;


int N = 64;




typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;




} vector;

typedef struct {
  scalar * x;

  scalar * y;




} vectorl;

typedef struct {
  vector x;

  vector y;




} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;

OMP(omp declare reduction (+ : coord :
      omp_out.x += omp_in.x,
      omp_out.y += omp_in.y,
      omp_out.z += omp_in.z))
#line 922 "/home/jiarongw/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 925

    norm += sq(n->x);
#line 925

    norm += sq(n->y);}
  norm = sqrt(norm);
  {
#line 928

    n->x /= norm;
#line 928

    n->y /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }






  enum { right, left, top, bottom };



int nboundary = 2*2;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;



#define strongif(x) if(x)
#define IF(x) if((x)||1)

typedef struct { int i[2]; } IJK;

void _stencil_fprintf (const char * file, int line,
         FILE * stream, const char *format, ...) {}
void _stencil_printf (const char * file, int line,
        const char *format, ...) {}
void _stencil_fputc (const char * file, int line,
       int c, FILE * stream) {}
void _stencil_fputs (const char * file, int line,
       const char * s, FILE * stream) {}
#define _stencil_qassert(...) ((void)0)

int _stencil_access (scalar s, IJK i, const char * file, int line);

#if 2 == 1
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1}}, file, line)])\

#line 989

#elif 2 == 2
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1, point.j + i2}}, file, line)])\

#line 994

#else
#define _stencil_val(file, line, s, i1, i2, i3)\
  (_attribute[is_constant(s) || s.i < 0 ? -1 : s.i].write\
   [_stencil_access(s, (IJK){{point.i + i1, point.j + i2, point.k + i3}},\
      file, line)])\

#line 1000

#endif




#define _stencil_fine(file, line, s, i1, i2, i3)\
  _stencil_val(file, line, s, i1, i2, i3)\

#line 1008

#define _stencil_coarse(file, line, s, i1, i2, i3)\
  _stencil_val(file, line, s, i1, i2, i3)\

#line 1011


#line 1 "/home/jiarongw/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries () {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/jiarongw/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 1014 "/home/jiarongw/basilisk/src/common.h"



typedef struct {

#line 1019 "/home/jiarongw/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;




  } d;
  vector v;
  int face;
  bool nodump, freed;
  int block;
  scalar * depends;

#line 18 "/home/jiarongw/basilisk/src/grid/stencils.h"

  double * write;
  int * read;
  int dirty;




#line 17 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

} _Attributes;
_Attributes * _attribute;
#define _call_vertical_diffusion 1
#define _call_evaluate_expression 1
#define _call_facet_normal 1
#define _call_mycs 1
#define _call_periodic_bc 1
#if _call_glvertex_normal3d
#  undef _call_interp
#  define _call_interp 1
#endif
#if _call_evaluate_expression
#  undef _call_evaluate_expression
#  define _call_evaluate_expression 1
#endif
#if _call_interp
#  undef _call_interpolate_linear
#  define _call_interpolate_linear 1
#endif
#if _call_facet_normal
#  undef _call_mycs
#  define _call_mycs 1
#endif
#if _call_multigrid_debug
#  undef _call_cartesian_debug
#  define _call_cartesian_debug 1
#endif
#if _call_refine_biquadratic
#  undef _call_biquadratic
#  define _call_biquadratic 1
#endif
#if _call_refine_bilinear
#  undef _call_bilinear
#  define _call_bilinear 1
#endif
#if _call_restriction_face
#  undef _call_face_average
#  define _call_face_average 1
#endif
#line 1017 "/home/jiarongw/basilisk/src/common.h"





























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  strongif (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  strongif (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    strongif (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    strongif (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  strongif (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  strongif (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  strongif (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  strongif (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 1137

      if (w.x.i != v.x.i)
 id = false;
#line 1137

      if (w.y.i != v.y.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    strongif (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1160
 {
      if (!(s->i >= 0)) qassert ("/home/jiarongw/basilisk/src/common.h", 1161, "s->i >= 0");
      v.x = *s++;
    }
#line 1160
 {
      if (!(s->i >= 0)) qassert ("/home/jiarongw/basilisk/src/common.h", 1161, "s->i >= 0");
      v.y = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  strongif (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1191
 {
      if (!(v->x.i >= 0)) qassert ("/home/jiarongw/basilisk/src/common.h", 1192, "v->x.i >= 0");
      t.x = *v++;
    }
#line 1191
 {
      if (!(v->y.i >= 0)) qassert ("/home/jiarongw/basilisk/src/common.h", 1192, "v->x.i >= 0");
      t.y = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;
scalar * baseblock = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1}};
vector unityf= {{_NVARMAX + 2},{_NVARMAX + 3}};
scalar unity= {_NVARMAX + 4};
scalar zeroc= {_NVARMAX + 5};



 vector fm = {{_NVARMAX + 2},{_NVARMAX + 3}};
 scalar cm = {(_NVARMAX + 4)};
#line 1296 "/home/jiarongw/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all ()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



void init_solver ()
{
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

void allocate_globals (int nvar)
{
  _attribute = (_Attributes *) pcalloc (nvar + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
  _attribute[0].write = (double *) pcalloc (1, sizeof(double),__func__,__FILE__,__LINE__);
  _attribute++;
  all = (scalar *) pmalloc (sizeof (scalar)*(nvar + 1),__func__,__FILE__,__LINE__);
  baseblock = (scalar *) pmalloc (sizeof (scalar)*(nvar + 1),__func__,__FILE__,__LINE__);
  for (int i = 0; i < nvar; i++)
    baseblock[i].i = all[i].i = i;
  baseblock[nvar].i = all[nvar].i = -1;
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults () {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, p.commands);
  }
}




#line 1 "/home/jiarongw/basilisk/src/grid/layers.h"






int nl = 1;





int _layer = 0;


#define foreach_block() for (_layer = 0; _layer < nl; _layer++)
#define end_foreach_block() _layer = 0





#define foreach_block_inner() for (point.l = 0; point.l < nl; point.l++)
#define end_foreach_block_inner() point.l = 0





#define foreach_blockf(_f)\
  for (point.l = 0; point.l < _attribute[_f.i].block; point.l++)\

#line 32

#define end_foreach_blockf() point.l = 0





#undef _index
#define _index(a,m)\
  (a.i + (_layer + point.l + m < _attribute[a.i].block ?\
   _layer + point.l + m : 0))\

#line 43


#undef val
#define val(a,k,p,m) data(k,p,m)[_index(a,m)]
#line 1495 "/home/jiarongw/basilisk/src/common.h"
#line 15 "field_new_PM_restore-cpp.c"
#line 1 "field_new_PM_restore.c"




#line 1 "grid/multigrid.h"
#line 1 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#line 16 "/home/jiarongw/basilisk/src/grid/multigrid.h"
typedef struct {
  Grid g;
  char ** d;
} Multigrid;

struct _Point {
  int i;

  int j;




  int level, n;
#ifdef foreach_block
  int l;
#define _BLOCK_INDEX , point.l
#else
#define _BLOCK_INDEX
#endif
};
static Point last_point;
#line 49 "/home/jiarongw/basilisk/src/grid/multigrid.h"
static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*2;
  return sq(n);
}
#line 62 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define data(k,l,m)\
  ((double *)&((Multigrid *)grid)->d[point.level][((point.i + k)*((1 << point.level) +\
       2*2) +\
      (point.j + l))*datasize]) 
#line 64

#line 89 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define allocated(k,l,m) (point.i+k >= 0 && point.i+k < (1 << point.level) + 2*2 &&\
         point.j+l >= 0 && point.j+l < (1 << point.level) + 2*2)\

#line 91


#define allocated_child(k,l,m) (level < depth() &&\
         point.i > 0 && point.i <= (1 << point.level) + 2 &&\
         point.j > 0 && point.j <= (1 << point.level) + 2)\

#line 96

#line 117 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define depth() (grid->depth)
#line 136 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define fine(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level+1][((2*point.i-2 +k)*2*((1 << point.level) +\
        2) +\
     (2*point.j-2 +l))*datasize])[_index(a,m)]\

#line 141

#define coarse(a,k,l,m)\
  ((double *)\
   &((Multigrid *)grid)->d[point.level-1][(((point.i+2)/2+k)*((1 << point.level)/2 +\
        2*2) +\
     (point.j+2)/2+l)*datasize])[_index(a,m)]\

#line 147

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
  struct { int x, y; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1\
  }; NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2; parent.j = (point.j + 2)/2;\

#line 157

#line 191 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define foreach_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 208

#define end_foreach_level()\
\
 }\
\
  }\
}\

#line 215


#define foreach()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
\
\
 {\
\
          POINT_VARIABLES\

#line 234

#define end_foreach()\
\
 }\
\
  }\
}\

#line 241


#define is_active(cell) (true)
#define is_leaf(cell) (level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (ferr, "grid depths do not match. Aborting.\n");\
  if (!(0)) qassert ("/home/jiarongw/basilisk/src/grid/multigrid.h", 249, "0");\
} while (0)\

#line 251

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/jiarongw/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/jiarongw/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
  Point root = {2,2,0};\
\
\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
    struct { int l, i, j, stage; } stack[20];\
\
\
\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
\
\
\
\
\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
    Point root = {2,2,0};\
\
\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 254 "/home/jiarongw/basilisk/src/grid/multigrid.h"

#define foreach_face_generic()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k <= point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j <= point.n + 2; point.j++)\
\
\
\
        {\
\
   POINT_VARIABLES\

#line 272

#define end_foreach_face_generic()\
\
 }\
\
  }\
}\

#line 279


#define foreach_vertex()\
foreach_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
\
\

#line 290

#define end_foreach_vertex() } end_foreach_face_generic()

#define is_coarse() (point.level < depth())
#line 320 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define is_face_x() (point.j < point.n + 2)
#define is_face_y() (point.i < point.n + 2)

#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2;\
  point.level++;\
  point.n *= 2;\
  for (int _k = 0; _k < 2; _k++)\
    for (int _l = 0; _l < 2; _l++) {\
      point.i = _i + _k; point.j = _j + _l;\
      POINT_VARIABLES;\

#line 331

#define end_foreach_child()\
  }\
  point.i = (_i + 2)/2; point.j = (_j + 2)/2;\
  point.level--;\
  point.n /= 2;\
}\

#line 338

#define foreach_child_break() _k = _l = 2
#line 381 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/jiarongw/basilisk/src/grid/neighbors.h"
#line 17 "/home/jiarongw/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j;\
  for (int _k = - _nn; _k <= _nn; _k++) {\
    point.i = _i + _k;\
    for (int _l = - _nn; _l <= _nn; _l++) {\
      point.j = _j + _l;\
      POINT_VARIABLES;\

#line 25

#define end_foreach_neighbor()\
    }\
  }\
  point.i = _i; point.j = _j;\
}\

#line 31

#define foreach_neighbor_break() _k = _l = _nn + 1
#line 387 "/home/jiarongw/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Point p;
  p.level = depth(); p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < sq(p.n + 2*2); i++)
      strongif (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10) {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     ((double *)(&((Multigrid *)grid)->d[p.level][i*datasize]))[s.i + b] = val;
      }
}
#line 427 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define foreach_boundary_dir(l,d)\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l < 0 ? depth() : l;\
  point.n = 1 << point.level;\
  int * _i = &point.j;\
  if (d == left) {\
    point.i = 2;\
    ig = -1;\
  }\
  else if (d == right) {\
    point.i = point.n + 2 - 1;\
    ig = 1;\
  }\
  else if (d == bottom) {\
    point.j = 2;\
    _i = &point.i;\
    jg = -1;\
  }\
  else if (d == top) {\
    point.j = point.n + 2 - 1;\
    _i = &point.i;\
    jg = 1;\
  }\
  int _l;\
  OMP(omp for schedule(static))\
  for (_l = 0; _l < point.n + 2*2; _l++) {\
    *_i = _l;\
    {\
      POINT_VARIABLES\

#line 458

#define end_foreach_boundary_dir()\
    }\
  }\
}\

#line 463


#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.level, point.n _BLOCK_INDEX})\

#line 467

#define is_boundary(point) (point.i < 2 || point.i >= point.n + 2 ||\
    point.j < 2 || point.j >= point.n + 2)\

#line 470

#line 532 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#define foreach_boundary(b)\
  if (default_scalar_bc[b] != periodic_bc)\
    foreach_boundary_dir (depth(), b)\
      if (!is_boundary(point)) {\

#line 536

#define end_foreach_boundary() } end_foreach_boundary_dir()

#define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, void * data);
#if _call_periodic_bc
static double _periodic_bc (Point point, Point neighbor, scalar s, void * data);
#endif

#line 541

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int bghost = 1; bghost <= 2; bghost++)
    for (int d = 0; d < 2*2; d++) {

      scalar * list = NULL, * listb = NULL;
      strongif (scalars) for (scalar s = *scalars, *_i11 = scalars; ((scalar *)&s)->i >= 0; s = *++_i11)
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 2) % 2];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }

      if (list) {
  { foreach_boundary_dir (l, d){

#line 568 "/home/jiarongw/basilisk/src/grid/multigrid.h"
 {
   scalar s, sb;
   scalar * _i0 = list; scalar * _i1 = listb; strongif (list) for (s = *list, sb = *listb; ((scalar *)&s)->i >= 0; s = *++_i0, sb = *++_i1) {
     if (_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) {

       if (bghost == 1)
   { foreach_block_inner()
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL); end_foreach_block_inner(); }
     }
     else

        { foreach_block_inner()
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL); end_foreach_block_inner(); }
   }
 } } end_foreach_boundary_dir(); }
 pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 703 "/home/jiarongw/basilisk/src/grid/multigrid.h"
void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++)
    pfree (m->d[l],__func__,__FILE__,__LINE__);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  N = 1 << depth();

  grid->n = grid->tn = 1 << 2*depth();

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);

  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#line 748 "/home/jiarongw/basilisk/src/grid/multigrid.h"
  m->d = (char **) pmalloc(sizeof(Point *)*(depth() + 1),__func__,__FILE__,__LINE__);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l)*datasize;
    m->d[l] = (char *) pmalloc (len,__func__,__FILE__,__LINE__);


    double * v = (double *) m->d[l];
    for (int i = 0; i < len/sizeof(double); i++)
      v[i] = undefined;
  }
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l);
    p->d[l] = (char *) prealloc (p->d[l], (len*(datasize + size))*sizeof(char),__func__,__FILE__,__LINE__);
    char * data = p->d[l] + (len - 1)*datasize;
    for (int i = len - 1; i > 0; i--, data -= datasize)
      memmove (data + i*size, data, datasize);
  }
  datasize += size;
}


int mpi_dims[2], mpi_coords[2];
#line 786 "/home/jiarongw/basilisk/src/grid/multigrid.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 790 "/home/jiarongw/basilisk/src/grid/multigrid.h"

  point.level = -1, point.n = 1 << depth();

  point.i = (p.x - X0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[0]*point.n;
  if (point.i < 2 || point.i >= point.n + 2)
    return point;

  point.j = (p.y - Y0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[1]*point.n;
  if (point.j < 2 || point.j >= point.n + 2)
    return point;
#line 821 "/home/jiarongw/basilisk/src/grid/multigrid.h"
  point.level = depth();
  return point;
}

#line 1 "grid/multigrid-common.h"
#line 1 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/jiarongw/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  if (!(Events)) qassert ("/home/jiarongw/basilisk/src/grid/events.h", 87, "Events");
  if (!(!event.last)) qassert ("/home/jiarongw/basilisk/src/grid/events.h", 88, "!event.last");
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      if (!(parent < 0)) qassert ("/home/jiarongw/basilisk/src/grid/events.h", 92, "parent < 0");
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/jiarongw/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    if (!(n < INT_MAX)) qassert ("/home/jiarongw/basilisk/src/grid/events.h", 245, "n < INT_MAX");
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])
#define diagonalize(a)
#define val_diagonal(a,k,l,m) ((k) == 0 && (l) == 0 && (m) == 0)

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level)/mpi_dims[0]);\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
\
  double x = (ig/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
\
\
  double z = 0.;\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
\
  ;\

#line 34


#line 1 "grid/fpe.h"
#line 1 "/home/jiarongw/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb ()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 37 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

#line 1 "grid/stencils.h"
#line 1 "/home/jiarongw/basilisk/src/grid/stencils.h"
#line 18 "/home/jiarongw/basilisk/src/grid/stencils.h"








#line 43 "/home/jiarongw/basilisk/src/grid/stencils.h"
typedef struct {
  const char * fname;
  int line;
  int first;
  const char * each;
  int face;
  bool vertex;
} ForeachData;
#line 59 "/home/jiarongw/basilisk/src/grid/stencils.h"
static inline int stencil_index (scalar s, IJK i)
{
  int len = 1, index = 0;
  for (int d = 0; d < 2; d++) {
    if (i.i[d] < 0 || i.i[d] >= 5)
      return -1;
    index += len*i.i[d], len *= 5;
  }
  return index;
}




static inline IJK stencil_ijk (int index)
{
  IJK i;
  int len = (5*5);
  for (int d = 2 - 1; d >= 0; d--) {
    len /= 5;
    i.i[d] = index/len;
    index -= len*i.i[d];
    i.i[d] -= 5/2;
  }
  return i;
}




static void write_stencil_index (IJK i, int shift)
{
  sysfprintf (qstderr(), "[%d", i.i[0] - shift);
  for (int d = 1; d < 2; d++)
    sysfprintf (qstderr(), ",%d", i.i[d] - shift);
  sysfprintf (qstderr(), "]");
}






int _stencil_access (scalar s, IJK i, const char * file, int line)
{
  if (is_constant(s) || s.i < 0)
    return 0;
  int index = stencil_index (s, i);
  if (index < 0) {
    sysfprintf (qstderr(), "%s:%d: error: stencil overflow: %s",
  file, line, _attribute[s.i].name);
    write_stencil_index (i, 5/2);
    sysfprintf (qstderr(), "\n");
    fflush (qstderr());
    abort();
  }
  _attribute[s.i].read[index]++;
  return index;
}
#line 127 "/home/jiarongw/basilisk/src/grid/stencils.h"
#define foreach_stencil() {\
  strongif (baseblock) for (scalar _s = *baseblock, *_i12 = baseblock; ((scalar *)&_s)->i >= 0; _s = *++_i12) {\
    for (int _i = 0; _i < (5*5); _i++) {\
      _attribute[_s.i].read[_i] = 0;\
      _attribute[_s.i].write[_i] = 1.7759437274 + _s.i + _i;\
    }\
  }\
  _foreach_data.face = 0;\
  int ig = 0, jg = 0, kg = 0;\
  NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
\
\
\
\
  for (int _i = 0; _i < 2; _i++)\
    ((int *)&point)[_i] = 5/2;\
  if (sizeof(Point) >= (2 + 2)*sizeof(int))\
    ((int *)&point)[2 + 1] = 1;\
  POINT_VARIABLES\

#line 147






#define end_foreach_stencil()\
  strongif (baseblock) for (scalar _s = *baseblock, *_i13 = baseblock; ((scalar *)&_s)->i >= 0; _s = *++_i13) {\
    for (int _i = 0; _i < (5*5); _i++)\
      _attribute[_s.i].write[_i] = (_attribute[_s.i].write[_i] != 1.7759437274 + _s.i + _i);\
  }\
  end_stencil (&_foreach_data);\
}\

#line 160


#define foreach_face_stencil foreach_stencil
#define end_foreach_face_stencil() end_foreach_stencil()

#define foreach_vertex_stencil foreach_stencil
#define end_foreach_vertex_stencil() end_foreach_stencil()

#define is_stencil_face_x() ((_foreach_data.face |= (1 << 0)))
#define is_stencil_face_y() ((_foreach_data.face |= (1 << 1)))
#define is_stencil_face_z() ((_foreach_data.face |= (1 << 2)))




void reduction_warning (const char * fname, int line, const char * var)
{
  fprintf (ferr,
  "%s:%d: warning: variable '%s' is modified by this foreach loop:\n"
  "%s:%d: warning: use a loop-local variable, a reduction operation\n"
  "%s:%d: warning: or a serial loop to get rid of this warning\n",
    fname, line, var, fname, line, fname, line);
}







static inline bool scalar_is_dirty (scalar s)
{
  if (_attribute[s.i].dirty)
    return true;
  scalar * depends = _attribute[s.i].depends;
  strongif (depends) for (scalar d = *depends, *_i14 = depends; ((scalar *)&d)->i >= 0; d = *++_i14)
    if (_attribute[d.i].dirty)
      return true;
  return false;
}




static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = _attribute[a.i].depends;
  strongif (depends) for (scalar s = *depends, *_i15 = depends; ((scalar *)&s)->i >= 0; s = *++_i15)
    if (s.i == b.i)
      return true;
  return false;
}







void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face) (vectorl);







void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  loop->vertex = !strcmp (loop->each, "foreach_vertex");




  strongif (baseblock) for (scalar s = *baseblock, *_i16 = baseblock; ((scalar *)&s)->i >= 0; s = *++_i16) {
    bool write = false, read = false;
    int max = 0;




    for (int n = 0; n < (5*5); n++)
      if (_attribute[s.i].write[n] || _attribute[s.i].read[n]) {
 IJK i = stencil_ijk (n);





 if (_attribute[s.i].write[n]) {
   for (int d = 0; d < 2; d++)
     if (i.i[d] != 0) {
       fprintf (ferr,
         "%s:%d: error: illegal write within this loop: %s",
         loop->fname, loop->line, _attribute[s.i].name);
       write_stencil_index (i, 0);
       fprintf (ferr, "\n");
       fflush (ferr);
       abort();
     }
   write = true;
 }





 if (_attribute[s.i].read[n]) {
   read = true;
   int d = 0;
   {
#line 274
 {
     if ((!_attribute[s.i].face || _attribute[s.i].v.x.i != s.i) && abs(i.i[d]) > max)
       max = abs(i.i[d]);
     d++;
   }
#line 274
 {
     if ((!_attribute[s.i].face || _attribute[s.i].v.y.i != s.i) && abs(i.i[d]) > max)
       max = abs(i.i[d]);
     d++;
   }}
 }
      }





    if (_layer == 0 || _attribute[s.i].block == 1)

    {





      if (read && scalar_is_dirty (s)) {





 if (_attribute[s.i].face) {
   if (max > 0)
     listc = list_append (listc, s);
   else if (!write) {
     scalar sn = _attribute[s.i].v.x.i >= 0 ? _attribute[s.i].v.x : s;
     {
#line 305

       if (_attribute[s.i].v.x.i == s.i) {




  if (_attribute[sn.i].boundary[left] || _attribute[sn.i].boundary[right])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.x = list_append (listf.x, s);
    flux = true;
  }
       }
#line 305

       if (_attribute[s.i].v.y.i == s.i) {




  if (_attribute[sn.i].boundary[bottom] || _attribute[sn.i].boundary[top])
    listc = list_append (listc, s);
  else if (_attribute[s.i].dirty != 2) {
    listf.y = list_append (listf.y, s);
    flux = true;
  }
       }}
   }
 }





 else if (max > 0)
   listc = list_append (listc, s);
      }





      if (write) {
 if (2 > 1 && !loop->vertex && loop->first) {
   bool vertex = true;
   {
#line 336

     if (_attribute[s.i].d.x != -1)
       vertex = false;
#line 336

     if (_attribute[s.i].d.y != -1)
       vertex = false;}
   if (vertex)
     fprintf (ferr,
       "%s:%d: warning: vertex scalar '%s' should be assigned with"
       " a foreach_vertex() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 if (_attribute[s.i].face) {
   if (loop->face == 0 && loop->first)
     fprintf (ferr,
       "%s:%d: warning: face vector '%s' should be assigned with"
       " a foreach_face() loop\n",
       loop->fname, loop->line, _attribute[s.i].name);
 }
 else if (loop->face) {
   if (_attribute[s.i].v.x.i < 0) {
     int d = 1, i = 0;
     {
#line 355
 {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.x.i = s.i;
  _attribute[s.i].boundary[left] = _attribute[s.i].boundary[right] = NULL;





       }
       d *= 2, i++;
     }
#line 355
 {
       if (loop->face == d) {
  _attribute[s.i].face = 2, _attribute[s.i].v.y.i = s.i;
  _attribute[s.i].boundary[bottom] = _attribute[s.i].boundary[top] = NULL;





       }
       d *= 2, i++;
     }}
     if (!_attribute[s.i].face && loop->first)
       fprintf (ferr,
         "%s:%d: warning: scalar '%s' should be assigned with "
         "a foreach_face(x|y|z) loop\n",
         loop->fname, loop->line, _attribute[s.i].name);
   }
   else {
     char * name = NULL;
     if (_attribute[s.i].name) {
       name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
       char * s = name + strlen(name) - 1;
       while (s != name && *s != '.') s--;
       if (s != name) *s = '\0';
     }
     init_face_vector (_attribute[s.i].v, name);




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }
 else if (loop->vertex) {
   bool vertex = true;
   {
#line 391

     if (_attribute[s.i].d.x != -1)
       vertex = false;
#line 391

     if (_attribute[s.i].d.y != -1)
       vertex = false;}
   if (!vertex) {
     char * name = NULL;
     if (_attribute[s.i].name) name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     init_vertex_scalar (s, name);
     {
#line 398

       _attribute[s.i].v.x.i = -1;
#line 398

       _attribute[s.i].v.y.i = -1;}




     pfree (name,__func__,__FILE__,__LINE__);
   }
 }





 dirty = list_append (dirty, s);
 strongif (baseblock) for (scalar d = *baseblock, *_i17 = baseblock; ((scalar *)&d)->i >= 0; d = *++_i17)
   if (scalar_depends_from (d, s))
     dirty = list_append (dirty, d);
      }
    }
  }




  if (flux) {
#line 436 "/home/jiarongw/basilisk/src/grid/stencils.h"
    boundary_face (listf);
    {
#line 437

      pfree (listf.x,__func__,__FILE__,__LINE__);
#line 437

      pfree (listf.y,__func__,__FILE__,__LINE__);}
  }




  if (listc) {






    boundary_internal (listc, loop->fname, loop->line);
    pfree (listc,__func__,__FILE__,__LINE__);
  }





  if (dirty) {






    strongif (dirty) for (scalar s = *dirty, *_i18 = dirty; ((scalar *)&s)->i >= 0; s = *++_i18)
      _attribute[s.i].dirty = true;
    pfree (dirty,__func__,__FILE__,__LINE__);
  }
}
#line 41 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"



static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    _attribute[sb.i].block = block;
    init_scalar (sb, bname);
    baseblock = list_append (baseblock, sb);
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    _attribute[sb.i].block = - n;
    init_scalar (sb, bname);
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  if (!(nvar + block <= _NVARMAX)) qassert ("/home/jiarongw/basilisk/src/grid/cartesian-common.h", 83, "nvar + block <= _NVARMAX");





  if (_attribute == NULL) {
    _attribute = (_Attributes *) pcalloc (nvar + block + 1, sizeof (_Attributes),__func__,__FILE__,__LINE__);
    _attribute[0].write = (double *) pcalloc (1, sizeof(double),__func__,__FILE__,__LINE__);
  }
  else
    _attribute = (_Attributes *)
      prealloc (_attribute - 1, (nvar + block + 1)*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  _attribute++;
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (new_scalar (name), name);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 122

    v.x = new_block_scalar (name, ext.x, block);
#line 122

    v.y = new_block_scalar (name, ext.y, block);}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 146

      vb.x.i = v.x.i + i;
#line 146

      vb.y.i = v.y.i + i;}
    init_vector (vb, NULL);
    {
#line 149

      _attribute[vb.x.i].block = - i;
#line 149

      _attribute[vb.y.i].block = - i;}
  }
  {
#line 152

    _attribute[v.x.i].block = block;
#line 152

    _attribute[v.y.i].block = block;}
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 162

      vb.x.i = v.x.i + i;
#line 162

      vb.y.i = v.y.i + i;}
    init_face_vector (vb, NULL);
    {
#line 165

      _attribute[vb.x.i].block = - i;
#line 165

      _attribute[vb.y.i].block = - i;}
  }
  {
#line 168

    _attribute[v.x.i].block = block;
#line 168

    _attribute[v.y.i].block = block;}
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 178
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 178
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 191
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 191
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;
#line 211 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 235

    init_const_scalar (v.x, name, *val++);
#line 235

    init_const_scalar (v.y, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 242

    v.x.i = _NVARMAX + i++;
#line 242

    v.y.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double * write = _attribute[a.i].write;
  int * read = _attribute[a.i].read;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  if (!(_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block)) qassert ("/home/jiarongw/basilisk/src/grid/cartesian-common.h", 256, "b.block > 0 && a.block == b.block");
  pfree (_attribute[a.i].depends,__func__,__FILE__,__LINE__);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].write = write;
  _attribute[a.i].read = read;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
  _attribute[a.i].depends = list_copy (_attribute[b.i].depends);
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  strongif (l) for (scalar s = *l, *_i19 = l; ((scalar *)&s)->i >= 0; s = *++_i19) {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  strongif (list) for (scalar s = *list, *_i20 = list; ((scalar *)&s)->i >= 0; s = *++_i20)
    {
#line 285

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 285

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  strongif (list) for (scalar f = *list, *_i21 = list; ((scalar *)&f)->i >= 0; f = *++_i21) {
    if (_attribute[f.i].block > 0) {
      pfree (_attribute[f.i].write,__func__,__FILE__,__LINE__); _attribute[f.i].write = NULL;
      pfree (_attribute[f.i].read,__func__,__FILE__,__LINE__); _attribute[f.i].read = NULL;
    }
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      _attribute[fb.i].read = NULL;
      _attribute[fb.i].write = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      pfree (_attribute[fb.i].depends,__func__,__FILE__,__LINE__); _attribute[fb.i].depends = NULL;
      _attribute[fb.i].freed = true;
    }
  }

  if (list == all) {
    all[0].i = -1;
    baseblock[0].i = -1;
    return;
  }

  trash (list);
  strongif (list) for (scalar f = *list, *_i22 = list; ((scalar *)&f)->i >= 0; f = *++_i22) {
    if (_attribute[f.i].block > 0) {
      scalar * s;
      for (s = all; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
      for (s = baseblock; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[1].i >= 0; s++)
   s[0] = s[1];
 s->i = -1;
      }
    }
  }
}

void free_solver ()
{
  if (!(_val_higher_dimension == 0.)) qassert ("/home/jiarongw/basilisk/src/grid/cartesian-common.h", 343, "_val_higher_dimension == 0.");

  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (baseblock,__func__,__FILE__,__LINE__); baseblock = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree ((_attribute - 1)->write,__func__,__FILE__,__LINE__);
  pfree (_attribute - 1,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_face) (vectorl);




void boundary_flux (vector * list) __attribute__ ((deprecated));

void boundary_flux (vector * list)
{
  vectorl list1 = {NULL};
  strongif (list) for (vector v = *list, *_i23 = list; ((scalar *)&v)->i >= 0; v = *++_i23)
    {
#line 395

      list1.x = list_append (list1.x, v.x);
#line 395

      list1.y = list_append (list1.y, v.y);}
  boundary_face (list1);
  {
#line 398

    pfree (list1.x,__func__,__FILE__,__LINE__);
#line 398

    pfree (list1.y,__func__,__FILE__,__LINE__);}
}

static scalar * list_add_depends (scalar * list, scalar s)
{
  strongif (list) for (scalar t = *list, *_i24 = list; ((scalar *)&t)->i >= 0; t = *++_i24)
    if (t.i == s.i)
      return list;
  scalar * list1 = list;
  strongif (_attribute[s.i].depends) for (scalar d = *_attribute[s.i].depends, *_i25 = _attribute[s.i].depends; ((scalar *)&d)->i >= 0; d = *++_i25)
    if (_attribute[d.i].dirty)
      list1 = list_add_depends (list1, d);
  return list_append (list1, s);
}


void boundary_internal (scalar * list, const char * fname, int line)
{ trace ("boundary_internal", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 416);
  if (list == NULL)
    { ; end_trace("boundary_internal", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 418);  return; }
  scalar * listc = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  strongif (list) for (scalar s = *list, *_i26 = list; ((scalar *)&s)->i >= 0; s = *++_i26)
    if (!is_constant(s) && _attribute[s.i].block > 0) {
      if (scalar_is_dirty (s)) {
 if (_attribute[s.i].face && _attribute[s.i].dirty != 2)
   {
#line 426

     if (_attribute[s.i].v.x.i == s.i)
       listf.x = list_add (listf.x, s), flux = true;
#line 426

     if (_attribute[s.i].v.y.i == s.i)
       listf.y = list_add (listf.y, s), flux = true;}
 if (!is_constant(cm) && _attribute[cm.i].dirty)
   listc = list_add_depends (listc, cm);
 if (_attribute[s.i].face != 2)
   listc = list_add_depends (listc, s);
      }




    }
  if (flux) {
    boundary_face (listf);
    {
#line 441

      pfree (listf.x,__func__,__FILE__,__LINE__);
#line 441

      pfree (listf.y,__func__,__FILE__,__LINE__);}
  }
  if (listc) {
    boundary_level (listc, -1);
    strongif (listc) for (scalar s = *listc, *_i27 = listc; ((scalar *)&s)->i >= 0; s = *++_i27)
      _attribute[s.i].dirty = false;
    pfree (listc,__func__,__FILE__,__LINE__);
  }
 end_trace("boundary_internal", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 450); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_face (vectorl list)
{
  {
#line 459

    strongif (list.x) for (scalar s = *list.x, *_i28 = list.x; ((scalar *)&s)->i >= 0; s = *++_i28)
      _attribute[s.i].dirty = 2;
#line 459

    strongif (list.y) for (scalar s = *list.y, *_i28 = list.y; ((scalar *)&s)->i >= 0; s = *++_i28)
      _attribute[s.i].dirty = 2;}
}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 465 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);

#if _call_symmetry
}
#define _IN_STENCIL 1

#line 464
static double _symmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 465 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return _stencil_val(__FILE__,__LINE__,s,0,0,0);

#undef _IN_STENCIL

#endif

#line 467
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 470 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);

#if _call_antisymmetry
}
#define _IN_STENCIL 1

#line 469
static double _antisymmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 470 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return -_stencil_val(__FILE__,__LINE__,s,0,0,0);

#undef _IN_STENCIL

#endif

#line 472
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

static double centered (double s0, double s1, double s2) {
  return (s2 - s0)/2.;
}

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  int block = _attribute[s.i].block;
  double * write = _attribute[s.i].write;
  int * read = _attribute[s.i].read;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[s.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[s.i].boundary_homogeneous;

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;
  if (block < 0) {
    scalar base = {s.i + block};
    _attribute[s.i].block = block;
    _attribute[s.i].write = _attribute[base.i].write;
    _attribute[s.i].read = _attribute[base.i].read;
  }
  else {
    _attribute[s.i].block = block > 0 ? block : 1;
    _attribute[s.i].write = write ? write : pmalloc((5*5)*sizeof(double),__func__,__FILE__,__LINE__);
    _attribute[s.i].read = read ? read : pmalloc((5*5)*sizeof(int),__func__,__FILE__,__LINE__);
  }

  _attribute[s.i].boundary = boundary ? boundary :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = boundary_homogeneous ? boundary_homogeneous :
    (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*2 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = centered;
  {
#line 523
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 523
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  {
#line 533

    _attribute[s.i].d.x = -1;
#line 533

    _attribute[s.i].d.y = -1;}
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 549
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 549
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*2 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 569
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 569
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 581
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 581
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }}






    for (int b = 0; b < nboundary; b++) {
      _attribute[t.x.x.i].boundary[b] = _attribute[t.y.x.i].boundary[b] =
 _attribute[t.x.x.i].boundary_homogeneous[b] = _attribute[t.y.y.i].boundary_homogeneous[b] =
 b < 2*2 ? default_scalar_bc[b] : symmetry;
      _attribute[t.x.y.i].boundary[b] = _attribute[t.y.y.i].boundary[b] =
 _attribute[t.x.y.i].boundary_homogeneous[b] = _attribute[t.y.x.i].boundary_homogeneous[b] =
 b < 2*2 ? default_vector_bc[b] : antisymmetry;
    }



  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/grid/cartesian-common.h", .line = 619,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 619 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
 {
    bool inside = true;
    coord o = {x,y,z};
    {
#line 622

      IF (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
#line 622

      IF (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;}
    IF (inside) {
      Delta /= 2.;



      _stencil_fprintf (__FILE__,__LINE__,p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 651 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
    }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 652
foreach(){

#line 619 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
 {
    bool inside = true;
    coord o = {x,y,z};
    {
#line 622

      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
#line 622

      if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;}
    if (inside) {
      Delta /= 2.;



      fprintf (p.fp, "%g %g\n%g %g\n%g %g\n%g %g\n%g %g\n\n",
        x - Delta, y - Delta,
        x - Delta, y + Delta,
        x + Delta, y + Delta,
        x + Delta, y - Delta,
        x - Delta, y - Delta);
#line 651 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
    }
  } } end_foreach(); }
  fflush (p.fp);
}
#line 663 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"




    "  plot '%s' w l lc 0, "
    "'%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 1 title columnhead(3+3*v)",





    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  strongif (all) for (scalar v = *all, *_i29 = all; ((scalar *)&v)->i >= 0; v = *++_i29)



    fprintf (fp, "x y %s ", _attribute[v.i].name);



  fputc ('\n', fp);
#line 730 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 strongif (all) for (scalar v = *all, *_i30 = all; ((scalar *)&v)->i >= 0; v = *++_i30) {
   fprintf (fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   if (allocated(k,l,0))
     fprintf (fp, "%g ", val(v,k,l,0));
   else
     fputs ("n/a ", fp);
 }
 fputc ('\n', fp);
      }
#line 760 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  strongif (all) for (scalar s = *all, *_i31 = all; ((scalar *)&s)->i >= 0; s = *++_i31) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);

#if _call_cartesian_debug
}
#define _IN_STENCIL 1

#line 696
static void _cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 697 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  IF (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  IF (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  strongif (all) for (scalar v = *all, *_i29 = all; ((scalar *)&v)->i >= 0; v = *++_i29)



    _stencil_fprintf (__FILE__,__LINE__,fp, "x y %s ", _attribute[v.i].name);



  _stencil_fputc (__FILE__,__LINE__,'\n', fp);
#line 730 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++) {
 strongif (all) for (scalar v = *all, *_i30 = all; ((scalar *)&v)->i >= 0; v = *++_i30) {
   _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g ",
     x + k*Delta + _attribute[v.i].d.x*Delta/2.,
     y + l*Delta + _attribute[v.i].d.y*Delta/2.);
   IF (allocated(k,l,0))
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g ", _stencil_val(__FILE__,__LINE__,v,k,l,0));
   
     _stencil_fputs (__FILE__,__LINE__,"n/a ", fp);
 }
 _stencil_fputc (__FILE__,__LINE__,'\n', fp);
      }
#line 760 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"
  fclose (fp);

  fp = fopen ("debug.plot", "w");
  _stencil_fprintf (__FILE__,__LINE__,fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  strongif (all) for (scalar s = *all, *_i31 = all; ((scalar *)&s)->i >= 0; s = *++_i31) {
    char * name = replace_ (_attribute[s.i].name);
    _stencil_fprintf (__FILE__,__LINE__,fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  _stencil_fprintf (__FILE__,__LINE__,ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);

#undef _IN_STENCIL

#endif

#line 781
}

void cartesian_methods ()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_face = cartesian_boundary_face;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 801 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
   (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y);
#line 829 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

#if _call_interpolate_linear
}
#define _IN_STENCIL 1

#line 800
static double _interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 801 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;







  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  int i = sign(x), j = sign(y);
  x = fabs(x); y = fabs(y);

  return ((_stencil_val(__FILE__,__LINE__,v,0,0,0)*(1. - x) + _stencil_val(__FILE__,__LINE__,v,i,0,0)*x)*(1. - y) +
   (_stencil_val(__FILE__,__LINE__,v,0,j,0)*(1. - x) + _stencil_val(__FILE__,__LINE__,v,i,j,0)*x)*y);
#line 829 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

#undef _IN_STENCIL

#endif

#line 829
}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 833);
  scalar v = p.v;
  boundary_internal ((scalar *)(((scalar []){v,{-1}})), "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 835);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 836 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 838);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 839);  return _ret; }
 end_trace("interpolate", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 840); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 844);
  boundary_internal ((scalar *)(list), "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 845);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 848 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      strongif (list) for (scalar s = *list, *_i32 = list; ((scalar *)&s)->i >= 0; s = *++_i32)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      strongif (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
 v[j++] = nodata;
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/jiarongw/basilisk/src/grid/cartesian-common.h", 867); }



typedef int bid;

bid new_bid ()
{
  int b = nboundary++;
  strongif (all) for (scalar s = *all, *_i34 = all; ((scalar *)&s)->i >= 0; s = *++_i34) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  strongif (all) for (scalar s = *all, *_i35 = all; ((scalar *)&s)->i >= 0; s = *++_i35) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 887

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 887

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 899 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return nodata;

#if _call_periodic_bc
}
#define _IN_STENCIL 1

#line 898
static double _periodic_bc (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 899 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return nodata;

#undef _IN_STENCIL

#endif

#line 901
}

static void periodic_boundary (int d)
{

  strongif (all) for (scalar s = *all, *_i36 = all; ((scalar *)&s)->i >= 0; s = *++_i36)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  strongif (all) for (scalar s = *all, *_i37 = all; ((scalar *)&s)->i >= 0; s = *++_i37)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{



    if (!(dir <= bottom)) qassert ("/home/jiarongw/basilisk/src/grid/cartesian-common.h", 924, "dir <= bottom");




  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 937 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return val(s,i,j,k);

#if _call_getvalue
}
#define _IN_STENCIL 1

#line 936
static double _getvalue (Point point, scalar s, int i, int j, int k)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 937 "/home/jiarongw/basilisk/src/grid/cartesian-common.h"

  return _stencil_val(__FILE__,__LINE__,s,i,j,k);

#undef _IN_STENCIL

#endif

#line 939
}
#line 4 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 27 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2);

#if _call_restriction_average
}
#define _IN_STENCIL 1

#line 26
static void _restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 27 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += _stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2);

#undef _IN_STENCIL

#endif

#line 32
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
#if _call_restriction_volume_average
}
#define _IN_STENCIL 1

#line 34
static void _restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 35 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*_stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*_stencil_val(__FILE__,__LINE__,s,0,0,0); end_foreach_child(); }
  _stencil_val(__FILE__,__LINE__,s,0,0,0) = sum/(1 << 2)/(val_cm(cm,0,0,0) + 1e-30);
 }
#undef _IN_STENCIL

#endif

#line 40
}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {




      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0))/2.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0))/2.;






  }
#line 44
 {




      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,1,0,0))/2.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,1,2,0))/2.;






  }}

#if _call_face_average
}
#define _IN_STENCIL 1

#line 42
static void _face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 43 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {




      _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = (_stencil_fine(__FILE__,__LINE__,v.x,0,0,0) + _stencil_fine(__FILE__,__LINE__,v.x,0,1,0))/2.;
      _stencil_val(__FILE__,__LINE__,v.x,1,0,0) = (_stencil_fine(__FILE__,__LINE__,v.x,2,0,0) + _stencil_fine(__FILE__,__LINE__,v.x,2,1,0))/2.;






  }
#line 44
 {




      _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = (_stencil_fine(__FILE__,__LINE__,v.y,0,0,0) + _stencil_fine(__FILE__,__LINE__,v.y,1,0,0))/2.;
      _stencil_val(__FILE__,__LINE__,v.y,0,1,0) = (_stencil_fine(__FILE__,__LINE__,v.y,0,2,0) + _stencil_fine(__FILE__,__LINE__,v.y,1,2,0))/2.;






  }}

#undef _IN_STENCIL

#endif

#line 58
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 61 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);

#if _call_restriction_face
}
#define _IN_STENCIL 1
#define face_average _face_average

#line 60
static void _restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 61 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);

#undef face_average
#undef _IN_STENCIL

#endif

#line 63
}

static inline void restriction_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);





  }

#if _call_restriction_vertex
}
#define _IN_STENCIL 1

#line 65
static void _restriction_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 66 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  for (int i = 0; i <= 1; i++) {
    _stencil_val(__FILE__,__LINE__,s,i,0,0) = _stencil_fine(__FILE__,__LINE__,s,2*i,0,0);

    _stencil_val(__FILE__,__LINE__,s,i,1,0) = _stencil_fine(__FILE__,__LINE__,s,2*i,2,0);





  }

#undef _IN_STENCIL

#endif

#line 77
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#if _call_no_restriction
}
#define _IN_STENCIL 1

#line 79
static void _no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 78 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#undef _IN_STENCIL

#endif

#line 79
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 81 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }

#if _call_no_data
}
#define _IN_STENCIL 1

#line 81
static void _no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 81 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = nodata; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 84
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level (l){

#line 90 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
       { foreach_child()
        val(w,0,0,0) = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
       { foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){w,{-1}}), l + 1);
  }

   { foreach_level(0){

#line 104 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

    val(w,0,0,0) = val(s,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
   { foreach_level(0){

#line 111 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

    val(s,0,0,0) = val(w,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){s,{-1}}), 0);
  for (int l = 0; l <= depth() - 1; l++) {
     { foreach_coarse_level (l){

#line 115 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
      _attribute[s.i].prolongation (point, s);
       { foreach_child()
        val(s,0,0,0) += val(w,0,0,0); end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 125 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




    return (9.*coarse(s,0,0,0) +
     3.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0)) +
     coarse(s,child.x,child.y,0))/16.;
#line 140 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#if _call_bilinear
}
#define _IN_STENCIL 1

#line 124
static double _bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 125 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




    return (9.*_stencil_coarse(__FILE__,__LINE__,s,0,0,0) +
     3.*(_stencil_coarse(__FILE__,__LINE__,s,child.x,0,0) + _stencil_coarse(__FILE__,__LINE__,s,0,child.y,0)) +
     _stencil_coarse(__FILE__,__LINE__,s,child.x,child.y,0))/16.;
#line 140 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

#undef _IN_STENCIL

#endif

#line 140
}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 143 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }

#if _call_refine_bilinear
}
#define _IN_STENCIL 1
#define bilinear _bilinear

#line 142
static void _refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 143 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = bilinear (point, s); end_foreach_child(); }

#undef bilinear
#undef _IN_STENCIL

#endif

#line 146
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 154 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return
    quadratic (quadratic (coarse(s,0,0,0),
     coarse(s,child.x,0,0),
     coarse(s,-child.x,0,0)),
        quadratic (coarse(s,0,child.y,0),
     coarse(s,child.x,child.y,0),
     coarse(s,-child.x,child.y,0)),
        quadratic (coarse(s,0,-child.y,0),
     coarse(s,child.x,-child.y,0),
     coarse(s,-child.x,-child.y,0)));





#if _call_biquadratic
}
#define _IN_STENCIL 1

#line 153
static double _biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 154 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return
    quadratic (quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,0,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,0,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,0,0)),
        quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,child.y,0)),
        quadratic (_stencil_coarse(__FILE__,__LINE__,s,0,-child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,child.x,-child.y,0),
     _stencil_coarse(__FILE__,__LINE__,s,-child.x,-child.y,0)));





#undef _IN_STENCIL

#endif

#line 172
}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 175 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return (36.*val(s,0,0,0) + 18.*(val(s,-1,0,0) + val(s,0,-1,0)) - 6.*(val(s,1,0,0) + val(s,0,1,0)) +
   9.*val(s,-1,-1,0) - 3.*(val(s,1,-1,0) + val(s,-1,1,0)) + val(s,1,1,0))/64.;





#if _call_biquadratic_vertex
}
#define _IN_STENCIL 1

#line 174
static double _biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 175 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"




  return (36.*_stencil_val(__FILE__,__LINE__,s,0,0,0) + 18.*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) + _stencil_val(__FILE__,__LINE__,s,0,-1,0)) - 6.*(_stencil_val(__FILE__,__LINE__,s,1,0,0) + _stencil_val(__FILE__,__LINE__,s,0,1,0)) +
   9.*_stencil_val(__FILE__,__LINE__,s,-1,-1,0) - 3.*(_stencil_val(__FILE__,__LINE__,s,1,-1,0) + _stencil_val(__FILE__,__LINE__,s,-1,1,0)) + _stencil_val(__FILE__,__LINE__,s,1,1,0))/64.;





#undef _IN_STENCIL

#endif

#line 185
}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 188 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }

#if _call_refine_biquadratic
}
#define _IN_STENCIL 1
#define biquadratic _biquadratic

#line 187
static void _refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 188 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = biquadratic (point, s); end_foreach_child(); }

#undef biquadratic
#undef _IN_STENCIL

#endif

#line 191
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 194 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  if (!(fabs(sum) < 1e-10)) qassert ("/home/jiarongw/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  if (!(fabs(sum) < 1e-10)) qassert ("/home/jiarongw/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
#if _call_refine_linear
}
#define _IN_STENCIL 1

#line 193
static void _refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 194 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 194

  coord g;
  IF (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0));}
  
    {
#line 200

      g.x = (_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,-1,0,0))/2.;
#line 200

      g.y = (_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,-1,0))/2.;}

  double sc = _stencil_val(__FILE__,__LINE__,s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
    {
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  IF (!(fabs(sum) < 1e-10)) _stencil_qassert (__FILE__,__LINE__,"/home/jiarongw/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 194

  coord g;
  IF (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0));}
  
    {
#line 200

      g.x = (_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,-1,0,0))/2.;
#line 200

      g.y = (_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,-1,0))/2.;}

  double sc = _stencil_val(__FILE__,__LINE__,s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 2);
   { foreach_child() {
    _stencil_val(__FILE__,__LINE__,s,0,0,0) = sc;
    {
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      _stencil_val(__FILE__,__LINE__,s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  IF (!(fabs(sum) < 1e-10)) _stencil_qassert (__FILE__,__LINE__,"/home/jiarongw/basilisk/src/grid/multigrid-common.h", 210, "fabs(sum) < 1e-10");
 }
#undef _IN_STENCIL

#endif

#line 211
}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 214 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }

#if _call_refine_reset
}
#define _IN_STENCIL 1

#line 213
static void _refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 214 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    _stencil_val(__FILE__,__LINE__,v,0,0,0) = 0.; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 217
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 220 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }

#if _call_refine_injection
}
#define _IN_STENCIL 1

#line 219
static void _refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 220 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  double val = _stencil_val(__FILE__,__LINE__,v,0,0,0);
   { foreach_child()
    _stencil_val(__FILE__,__LINE__,v,0,0,0) = val; end_foreach_child(); }

#undef _IN_STENCIL

#endif

#line 224
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 244

    _attribute[v.y.i].restriction = no_restriction;
#line 244

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 251 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   strongif (all) for (scalar v = *all, *_i38 = all; ((scalar *)&v)->i >= 0; v = *++_i38)
     fprintf (fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       coarse(v,k*child.x,l*child.y,0));
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   strongif (all) for (scalar v = *all, *_i39 = all; ((scalar *)&v)->i >= 0; v = *++_i39) {
     fprintf (fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     if (allocated_child(k,l,0))
       fprintf (fp, "%g ", fine(v,k,l,0));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }
      fprintf (ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      fprintf (plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);

#if _call_multigrid_debug
}
#define _IN_STENCIL 1
#define cartesian_debug _cartesian_debug

#line 250
static void _multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 251 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  IF (point.level > 0) {
    char name[80] = "coarse";
    IF (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 271 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++) {
   strongif (all) for (scalar v = *all, *_i38 = all; ((scalar *)&v)->i >= 0; v = *++_i38)
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g %g ",
       xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
       yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
       _stencil_coarse(__FILE__,__LINE__,v,k*child.x,l*child.y,0));
   _stencil_fputc (__FILE__,__LINE__,'\n', fp);
 }
      _stencil_fprintf (__FILE__,__LINE__,ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
      _stencil_fprintf (__FILE__,__LINE__,plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 3 t ''", name);
#line 302 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }

  IF (is_coarse()) {
    char name[80] = "fine";
    IF (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 324 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++) {
   strongif (all) for (scalar v = *all, *_i39 = all; ((scalar *)&v)->i >= 0; v = *++_i39) {
     _stencil_fprintf (__FILE__,__LINE__,fp, "%g %g ",
       xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
       yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.);
     IF (allocated_child(k,l,0))
       _stencil_fprintf (__FILE__,__LINE__,fp, "%g ", _stencil_fine(__FILE__,__LINE__,v,k,l,0));
     
       _stencil_fputs (__FILE__,__LINE__,"n/a ", fp);
   }
   _stencil_fputc (__FILE__,__LINE__,'\n', fp);
 }
      _stencil_fprintf (__FILE__,__LINE__,ferr, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
      _stencil_fprintf (__FILE__,__LINE__,plot, ", '%s' u 1+3*v:2+3*v:3+3*v w labels tc lt 2 t ''", name);
#line 362 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);

#undef cartesian_debug
#undef _IN_STENCIL

#endif

#line 366
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  strongif (list) for (scalar s = *list, *_i40 = list; ((scalar *)&s)->i >= 0; s = *++_i40)
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 380

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 380

     list2 = list_add (list2, _attribute[s.i].v.y);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 389 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
 strongif (listdef) for (scalar s = *listdef, *_i41 = listdef; ((scalar *)&s)->i >= 0; s = *++_i41)
    { foreach_block_inner()
     restriction_average (point, s); end_foreach_block_inner(); }
 strongif (listc) for (scalar s = *listc, *_i42 = listc; ((scalar *)&s)->i >= 0; s = *++_i42) {
    { foreach_block_inner()
     _attribute[s.i].restriction (point, s); end_foreach_block_inner(); }
 }
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods ()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/grid/multigrid-common.h", .line = 428,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 428 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

    _stencil_val(__FILE__,__LINE__,size,0,0,0) = 1; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 429
foreach(){

#line 428 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 437 "/home/jiarongw/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 826 "/home/jiarongw/basilisk/src/grid/multigrid.h"

struct Dimensions {
  int nx, ny, nz;
};

void dimensions (struct Dimensions p)
{

  for (int i = 0; i < 2; i++)
    mpi_dims[i] = (&p.nx)[i];

}



#if 2 == 1

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\

#line 847

#define end_foreach_slice_x() }

#elif 2 == 2

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\

#line 857

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\

#line 865

#define end_foreach_slice_y() }

#elif 2 == 3

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 876

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 885

#define end_foreach_slice_y() }

#define foreach_slice_z(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = start; point.k < end; point.k++)\

#line 894

#define end_foreach_slice_z() }

#endif

#line 1 "grid/multigrid-mpi.h"
#line 1 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"


typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;

#define BUF void *


#line 10

static BUF snd_x (int i, int dst, int tag, int level, scalar * list,
    MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  strongif (list) for (scalar s = *list, *_i43 = list; ((scalar *)&s)->i >= 0; s = *++_i43)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 2 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   { foreach_slice_x (i, i + 2, level){

#line 21 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"

    strongif (list) for (scalar s = *list, *_i44 = list; ((scalar *)&s)->i >= 0; s = *++_i44) {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_x(); }
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}
#line 10

static BUF snd_y (int i, int dst, int tag, int level, scalar * list,
    MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  strongif (list) for (scalar s = *list, *_i43 = list; ((scalar *)&s)->i >= 0; s = *++_i43)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 2 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   { foreach_slice_y (i, i + 2, level){

#line 21 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"

    strongif (list) for (scalar s = *list, *_i44 = list; ((scalar *)&s)->i >= 0; s = *++_i44) {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_y(); }
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}


#line 30

static void rcv_x (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  strongif (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 2 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   { foreach_slice_x (i, i + 2, level){

#line 42 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"

    strongif (list) for (scalar s = *list, *_i46 = list; ((scalar *)&s)->i >= 0; s = *++_i46) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_x(); }
  pfree (buf,__func__,__FILE__,__LINE__);
}
#line 30

static void rcv_y (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  strongif (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 2 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   { foreach_slice_y (i, i + 2, level){

#line 42 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"

    strongif (list) for (scalar s = *list, *_i46 = list; ((scalar *)&s)->i >= 0; s = *++_i46) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_y(); }
  pfree (buf,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{ trace ("mpi_boundary_level", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 52);
  scalar * list1 = NULL;
  strongif (list) for (scalar s = *list, *_i47 = list; ((scalar *)&s)->i >= 0; s = *++_i47)
    if (!is_constant(s) && _attribute[s.i].block > 0)
      list1 = list_add (list1, s);
  if (!list1)
    { ; end_trace("mpi_boundary_level", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 58);  return; }

  prof_start ("mpi_boundary_level");

  if (level < 0) level = depth();
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
  {
#line 65
 {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_x (npl - 2*2, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list1);
    rcv_x (npl - 2, right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }
#line 65
 {
    int bottom, top;
    MPI_Cart_shift (mpi->cartcomm, dir.y, 1, &bottom, &top);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_y (npl - 2*2, top, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_y (2, bottom, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_y (0, bottom, 0, level, list1);
    rcv_y (npl - 2, top, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }}

  pfree (list1,__func__,__FILE__,__LINE__);

  prof_stop();
 end_trace("mpi_boundary_level", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 85); }

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  pfree (m,__func__,__FILE__,__LINE__);
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  MPI_Dims_create (npe(), 2, mpi_dims);
  MPI_Cart_create (MPI_COMM_WORLD, 2,
     mpi_dims, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), 2, mpi_coords);


  struct { int x, y, z; } dir = {0,1,2};
  {
#line 104
 {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  }
#line 104
 {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.y, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (bottom);
    if (r != MPI_PROC_NULL)
      periodic_boundary (top);
  }}


  N /= mpi_dims[0];
  int r = 0;
  while (N > 1)
    N /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = mpi_dims[0]*(1 << r);
  grid->n = 1 << 2*depth();
  grid->tn = npe()*grid->n;


  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}


double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 134);
  long i;
  if (leaves)
    i = pid()*(1 << 2*depth());
  else
    i = pid()*((1 << 2*(depth() + 1)) - 1)/((1 << 2) - 1);
   { foreach_cell(){

#line 140 "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h"
 {
    if (!leaves || is_leaf(cell))
      val(index,0,0,0) = i++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary_internal ((scalar *)(((scalar []){index,{-1}})), "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 146);
  { double _ret =  pid() == 0 ? i*npe() - 1 : -1; end_trace("z_indexing", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 147);  return _ret; }
 end_trace("z_indexing", "/home/jiarongw/basilisk/src/grid/multigrid-mpi.h", 148); }
#line 900 "/home/jiarongw/basilisk/src/grid/multigrid.h"
#line 6 "field_new_PM_restore.c"
#line 1 "view.h"
#line 1 "/home/jiarongw/basilisk/src/view.h"
#line 67 "/home/jiarongw/basilisk/src/view.h"
#include <gl/framebuffer.h>
#include <gl/trackball.h>
#include <gl/utils.h>


#line 1 "utils.h"
#line 1 "/home/jiarongw/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf () {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _n = n;
{ double n = _n; NOT_UNUSED(n);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 69,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 69 "/home/jiarongw/basilisk/src/utils.h"
 n++; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 69

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:n)) {

#line 69
foreach(){

#line 69 "/home/jiarongw/basilisk/src/utils.h"
 n++; } end_foreach();mpi_all_reduce_array (&n, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Multigrid"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _max = max;
 double _avg = avg;
 double _rms = rms;
 double _volume = volume;
{ double max = _max; NOT_UNUSED(max);
 double avg = _avg; NOT_UNUSED(avg);
 double rms = _rms; NOT_UNUSED(rms);
 double volume = _volume; NOT_UNUSED(volume);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 135,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 135
foreach_stencil(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach_stencil(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 143

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)  reduction(+:avg) 
   reduction(+:rms)  reduction(+:volume)) {

#line 135

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/jiarongw/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > max) max = v;
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      avg += (sq(Delta)*val_cm(cm,0,0,0))*v;
      rms += (sq(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }mpi_all_reduce_array (&max, double, MPI_MAX, 1);
mpi_all_reduce_array (&avg, double, MPI_SUM, 1);
mpi_all_reduce_array (&rms, double, MPI_SUM, 1);
mpi_all_reduce_array (&volume, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _sum = sum;
 double _sum2 = sum2;
 double _volume = volume;
 double _max = max;
 double _min = min;
{ double sum = _sum; NOT_UNUSED(sum);
 double sum2 = _sum2; NOT_UNUSED(sum2);
 double volume = _volume; NOT_UNUSED(volume);
 double max = _max; NOT_UNUSED(max);
 double min = _min; NOT_UNUSED(min);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 163,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 163
foreach_stencil(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && _stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) > max) max = _stencil_val(__FILE__,__LINE__,f,0,0,0);
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) < min) min = _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach_stencil(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && _stencil_val(__FILE__,__LINE__,f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(_stencil_val(__FILE__,__LINE__,f,0,0,0));
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) > max) max = _stencil_val(__FILE__,__LINE__,f,0,0,0);
      IF (_stencil_val(__FILE__,__LINE__,f,0,0,0) < min) min = _stencil_val(__FILE__,__LINE__,f,0,0,0);
    } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 171

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)  reduction(+:sum2)  reduction(+:volume) 
   reduction(max:max)  reduction(min:min)) {

#line 163

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/jiarongw/basilisk/src/utils.h"

    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      volume += (sq(Delta)*val_cm(cm,0,0,0));
      sum += (sq(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (sq(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > max) max = val(f,0,0,0);
      if (val(f,0,0,0) < min) min = val(f,0,0,0);
    } } end_foreach(); }mpi_all_reduce_array (&sum, double, MPI_SUM, 1);
mpi_all_reduce_array (&sum2, double, MPI_SUM, 1);
mpi_all_reduce_array (&volume, double, MPI_SUM, 1);
mpi_all_reduce_array (&max, double, MPI_MAX, 1);
mpi_all_reduce_array (&min, double, MPI_MIN, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/jiarongw/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/jiarongw/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/jiarongw/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  if (!(list_len(f) == vectors_len(g))) qassert ("/home/jiarongw/basilisk/src/utils.h", 239, "list_len(f) == vectors_len(g)");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 240,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 240 "/home/jiarongw/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i2 = f; vector * _i3 = g; strongif (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i2, v = *++_i3)
      {
#line 243
 {





   _stencil_val(__FILE__,__LINE__,v.x,0,0,0) = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,0,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,1,0,0))/Delta;
      }
#line 243
 {





   _stencil_val(__FILE__,__LINE__,v.y,0,0,0) = _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,-1,0), _stencil_val(__FILE__,__LINE__,s,0,0,0), _stencil_val(__FILE__,__LINE__,s,0,1,0))/Delta;
      }}
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 251
foreach(){

#line 240 "/home/jiarongw/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i2 = f; vector * _i3 = g; strongif (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i2, v = *++_i3)
      {
#line 243
 {





   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
      }
#line 243
 {





   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
      }}
  } } end_foreach(); }
}
#line 269 "/home/jiarongw/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 271,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 271
foreach_stencil(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 271
foreach_stencil(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach_stencil(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach_stencil(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    _stencil_val(__FILE__,__LINE__,omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*_stencil_val(__FILE__,__LINE__,u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*_stencil_val(__FILE__,__LINE__,u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*_stencil_val(__FILE__,__LINE__,u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*_stencil_val(__FILE__,__LINE__,u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*_stencil_val(__FILE__,__LINE__,u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*_stencil_val(__FILE__,__LINE__,u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 275

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 271
foreach(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 271
foreach(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 271
foreach(){

#line 271 "/home/jiarongw/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); } }
}





double change (scalar s, scalar sn)
{
  double max = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _max = max;
{ double max = _max; NOT_UNUSED(max);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/utils.h", .line = 285,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 285
foreach_stencil(){

#line 285 "/home/jiarongw/basilisk/src/utils.h"
 {
    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (_stencil_val(__FILE__,__LINE__,s,0,0,0) - _stencil_val(__FILE__,__LINE__,sn,0,0,0));
      IF (ds > max)
 max = ds;
    }
    _stencil_val(__FILE__,__LINE__,sn,0,0,0) = _stencil_val(__FILE__,__LINE__,s,0,0,0);
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 285
foreach_stencil(){

#line 285 "/home/jiarongw/basilisk/src/utils.h"
 {
    IF ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (_stencil_val(__FILE__,__LINE__,s,0,0,0) - _stencil_val(__FILE__,__LINE__,sn,0,0,0));
      IF (ds > max)
 max = ds;
    }
    _stencil_val(__FILE__,__LINE__,sn,0,0,0) = _stencil_val(__FILE__,__LINE__,s,0,0,0);
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 292

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:max)) {

#line 285

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 285
foreach(){

#line 285 "/home/jiarongw/basilisk/src/utils.h"
 {
    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 285
foreach(){

#line 285 "/home/jiarongw/basilisk/src/utils.h"
 {
    if ((sq(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > max)
 max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }mpi_all_reduce_array (&max, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 292
 }
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    strongif (all) for (scalar s = *all, *_i48 = all; ((scalar *)&s)->i >= 0; s = *++_i48)
      if (!strcmp (_attribute[s.i].name, name))
 return s;
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    strongif (all) for (scalar s = *all, *_i49 = all; ((scalar *)&s)->i >= 0; s = *++_i49)
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;
  }
  return (vector){{-1}};
}







#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  if (!(norm > 0.)) qassert ("/home/jiarongw/basilisk/src/utils.h", 331, "norm > 0.");\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
      {\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }}\
      if (_n == 2) {\

#line 354

#define end_foreach_segment() } } end_foreach(); }




void fields_stats ()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  strongif (all) for (scalar s = *all, *_i50 = all; ((scalar *)&s)->i >= 0; s = *++_i50)
    fprintf (ferr, " %s", _attribute[s.i].name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  strongif (all) for (scalar s = *all, *_i51 = all; ((scalar *)&s)->i >= 0; s = *++_i51) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

#line 1 "output.h"
#line 1 "/home/jiarongw/basilisk/src/output.h"
#line 37 "/home/jiarongw/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/jiarongw/basilisk/src/output.h", 47);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  boundary_internal ((scalar *)(p.list), "/home/jiarongw/basilisk/src/output.h", 58);
  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i52 = p.list; ((scalar *)&s)->i >= 0; s = *++_i52)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});
      }
      else {
 Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 73 "/home/jiarongw/basilisk/src/output.h"

 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i53 = p.list; ((scalar *)&s)->i >= 0; s = *++_i53)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    strongif (p.list) for (scalar s = *p.list, *_i54 = p.list; ((scalar *)&s)->i >= 0; s = *++_i54)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 strongif (p.list) for (scalar s = *p.list, *_i55 = p.list; ((scalar *)&s)->i >= 0; s = *++_i55)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/jiarongw/basilisk/src/output.h", 113); }
#line 141 "/home/jiarongw/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/jiarongw/basilisk/src/output.h", 150);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.linear) {
    scalar f = p.f;
    boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/jiarongw/basilisk/src/output.h", 155);
  }
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 172 "/home/jiarongw/basilisk/src/output.h"

 if (!(point.level >= 0)) qassert ("/home/jiarongw/basilisk/src/output.h", 173, "point.level >= 0");
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/jiarongw/basilisk/src/output.h", 180); }
#line 189 "/home/jiarongw/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  if (!(i < 127 - 1)) qassert ("/home/jiarongw/basilisk/src/output.h", 321, "i < NCMAP - 1");
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 340 "/home/jiarongw/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup ()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  if (!(pid() == 0)) qassert ("/home/jiarongw/basilisk/src/output.h", 422, "pid() == 0");
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  if (!(pid() == 0)) qassert ("/home/jiarongw/basilisk/src/output.h", 496, "pid() == 0");
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 571 "/home/jiarongw/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/jiarongw/basilisk/src/output.h", 586);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)(((scalar []){f,mask,{-1}})), "/home/jiarongw/basilisk/src/output.h", 608);
    else
      boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/jiarongw/basilisk/src/output.h", 610);
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 636 "/home/jiarongw/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 646 "/home/jiarongw/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/jiarongw/basilisk/src/output.h", 678); }
#line 710 "/home/jiarongw/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/jiarongw/basilisk/src/output.h", 721);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }
  if (p.linear) {
    scalar f = p.f, mask = p.mask;
    if (mask.i)
      boundary_internal ((scalar *)(((scalar []){f,mask,{-1}})), "/home/jiarongw/basilisk/src/output.h", 733);
    else
      boundary_internal ((scalar *)(((scalar []){f,{-1}})), "/home/jiarongw/basilisk/src/output.h", 735);
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 764 "/home/jiarongw/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 774 "/home/jiarongw/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/jiarongw/basilisk/src/output.h", 786); }
#line 813 "/home/jiarongw/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/jiarongw/basilisk/src/output.h", 843);
  char * fname = p.file;

#if 1

  not_mpi_compatible();

  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);




  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  strongif (list) for (scalar s = *list, *_i56 = list; ((scalar *)&s)->i >= 0; s = *++_i56)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 921 "/home/jiarongw/basilisk/src/output.h"
 {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :



      child.x == -1 && child.y == -1 ? 0 :
 child.x == -1 && child.y == 1 ? 1 :
 child.x == 1 && child.y == -1 ? 2 :
 3;
#line 951 "/home/jiarongw/basilisk/src/output.h"
      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      strongif (list) for (scalar s = *list, *_i57 = list; ((scalar *)&s)->i >= 0; s = *++_i57)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != nodata ? - val(s,0,0,0) : (double) DBL_MAX;
     }





   }
   else
     a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if 1
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/jiarongw/basilisk/src/output.h", 1019); }
#line 1043 "/home/jiarongw/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  strongif (lista) for (scalar s = *lista, *_i58 = lista; ((scalar *)&s)->i >= 0; s = *++_i58)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  strongif (list) for (scalar s = *list, *_i59 = list; ((scalar *)&s)->i >= 0; s = *++_i59) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1

void dump (struct Dump p)
{ trace ("dump", "/home/jiarongw/basilisk/src/output.h", 1097);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  if (!(fp)) qassert ("/home/jiarongw/basilisk/src/output.h", 1112, "fp");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 1123 "/home/jiarongw/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    strongif (list) for (scalar s = *list, *_i60 = list; ((scalar *)&s)->i >= 0; s = *++_i60)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/jiarongw/basilisk/src/output.h", 1145); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/jiarongw/basilisk/src/output.h", 1149);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };


  for (int i = 0; i < 2; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);


  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  strongif (list) for (scalar s = *list, *_i61 = list; ((scalar *)&s)->i >= 0; s = *++_i61)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

   { foreach_cell(){

#line 1195 "/home/jiarongw/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      strongif (list) for (scalar s = *list, *_i62 = list; ((scalar *)&s)->i >= 0; s = *++_i62)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/jiarongw/basilisk/src/output.h", 1219); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/jiarongw/basilisk/src/output.h", 1224);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1228);  return _ret; }
  if (!(fp)) qassert ("/home/jiarongw/basilisk/src/output.h", 1229, "fp");

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }
#line 1246 "/home/jiarongw/basilisk/src/output.h"
  if (header.npe != npe()) {
    fprintf (ferr,
      "restore(): error: the number of processes don't match:"
      " %d != %d\n",
      header.npe, npe());
    exit (1);
  }
  dimensions ((struct Dimensions){header.n.x, header.n.y, header.n.z});
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);





  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 strongif (list) for (scalar s = *list, *_i63 = list; ((scalar *)&s)->i >= 0; s = *++_i63)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }


  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << 2*(header.depth + 1)) - 1)/
    ((1 << 2) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }


  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y},{{-1},{-1}}});



   { foreach_cell(){

#line 1343 "/home/jiarongw/basilisk/src/output.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    strongif (list) for (scalar s = *list, *_i64 = list; ((scalar *)&s)->i >= 0; s = *++_i64) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  strongif (all) for (scalar s = *all, *_i65 = all; ((scalar *)&s)->i >= 0; s = *++_i65)
    _attribute[s.i].dirty = true;


  scalar * other = NULL;
  strongif (all) for (scalar s = *all, *_i66 = all; ((scalar *)&s)->i >= 0; s = *++_i66)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1389);  return _ret; }
 end_trace("restore", "/home/jiarongw/basilisk/src/output.h", 1390); }
#line 376 "/home/jiarongw/basilisk/src/utils.h"
#line 73 "/home/jiarongw/basilisk/src/view.h"
#line 1 "input.h"
#line 1 "/home/jiarongw/basilisk/src/input.h"
#line 16 "/home/jiarongw/basilisk/src/input.h"
struct InputPGM {

  scalar s;
  FILE * fp;

  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (ferr, "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((width*height)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  int _width = width;
 int _height = height;
 char _a[10] = {0};
 int _maxval = maxval;
{ int width = _width; NOT_UNUSED(width);
 int height = _height; NOT_UNUSED(height);
 char * a = _a; NOT_UNUSED(a);
 int maxval = _maxval; NOT_UNUSED(maxval);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/input.h", .line = 73,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 73 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      IF (i >= 0 && i < width && j >= 0 && j < height)
 _stencil_val(__FILE__,__LINE__,s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      
 _stencil_val(__FILE__,__LINE__,s,0,0,0) = 0.;
    } } end_foreach_stencil(); if (_first_call) {
 if (width != _width)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 73, "width");
 }
 if (_first_call) {
 if (height != _height)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 73, "height");
 }
 if (_first_call) {
 for (int i = 0; i < (10*sizeof(char)); i++)
   if (((char *)_a)[i] != 0) {
     reduction_warning ("/home/jiarongw/basilisk/src/input.h", 73, "a");
     break; }
 }
 if (_first_call) {
 if (maxval != _maxval)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 73, "maxval");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 79
foreach(){

#line 73 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((width*height)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  int _width = width;
 int _height = height;
 short _a[10] = {0};
 int _maxval = maxval;
{ int width = _width; NOT_UNUSED(width);
 int height = _height; NOT_UNUSED(height);
 short * a = _a; NOT_UNUSED(a);
 int maxval = _maxval; NOT_UNUSED(maxval);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/input.h", .line = 96,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 96 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      IF (i >= 0 && i < width && j >= 0 && j < height)
 _stencil_val(__FILE__,__LINE__,s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      
 _stencil_val(__FILE__,__LINE__,s,0,0,0) = 0.;
    } } end_foreach_stencil(); if (_first_call) {
 if (width != _width)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 96, "width");
 }
 if (_first_call) {
 if (height != _height)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 96, "height");
 }
 if (_first_call) {
 for (int i = 0; i < (10*sizeof(short)); i++)
   if (((char *)_a)[i] != 0) {
     reduction_warning ("/home/jiarongw/basilisk/src/input.h", 96, "a");
     break; }
 }
 if (_first_call) {
 if (maxval != _maxval)
   reduction_warning ("/home/jiarongw/basilisk/src/input.h", 96, "maxval");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 102
foreach(){

#line 96 "/home/jiarongw/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (ferr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 166 "/home/jiarongw/basilisk/src/input.h"

void input_gfs (struct OutputGfs p)
{ trace ("input_gfs", "/home/jiarongw/basilisk/src/input.h", 168);
  not_mpi_compatible();

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  bool input_all = (p.list == all);
  if (p.list == NULL) p.list = all;





  next_char (p.fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (ferr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    strongif (p.list) for (scalar s = *p.list, *_i67 = p.list; ((scalar *)&s)->i >= 0; s = *++_i67)
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (ferr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = ((scalar []){cm,fm.x,fm.y,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

   { foreach_cell(){

#line 273 "/home/jiarongw/basilisk/src/input.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (ferr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    strongif (input) for (scalar s = *input, *_i68 = input; ((scalar *)&s)->i >= 0; s = *++_i68) {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
 fprintf (ferr, "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }





 }
 else
   val(s,0,0,0) = a;
      }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  strongif (listm) for (scalar s = *listm, *_i69 = listm; ((scalar *)&s)->i >= 0; s = *++_i69)
    if (!is_constant(s))
      _attribute[s.i].dirty = true;
  strongif (input) for (scalar s = *input, *_i70 = input; ((scalar *)&s)->i >= 0; s = *++_i70)
    if (!is_constant(s))
      _attribute[s.i].dirty = true;

  pfree (input,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
 end_trace("input_gfs", "/home/jiarongw/basilisk/src/input.h", 332); }
#line 371 "/home/jiarongw/basilisk/src/input.h"
struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear, periodic, zero;
  int smooth;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  if (fscanf (p.fp, "%s %d", waste, &nx) != 2) {
    fprintf (ferr, "input_grd(): error reading 'nx'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %d", waste, &ny) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ny'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &XG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'XG0'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &YG0) != 2) {
    fprintf (ferr, "input_grd(): error reading 'YG0'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &DeltaGRD) != 2) {
    fprintf (ferr, "input_grd(): error reading 'DeltaGRD'\n");
    if (opened) fclose (p.fp);
    return;
  }
  if (fscanf (p.fp, "%s %lf", waste, &ndv) != 2) {
    fprintf (ferr, "input_grd(): error reading 'ndv'\n");
    if (opened) fclose (p.fp);
    return;
  }


  if (!p.nodatavalue)
    p.nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      if (fscanf (p.fp, "%lf ", &value[j + i*nx]) != 1) {
 fprintf (ferr, "input_grd(): error reading value %d,%d\n", i, j);
 if (opened) fclose (p.fp);
 pfree (value,__func__,__FILE__,__LINE__);
 return;
      }
      if (p.zero && value[j + i*nx] == ndv)
 value[j + i*nx] = 0.;
    }


  if (p.smooth > 0) {
    double * smoothed = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
    for (int s = 0; s < p.smooth; s++) {
      for (int i = 0; i < ny; i++)
 for (int j = 0 ; j < nx; j++) {
   int n = 0;
   smoothed[j + i*nx] = 0.;
   for (int k = -1; k <= 1; k++)
     for (int l = -1; l <= 1; l++)
       if ((l != 0 || k != 0) &&
    i + k >= 0 && i + k < ny &&
    j + l >= 0 && j + l < nx &&
    value[j + l + (i + k)*nx] != ndv)
  smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
   if (n == 0)
     smoothed[j + i*nx] = p.zero ? 0. : ndv;
   else
     smoothed[j + i*nx] /= n;
 }
      swap (double *, value, smoothed);
    }
    pfree (smoothed,__func__,__FILE__,__LINE__);
  }

  bool warning = false;
   { 
#define periodic_bc _periodic_bc
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _XG0 = XG0;
 int _nx = nx;
 double _DeltaGRD = DeltaGRD;
 double _YG0 = YG0;
 int _ny = ny;
 double _value[10] = {0};
 double _ndv = ndv;
 bool _warning = warning;
{ double XG0 = _XG0; NOT_UNUSED(XG0);
 int nx = _nx; NOT_UNUSED(nx);
 double DeltaGRD = _DeltaGRD; NOT_UNUSED(DeltaGRD);
 double YG0 = _YG0; NOT_UNUSED(YG0);
 int ny = _ny; NOT_UNUSED(ny);
 double * value = _value; NOT_UNUSED(value);
 double ndv = _ndv; NOT_UNUSED(ndv);
 bool warning = _warning; NOT_UNUSED(warning);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/input.h", .line = 478,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 478 "/home/jiarongw/basilisk/src/input.h"
 {
    IF (p.periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      IF (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      IF (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    IF (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      IF (p.linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      
 val = value[j + i*nx];
      IF (val == ndv)
 _stencil_val(__FILE__,__LINE__,input,0,0,0) = nodata;
      
 _stencil_val(__FILE__,__LINE__,input,0,0,0) = val;
    }
     {
      _stencil_val(__FILE__,__LINE__,input,0,0,0) = nodata;
      warning = true;
    }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef periodic_bc
#line 517

#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
#line 478
foreach (){

#line 478 "/home/jiarongw/basilisk/src/input.h"
 {
    if (p.periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      else if (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (p.linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      else
 val = value[j + i*nx];
      if (val == ndv)
 val(input,0,0,0) = nodata;
      else
 val(input,0,0,0) = val;
    }
    else {
      val(input,0,0,0) = nodata;
      warning = true;
    }
  } } end_foreach();
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif
#line 517
 }
  pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (ferr,
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (opened)
    fclose (p.fp);
}
#line 74 "/home/jiarongw/basilisk/src/view.h"







typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {


      cexpr tmp = *c;
      while ((c + 1)->expr)
 *c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
     const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = prealloc (cache, sizeof(cexpr)*(len + 2),__func__,__FILE__,__LINE__);
    c = &cache[len];
  }
  else {

    c = cache;
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar []){s,{-1}}));

    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = pstrdup (expr,__func__,__FILE__,__LINE__);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar []){s,{-1}}));
    c++;
  }
  pfree (cache,__func__,__FILE__,__LINE__);
}






struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;
  bool reversed;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  void (* map) (coord *);

  int ni;

  bool active;

  cexpr * cache;
  int maxlen;
};

typedef struct _bview bview;




bview * bview_new()
{
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);


  p->bg[0] = 1; p->bg[1] = 1; p->bg[2] = 1;



  p->res = 1.;
  p->lc = 0.001;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view ()
{
  if (!(_view)) qassert ("/home/jiarongw/basilisk/src/view.h", 237, "_view");
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw () {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  if (view->far <= view->near) {
    double max = 2.;
    gluPerspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else {
    gluPerspective (view->fov, view->width/(float)view->height,
      view->near, view->far);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else


    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  return view;
}







typedef void * pointer;
#line 329 "/home/jiarongw/basilisk/src/view.h"
typedef struct {
  GLubyte a[4];
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
         MPI_Datatype * dptr)
{
  RGBA * in = pin, * out = pout;
  for (int i = 0; i < *len; i++,in++,out++)
    if (out->a[3] == 0)
      *out = *in;
}


static pointer compose_image (bview * view)
{ trace ("compose_image", "/home/jiarongw/basilisk/src/view.h", 344);
  unsigned char * image = framebuffer_image (view->fb);
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_contiguous (4, MPI_BYTE, &rgba);
    MPI_Type_commit (&rgba);
    int size = view->width*view->height;
    if (pid() == 0)
      MPI_Reduce (MPI_IN_PLACE, image, size, rgba, op, 0, MPI_COMM_WORLD);
    else
      MPI_Reduce (image, image, size, rgba, op, 0, MPI_COMM_WORLD);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  { pointer _ret =  image; end_trace("compose_image", "/home/jiarongw/basilisk/src/view.h", 360);  return _ret; }
 end_trace("compose_image", "/home/jiarongw/basilisk/src/view.h", 361); }
#line 419 "/home/jiarongw/basilisk/src/view.h"
#line 1 "vertexbuffer.h"
#line 1 "/home/jiarongw/basilisk/src/vertexbuffer.h"
#line 14 "/home/jiarongw/basilisk/src/vertexbuffer.h"
struct {

  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;
  int dim;
  int vertex, nvertex;
  bool visible;


  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false,
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup ()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free ()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (GL_MODELVIEW_MATRIX, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    swap (float, q[1], q[4]);
    swap (float, q[2], q[8]);
    swap (float, q[6], q[9]);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case GL_LINE_LOOP:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case GL_LINES:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case GL_LINE_STRIP:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;
    case GL_QUADS:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case GL_POLYGON:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case GL_TRIANGLE_FAN:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (ferr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  glBegin (state);
}

static void vertex_buffer_glEnd ()
{
  glEnd();
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case GL_LINE_LOOP:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;

    case GL_LINES:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case GL_LINE_STRIP:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case GL_QUADS:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
 for (int j = 1; j <= 2; j++) {
   vertex_buffer_push_index (i);
   vertex_buffer_push_index (i + j);
   vertex_buffer_push_index (i + j + 1);
 }
      type = 1;
      break;

    case GL_POLYGON:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
    j++) {
 vertex_buffer_push_index (VertexBuffer.polygon);
 vertex_buffer_push_index (VertexBuffer.polygon + j);
 vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;

    case GL_TRIANGLE_FAN:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (VertexBuffer.fan);
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;

    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0) {

      if (!(VertexBuffer.type == type)) qassert ("/home/jiarongw/basilisk/src/vertexbuffer.h", 179, "VertexBuffer.type == type");
    }
    else
      VertexBuffer.type = type;
  }
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  glColor3f (r, g, b);
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b};
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  glNormal3d (nx, ny, nz);
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  glVertex3d (x, y, z);

  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  glVertex3d (x, y, 0.);

  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
}
#line 420 "/home/jiarongw/basilisk/src/view.h"






#line 1 "draw.h"
#line 1 "/home/jiarongw/basilisk/src/draw.h"




#line 1 "fractions.h"
#line 1 "/home/jiarongw/basilisk/src/fractions.h"
#line 12 "/home/jiarongw/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/jiarongw/basilisk/src/geometry.h"
#line 28 "/home/jiarongw/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}
#line 133 "/home/jiarongw/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}
#line 237 "/home/jiarongw/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }}
  return line_area(n1.x, n1.y, alpha);
}
#line 262 "/home/jiarongw/basilisk/src/geometry.h"
int facets (coord n, double alpha, coord p[2])
{
  int i = 0;
  for (double s = -0.5; s <= 0.5; s += 1.)
    {
#line 266

      if (fabs (n.y) > 1e-4 && i < 2) {
 double a = (alpha - s*n.x)/n.y;
 if (a >= -0.5 && a <= 0.5) {
   p[i].x = s;
   p[i++].y = a;
 }
      }
#line 266

      if (fabs (n.x) > 1e-4 && i < 2) {
 double a = (alpha - s*n.y)/n.x;
 if (a >= -0.5 && a <= 0.5) {
   p[i].y = s;
   p[i++].x = a;
 }
      }}
  return i;
}
#line 352 "/home/jiarongw/basilisk/src/geometry.h"
double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 357

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 357

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 368

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 368

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 394
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 394
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}
#line 482 "/home/jiarongw/basilisk/src/geometry.h"
void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 487

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 487

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  {
#line 504

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
#line 504

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }}

  p->x = p->y = cube(alpha);

  {
#line 513
 {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  }
#line 513
 {
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }}

  {
#line 521
 {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  }
#line 521
 {
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }}
}
#line 13 "/home/jiarongw/basilisk/src/fractions.h"






#line 1 "myc2d.h"
#line 1 "/home/jiarongw/basilisk/src/myc2d.h"





coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 7 "/home/jiarongw/basilisk/src/myc2d.h"

  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = val(c,-1,1,0) + val(c,0,1,0) + val(c,1,1,0);
  c_b = val(c,-1,-1,0) + val(c,0,-1,0) + val(c,1,-1,0);
  c_r = val(c,1,-1,0) + val(c,1,0,0) + val(c,1,1,0);
  c_l = val(c,-1,-1,0) + val(c,-1,0,0) + val(c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  if (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
  else {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = val(c,-1,-1,0) + 2.0*val(c,-1,0,0) + val(c,-1,1,0);
  mm2 = val(c,1,-1,0) + 2.0*val(c,1,0,0) + val(c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = val(c,-1,-1,0) + 2.0*val(c,0,-1,0) + val(c,1,-1,0);
  mm2 = val(c,-1,1,0) + 2.0*val(c,0,1,0) + val(c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  if (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    if (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
  else {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    if (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;

#if _call_mycs
}
#define _IN_STENCIL 1

#line 6
static coord _mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 7 "/home/jiarongw/basilisk/src/myc2d.h"

  int ix;
  double c_t,c_b,c_r,c_l;
  double mx0,my0,mx1,my1,mm1,mm2;


  c_t = _stencil_val(__FILE__,__LINE__,c,-1,1,0) + _stencil_val(__FILE__,__LINE__,c,0,1,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  c_b = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + _stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,-1,0);
  c_r = _stencil_val(__FILE__,__LINE__,c,1,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,0,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  c_l = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + _stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,1,0);



  mx0 = 0.5*(c_l-c_r);
  my0 = 0.5*(c_b-c_t);


  IF (fabs(mx0) <= fabs(my0)) {
    my0 = my0 > 0. ? 1. : -1.;
    ix = 1;
  }
   {
    mx0 = mx0 > 0. ? 1. : -1.;
    ix = 0;
  }


  mm1 = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,1,0);
  mm2 = _stencil_val(__FILE__,__LINE__,c,1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,1,0,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  mx1 = mm1 - mm2 + 1.e-30;
  mm1 = _stencil_val(__FILE__,__LINE__,c,-1,-1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,1,-1,0);
  mm2 = _stencil_val(__FILE__,__LINE__,c,-1,1,0) + 2.0*_stencil_val(__FILE__,__LINE__,c,0,1,0) + _stencil_val(__FILE__,__LINE__,c,1,1,0);
  my1 = mm1 - mm2 + 1.e-30;


  IF (ix) {
    mm1 = fabs(my1);
    mm1 = fabs(mx1)/mm1;
    IF (mm1 > fabs(mx0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }
   {
    mm1 = fabs(mx1);
    mm1 = fabs(my1)/mm1;
    IF (mm1 > fabs(my0)) {
      mx0 = mx1;
      my0 = my1;
    }
  }



  mm1 = fabs(mx0) + fabs(my0);
  coord n = {mx0/mm1, my0/mm1};

  return n;

#undef _IN_STENCIL

#endif

#line 65
}
#line 20 "/home/jiarongw/basilisk/src/fractions.h"
#line 121 "/home/jiarongw/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
  double val;
};


void fractions (struct Fractions a)
{ trace ("fractions", "/home/jiarongw/basilisk/src/fractions.h", 130);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
  double val = a.val;
#line 145 "/home/jiarongw/basilisk/src/fractions.h"
  vector p;
  p.x = s.y; p.y = s.x;
#line 155 "/home/jiarongw/basilisk/src/fractions.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _val = val;
{ double val = _val; NOT_UNUSED(val);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/fractions.h", .line = 155,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 155
{

#line 155 "/home/jiarongw/basilisk/src/fractions.h"
 {





    IF ((_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)*(_stencil_val(__FILE__,__LINE__,Phi,1,0,0) - val) < 0.) {






      _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)/(_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - _stencil_val(__FILE__,__LINE__,Phi,1,0,0));
      IF (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) < val)
 _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = 1. - _stencil_val(__FILE__,__LINE__,p.x,0,0,0);
    }
#line 180 "/home/jiarongw/basilisk/src/fractions.h"
    
      _stencil_val(__FILE__,__LINE__,p.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) > val || _stencil_val(__FILE__,__LINE__,Phi,1,0,0) > val);
  } }  }}  { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 155
{

#line 155 "/home/jiarongw/basilisk/src/fractions.h"
 {





    IF ((_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)*(_stencil_val(__FILE__,__LINE__,Phi,0,1,0) - val) < 0.) {






      _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - val)/(_stencil_val(__FILE__,__LINE__,Phi,0,0,0) - _stencil_val(__FILE__,__LINE__,Phi,0,1,0));
      IF (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) < val)
 _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = 1. - _stencil_val(__FILE__,__LINE__,p.y,0,0,0);
    }
#line 180 "/home/jiarongw/basilisk/src/fractions.h"
    
      _stencil_val(__FILE__,__LINE__,p.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,Phi,0,0,0) > val || _stencil_val(__FILE__,__LINE__,Phi,0,1,0) > val);
  } }  }}  end_foreach_face_stencil()
#line 182
 if (_first_call) {
 if (val != _val)
   reduction_warning ("/home/jiarongw/basilisk/src/fractions.h", 155, "val");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 182
foreach_face_generic() { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 155
{

#line 155 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 180 "/home/jiarongw/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  } }  }}  { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 155
{

#line 155 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 180 "/home/jiarongw/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  } }  }}  end_foreach_face_generic()
#line 182
 end_foreach_face(); }
#line 205 "/home/jiarongw/basilisk/src/fractions.h"
  scalar s_z = c;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _val = val;
{ double val = _val; NOT_UNUSED(val);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/fractions.h", .line = 206,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 206 "/home/jiarongw/basilisk/src/fractions.h"


  {
#line 240 "/home/jiarongw/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 242
 {
      n.x = _stencil_val(__FILE__,__LINE__,p.y,0,0,0) - _stencil_val(__FILE__,__LINE__,p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 242
 {
      n.y = _stencil_val(__FILE__,__LINE__,p.x,0,0,0) - _stencil_val(__FILE__,__LINE__,p.x,0,1,0);
      nn += fabs(n.y);
    }}





    IF (nn == 0.)
      _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = _stencil_val(__FILE__,__LINE__,p.x,0,0,0);
     {





      {
#line 259

 n.x /= nn;
#line 259

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 269

   IF (_stencil_val(__FILE__,__LINE__,p.x,0,i,0) > 0. && _stencil_val(__FILE__,__LINE__,p.x,0,i,0) < 1.) {
     double a = sign(_stencil_val(__FILE__,__LINE__,Phi,0,i,0) - val)*(_stencil_val(__FILE__,__LINE__,p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 269

   IF (_stencil_val(__FILE__,__LINE__,p.y,i,0,0) > 0. && _stencil_val(__FILE__,__LINE__,p.y,i,0,0) < 1.) {
     double a = sign(_stencil_val(__FILE__,__LINE__,Phi,i,0,0) - val)*(_stencil_val(__FILE__,__LINE__,p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/home/jiarongw/basilisk/src/fractions.h"
      IF (ni == 0)
 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = max (_stencil_val(__FILE__,__LINE__,p.x,0,0,0), _stencil_val(__FILE__,__LINE__,p.y,0,0,0));
      IF (ni != 4)
 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
       {



 _stencil_val(__FILE__,__LINE__,s_z,0,0,0) = 0.;

      }
    }
  } } end_foreach_stencil(); if (_first_call) {
 if (val != _val)
   reduction_warning ("/home/jiarongw/basilisk/src/fractions.h", 206, "val");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 295
foreach(){

#line 206 "/home/jiarongw/basilisk/src/fractions.h"


  {
#line 240 "/home/jiarongw/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 242
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 242
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 259

 n.x /= nn;
#line 259

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 269

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 269

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 283 "/home/jiarongw/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {



 val(s_z,0,0,0) = 0.;

      }
    }
  } } end_foreach(); }
#line 347 "/home/jiarongw/basilisk/src/fractions.h"
 { strongif (!(a.s).x.i) delete (((scalar []){s.x,s.y,{-1}})); }  end_trace("fractions", "/home/jiarongw/basilisk/src/fractions.h", 347); }
#line 391 "/home/jiarongw/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 392 "/home/jiarongw/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  if (!(2 == 2)) qassert ("/home/jiarongw/basilisk/src/fractions.h", 395, "dimension == 2");
  {
#line 396
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 396
 {
    n.y = (val(c,1,-1,0) + 2.*val(c,0,-1,0) + val(c,-1,-1,0) -
    val(c,1,+1,0) - 2.*val(c,0,+1,0) - val(c,-1,+1,0));
    nn += fabs(n.y);
  }}

  if (nn > 0.)
    {
#line 403

      n.x /= nn;
#line 403

      n.y /= nn;}
  else
    n.x = 1.;
  return n;

#if _call_youngs_normal
}
#define _IN_STENCIL 1

#line 391
static coord _youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 392 "/home/jiarongw/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  IF (!(2 == 2)) _stencil_qassert (__FILE__,__LINE__,"/home/jiarongw/basilisk/src/fractions.h", 395, "dimension == 2");
  {
#line 396
 {
    n.x = (_stencil_val(__FILE__,__LINE__,c,-1,1,0) + 2.*_stencil_val(__FILE__,__LINE__,c,-1,0,0) + _stencil_val(__FILE__,__LINE__,c,-1,-1,0) -
    _stencil_val(__FILE__,__LINE__,c,+1,1,0) - 2.*_stencil_val(__FILE__,__LINE__,c,+1,0,0) - _stencil_val(__FILE__,__LINE__,c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 396
 {
    n.y = (_stencil_val(__FILE__,__LINE__,c,1,-1,0) + 2.*_stencil_val(__FILE__,__LINE__,c,0,-1,0) + _stencil_val(__FILE__,__LINE__,c,-1,-1,0) -
    _stencil_val(__FILE__,__LINE__,c,1,+1,0) - 2.*_stencil_val(__FILE__,__LINE__,c,0,+1,0) - _stencil_val(__FILE__,__LINE__,c,-1,+1,0));
    nn += fabs(n.y);
  }}

  IF (nn > 0.)
    {
#line 403

      n.x /= nn;
#line 403

      n.y /= nn;}
  
    n.x = 1.;
  return n;

#undef _IN_STENCIL

#endif

#line 408
}





coord facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 415 "/home/jiarongw/basilisk/src/fractions.h"

  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    {
#line 419
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 419
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }}
    if (nn > 0.)
      {
#line 424

 n.x /= nn;
#line 424

 n.y /= nn;}
    else
      {
#line 427

 n.x = 1./2;
#line 427

 n.y = 1./2;}
    return n;
  }
  return mycs (point, c);

#if _call_facet_normal
}
#define _IN_STENCIL 1
#define mycs _mycs

#line 414
static coord _facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 415 "/home/jiarongw/basilisk/src/fractions.h"

  IF (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    {
#line 419
 {
      n.x = _stencil_val(__FILE__,__LINE__,s.x,0,0,0) - _stencil_val(__FILE__,__LINE__,s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 419
 {
      n.y = _stencil_val(__FILE__,__LINE__,s.y,0,0,0) - _stencil_val(__FILE__,__LINE__,s.y,0,1,0);
      nn += fabs(n.y);
    }}
    IF (nn > 0.)
      {
#line 424

 n.x /= nn;
#line 424

 n.y /= nn;}
    
      {
#line 427

 n.x = 1./2;
#line 427

 n.y = 1./2;}
    return n;
  }
  return mycs (point, c);

#undef mycs
#undef _IN_STENCIL

#endif

#line 432
}
#line 441 "/home/jiarongw/basilisk/src/fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "/home/jiarongw/basilisk/src/fractions.h", 443);
   { 
#define mycs _mycs
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/fractions.h", .line = 444,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 444 "/home/jiarongw/basilisk/src/fractions.h"
 {





    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) <= 0. || _stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1.) {
      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) = 0.;
      {
#line 452

 _stencil_val(__FILE__,__LINE__,n.x,0,0,0) = 0.;
#line 452

 _stencil_val(__FILE__,__LINE__,n.y,0,0,0) = 0.;}
    }
     {






      coord m = mycs (point, c);
      {
#line 463

 _stencil_val(__FILE__,__LINE__,n.x,0,0,0) = m.x;
#line 463

 _stencil_val(__FILE__,__LINE__,n.y,0,0,0) = m.y;}
      _stencil_val(__FILE__,__LINE__,alpha,0,0,0) = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), m);
    }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef mycs
#line 467
foreach(){

#line 444 "/home/jiarongw/basilisk/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 452

 val(n.x,0,0,0) = 0.;
#line 452

 val(n.y,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);
      {
#line 463

 val(n.x,0,0,0) = m.x;
#line 463

 val(n.y,0,0,0) = m.y;}
      val(alpha,0,0,0) = line_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 485 "/home/jiarongw/basilisk/src/fractions.h"
 end_trace("reconstruction", "/home/jiarongw/basilisk/src/fractions.h", 485); }
#line 505 "/home/jiarongw/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "/home/jiarongw/basilisk/src/fractions.h", 513);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = fout;
  if (!s.x.i) s.x.i = -1;

   { 
#define facet_normal _facet_normal
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/fractions.h", .line = 519,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 519 "/home/jiarongw/basilisk/src/fractions.h"

    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 1e-6 && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), n);

      coord segment[2];
      IF (facets (n, alpha, segment) == 2)
 _stencil_fprintf (__FILE__,__LINE__,p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/home/jiarongw/basilisk/src/fractions.h"
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef facet_normal
#line 538
foreach(){

#line 519 "/home/jiarongw/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = line_alpha (val(c,0,0,0), n);

      coord segment[2];
      if (facets (n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
#line 538 "/home/jiarongw/basilisk/src/fractions.h"
    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "/home/jiarongw/basilisk/src/fractions.h", 541); }








double interface_area (scalar c)
{ trace ("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 551);
  double area = 0.;
   { 
#define mycs _mycs
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _area = area;
{ double area = _area; NOT_UNUSED(area);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/fractions.h", .line = 553,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 553 "/home/jiarongw/basilisk/src/fractions.h"

    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > 1e-6 && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (_stencil_val(__FILE__,__LINE__,c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef mycs
#line 558

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:area)) {

#line 553
foreach (){

#line 553 "/home/jiarongw/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = line_alpha (val(c,0,0,0), n);
      area += pow(Delta, 2 - 1)*line_length_center(n,alpha,&p);
    } } end_foreach();mpi_all_reduce_array (&area, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 558
 }
  { double _ret =  area; end_trace("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 559);  return _ret; }
 end_trace("interface_area", "/home/jiarongw/basilisk/src/fractions.h", 560); }
#line 6 "/home/jiarongw/basilisk/src/draw.h"
#line 1 "gl/font.h"
#line 1 "/home/jiarongw/basilisk/src/gl/font.h"
#line 27 "/home/jiarongw/basilisk/src/gl/font.h"
#include <stdio.h>
#line 1 "gl/og_font.h"
#line 1 "/home/jiarongw/basilisk/src/gl/og_font.h"




typedef struct tagSOG_StrokeVertex SOG_StrokeVertex;
struct tagSOG_StrokeVertex
{
    GLfloat X, Y;
};

typedef struct tagSOG_StrokeStrip SOG_StrokeStrip;
struct tagSOG_StrokeStrip
{
    int Number;
    const SOG_StrokeVertex *Vertices;
};

typedef struct tagSOG_StrokeChar SOG_StrokeChar;
struct tagSOG_StrokeChar
{
    GLfloat Right;
    int Number;
    const SOG_StrokeStrip* Strips;
};

typedef struct tagSOG_StrokeFont SOG_StrokeFont;
struct tagSOG_StrokeFont
{
    char *Name;
    int Quantity;
    GLfloat Height;
    const SOG_StrokeChar **Characters;
};
#line 29 "/home/jiarongw/basilisk/src/gl/font.h"
#line 39 "/home/jiarongw/basilisk/src/gl/font.h"
extern SOG_StrokeFont ogStrokeMonoRoman;
#line 48 "/home/jiarongw/basilisk/src/gl/font.h"
static SOG_StrokeFont *oghStrokeByID( void *font )
{


    if( font == ((void *)0x0001) )
        return &ogStrokeMonoRoman;

    fprintf (ferr, "stroke font %p not found", font );
    return 0;
}
#line 83 "/home/jiarongw/basilisk/src/gl/font.h"
void gl_StrokeCharacter ( int character )
{
    void *fontID = ((void *)0x0001);
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            vertex_buffer_glBegin( GL_LINE_STRIP );
            for( j = 0; j < strip->Number; j++ )
                vertex_buffer_glVertex2d( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            vertex_buffer_glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}
#line 147 "/home/jiarongw/basilisk/src/gl/font.h"
void gl_StrokeString ( const char *string )
{
    void *fontID = ((void *)0x0001);
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )





        while(( c = *string++ ))
       if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            vertex_buffer_glBegin( GL_LINE_STRIP );

                            for( j = 0; j < strip->Number; j++ )
                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            vertex_buffer_glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
     }
}
#line 226 "/home/jiarongw/basilisk/src/gl/font.h"
float gl_StrokeWidth ( int character )
{
    void *fontID = ((void *)0x0001);
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}
#line 269 "/home/jiarongw/basilisk/src/gl/font.h"
float gl_StrokeLength ( const char *string )
{
    void *fontID = ((void *)0x0001);
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}
#line 321 "/home/jiarongw/basilisk/src/gl/font.h"
GLfloat gl_StrokeHeight ()
{
    void *fontID = ((void *)0x0001);
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}
#line 7 "/home/jiarongw/basilisk/src/draw.h"




void clear ()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 49 "/home/jiarongw/basilisk/src/draw.h"
struct _view_set {
  float tx, ty;
  float fov;
  float quat[4];
  float sx, sy, sz;
  unsigned width, height, samples;
  float bg[3];
  float theta, phi, psi;
  bool relative;
  float tz, near, far;
  float res;
  char * camera;
  void (* map) (coord *);
  int cache;
  float p1x, p1y, p2x, p2y;
  bview * view;
};

void view (struct _view_set p)
{
  bview * v = p.view ? p.view : get_view();
  if (p.fov) {
    if (p.relative)
      v->fov += (0.1 + 3.*v->fov)*p.fov;
    else
      v->fov = p.fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (p.quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = p.quat[j];
      break;
    }
  if (p.tx) v->tx = p.relative ? v->tx + p.tx*0.02*(0.01 + 3.*v->fov) : p.tx;
  if (p.ty) v->ty = p.relative ? v->ty + p.ty*0.02*(0.01 + 3.*v->fov) : p.ty;
  if (p.sx) v->sx = p.sx;
  if (p.sy) v->sy = p.sy;
  if (p.sz) v->sz = p.sz;
  if (p.bg[0] || p.bg[1] || p.bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = p.bg[i];

  if (p.camera) {
    v->gfsview = false;
    if (strlen(p.camera) >= 4 &&
 !strcmp (&p.camera[strlen(p.camera) - 4], ".gfv")) {
      FILE * fp = fopen (p.camera, "r");
      if (!fp) {
 perror (p.camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (ferr, "%s: not a valid gfv file\n", p.camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (p.camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (p.camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, pi, v->quat);
    else if (!strcmp (p.camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, pi/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - pi/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (ferr, "view(): unknown camera '%s'\n", p.camera);
      exit (1);
    }
  }
  else if (p.theta || p.phi || p.psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - p.phi, q);
    if (p.relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (p.map)
    v->map = p.map;

  if (p.p1x || p.p1y || p.p2x || p.p2y) {
    float q[4];
    gl_trackball(q, p.p1x, p.p1y, p.p2x, p.p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (p.far > p.near) {
    v->tz = p.tz;
    v->far = p.far;
    v->near = p.near;
  }

  if (p.res)
    v->res = p.res;

  if ((p.width && p.width != v->width) ||
      (p.height && p.height != v->height) ||
      (p.samples && p.samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (p.width) v->width = p.width;
    if (p.height) v->height = p.height;
    if (p.samples) v->samples = p.samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);
    v->fb = framebuffer_new (v->width, v->height);
    init_gl();
  }

  if (p.cache > 0) {
    v->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
    v->maxlen = p.cache;
  }

  clear();
}







struct _translate {
  float x, y, z;
};

void begin_translate (struct _translate p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glTranslatef (p.x, p.y, p.z);
  gl_get_frustum (&view->frustum);
}

void end_translate ()
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 240 "/home/jiarongw/basilisk/src/draw.h"
struct _mirror {
  coord n;
  double alpha;
};

void begin_mirror (struct _mirror p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  normalize (&p.n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*p.n.x*p.n.x;
  s[1] = - 2.*p.n.x*p.n.y; s[2] = - 2.*p.n.x*p.n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*p.n.y*p.n.y; s[6] = - 2.*p.n.y*p.n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*p.n.z*p.n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = 1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = 1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = 1.; t[11] = 0.;
  t[12] = - 2.*p.n.x*p.alpha;
  t[13] = - 2.*p.n.y*p.alpha;
  t[14] = - 2.*p.n.z*p.alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
  view->reversed = !view->reversed;
}

void end_mirror () {
  end_translate();
  bview * view = draw();
  view->reversed = !view->reversed;
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
  double _r = Delta*0.71;\
\
\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (VertexBuffer.visible &&\
      !sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
       < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 322

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 329

#line 380 "/home/jiarongw/basilisk/src/draw.h"
static bool _reversed = false;

static void begin_draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  vertex_buffer_glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_lines ()
{
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

static inline double interp (Point point, coord p, scalar col) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 401 "/home/jiarongw/basilisk/src/draw.h"

  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);

#if _call_interp
}
#define _IN_STENCIL 1
#define interpolate_linear _interpolate_linear

#line 401
static double _interp (Point point, coord p, scalar col) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 401 "/home/jiarongw/basilisk/src/draw.h"

  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);

#undef interpolate_linear
#undef _IN_STENCIL

#endif

#line 404
}

static double evaluate_expression (Point point, Node * n)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 407 "/home/jiarongw/basilisk/src/draw.h"

  if (!(n)) qassert ("/home/jiarongw/basilisk/src/draw.h", 408, "n");
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      strongif (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return val(s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    fprintf (ferr, "unknown operation type '%c'\n", n->type);
    if (!(false)) qassert ("/home/jiarongw/basilisk/src/draw.h", 456, "false");
  }
  return undefined;

#if _call_evaluate_expression
}
#define _IN_STENCIL 1
#define evaluate_expression _evaluate_expression

#line 406
static double _evaluate_expression (Point point, Node * n)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 407 "/home/jiarongw/basilisk/src/draw.h"

  IF (!(n)) _stencil_qassert (__FILE__,__LINE__,"/home/jiarongw/basilisk/src/draw.h", 408, "n");
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      strongif (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return _stencil_val(__FILE__,__LINE__,s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    _stencil_fprintf (__FILE__,__LINE__,ferr, "unknown operation type '%c'\n", n->type);
    IF (!(false)) _stencil_qassert (__FILE__,__LINE__,"/home/jiarongw/basilisk/src/draw.h", 456, "false");
  }
  return undefined;

#undef evaluate_expression
#undef _IN_STENCIL

#endif

#line 459
}

static bool assemble_node (Node * n)
{
  if (n->type == 'v') {
    char * id = n->d.id;
    scalar s = lookup_field (id);
    if (s.i >= 0)
      n->s = s.i;
    else {
      n->s = -1;
      if (!strcmp (id, "Delta"))
 reset_node_type (n, 'D');
      else if (!strcmp (id, "x"))
 reset_node_type (n, 'x');
      else if (!strcmp (id, "y"))
 reset_node_type (n, 'y');
      else if (!strcmp (id, "z"))
 reset_node_type (n, 'z');
      else {
 typedef struct { char * name; double val; } Constant;
 static Constant constants[] = {
   {"pi", pi },
   {"nodata", nodata },
   {"HUGE", HUGE },
   { NULL },
 };
 Constant * p = constants;
 while (p->name) {
   if (!strcmp (p->name, id)) {
     reset_node_type (n, '1');
     n->d.value = p->val;
     break;
   }
   p++;
 }
 if (n->type == 'v') {
   fprintf (ferr, "unknown identifier '%s'\n", id);
   return false;
 }
      }
    }
  }
  for (int i = 0; i < 3; i++)
    if (n->e[i] && !assemble_node (n->e[i]))
      return false;
  return true;
}

static scalar compile_expression (char * expr, bool * isexpr)
{
  *isexpr = false;
  if (!expr)
    return (scalar){-1};

  bview * view = get_view();
  scalar s;
  if (view->cache && (s = get_cexpr (view->cache, expr)).i >= 0)
    return s;

  Node * node = parse_node (expr);
  if (node == NULL) {
    fprintf (ferr, "'%s': syntax error\n", expr);
    return (scalar){-1};
  }
  if (!assemble_node (node)) {
    free_node (node);
    return (scalar){-1};
  }
  if (node->type == 'v' && node->e[0] == NULL) {
    scalar s = {node->s};
    free_node (node);
    return s;
  }
  s = new_scalar("s");
  pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
  _attribute[s.i].name = pstrdup (expr,__func__,__FILE__,__LINE__);
   { 
#define evaluate_expression _evaluate_expression
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/draw.h", .line = 536,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 536 "/home/jiarongw/basilisk/src/draw.h"

    _stencil_val(__FILE__,__LINE__,s,0,0,0) = evaluate_expression (point, node); } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef evaluate_expression
#line 537
foreach(){

#line 536 "/home/jiarongw/basilisk/src/draw.h"

    val(s,0,0,0) = evaluate_expression (point, node); } end_foreach(); }
  restriction (((scalar []){s,{-1}}));
  free_node (node);

  if (view->cache)
    view->cache = add_cexpr (view->cache, view->maxlen, expr, s);
  else
    *isexpr = true;
  return s;
}
#line 605 "/home/jiarongw/basilisk/src/draw.h"
static void begin_colorized (float fc[3], bool constant_color,
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      color j = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = j.r/255.;
      texture[3*i + 1] = j.g/255.;
      texture[3*i + 2] = j.b/255.;
    }
    glTexImage1D (GL_TEXTURE_1D, 0, GL_RGB, 256,0, GL_RGB, GL_FLOAT, texture);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glEnable (GL_TEXTURE_1D);
  }
  if (constant_color)
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized () {
  glDisable (GL_TEXTURE_1D);
}
#line 662 "/home/jiarongw/basilisk/src/draw.h"
struct _draw_vof {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};







static bool cfilter (Point point, scalar c, double cmin)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 684 "/home/jiarongw/basilisk/src/draw.h"

  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    {
#line 687

      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
#line 687

      if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;}
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    {
#line 693

      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
#line 693

      if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;}
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
   { foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 2))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  } end_foreach_neighbor(); }
  return max - min > 0.5;

#if _call_cfilter
}
#define _IN_STENCIL 1

#line 683
static bool _cfilter (Point point, scalar c, double cmin)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 684 "/home/jiarongw/basilisk/src/draw.h"

  double cmin1 = 4.*cmin;
  IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) <= cmin) {
    {
#line 687

      IF (_stencil_val(__FILE__,__LINE__,c,1,0,0) >= 1. - cmin1 || _stencil_val(__FILE__,__LINE__,c,-1,0,0) >= 1. - cmin1)
 return true;
#line 687

      IF (_stencil_val(__FILE__,__LINE__,c,0,1,0) >= 1. - cmin1 || _stencil_val(__FILE__,__LINE__,c,0,-1,0) >= 1. - cmin1)
 return true;}
    return false;
  }
  IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) >= 1. - cmin) {
    {
#line 693

      IF (_stencil_val(__FILE__,__LINE__,c,1,0,0) <= cmin1 || _stencil_val(__FILE__,__LINE__,c,-1,0,0) <= cmin1)
 return true;
#line 693

      IF (_stencil_val(__FILE__,__LINE__,c,0,1,0) <= cmin1 || _stencil_val(__FILE__,__LINE__,c,0,-1,0) <= cmin1)
 return true;}
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
   { foreach_neighbor(1) {
    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > cmin && _stencil_val(__FILE__,__LINE__,c,0,0,0) < 1. - cmin && ++n >= (1 << 2))
      return true;
    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) > max) max = _stencil_val(__FILE__,__LINE__,c,0,0,0);
    IF (_stencil_val(__FILE__,__LINE__,c,0,0,0) < min) min = _stencil_val(__FILE__,__LINE__,c,0,0,0);
  } end_foreach_neighbor(); }
  return max - min > 0.5;

#undef _IN_STENCIL

#endif

#line 707
}

static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    vertex_buffer_glVertex3d (p.x, p.y, p.z);
  }
  else
    vertex_buffer_glVertex3d (x, y, z);
}


static void glvertex2d (bview * view, double x, double y) {
  if (view->map) {
    coord p = {x, y, 0.};
    view->map (&p);
    vertex_buffer_glVertex2d (p.x, p.y);
  }
  else
    vertex_buffer_glVertex2d (x, y);
}

static void glvertex_normal3d (bview * view, Point point, vector n,
          double xp, double yp, double zp)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 732 "/home/jiarongw/basilisk/src/draw.h"

  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  {
#line 734

    np.x = - interp (point, v, n.x);
#line 734

    np.y = - interp (point, v, n.y);}
  vertex_buffer_glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);

#if _call_glvertex_normal3d
}
#define _IN_STENCIL 1
#define interp _interp

#line 730
static void _glvertex_normal3d (bview * view, Point point, vector n,
          double xp, double yp, double zp)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 732 "/home/jiarongw/basilisk/src/draw.h"

  coord v = {(xp - x)/Delta, (yp - y)/Delta}, np;
  {
#line 734

    np.x = - interp (point, v, n.x);
#line 734

    np.y = - interp (point, v, n.y);}
  vertex_buffer_glNormal3d (np.x, np.y, 1.);
  glvertex3d (view, xp, yp, zp);

#undef interp
#undef _IN_STENCIL

#endif

#line 738
}



bool draw_vof (struct _draw_vof p)
{ trace ("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 743);
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (ferr, "draw_vof(): no field named '%s'\n", p.c);
    { bool _ret =  false; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 747);  return _ret; }
  }
  vector s = lookup_vector (p.s);

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 751);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.; if (p.linear && col.i >= 0) boundary_internal ((scalar *)(((scalar []){col,{-1}})), "/home/jiarongw/basilisk/src/draw.h", 751);;

  double cmin = 1e-3;
#line 763 "/home/jiarongw/basilisk/src/draw.h"
  boundary_internal ((scalar *)(((scalar []){c,{-1}})), "/home/jiarongw/basilisk/src/draw.h", 763);

  bview * view = draw();

  if (p.filled) {
    vertex_buffer_glColor3f (p.fc[0], p.fc[1], p.fc[2]);
    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
     { foreach_visible (view){

#line 770 "/home/jiarongw/basilisk/src/draw.h"
 {
      if ((p.filled > 0 && val(c,0,0,0) >= 1.) || (p.filled < 0 && val(c,0,0,0) <= 0.)) {
 vertex_buffer_glBegin (GL_QUADS);
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 vertex_buffer_glEnd();
 view->ni++;
      }
      else if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
 coord n = facet_normal (point, c, s), s = {1.,1.};
 if (p.filled < 0)
   {
#line 783

     n.x = - n.x;
#line 783

     n.y = - n.y;}
 double alpha = line_alpha (p.filled < 0. ? 1. - val(c,0,0,0) : val(c,0,0,0), n);
 alpha += (n.x + n.y)/2.;
 {
#line 787

   if (n.x < 0.) alpha -= n.x, n.x = - n.x, s.x = - 1.;
#line 787

   if (n.y < 0.) alpha -= n.y, n.y = - n.y, s.y = - 1.;}
 coord v[5];
 int nv = 0;
 if (alpha >= 0. && alpha <= n.x) {
   v[nv].x = alpha/n.x, v[nv++].y = 0.;
   if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   else if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   v[nv].x = 0., v[nv++].y = 0.;
 }
 else if (alpha >= n.x && alpha - n.x <= n.y) {
   v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;
   if (alpha >= n.y && alpha - n.y <= n.x) {
     v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
     v[nv].x = 0., v[nv++].y = 1.;
   }
   else if (alpha <= n.y)
     v[nv].x = 0., v[nv++].y = alpha/n.y;
   v[nv].x = 0., v[nv++].y = 0.;
   v[nv].x = 1., v[nv++].y = 0.;
 }
 vertex_buffer_glBegin (GL_POLYGON);
 if (s.x*s.y < 0.)
   for (int i = nv - 1; i >= 0; i--)
     glvertex2d (view, x + s.x*(v[i].x - 0.5)*Delta,
   y + s.y*(v[i].y - 0.5)*Delta);
 else
   for (int i = 0; i < nv; i++)
     glvertex2d (view, x + s.x*(v[i].x - 0.5)*Delta,
   y + s.y*(v[i].y - 0.5)*Delta);
 vertex_buffer_glEnd ();
 view->ni++;
      }
    } } end_foreach_visible(); }
  }
  else
    { begin_draw_lines (view, p.lc, p.lw); {
      vertex_buffer_glBegin (GL_LINES);
       { foreach_visible (view){

#line 829 "/home/jiarongw/basilisk/src/draw.h"

 if (cfilter (point, c, cmin)) {
   coord n = facet_normal (point, c, s);
   double alpha = line_alpha (val(c,0,0,0), n);
   coord segment[2];
   if (facets (n, alpha, segment) == 2) {
     glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);
     glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);
     view->ni++;
   }
 } } end_foreach_visible(); }
      vertex_buffer_glEnd ();
    } end_draw_lines(); }
#line 899 "/home/jiarongw/basilisk/src/draw.h"
  if (p.expr) delete(((scalar []){col,{-1}}));
  { bool _ret =  true; end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 900);  return _ret; }
 end_trace("draw_vof", "/home/jiarongw/basilisk/src/draw.h", 901); }
#line 912 "/home/jiarongw/basilisk/src/draw.h"
struct _isoline {
  char * phi;
  double val;
  int n;


  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};


bool isoline (struct _isoline p)
{ trace ("isoline", "/home/jiarongw/basilisk/src/draw.h", 934);

  if (!p.color) p.color = p.phi;
  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("isoline", "/home/jiarongw/basilisk/src/draw.h", 937);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.; if (p.linear && col.i >= 0) boundary_internal ((scalar *)(((scalar []){col,{-1}})), "/home/jiarongw/basilisk/src/draw.h", 937);;
  scalar phi = col, fiso= new_scalar("fiso");
  vector siso= new_face_vector("siso");
  p.c = "fiso", p.s = "siso";
  struct _draw_vof a = *((struct _draw_vof *)&p.c);
  if (p.n < 2) {
    fractions ((struct Fractions){phi, fiso, siso, p.val});
    draw_vof (a);
  }
  else if (p.max > p.min) {
    double dv = (p.max - p.min)/(p.n - 1);
    for (p.val = p.min; p.val <= p.max; p.val += dv) {
      fractions ((struct Fractions){phi, fiso, siso, p.val});
      draw_vof (a);
    }
  }
  if (p.expr) delete(((scalar []){col,{-1}}));



  { bool _ret =  true; delete (((scalar []){siso.x,siso.y,fiso,{-1}}));  end_trace("isoline", "/home/jiarongw/basilisk/src/draw.h", 957);  return _ret; }
 delete (((scalar []){siso.x,siso.y,fiso,{-1}}));  end_trace("isoline", "/home/jiarongw/basilisk/src/draw.h", 958); }
#line 971 "/home/jiarongw/basilisk/src/draw.h"
struct _cells {
  coord n;
  double alpha;
  float lc[3], lw;
};


bool cells (struct _cells p)
{ trace ("cells", "/home/jiarongw/basilisk/src/draw.h", 979);
  bview * view = draw();
  { begin_draw_lines (view, p.lc, p.lw); {

     { foreach_visible (view){

#line 983 "/home/jiarongw/basilisk/src/draw.h"
 {
      vertex_buffer_glBegin (GL_LINE_LOOP);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd();
      view->ni++;
    } } end_foreach_visible(); }
#line 1005 "/home/jiarongw/basilisk/src/draw.h"
  } end_draw_lines(); }
  { bool _ret =  true; end_trace("cells", "/home/jiarongw/basilisk/src/draw.h", 1006);  return _ret; }
 end_trace("cells", "/home/jiarongw/basilisk/src/draw.h", 1007); }






struct _vectors {
  char * u;
  double scale;
  float lc[3], lw;
};


bool vectors (struct _vectors p)
{ trace ("vectors", "/home/jiarongw/basilisk/src/draw.h", 1022);

  vector u;
  struct { char x, y, z; } index = {'x', 'y', 'z'};
  {
#line 1026
 {
    char name[80];
    sprintf (name, "%s.%c", p.u, index.x);
    u.x = lookup_field (name);
  }
#line 1026
 {
    char name[80];
    sprintf (name, "%s.%c", p.u, index.y);
    u.y = lookup_field (name);
  }}
  bview * view = draw();
  float res = view->res;
  if (view->res < 15*view->samples)
    view->res = 15*view->samples;
  { begin_draw_lines (view, p.lc, p.lw); {
    double scale = (p.scale ? p.scale : 1.)*view->res/view->samples;
    vertex_buffer_glBegin (GL_LINES);
     { foreach_visible (view){

#line 1038 "/home/jiarongw/basilisk/src/draw.h"

      if (val(u.x,0,0,0) != nodata) {
 coord f = { scale*val(u.x,0,0,0), scale*val(u.y,0,0,0) };
 glvertex2d (view, x + f.x - (f.x - f.y/2.)/5.,
      y + f.y - (f.x/2. + f.y)/5.);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x, y + f.y);
 glvertex2d (view, x + f.x - (f.x + f.y/2.)/5.,
      y + f.y + (f.x/2. - f.y)/5.);
 glvertex2d (view, x, y);
 glvertex2d (view, x + f.x, y + f.y);
 view->ni++;
      } } end_foreach_visible(); }
    vertex_buffer_glEnd();
  } end_draw_lines(); }
  view->res = res;



  { bool _ret =  true; end_trace("vectors", "/home/jiarongw/basilisk/src/draw.h", 1057);  return _ret; }
 end_trace("vectors", "/home/jiarongw/basilisk/src/draw.h", 1058); }
#line 1078 "/home/jiarongw/basilisk/src/draw.h"
struct _squares {
  char * color;
  char * z;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];
  bool expr;

  coord n;
  double alpha;
};


bool squares (struct _squares p)
{ trace ("squares", "/home/jiarongw/basilisk/src/draw.h", 1093);

  scalar Z = {-1};
  vector n;
  bool zexpr = false;
  if (p.z) {
    Z = compile_expression (p.z, &zexpr);
    if (Z.i < 0)
      { bool _ret =  false; end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 1101);  return _ret; }
    n = new_vector("n");
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/draw.h", .line = 1103,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 1103 "/home/jiarongw/basilisk/src/draw.h"

      {
#line 1104

        _stencil_val(__FILE__,__LINE__,n.x,0,0,0) = (_stencil_val(__FILE__,__LINE__,Z,1,0,0) - _stencil_val(__FILE__,__LINE__,Z,-1,0,0))/(2.*Delta_x);
#line 1104

        _stencil_val(__FILE__,__LINE__,n.y,0,0,0) = (_stencil_val(__FILE__,__LINE__,Z,0,1,0) - _stencil_val(__FILE__,__LINE__,Z,0,-1,0))/(2.*Delta_y);}; } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 1105
foreach(){

#line 1103 "/home/jiarongw/basilisk/src/draw.h"

      {
#line 1104

        val(n.x,0,0,0) = (val(Z,1,0,0) - val(Z,-1,0,0))/(2.*Delta_x);
#line 1104

        val(n.y,0,0,0) = (val(Z,0,1,0) - val(Z,0,-1,0))/(2.*Delta_y);}; } end_foreach(); }
    boundary_internal ((scalar *)((scalar *)((vector []){{n.x,n.y},{{-1},{-1}}})), "/home/jiarongw/basilisk/src/draw.h", 1106);
  }

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 1109);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((2 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.; if (p.linear && col.i >= 0) boundary_internal ((scalar *)(((scalar []){col,{-1}})), "/home/jiarongw/basilisk/src/draw.h", 1109);;
  scalar f = col;

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (p.linear) {
    { begin_colorized (p.fc, !VertexBuffer.color || !p.color, cmap, !VertexBuffer.color && p.color && p.linear && col.i >= 0); {

      if (Z.i < 0) {
 vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
  { foreach_visible (view){

#line 1119 "/home/jiarongw/basilisk/src/draw.h"

   if (val(f,0,0,0) != nodata) {
     vertex_buffer_glBegin (GL_TRIANGLE_FAN);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } }



      ;
     glvertex2d (view, x, y);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
     vertex_buffer_glEnd();
     view->ni++;
   } } end_foreach_visible(); }
      }
      else
  { foreach_leaf(){

#line 1143 "/home/jiarongw/basilisk/src/draw.h"

   if (val(f,0,0,0) != nodata) {
     vertex_buffer_glBegin (GL_TRIANGLE_FAN);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (4.*val(f,0,0,0) + 2.*(val(f,1,0,0) + val(f,-1,0,0) + val(f,0,1,0) + val(f,0,-1,0)) + val(f,-1,-1,0) + val(f,1,1,0) + val(f,-1,1,0) + val(f,1,-1,0))/16.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } }

                                                    ;
     glvertex_normal3d (view, point, n, x, y, val(Z,0,0,0));
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex_normal3d (view, point, n, x + Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,-1,0) + val(Z,0,-1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,1,0,0) + val(f,1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex_normal3d (view, point, n, x + Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,1,0,0) + val(Z,1,1,0) + val(Z,0,1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,1,0) + val(f,0,1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y + Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,1,0) + val(Z,0,1,0))/4.);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4., p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = (val(f,0,0,0) + val(f,-1,0,0) + val(f,-1,-1,0) + val(f,0,-1,0))/4.; glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex_normal3d (view, point, n, x - Delta_x/2., y - Delta_y/2.,
          (val(Z,0,0,0) + val(Z,-1,0,0) + val(Z,-1,-1,0) + val(Z,0,-1,0))/4.);
     vertex_buffer_glEnd();
     view->ni++;
   } } end_foreach_leaf(); }
#line 1194 "/home/jiarongw/basilisk/src/draw.h"
    } end_colorized(); }
  }
  else {

    vertex_buffer_glNormal3d (0, 0, view->reversed ? -1 : 1);
    vertex_buffer_glBegin (GL_QUADS);
     { foreach_visible (view){

#line 1200 "/home/jiarongw/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
 view->ni++;
      } } end_foreach_visible(); }
    vertex_buffer_glEnd();
#line 1230 "/home/jiarongw/basilisk/src/draw.h"
  }
  if (p.expr) delete (((scalar []){col,{-1}}));

  if (zexpr) delete (((scalar []){Z,{-1}}));
  if (p.z) delete ((scalar *)((vector []){{n.x,n.y},{{-1},{-1}}}));

  { bool _ret =  true; end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 1236);  return _ret; }
 end_trace("squares", "/home/jiarongw/basilisk/src/draw.h", 1237); }
#line 1248 "/home/jiarongw/basilisk/src/draw.h"
struct _box {
  bool notics;
  float lc[3], lw;
};


bool box (struct _box p)
{ trace ("box", "/home/jiarongw/basilisk/src/draw.h", 1255);
  bview * view = draw();
  { begin_draw_lines (view, p.lc, p.lw); {

    float height = 0.5*gl_StrokeHeight();
    float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
    float Z1 = 2 == 2 ? 0. : Z0;
    char label[80];

    glMatrixMode (GL_MODELVIEW);

    if (!p.notics) {
      int nt = 8;
      for (int i = 0; i <= nt; i++) {
 glPushMatrix();
 glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
 glRotatef (-90, 0, 0, 1);
 glScalef (scale, scale, 1.);
 sprintf (label, "%g", X0 + i*L0/nt);
 gl_StrokeString (label);
 glPopMatrix();

 glPushMatrix();
 sprintf (label, "%g", Y0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0 + i*L0/nt - height/2.*scale, Z1);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();
#line 1297 "/home/jiarongw/basilisk/src/draw.h"
      }

      glPushMatrix();
      sprintf (label, "%g", X0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("X");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Y0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0 + L0/2. - height*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Y");
      glPopMatrix();
#line 1328 "/home/jiarongw/basilisk/src/draw.h"
    }


     { foreach_level (0){

#line 1331 "/home/jiarongw/basilisk/src/draw.h"
 {
      vertex_buffer_glBegin (GL_LINE_LOOP);
      glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
      glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
      glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
      vertex_buffer_glEnd ();
      view->ni++;
    } } end_foreach_level(); }
#line 1360 "/home/jiarongw/basilisk/src/draw.h"
  } end_draw_lines(); }
  { bool _ret =  true; end_trace("box", "/home/jiarongw/basilisk/src/draw.h", 1361);  return _ret; }
 end_trace("box", "/home/jiarongw/basilisk/src/draw.h", 1362); }
#line 1375 "/home/jiarongw/basilisk/src/draw.h"
struct _isosurface {
  char * f;
  double v;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};


bool isosurface (struct _isosurface p)
{ trace ("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1389);
#line 1450 "/home/jiarongw/basilisk/src/draw.h"
  { bool _ret =  true; end_trace("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1450);  return _ret; }
 end_trace("isosurface", "/home/jiarongw/basilisk/src/draw.h", 1451); }
#line 1461 "/home/jiarongw/basilisk/src/draw.h"
struct _travelling {
  double start, end;
  float tx, ty, quat[4], fov;
};




void travelling (struct _travelling p)
{
  static float tx, ty, quat[4], fov;
  static double told = -1.;
  if (told < p.start && t >= p.start) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  if (t >= p.start && t <= p.end)
    view ((struct _view_set){.tx = (!p.tx ? tx : ((t - p.start)*(p.tx) + (p.end - t)*(tx))/(p.end - p.start)), .ty = (!p.ty ? ty : ((t - p.start)*(p.ty) + (p.end - t)*(ty))/(p.end - p.start)),
   .fov = (!p.fov ? fov : ((t - p.start)*(p.fov) + (p.end - t)*(fov))/(p.end - p.start)),
   .quat = {(!p.quat[0] ? quat[0] : ((t - p.start)*(p.quat[0]) + (p.end - t)*(quat[0]))/(p.end - p.start)), (!p.quat[1] ? quat[1] : ((t - p.start)*(p.quat[1]) + (p.end - t)*(quat[1]))/(p.end - p.start)),
           (!p.quat[2] ? quat[2] : ((t - p.start)*(p.quat[2]) + (p.end - t)*(quat[2]))/(p.end - p.start)), (!p.quat[3] ? quat[3] : ((t - p.start)*(p.quat[3]) + (p.end - t)*(quat[3]))/(p.end - p.start))}});
  if (told < p.end && t >= p.end) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  told = t;
}
#line 1508 "/home/jiarongw/basilisk/src/draw.h"
struct _draw_string {
  char * str;
  int pos;
  float size;
  float lc[3], lw;
};


bool draw_string (struct _draw_string p)
{ trace ("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1517);
  bview * view = draw();

  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  vertex_buffer_glColor3f (p.lc[0], p.lc[1], p.lc[2]);
  glLineWidth (view->samples*(p.lw > 0. ? p.lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!p.size)
    p.size = 40;
  float hscale = 2./(p.size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (p.pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (p.pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (p.pos == 2)
    glTranslatef (1. - strlen(p.str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(p.str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (p.str);

  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();

  { bool _ret =  true; end_trace("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1552);  return _ret; }
 end_trace("draw_string", "/home/jiarongw/basilisk/src/draw.h", 1553); }




struct _labels {
  char * f;
  float lc[3], lw;
};


bool labels (struct _labels p)
{ trace ("labels", "/home/jiarongw/basilisk/src/draw.h", 1565);

  bool expr = false;
  scalar f = compile_expression (p.f, &expr);
  if (f.i < 0)
    { bool _ret =  false; end_trace("labels", "/home/jiarongw/basilisk/src/draw.h", 1570);  return _ret; }
  bview * view = draw();
  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  float res = view->res;
  if (view->res < 150*view->samples)
    view->res = 150*view->samples;
  { begin_draw_lines (view, p.lc, p.lw); {
    glMatrixMode (GL_MODELVIEW);
     { foreach_visible (view){

#line 1578 "/home/jiarongw/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 glPushMatrix();
 char s[80];
 sprintf (s, "%g", val(f,0,0,0));
 float scale = 0.8*Delta_x/(strlen(s)*width);
 glTranslatef (x - 0.4*Delta_x, y - scale*height/3., 0.);
 glScalef (scale, scale, 1.);
 gl_StrokeString (s);
 glPopMatrix();
      } } end_foreach_visible(); }
  } end_draw_lines(); }
  view->res = res;
  if (expr) delete (((scalar []){f,{-1}}));
  { bool _ret =  true; end_trace("labels", "/home/jiarongw/basilisk/src/draw.h", 1592);  return _ret; }




 end_trace("labels", "/home/jiarongw/basilisk/src/draw.h", 1597); }







#line 1 "draw_json.h"
#line 1 "/home/jiarongw/basilisk/src/draw_json.h"



int _view_set_json (void * q, char * s, int len) {
  struct _view_set * p = (struct _view_set *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"view_set\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sy\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sy);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"width\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->width);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"height\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->height);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"samples\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->samples);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"bg\": { \"type\": \"pfloat\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->bg[0], p->bg[1], p->bg[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"theta\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->theta);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"phi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"psi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->psi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"relative\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->relative);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"near\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->near);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"far\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->far);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"res\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->res);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"camera\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->camera);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"cache\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->cache);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _translate_json (void * q, char * s, int len) {
  struct _translate * p = (struct _translate *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"translate\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _mirror_json (void * q, char * s, int len) {
  struct _mirror * p = (struct _mirror *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"mirror\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_vof_json (void * q, char * s, int len) {
  struct _draw_vof * p = (struct _draw_vof *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_vof\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isoline_json (void * q, char * s, int len) {
  struct _isoline * p = (struct _isoline *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isoline\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"phi\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"val\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->val);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->n);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _cells_json (void * q, char * s, int len) {
  struct _cells * p = (struct _cells *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"cells\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _vectors_json (void * q, char * s, int len) {
  struct _vectors * p = (struct _vectors *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"vectors\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"u\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->u);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"scale\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->scale);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _squares_json (void * q, char * s, int len) {
  struct _squares * p = (struct _squares *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"squares\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _box_json (void * q, char * s, int len) {
  struct _box * p = (struct _box *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"box\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"notics\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->notics);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isosurface_json (void * q, char * s, int len) {
  struct _isosurface * p = (struct _isosurface *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isosurface\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"v\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->v);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _travelling_json (void * q, char * s, int len) {
  struct _travelling * p = (struct _travelling *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"travelling\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"start\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->start);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"end\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->end);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_string_json (void * q, char * s, int len) {
  struct _draw_string * p = (struct _draw_string *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_string\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"str\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->str);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"pos\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->pos);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"size\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->size);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _labels_json (void * q, char * s, int len) {
  struct _labels * p = (struct _labels *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"labels\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
#line 1606 "/home/jiarongw/basilisk/src/draw.h"

struct {
  int (* json) (void * q, char * s, int len);
} bview_interface[] = {
  { _draw_vof_json },
  { _squares_json },
  { _cells_json },
  { _box_json },

  { _isoline_json },
  { _labels_json },
  { _vectors_json },



  { NULL }
};
#line 427 "/home/jiarongw/basilisk/src/view.h"
#line 442 "/home/jiarongw/basilisk/src/view.h"
struct _load {
  FILE * fp;
  char * file;
  Array * buf;
};

bool load (struct _load p);
#line 496 "/home/jiarongw/basilisk/src/view.h"
struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw;
  int sort, options;

  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}


bool save (struct _save p)
{ trace ("save", "/home/jiarongw/basilisk/src/view.h", 516);
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
 dot = s;
 s = strchr (s + 1, '.');
      }
      if (dot)
 p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();

  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      if (!fp) {
 perror (p.file);
 { bool _ret =  false; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 542);  return _ret; }
      }
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    { bool _ret =  true; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 547);  return _ret; }
  }

  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    { bool _ret =  false; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 552);  return _ret; }
  }
  if (!p.fp)
    p.fp = fout;

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (p.format, "bv")) {

    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);
#line 584 "/home/jiarongw/basilisk/src/view.h"
  }

  else if (!strcmp (p.format, "gnu") ||
    !strcmp (p.format, "obj") ||
    !strcmp (p.format, "kml") ||
    !strcmp (p.format, "ps") ||
    !strcmp (p.format, "eps") ||
    !strcmp (p.format, "tex") ||
    !strcmp (p.format, "pdf") ||
    !strcmp (p.format, "svg") ||
    !strcmp (p.format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    { bool _ret =  false; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 604);  return _ret; }
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

  { bool _ret =  true; end_trace("save", "/home/jiarongw/basilisk/src/view.h", 611);  return _ret; }
 end_trace("save", "/home/jiarongw/basilisk/src/view.h", 612); }







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}






#line 1 "draw_get.h"
#line 1 "/home/jiarongw/basilisk/src/draw_get.h"

#line 1 "parse.h"
#line 1 "/home/jiarongw/basilisk/src/parse.h"




enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}

static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
      "bool", "float", "double", "colormap" };
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (ferr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (ferr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0;
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
 nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((colormap *)p->val) = randomap;
    else {
      fprintf (ferr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (ferr, "expecting a %s for '%s' got %s\n",
        name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: if (!(false)) qassert ("/home/jiarongw/basilisk/src/parse.h", 77, "false");
      }
    }
    else {
      if (val[0] != '{') {
 fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
 return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
 char * s = strchr (val, ',');
 if (!s)
   s = strchr (val, '}');
 if (!s) {
   fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
   return false;
 }
 c = *s;
 *s++ = '\0';
 switch (p->type) {
 case pint: ((int *)p->val)[i++] = atoi (val); break;
 case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
 case pbool: ((bool *)p->val)[i++] = atobool (val); break;
 case pfloat: ((float *)p->val)[i++] = atof (val); break;
 case pdouble: ((double *)p->val)[i++] = atof (val); break;
 default: if (!(false)) qassert ("/home/jiarongw/basilisk/src/parse.h", 104, "false");
 }
 val = s;
      }
      if (c != '}') {
 fprintf (ferr, "expecting '}' for '%s' got %s\n", p->key, val);
 return false;
      }
    }
    break;

  default:
    if (!(false)) qassert ("/home/jiarongw/basilisk/src/parse.h", 116, "false");
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
 fprintf (ferr, "unknown key '%s'\n", key);
 return false;
      }
      if (!args (p, s))
 return false;
    }
    else {
      if (i < 0) {
 fprintf (ferr, "anonymous value '%s' after keys\n", key);
 return false;
      }
      if (i >= n) {
 fprintf (ferr, "too many parameters: '%s' %d %d\n", key, i, n);
 return false;
      }
      if (!args (&params[i], key))
 return false;
      i++;
    }
    s = next;
  }
  return true;
}
#line 3 "/home/jiarongw/basilisk/src/draw_get.h"

bool _view_set_get (struct _view_set * p) {
  Params params[] = {
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"fov", pfloat, &p->fov},
    {"quat", pfloat, p->quat, 4},
    {"sx", pfloat, &p->sx},
    {"sy", pfloat, &p->sy},
    {"sz", pfloat, &p->sz},
    {"width", punsigned, &p->width},
    {"height", punsigned, &p->height},
    {"samples", punsigned, &p->samples},
    {"bg", pfloat, p->bg, 3},
    {"theta", pfloat, &p->theta},
    {"phi", pfloat, &p->phi},
    {"psi", pfloat, &p->psi},
    {"relative", pbool, &p->relative},
    {"tz", pfloat, &p->tz},
    {"near", pfloat, &p->near},
    {"far", pfloat, &p->far},
    {"res", pfloat, &p->res},
    {"camera", pstring, &p->camera},
    {"cache", pint, &p->cache},
    {"p2y", pfloat, &p->p2y},
    {"p1x", pfloat, &p->p1x},
    {"p1y", pfloat, &p->p1y},
    {"p2x", pfloat, &p->p2x},
    {NULL}
  };
  return parse_params (params);
}

bool _translate_get (struct _translate * p) {
  Params params[] = {
    {"x", pfloat, &p->x},
    {"y", pfloat, &p->y},
    {"z", pfloat, &p->z},
    {NULL}
  };
  return parse_params (params);
}

bool _mirror_get (struct _mirror * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_vof_get (struct _draw_vof * p) {
  Params params[] = {
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _isoline_get (struct _isoline * p) {
  Params params[] = {
    {"phi", pstring, &p->phi},
    {"val", pdouble, &p->val},
    {"n", pint, &p->n},
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _cells_get (struct _cells * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _vectors_get (struct _vectors * p) {
  Params params[] = {
    {"u", pstring, &p->u},
    {"scale", pdouble, &p->scale},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _squares_get (struct _squares * p) {
  Params params[] = {
    {"color", pstring, &p->color},
    {"z", pstring, &p->z},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"expr", pbool, &p->expr},
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _box_get (struct _box * p) {
  Params params[] = {
    {"notics", pbool, &p->notics},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _isosurface_get (struct _isosurface * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"v", pdouble, &p->v},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _travelling_get (struct _travelling * p) {
  Params params[] = {
    {"start", pdouble, &p->start},
    {"end", pdouble, &p->end},
    {"fov", pfloat, &p->fov},
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"quat", pfloat, p->quat, 4},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_string_get (struct _draw_string * p) {
  Params params[] = {
    {"str", pstring, &p->str},
    {"pos", pint, &p->pos},
    {"size", pfloat, &p->size},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _labels_get (struct _labels * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}
#line 642 "/home/jiarongw/basilisk/src/view.h"

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * buf = pstrdup (line,__func__,__FILE__,__LINE__);
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
 free_cexpr (view->cache);
 view->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
      }
      if (!restore ((struct Dump){.file = file, .list = all}))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ((struct Dump){.file = file});
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ((struct OutputGfs){.file = file, .list = all});
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ((struct _save){.file = file});
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load ((struct _load){.file = file});
  }

  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    if (!_cells_get (&p) || !cells (p))
      return false;
  }

  else if (!strcmp (s, "vectors")) {
    struct _vectors p = {0};
    if (!_vectors_get (&p) || !vectors (p))
      return false;
  }

  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    if (!_draw_vof_get (&p) || !draw_vof (p))
      return false;
  }

  else if (!strcmp (s, "isoline")) {
    struct _isoline p = {0};
    if (!_isoline_get (&p) || !isoline (p))
      return false;
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
  }

  else if (!strcmp (s, "end_translate"))
    end_translate();

  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
  }

  else if (!strcmp (s, "end_mirror"))
    end_mirror();

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    if (!_isosurface_get (&p) || !isosurface (p))
      return false;
  }

  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    if (!_draw_string_get (&p) || !draw_string (p))
      return false;
  }

  else if (!strcmp (s, "labels")) {
    struct _labels p = {0};
    if (!_labels_get (&p) || !labels (p))
      return false;
  }

  else if (!strcmp (s, "clear"))
    clear();

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    if (!_box_get (&p) || !box (p))
      return false;
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
  }

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  pfree (buf,__func__,__FILE__,__LINE__);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  if (p.fp) {
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line));
  }
  else if (p.buf) {
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line);
      }
    }
  }
  return true;
}
#line 7 "field_new_PM_restore.c"
#line 1 "layered/hydro.h"
#line 1 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 46 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 1 "./utils.h"
#line 47 "/home/jiarongw/basilisk/src/layered/hydro.h"

scalar zb= {0}, eta, h;
vector u;
double G = 1., dry = 1e-12, CFL_H = 1e40;
double (* gradient) (double, double, double) = minmod2;

scalar * tracers = NULL;
bool linearised = false;
#line 87 "/home/jiarongw/basilisk/src/layered/hydro.h"
static int defaults0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults0 (const int i, const double t, Event * _ev) { trace ("defaults0", "/home/jiarongw/basilisk/src/layered/hydro.h", 87); 
{
  if (!(nl > 0)) qassert ("/home/jiarongw/basilisk/src/layered/hydro.h", 89, "nl > 0");
  h = new_block_scalar("h", "", nl);
  _attribute[h.i].gradient = gradient;




  eta = new_scalar("eta");
  reset (((scalar []){h,zb,{-1}}), 0.);




  _attribute[zb.i].gradient = gradient;
  _attribute[eta.i].gradient = gradient;







 end_trace("defaults0", "/home/jiarongw/basilisk/src/layered/hydro.h", 111); } return 0; } 

#line 1 "./run.h"
#line 1 "/home/jiarongw/basilisk/src/run.h"
#line 9 "/home/jiarongw/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 12 "/home/jiarongw/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/jiarongw/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/jiarongw/basilisk/src/run.h", 37); }




static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/jiarongw/basilisk/src/run.h", 42);  {
  display ((struct _display){"box();"});
 end_trace("defaults", "/home/jiarongw/basilisk/src/run.h", 44); } return 0; } 





static int cleanup_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup (const int i, const double t, Event * _ev) { trace ("cleanup", "/home/jiarongw/basilisk/src/run.h", 50);  {
  display ((struct _display){"", true});
 end_trace("cleanup", "/home/jiarongw/basilisk/src/run.h", 52); } return 0; } 
#line 114 "/home/jiarongw/basilisk/src/layered/hydro.h"





static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/jiarongw/basilisk/src/layered/hydro.h", 119); 
{






  CFL = 1./(2.*2);
  if (CFL_H == 1e40)
    CFL_H = 0.5;

  u = new_block_vector("u", nl);
  reset (((vector []){{u.x,u.y},{{-1},{-1}}}), 0.);

  if (!linearised)
    {
#line 135

      tracers = list_append (tracers, u.x);
#line 135

      tracers = list_append (tracers, u.y);}





  strongif (tracers) for (scalar s = *tracers, *_i71 = tracers; ((scalar *)&s)->i >= 0; s = *++_i71) {
    _attribute[s.i].gradient = gradient;




  }




  display ((struct _display){"squares (color = 'eta > zb ? eta : nodata', spread = -1);"});
 end_trace("defaults_0", "/home/jiarongw/basilisk/src/layered/hydro.h", 154); } return 0; } 




static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/jiarongw/basilisk/src/layered/hydro.h", 159); 
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 161,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 161 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    _stencil_val(__FILE__,__LINE__,eta,0,0,0) = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner()
      _stencil_val(__FILE__,__LINE__,eta,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block_inner(); }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 165
foreach(){

#line 161 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    val(eta,0,0,0) = val(zb,0,0,0);
     { foreach_block_inner()
      val(eta,0,0,0) += val(h,0,0,0); end_foreach_block_inner(); }
  } } end_foreach(); }
 end_trace("init", "/home/jiarongw/basilisk/src/layered/hydro.h", 166); } return 0; } 





double dtmax;
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/jiarongw/basilisk/src/layered/hydro.h", 173);  dtmax = DT; end_trace("set_dtmax", "/home/jiarongw/basilisk/src/layered/hydro.h", 173);  return 0; } 
#line 193 "/home/jiarongw/basilisk/src/layered/hydro.h"
static bool hydrostatic = true;
vector hu, hf, ha;

static int face_fields_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int face_fields (const int i, const double t, Event * _ev) { trace ("face_fields", "/home/jiarongw/basilisk/src/layered/hydro.h", 196); 
{
  hu = new_block_face_vector("hu", nl);
  hf = new_block_face_vector("hf", nl);
  ha = new_block_face_vector("ha", nl);
#line 209 "/home/jiarongw/basilisk/src/layered/hydro.h"
  static double pdt = 1e-6;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _G = G;
 double _dry = dry;
 double _pdt = pdt;
 double _CFL = CFL;
 bool _hydrostatic = hydrostatic;
 double _CFL_H = CFL_H;
 double _dtmax = dtmax;
{ double G = _G; NOT_UNUSED(G);
 double dry = _dry; NOT_UNUSED(dry);
 double pdt = _pdt; NOT_UNUSED(pdt);
 double CFL = _CFL; NOT_UNUSED(CFL);
 bool hydrostatic = _hydrostatic; NOT_UNUSED(hydrostatic);
 double CFL_H = _CFL_H; NOT_UNUSED(CFL_H);
 double dtmax = _dtmax; NOT_UNUSED(dtmax);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 210,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 210
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,i-1,0,0), _stencil_val(__FILE__,__LINE__,h,i,0,0), _stencil_val(__FILE__,__LINE__,h,i+1,0,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,i+1,0,0) - _stencil_val(__FILE__,__LINE__,h,i-1,0,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,0,i-1,0), _stencil_val(__FILE__,__LINE__,h,0,i,0), _stencil_val(__FILE__,__LINE__,h,0,i+1,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,0,i+1,0) - _stencil_val(__FILE__,__LINE__,h,0,i-1,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_stencil()
#line 260
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 210
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,i-1,0,0), _stencil_val(__FILE__,__LINE__,h,i,0,0), _stencil_val(__FILE__,__LINE__,h,i+1,0,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,i+1,0,0) - _stencil_val(__FILE__,__LINE__,h,i-1,0,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,0,i-1,0), _stencil_val(__FILE__,__LINE__,h,0,i,0), _stencil_val(__FILE__,__LINE__,h,0,i+1,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,0,i+1,0) - _stencil_val(__FILE__,__LINE__,h,0,i-1,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_stencil()
#line 260
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 210
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,i-1,0,0), _stencil_val(__FILE__,__LINE__,h,i,0,0), _stencil_val(__FILE__,__LINE__,h,i+1,0,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,i+1,0,0) - _stencil_val(__FILE__,__LINE__,h,i-1,0,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,0,i-1,0), _stencil_val(__FILE__,__LINE__,h,0,i,0), _stencil_val(__FILE__,__LINE__,h,0,i+1,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,0,i+1,0) - _stencil_val(__FILE__,__LINE__,h,0,i-1,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_stencil()
#line 260
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 210
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,i-1,0,0), _stencil_val(__FILE__,__LINE__,h,i,0,0), _stencil_val(__FILE__,__LINE__,h,i+1,0,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,i+1,0,0) - _stencil_val(__FILE__,__LINE__,h,i-1,0,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (_stencil_val(__FILE__,__LINE__,h,0,i-1,0), _stencil_val(__FILE__,__LINE__,h,0,i,0), _stencil_val(__FILE__,__LINE__,h,0,i+1,0))/Delta :
 (_stencil_val(__FILE__,__LINE__,h,0,i+1,0) - _stencil_val(__FILE__,__LINE__,h,0,i-1,0))/(2.*Delta);
      hff = _stencil_val(__FILE__,__LINE__,h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      _stencil_val(__FILE__,__LINE__,hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      IF (fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) > um)
 um = fabs(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0));

      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) *= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) = _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    IF (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      IF (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 IF (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_stencil()
#line 260
 } if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "G");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "dry");
 }
 if (_first_call) {
 if (pdt != _pdt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "pdt");
 }
 if (_first_call) {
 if (CFL != _CFL)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "CFL");
 }
 if (_first_call) {
 if (hydrostatic != _hydrostatic)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "hydrostatic");
 }
 if (_first_call) {
 if (CFL_H != _CFL_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 210, "CFL_H");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 260

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction (min:dtmax)) {

#line 210

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 210
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,i-1,0,0), val(h,i,0,0), val(h,i+1,0,0))/Delta :
 (val(h,i+1,0,0) - val(h,i-1,0,0))/(2.*Delta);
      hff = val(h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      if (fabs(val(hu.x,0,0,0)) > um)
 um = fabs(val(hu.x,0,0,0));

      val(hu.x,0,0,0) *= val(hf.x,0,0,0);
      val(ha.x,0,0,0) = val(hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,0,i-1,0), val(h,0,i,0), val(h,0,i+1,0))/Delta :
 (val(h,0,i+1,0) - val(h,0,i-1,0))/(2.*Delta);
      hff = val(h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      if (fabs(val(hu.y,0,0,0)) > um)
 um = fabs(val(hu.y,0,0,0));

      val(hu.y,0,0,0) *= val(hf.y,0,0,0);
      val(ha.y,0,0,0) = val(hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_generic()
#line 260
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 210
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,i-1,0,0), val(h,i,0,0), val(h,i+1,0,0))/Delta :
 (val(h,i+1,0,0) - val(h,i-1,0,0))/(2.*Delta);
      hff = val(h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      if (fabs(val(hu.x,0,0,0)) > um)
 um = fabs(val(hu.x,0,0,0));

      val(hu.x,0,0,0) *= val(hf.x,0,0,0);
      val(ha.x,0,0,0) = val(hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,0,i-1,0), val(h,0,i,0), val(h,0,i+1,0))/Delta :
 (val(h,0,i+1,0) - val(h,0,i-1,0))/(2.*Delta);
      hff = val(h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      if (fabs(val(hu.y,0,0,0)) > um)
 um = fabs(val(hu.y,0,0,0));

      val(hu.y,0,0,0) *= val(hf.y,0,0,0);
      val(ha.y,0,0,0) = val(hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_generic()
#line 260
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 210
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,i-1,0,0), val(h,i,0,0), val(h,i+1,0,0))/Delta :
 (val(h,i+1,0,0) - val(h,i-1,0,0))/(2.*Delta);
      hff = val(h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      if (fabs(val(hu.x,0,0,0)) > um)
 um = fabs(val(hu.x,0,0,0));

      val(hu.x,0,0,0) *= val(hf.x,0,0,0);
      val(ha.x,0,0,0) = val(hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,0,i-1,0), val(h,0,i,0), val(h,0,i+1,0))/Delta :
 (val(h,0,i+1,0) - val(h,0,i-1,0))/(2.*Delta);
      hff = val(h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      if (fabs(val(hu.y,0,0,0)) > um)
 um = fabs(val(hu.y,0,0,0));

      val(hu.y,0,0,0) *= val(hf.y,0,0,0);
      val(ha.y,0,0,0) = val(hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_generic()
#line 260
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 210
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.x,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.x,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,i-1,0,0), val(h,i,0,0), val(h,i+1,0,0))/Delta :
 (val(h,i+1,0,0) - val(h,i-1,0,0))/(2.*Delta);
      hff = val(h,i,0,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*hff;





      if (fabs(val(hu.x,0,0,0)) > um)
 um = fabs(val(hu.x,0,0,0));

      val(hu.x,0,0,0) *= val(hf.x,0,0,0);
      val(ha.x,0,0,0) = val(hf.x,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,-1,0,0))*Delta/(c*val_fm_x(fm.x,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 210
{

#line 210 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double ax = (G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double H = 0., um = 0.;
     { foreach_block_inner() {





      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      val(hu.y,0,0,0) = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;





      double hff, un = pdt*(val(hu.y,0,0,0) + pdt*ax)/Delta, a = sign(un);
      int i = - (a + 1.)/2.;
      double g = _attribute[h.i].gradient ? _attribute[h.i].gradient (val(h,0,i-1,0), val(h,0,i,0), val(h,0,i+1,0))/Delta :
 (val(h,0,i+1,0) - val(h,0,i-1,0))/(2.*Delta);
      hff = val(h,0,i,0) + a*(1. - a*un)*g*Delta/2.;
      val(hf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*hff;





      if (fabs(val(hu.y,0,0,0)) > um)
 um = fabs(val(hu.y,0,0,0));

      val(hu.y,0,0,0) *= val(hf.y,0,0,0);
      val(ha.y,0,0,0) = val(hf.y,0,0,0)*ax;

      H += hff;
    } end_foreach_block_inner(); }






    if (H > dry) {
      double c = um/CFL + sqrt(G*(hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      if (c > 0.) {
 double dt = min(val_cm(cm,0,0,0), val_cm(cm,0,-1,0))*Delta/(c*val_fm_y(fm.y,0,0,0));
 if (dt < dtmax)
   dtmax = dt;
      }
    }
  } }  }}  end_foreach_face_generic()
#line 260
 end_foreach_face(); }mpi_all_reduce_array (&dtmax, double, MPI_MIN, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 260
 }





  pdt = dt = dtnext (dtmax);
 end_trace("face_fields", "/home/jiarongw/basilisk/src/layered/hydro.h", 267); } return 0; } 







void advect (scalar * tracers, vector hu, vector hf, double dt)
{







   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _CFL = CFL;
 int _nl = nl;
{ double dt = _dt; NOT_UNUSED(dt);
 double CFL = _CFL; NOT_UNUSED(CFL);
 int nl = _nl; NOT_UNUSED(nl);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 284,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 284
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      IF (hul*dt/(Delta*val_cm(cm,-1,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,-1,0,0))
 hul = CFL*_stencil_val(__FILE__,__LINE__,h,-1,0,0)*Delta*val_cm(cm,-1,0,0)/dt;
      IF (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0))
 hul = - CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      IF (hul != _stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) {
 IF (point.l < nl - 1)
   _stencil_val(__FILE__,__LINE__,hu.x,0,0,1) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) - hul;
 IF (nl > 1)
   _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      IF (hul*dt/(Delta*val_cm(cm,0,-1,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,-1,0))
 hul = CFL*_stencil_val(__FILE__,__LINE__,h,0,-1,0)*Delta*val_cm(cm,0,-1,0)/dt;
      IF (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0))
 hul = - CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      IF (hul != _stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) {
 IF (point.l < nl - 1)
   _stencil_val(__FILE__,__LINE__,hu.y,0,0,1) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) - hul;
 IF (nl > 1)
   _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  end_foreach_face_stencil()
#line 307
 }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 284
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      IF (hul*dt/(Delta*val_cm(cm,-1,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,-1,0,0))
 hul = CFL*_stencil_val(__FILE__,__LINE__,h,-1,0,0)*Delta*val_cm(cm,-1,0,0)/dt;
      IF (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0))
 hul = - CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      IF (hul != _stencil_val(__FILE__,__LINE__,hu.x,0,0,0)) {
 IF (point.l < nl - 1)
   _stencil_val(__FILE__,__LINE__,hu.x,0,0,1) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) - hul;
 IF (nl > 1)
   _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      IF (hul*dt/(Delta*val_cm(cm,0,-1,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,-1,0))
 hul = CFL*_stencil_val(__FILE__,__LINE__,h,0,-1,0)*Delta*val_cm(cm,0,-1,0)/dt;
      IF (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0))
 hul = - CFL*_stencil_val(__FILE__,__LINE__,h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      IF (hul != _stencil_val(__FILE__,__LINE__,hu.y,0,0,0)) {
 IF (point.l < nl - 1)
   _stencil_val(__FILE__,__LINE__,hu.y,0,0,1) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) - hul;
 IF (nl > 1)
   _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  end_foreach_face_stencil()
#line 307
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 284, "dt");
 }
 if (_first_call) {
 if (CFL != _CFL)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 284, "CFL");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 284, "nl");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 307

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 284
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = val(hu.x,0,0,0);
      if (hul*dt/(Delta*val_cm(cm,-1,0,0)) > CFL*val(h,-1,0,0))
 hul = CFL*val(h,-1,0,0)*Delta*val_cm(cm,-1,0,0)/dt;
      else if (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*val(h,0,0,0))
 hul = - CFL*val(h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      if (hul != val(hu.x,0,0,0)) {
 if (point.l < nl - 1)
   val(hu.x,0,0,1) += val(hu.x,0,0,0) - hul;
 else if (nl > 1)
   fprintf (ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 val(hu.x,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = val(hu.y,0,0,0);
      if (hul*dt/(Delta*val_cm(cm,0,-1,0)) > CFL*val(h,0,-1,0))
 hul = CFL*val(h,0,-1,0)*Delta*val_cm(cm,0,-1,0)/dt;
      else if (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*val(h,0,0,0))
 hul = - CFL*val(h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      if (hul != val(hu.y,0,0,0)) {
 if (point.l < nl - 1)
   val(hu.y,0,0,1) += val(hu.y,0,0,0) - hul;
 else if (nl > 1)
   fprintf (ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 val(hu.y,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  end_foreach_face_generic()
#line 307
 end_foreach_face(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 284
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = val(hu.x,0,0,0);
      if (hul*dt/(Delta*val_cm(cm,-1,0,0)) > CFL*val(h,-1,0,0))
 hul = CFL*val(h,-1,0,0)*Delta*val_cm(cm,-1,0,0)/dt;
      else if (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*val(h,0,0,0))
 hul = - CFL*val(h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      if (hul != val(hu.x,0,0,0)) {
 if (point.l < nl - 1)
   val(hu.x,0,0,1) += val(hu.x,0,0,0) - hul;
 else if (nl > 1)
   fprintf (ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 val(hu.x,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 284
{

#line 284 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      double hul = val(hu.y,0,0,0);
      if (hul*dt/(Delta*val_cm(cm,0,-1,0)) > CFL*val(h,0,-1,0))
 hul = CFL*val(h,0,-1,0)*Delta*val_cm(cm,0,-1,0)/dt;
      else if (- hul*dt/(Delta*val_cm(cm,0,0,0)) > CFL*val(h,0,0,0))
 hul = - CFL*val(h,0,0,0)*Delta*val_cm(cm,0,0,0)/dt;
#line 299 "/home/jiarongw/basilisk/src/layered/hydro.h"
      if (hul != val(hu.y,0,0,0)) {
 if (point.l < nl - 1)
   val(hu.y,0,0,1) += val(hu.y,0,0,0) - hul;
 else if (nl > 1)
   fprintf (ferr, "warning: could not conserve barotropic flux "
     "at %g,%g,%d\n", x, y, point.l);
 val(hu.y,0,0,0) = hul;
      }
    } end_foreach_block_inner(); } }  }}  end_foreach_face_generic()
#line 307
 end_foreach_face(); } }

  vector flux= new_face_vector("flux");
   { foreach_block() {





    strongif (tracers) for (scalar s = *tracers, *_i72 = tracers; ((scalar *)&s)->i >= 0; s = *++_i72) {
       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _dry = dry;
{ double dt = _dt; NOT_UNUSED(dt);
 double dry = _dry; NOT_UNUSED(dry);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 317,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 double un = dt*_stencil_val(__FILE__,__LINE__,hu.x,0,0,0)/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*Delta + dry), a = sign(un);
 int i = -(a + 1.)/2.;
 double g = _attribute[s.i].gradient ?
   _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,i-1,0,0), _stencil_val(__FILE__,__LINE__,s,i,0,0), _stencil_val(__FILE__,__LINE__,s,i+1,0,0))/Delta :
   (_stencil_val(__FILE__,__LINE__,s,i+1,0,0) - _stencil_val(__FILE__,__LINE__,s,i-1,0,0))/(2.*Delta);
 double s2 = _stencil_val(__FILE__,__LINE__,s,i,0,0) + a*(1. - a*un)*g*Delta/2.;


 IF (_stencil_val(__FILE__,__LINE__,hf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,i,1,0) > dry) {
   double vn = (_stencil_val(__FILE__,__LINE__,hu.y,i,0,0) + _stencil_val(__FILE__,__LINE__,hu.y,i,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,i,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,i,1,0));
   double syy = (_attribute[s.i].gradient ? _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,i,-1,0), _stencil_val(__FILE__,__LINE__,s,i,0,0), _stencil_val(__FILE__,__LINE__,s,i,1,0)) :
   vn < 0. ? _stencil_val(__FILE__,__LINE__,s,i,1,0) - _stencil_val(__FILE__,__LINE__,s,i,0,0) : _stencil_val(__FILE__,__LINE__,s,i,0,0) - _stencil_val(__FILE__,__LINE__,s,i,-1,0));
   s2 -= dt*vn*syy/(2.*Delta);
 }


 _stencil_val(__FILE__,__LINE__,flux.x,0,0,0) = s2*_stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 double un = dt*_stencil_val(__FILE__,__LINE__,hu.y,0,0,0)/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*Delta + dry), a = sign(un);
 int i = -(a + 1.)/2.;
 double g = _attribute[s.i].gradient ?
   _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,0,i-1,0), _stencil_val(__FILE__,__LINE__,s,0,i,0), _stencil_val(__FILE__,__LINE__,s,0,i+1,0))/Delta :
   (_stencil_val(__FILE__,__LINE__,s,0,i+1,0) - _stencil_val(__FILE__,__LINE__,s,0,i-1,0))/(2.*Delta);
 double s2 = _stencil_val(__FILE__,__LINE__,s,0,i,0) + a*(1. - a*un)*g*Delta/2.;


 IF (_stencil_val(__FILE__,__LINE__,hf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,i,0) > dry) {
   double vn = (_stencil_val(__FILE__,__LINE__,hu.x,0,i,0) + _stencil_val(__FILE__,__LINE__,hu.x,1,i,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,i,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,i,0));
   double syy = (_attribute[s.i].gradient ? _attribute[s.i].gradient (_stencil_val(__FILE__,__LINE__,s,-1,i,0), _stencil_val(__FILE__,__LINE__,s,0,i,0), _stencil_val(__FILE__,__LINE__,s,1,i,0)) :
   vn < 0. ? _stencil_val(__FILE__,__LINE__,s,1,i,0) - _stencil_val(__FILE__,__LINE__,s,0,i,0) : _stencil_val(__FILE__,__LINE__,s,0,i,0) - _stencil_val(__FILE__,__LINE__,s,-1,i,0));
   s2 -= dt*vn*syy/(2.*Delta);
 }


 _stencil_val(__FILE__,__LINE__,flux.y,0,0,0) = s2*_stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      } }  }}  end_foreach_face_stencil()
#line 335
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 317, "dt");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 317, "dry");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 335
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 double un = dt*val(hu.x,0,0,0)/(val(hf.x,0,0,0)*Delta + dry), a = sign(un);
 int i = -(a + 1.)/2.;
 double g = _attribute[s.i].gradient ?
   _attribute[s.i].gradient (val(s,i-1,0,0), val(s,i,0,0), val(s,i+1,0,0))/Delta :
   (val(s,i+1,0,0) - val(s,i-1,0,0))/(2.*Delta);
 double s2 = val(s,i,0,0) + a*(1. - a*un)*g*Delta/2.;


 if (val(hf.y,i,0,0) + val(hf.y,i,1,0) > dry) {
   double vn = (val(hu.y,i,0,0) + val(hu.y,i,1,0))/(val(hf.y,i,0,0) + val(hf.y,i,1,0));
   double syy = (_attribute[s.i].gradient ? _attribute[s.i].gradient (val(s,i,-1,0), val(s,i,0,0), val(s,i,1,0)) :
   vn < 0. ? val(s,i,1,0) - val(s,i,0,0) : val(s,i,0,0) - val(s,i,-1,0));
   s2 -= dt*vn*syy/(2.*Delta);
 }


 val(flux.x,0,0,0) = s2*val(hu.x,0,0,0);
      } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 317
{

#line 317 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 double un = dt*val(hu.y,0,0,0)/(val(hf.y,0,0,0)*Delta + dry), a = sign(un);
 int i = -(a + 1.)/2.;
 double g = _attribute[s.i].gradient ?
   _attribute[s.i].gradient (val(s,0,i-1,0), val(s,0,i,0), val(s,0,i+1,0))/Delta :
   (val(s,0,i+1,0) - val(s,0,i-1,0))/(2.*Delta);
 double s2 = val(s,0,i,0) + a*(1. - a*un)*g*Delta/2.;


 if (val(hf.x,0,i,0) + val(hf.x,1,i,0) > dry) {
   double vn = (val(hu.x,0,i,0) + val(hu.x,1,i,0))/(val(hf.x,0,i,0) + val(hf.x,1,i,0));
   double syy = (_attribute[s.i].gradient ? _attribute[s.i].gradient (val(s,-1,i,0), val(s,0,i,0), val(s,1,i,0)) :
   vn < 0. ? val(s,1,i,0) - val(s,0,i,0) : val(s,0,i,0) - val(s,-1,i,0));
   s2 -= dt*vn*syy/(2.*Delta);
 }


 val(flux.y,0,0,0) = s2*val(hu.y,0,0,0);
      } }  }}  end_foreach_face_generic()
#line 335
 end_foreach_face(); }





       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 341,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 341
foreach_stencil(){

#line 341 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 _stencil_val(__FILE__,__LINE__,s,0,0,0) *= _stencil_val(__FILE__,__LINE__,h,0,0,0);
 {
#line 343

   _stencil_val(__FILE__,__LINE__,s,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux.x,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 343

   _stencil_val(__FILE__,__LINE__,s,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux.y,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 341
foreach_stencil(){

#line 341 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 _stencil_val(__FILE__,__LINE__,s,0,0,0) *= _stencil_val(__FILE__,__LINE__,h,0,0,0);
 {
#line 343

   _stencil_val(__FILE__,__LINE__,s,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux.x,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 343

   _stencil_val(__FILE__,__LINE__,s,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,flux.y,0,0,0) - _stencil_val(__FILE__,__LINE__,flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 341, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 345

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 341
foreach(){

#line 341 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 val(s,0,0,0) *= val(h,0,0,0);
 {
#line 343

   val(s,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 343

   val(s,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 341
foreach(){

#line 341 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
 val(s,0,0,0) *= val(h,0,0,0);
 {
#line 343

   val(s,0,0,0) += dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 343

   val(s,0,0,0) += dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      } } end_foreach(); } }
    }
#line 358 "/home/jiarongw/basilisk/src/layered/hydro.h"
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 int __layer = _layer;
 double _t = t;
 double _dry = dry;
{ double dt = _dt; NOT_UNUSED(dt);
 int _layer = __layer; NOT_UNUSED(_layer);
 double t = _t; NOT_UNUSED(t);
 double dry = _dry; NOT_UNUSED(dry);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 358,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 358
foreach_stencil(){

#line 358 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
      double h1 = _stencil_val(__FILE__,__LINE__,h,0,0,0);
      {
#line 360

 h1 += dt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 360

 h1 += dt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      IF (h1 < - 1e-12)
 _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: h1 = %g < - 1e-12 at %g,%g,%d,%g\n",
   h1, x, y, _layer, t);
      _stencil_val(__FILE__,__LINE__,h,0,0,0) = fmax(h1, 0.);
      IF (h1 < dry) {
 strongif (tracers) for (scalar f = *tracers, *_i73 = tracers; ((scalar *)&f)->i >= 0; f = *++_i73)
   _stencil_val(__FILE__,__LINE__,f,0,0,0) = 0.;
      }
      
 strongif (tracers) for (scalar f = *tracers, *_i74 = tracers; ((scalar *)&f)->i >= 0; f = *++_i74)
   _stencil_val(__FILE__,__LINE__,f,0,0,0) /= h1;
    } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 358
foreach_stencil(){

#line 358 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
      double h1 = _stencil_val(__FILE__,__LINE__,h,0,0,0);
      {
#line 360

 h1 += dt*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 360

 h1 += dt*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      IF (h1 < - 1e-12)
 _stencil_fprintf (__FILE__,__LINE__,ferr, "warning: h1 = %g < - 1e-12 at %g,%g,%d,%g\n",
   h1, x, y, _layer, t);
      _stencil_val(__FILE__,__LINE__,h,0,0,0) = fmax(h1, 0.);
      IF (h1 < dry) {
 strongif (tracers) for (scalar f = *tracers, *_i73 = tracers; ((scalar *)&f)->i >= 0; f = *++_i73)
   _stencil_val(__FILE__,__LINE__,f,0,0,0) = 0.;
      }
      
 strongif (tracers) for (scalar f = *tracers, *_i74 = tracers; ((scalar *)&f)->i >= 0; f = *++_i74)
   _stencil_val(__FILE__,__LINE__,f,0,0,0) /= h1;
    } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 358, "dt");
 }
 if (_first_call) {
 if (_layer != __layer)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 358, "_layer");
 }
 if (_first_call) {
 if (t != _t)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 358, "t");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 358, "dry");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 373

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 358
foreach(){

#line 358 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
      double h1 = val(h,0,0,0);
      {
#line 360

 h1 += dt*(val(hu.x,0,0,0) - val(hu.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 360

 h1 += dt*(val(hu.y,0,0,0) - val(hu.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      if (h1 < - 1e-12)
 fprintf (ferr, "warning: h1 = %g < - 1e-12 at %g,%g,%d,%g\n",
   h1, x, y, _layer, t);
      val(h,0,0,0) = fmax(h1, 0.);
      if (h1 < dry) {
 strongif (tracers) for (scalar f = *tracers, *_i73 = tracers; ((scalar *)&f)->i >= 0; f = *++_i73)
   val(f,0,0,0) = 0.;
      }
      else
 strongif (tracers) for (scalar f = *tracers, *_i74 = tracers; ((scalar *)&f)->i >= 0; f = *++_i74)
   val(f,0,0,0) /= h1;
    } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 358
foreach(){

#line 358 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
      double h1 = val(h,0,0,0);
      {
#line 360

 h1 += dt*(val(hu.x,0,0,0) - val(hu.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 360

 h1 += dt*(val(hu.y,0,0,0) - val(hu.y,0,1,0))/(Delta*val_cm(cm,0,0,0));}
      if (h1 < - 1e-12)
 fprintf (ferr, "warning: h1 = %g < - 1e-12 at %g,%g,%d,%g\n",
   h1, x, y, _layer, t);
      val(h,0,0,0) = fmax(h1, 0.);
      if (h1 < dry) {
 strongif (tracers) for (scalar f = *tracers, *_i73 = tracers; ((scalar *)&f)->i >= 0; f = *++_i73)
   val(f,0,0,0) = 0.;
      }
      else
 strongif (tracers) for (scalar f = *tracers, *_i74 = tracers; ((scalar *)&f)->i >= 0; f = *++_i74)
   val(f,0,0,0) /= h1;
    } } end_foreach(); } }
  } end_foreach_block(); }
 delete (((scalar []){flux.x,flux.y,{-1}})); }






static int half_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int half_advection (const int i, const double t, Event * _ev) { trace ("half_advection", "/home/jiarongw/basilisk/src/layered/hydro.h", 382); ; end_trace("half_advection", "/home/jiarongw/basilisk/src/layered/hydro.h", 382);  return 0; } 




#line 1 "layered/diffusion.h"
#line 1 "/home/jiarongw/basilisk/src/layered/diffusion.h"
#line 31 "/home/jiarongw/basilisk/src/layered/diffusion.h"
void vertical_diffusion (Point point, scalar h, scalar s, double dt, double D,
    double dst, double s_b, double lambda_b)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 33 "/home/jiarongw/basilisk/src/layered/diffusion.h"

  double a[nl], b[nl], c[nl], rhs[nl];




   { foreach_block()
    rhs[_layer] = val(s,0,0,0)*val(h,0,0,0); end_foreach_block(); }
#line 56 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(val(h,0,0,l-1) + val(h,0,0,l));
    c[l] = - 2.*D*dt/(val(h,0,0,l) + val(h,0,0,l+1));
    b[l] = val(h,0,0,l) - a[l] - c[l];
  }
#line 79 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  a[nl-1] = - 2.*D*dt/(val(h,0,0,nl-2) + val(h,0,0,nl-1));
  b[nl-1] = val(h,0,0,nl-1) - a[nl-1];
  rhs[nl-1] += D*dt*dst;
#line 99 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  double den = val(h,0,0,0)*sq(val(h,0,0,0) + val(h,0,0,1))
    + 2.*lambda_b*(3.*val(h,0,0,0)*val(h,0,0,1) + 2.*sq(val(h,0,0,0)) + sq(val(h,0,0,1)));
  b[0] = val(h,0,0,0) + 2.*dt*D*(1./(val(h,0,0,0) + val(h,0,0,1)) +
     (sq(val(h,0,0,1)) + 3.*val(h,0,0,0)*val(h,0,0,1) + 3.*sq(val(h,0,0,0)))/den);
  c[0] = - 2.*dt*D*(1./(val(h,0,0,0) + val(h,0,0,1)) + sq(val(h,0,0,0))/den);
  rhs[0] += 2.*dt*D*s_b*(sq(val(h,0,0,1)) + 3.*val(h,0,0,0)*val(h,0,0,1) + 2.*sq(val(h,0,0,0)))/den;

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*val(h,0,0,0) - D*dt) * dst;
  }





  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  val(s,0,0,nl-1) = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    val(s,0,0,l) = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];

#if _call_vertical_diffusion
}
#define _IN_STENCIL 1

#line 31
static void _vertical_diffusion (Point point, scalar h, scalar s, double dt, double D,
    double dst, double s_b, double lambda_b)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 33 "/home/jiarongw/basilisk/src/layered/diffusion.h"

  double a[nl], b[nl], c[nl], rhs[nl];




   { foreach_block()
    rhs[_layer] = _stencil_val(__FILE__,__LINE__,s,0,0,0)*_stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block(); }
#line 56 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(_stencil_val(__FILE__,__LINE__,h,0,0,l-1) + _stencil_val(__FILE__,__LINE__,h,0,0,l));
    c[l] = - 2.*D*dt/(_stencil_val(__FILE__,__LINE__,h,0,0,l) + _stencil_val(__FILE__,__LINE__,h,0,0,l+1));
    b[l] = _stencil_val(__FILE__,__LINE__,h,0,0,l) - a[l] - c[l];
  }
#line 79 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  a[nl-1] = - 2.*D*dt/(_stencil_val(__FILE__,__LINE__,h,0,0,nl-2) + _stencil_val(__FILE__,__LINE__,h,0,0,nl-1));
  b[nl-1] = _stencil_val(__FILE__,__LINE__,h,0,0,nl-1) - a[nl-1];
  rhs[nl-1] += D*dt*dst;
#line 99 "/home/jiarongw/basilisk/src/layered/diffusion.h"
  double den = _stencil_val(__FILE__,__LINE__,h,0,0,0)*sq(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0,1))
    + 2.*lambda_b*(3.*_stencil_val(__FILE__,__LINE__,h,0,0,0)*_stencil_val(__FILE__,__LINE__,h,0,0,1) + 2.*sq(_stencil_val(__FILE__,__LINE__,h,0,0,0)) + sq(_stencil_val(__FILE__,__LINE__,h,0,0,1)));
  b[0] = _stencil_val(__FILE__,__LINE__,h,0,0,0) + 2.*dt*D*(1./(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0,1)) +
     (sq(_stencil_val(__FILE__,__LINE__,h,0,0,1)) + 3.*_stencil_val(__FILE__,__LINE__,h,0,0,0)*_stencil_val(__FILE__,__LINE__,h,0,0,1) + 3.*sq(_stencil_val(__FILE__,__LINE__,h,0,0,0)))/den);
  c[0] = - 2.*dt*D*(1./(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0,1)) + sq(_stencil_val(__FILE__,__LINE__,h,0,0,0))/den);
  rhs[0] += 2.*dt*D*s_b*(sq(_stencil_val(__FILE__,__LINE__,h,0,0,1)) + 3.*_stencil_val(__FILE__,__LINE__,h,0,0,0)*_stencil_val(__FILE__,__LINE__,h,0,0,1) + 2.*sq(_stencil_val(__FILE__,__LINE__,h,0,0,0)))/den;

  IF (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*_stencil_val(__FILE__,__LINE__,h,0,0,0) - D*dt) * dst;
  }





  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  _stencil_val(__FILE__,__LINE__,s,0,0,nl-1) = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    _stencil_val(__FILE__,__LINE__,s,0,0,l) = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];

#undef _IN_STENCIL

#endif

#line 123
}
#line 133 "/home/jiarongw/basilisk/src/layered/diffusion.h"
double nu = 0.;
 vector lambda_b = {{_NVARMAX + 0},{_NVARMAX + 1}}, dut = {{_NVARMAX + 0},{_NVARMAX + 1}}, u_b = {{_NVARMAX + 0},{_NVARMAX + 1}};
#line 144 "/home/jiarongw/basilisk/src/layered/diffusion.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/jiarongw/basilisk/src/layered/diffusion.h", 144); 
{
  if (nu > 0.) {
     { 
#define vertical_diffusion _vertical_diffusion
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _dry = dry;
 double _nu = nu;
{ double dt = _dt; NOT_UNUSED(dt);
 double dry = _dry; NOT_UNUSED(dry);
 double nu = _nu; NOT_UNUSED(nu);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/diffusion.h", .line = 147,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(dut.x) && !is_constant(u_b.x) && !is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (is_constant(dut.x) && !is_constant(u_b.x) && !is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (!is_constant(dut.x) && is_constant(u_b.x) && !is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (is_constant(dut.x) && is_constant(u_b.x) && !is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (!is_constant(dut.x) && !is_constant(u_b.x) && is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (is_constant(dut.x) && !is_constant(u_b.x) && is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (!is_constant(dut.x) && is_constant(u_b.x) && is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); }
strongif (is_constant(dut.x) && is_constant(u_b.x) && is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach_stencil(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 149

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   _stencil_val(__FILE__,__LINE__,u.x,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 155

   _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/diffusion.h", 147, "dt");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/diffusion.h", 147, "dry");
 }
 if (_first_call) {
 if (nu != _nu)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/diffusion.h", 147, "nu");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef vertical_diffusion
#line 157

strongif (!is_constant(dut.x) && !is_constant(u_b.x) && !is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) val(a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) coarse(a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) val(a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) coarse(a,i,j,k)
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (is_constant(dut.x) && !is_constant(u_b.x) && !is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) coarse(a,i,j,k)
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (!is_constant(dut.x) && is_constant(u_b.x) && !is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) val(a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) coarse(a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) val(a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) coarse(a,i,j,k)
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (is_constant(dut.x) && is_constant(u_b.x) && !is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) coarse(a,i,j,k)
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (!is_constant(dut.x) && !is_constant(u_b.x) && is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) val(a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) coarse(a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) val(a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (is_constant(dut.x) && !is_constant(u_b.x) && is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) coarse(a,i,j,k)
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) val(a,i,j,k)
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (!is_constant(dut.x) && is_constant(u_b.x) && is_constant(lambda_b.x)) {
#undef val_dut_x
#define val_dut_x(a,i,j,k) val(a,i,j,k)
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) coarse(a,i,j,k)
#undef val_dut_y
#define val_dut_y(a,i,j,k) val(a,i,j,k)
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); }
strongif (is_constant(dut.x) && is_constant(u_b.x) && is_constant(lambda_b.x)) {
const struct { double x, y; } _const_dut = {_constant[dut.x.i -_NVARMAX], _constant[dut.y.i - _NVARMAX]};
NOT_UNUSED(_const_dut);
#undef val_dut_x
#define val_dut_x(a,i,j,k) _const_dut.x
#undef fine_dut_x
#define fine_dut_x(a,i,j,k) _const_dut.x
#undef coarse_dut_x
#define coarse_dut_x(a,i,j,k) _const_dut.x
#undef val_dut_y
#define val_dut_y(a,i,j,k) _const_dut.y
#undef fine_dut_y
#define fine_dut_y(a,i,j,k) _const_dut.y
#undef coarse_dut_y
#define coarse_dut_y(a,i,j,k) _const_dut.y
const struct { double x, y; } _const_u_b = {_constant[u_b.x.i -_NVARMAX], _constant[u_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_u_b);
#undef val_u_b_x
#define val_u_b_x(a,i,j,k) _const_u_b.x
#undef fine_u_b_x
#define fine_u_b_x(a,i,j,k) _const_u_b.x
#undef coarse_u_b_x
#define coarse_u_b_x(a,i,j,k) _const_u_b.x
#undef val_u_b_y
#define val_u_b_y(a,i,j,k) _const_u_b.y
#undef fine_u_b_y
#define fine_u_b_y(a,i,j,k) _const_u_b.y
#undef coarse_u_b_y
#define coarse_u_b_y(a,i,j,k) _const_u_b.y
const struct { double x, y; } _const_lambda_b = {_constant[lambda_b.x.i -_NVARMAX], _constant[lambda_b.y.i - _NVARMAX]};
NOT_UNUSED(_const_lambda_b);
#undef val_lambda_b_x
#define val_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef fine_lambda_b_x
#define fine_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef coarse_lambda_b_x
#define coarse_lambda_b_x(a,i,j,k) _const_lambda_b.x
#undef val_lambda_b_y
#define val_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef fine_lambda_b_y
#define fine_lambda_b_y(a,i,j,k) _const_lambda_b.y
#undef coarse_lambda_b_y
#define coarse_lambda_b_y(a,i,j,k) _const_lambda_b.y
#line 147
foreach(){

#line 147 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
       { foreach_block_inner()
 {
#line 149

   val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 149

   val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
      {
#line 151

 vertical_diffusion (point, h, u.x, dt, nu,
       val_dut_x(dut.x,0,0,0), val_u_b_x(u_b.x,0,0,0), val_lambda_b_x(lambda_b.x,0,0,0));
#line 151

 vertical_diffusion (point, h, u.y, dt, nu,
       val_dut_y(dut.y,0,0,0), val_u_b_y(u_b.y,0,0,0), val_lambda_b_y(lambda_b.y,0,0,0));}
       { foreach_block_inner()
 {
#line 155

   val(u.x,0,0,0) -= dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 155

   val(u.y,0,0,0) -= dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);} end_foreach_block_inner(); }
    } } end_foreach(); } }
  }
 end_trace("viscous_term", "/home/jiarongw/basilisk/src/layered/diffusion.h", 159); } return 0; } 
#line 173 "/home/jiarongw/basilisk/src/layered/diffusion.h"
void horizontal_diffusion (scalar * list, double D, double dt)
{
  if (D > 0.) {
    scalar * d2sl = list_clone (list);
     { foreach_block() {
       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/diffusion.h", .line = 178,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 178
foreach_stencil(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,-1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));}
   _stencil_val(__FILE__,__LINE__,d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 178
foreach_stencil(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,-1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));}
   _stencil_val(__FILE__,__LINE__,d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 178
foreach_stencil(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,-1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));}
   _stencil_val(__FILE__,__LINE__,d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 178
foreach_stencil(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,-1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,1,0,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));
#line 182

     a += (_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,-1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)) +
    _stencil_val(__FILE__,__LINE__,hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(_stencil_val(__FILE__,__LINE__,s,0,1,0) - _stencil_val(__FILE__,__LINE__,s,0,0,0)));}
   _stencil_val(__FILE__,__LINE__,d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 187

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 178
foreach(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (val(hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(val(s,-1,0,0) - val(s,0,0,0)) +
    val(hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(val(s,1,0,0) - val(s,0,0,0)));
#line 182

     a += (val(hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(val(s,0,-1,0) - val(s,0,0,0)) +
    val(hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(val(s,0,1,0) - val(s,0,0,0)));}
   val(d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 178
foreach(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (val(hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(val(s,-1,0,0) - val(s,0,0,0)) +
    val(hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(val(s,1,0,0) - val(s,0,0,0)));
#line 182

     a += (val(hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(val(s,0,-1,0) - val(s,0,0,0)) +
    val(hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(val(s,0,1,0) - val(s,0,0,0)));}
   val(d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 178
foreach(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (val(hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(val(s,-1,0,0) - val(s,0,0,0)) +
    val(hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(val(s,1,0,0) - val(s,0,0,0)));
#line 182

     a += (val(hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(val(s,0,-1,0) - val(s,0,0,0)) +
    val(hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(val(s,0,1,0) - val(s,0,0,0)));}
   val(d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 178
foreach(){

#line 178 "/home/jiarongw/basilisk/src/layered/diffusion.h"
 {
 scalar s, d2s;
 scalar * _i4 = list; scalar * _i5 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i4, d2s = *++_i5) {
   double a = 0.;
   {
#line 182

     a += (val(hf.x,0,0,0)*val_fm_x(fm.x,0,0,0)/(val_cm(cm,-1,0,0) + val_cm(cm,0,0,0))*(val(s,-1,0,0) - val(s,0,0,0)) +
    val(hf.x,1,0,0)*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,0,0,0))*(val(s,1,0,0) - val(s,0,0,0)));
#line 182

     a += (val(hf.y,0,0,0)*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,-1,0) + val_cm(cm,0,0,0))*(val(s,0,-1,0) - val(s,0,0,0)) +
    val(hf.y,0,1,0)*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,0,0))*(val(s,0,1,0) - val(s,0,0,0)));}
   val(d2s,0,0,0) = 2.*a/(val_cm(cm,0,0,0)*sq(Delta));
        }
      } } end_foreach(); } }
       { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dry = dry;
 double _dt = dt;
{ double dry = _dry; NOT_UNUSED(dry);
 double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/diffusion.h", .line = 188,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 188 "/home/jiarongw/basilisk/src/layered/diffusion.h"

 IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) > dry) {
   scalar s, d2s;
   scalar * _i6 = list; scalar * _i7 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i6, d2s = *++_i7)
     _stencil_val(__FILE__,__LINE__,s,0,0,0) += dt*D*_stencil_val(__FILE__,__LINE__,d2s,0,0,0)/_stencil_val(__FILE__,__LINE__,h,0,0,0);
 } } end_foreach_stencil(); if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/diffusion.h", 188, "dry");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/diffusion.h", 188, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 193
foreach(){

#line 188 "/home/jiarongw/basilisk/src/layered/diffusion.h"

 if (val(h,0,0,0) > dry) {
   scalar s, d2s;
   scalar * _i6 = list; scalar * _i7 = d2sl; strongif (list) for (s = *list, d2s = *d2sl; ((scalar *)&s)->i >= 0; s = *++_i6, d2s = *++_i7)
     val(s,0,0,0) += dt*D*val(d2s,0,0,0)/val(h,0,0,0);
 } } end_foreach(); }
    } end_foreach_block(); }
    delete (d2sl);
    pfree (d2sl,__func__,__FILE__,__LINE__);
  }
}
#line 388 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 396 "/home/jiarongw/basilisk/src/layered/hydro.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/jiarongw/basilisk/src/layered/hydro.h", 396); ; end_trace("acceleration", "/home/jiarongw/basilisk/src/layered/hydro.h", 396);  return 0; } 

static int pressure_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int pressure (const int i, const double t, Event * _ev) { trace ("pressure", "/home/jiarongw/basilisk/src/layered/hydro.h", 398); 
{




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 404,
    .each = "foreach_face", .first = _first_call
  };
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 404
{

#line 404 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner()
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0); end_foreach_block_inner(); }; }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 404
{

#line 404 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner()
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0); end_foreach_block_inner(); }; }  }}  end_foreach_face_stencil()
#line 406
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 404, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 406
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 404
{

#line 404 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner()
      val(hu.x,0,0,0) += dt*val(ha.x,0,0,0); end_foreach_block_inner(); }; }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 404
{

#line 404 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner()
      val(hu.y,0,0,0) += dt*val(ha.y,0,0,0); end_foreach_block_inner(); }; }  }}  end_foreach_face_generic()
#line 406
 end_foreach_face(); }




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _dry = dry;
{ double dt = _dt; NOT_UNUSED(dt);
 double dry = _dry; NOT_UNUSED(dry);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 411,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 411
foreach_stencil(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 413

 _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = _stencil_val(__FILE__,__LINE__,u.x,0,0,0), uy = _stencil_val(__FILE__,__LINE__,u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*fG*uy;
      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 411
foreach_stencil(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 413

 _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = _stencil_val(__FILE__,__LINE__,u.x,0,0,0), uy = _stencil_val(__FILE__,__LINE__,u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*fG*uy;
      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 411
foreach_stencil(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 413

 _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = _stencil_val(__FILE__,__LINE__,u.x,0,0,0), uy = _stencil_val(__FILE__,__LINE__,u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*fG*uy;
      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 411
foreach_stencil(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry);
#line 413

 _stencil_val(__FILE__,__LINE__,u.y,0,0,0) += dt*(_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) + _stencil_val(__FILE__,__LINE__,ha.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = _stencil_val(__FILE__,__LINE__,u.x,0,0,0), uy = _stencil_val(__FILE__,__LINE__,u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      _stencil_val(__FILE__,__LINE__,u.x,0,0,0) += dt*fG*uy;
      _stencil_val(__FILE__,__LINE__,u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 411, "dt");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 411, "dry");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 424

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 411
foreach(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 413

 val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = val(u.x,0,0,0), uy = val(u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      val(u.x,0,0,0) += dt*fG*uy;
      val(u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 411
foreach(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 413

 val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = val(u.x,0,0,0), uy = val(u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      val(u.x,0,0,0) += dt*fG*uy;
      val(u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 411
foreach(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 413

 val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = val(u.x,0,0,0), uy = val(u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      val(u.x,0,0,0) += dt*fG*uy;
      val(u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 411
foreach(){

#line 411 "/home/jiarongw/basilisk/src/layered/hydro.h"

     { foreach_block_inner() {
      {
#line 413

 val(u.x,0,0,0) += dt*(val(ha.x,0,0,0) + val(ha.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry);
#line 413

 val(u.y,0,0,0) += dt*(val(ha.y,0,0,0) + val(ha.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry);}


      double dmdl = (val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double dmdt = (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))/(val_cm(cm,0,0,0)*Delta);
      double ux = val(u.x,0,0,0), uy = val(u.y,0,0,0);
      double fG = uy*dmdl - ux*dmdt;
      val(u.x,0,0,0) += dt*fG*uy;
      val(u.y,0,0,0) -= dt*fG*ux;

    } end_foreach_block_inner(); } } end_foreach(); } }
  delete ((scalar *)((vector []){{ha.x,ha.y},{{-1},{-1}}}));





  advect (tracers, hu, hf, dt);
 end_trace("pressure", "/home/jiarongw/basilisk/src/layered/hydro.h", 432); } return 0; } 





static int update_eta_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int update_eta (const int i, const double t, Event * _ev) { trace ("update_eta", "/home/jiarongw/basilisk/src/layered/hydro.h", 438); 
{
  delete ((scalar *)((vector []){{hu.x,hu.y},{hf.x,hf.y},{{-1},{-1}}}));
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 441,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 441 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double etap = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner()
      etap += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,eta,0,0,0) = etap;
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 446
foreach(){

#line 441 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double etap = val(zb,0,0,0);
     { foreach_block_inner()
      etap += val(h,0,0,0); end_foreach_block_inner(); }
    val(eta,0,0,0) = etap;
  } } end_foreach(); }
 end_trace("update_eta", "/home/jiarongw/basilisk/src/layered/hydro.h", 447); } return 0; } 






static int remap_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int remap (const int i, const double t, Event * _ev) { trace ("remap", "/home/jiarongw/basilisk/src/layered/hydro.h", 454); ; end_trace("remap", "/home/jiarongw/basilisk/src/layered/hydro.h", 454);  return 0; } 
#line 466 "/home/jiarongw/basilisk/src/layered/hydro.h"
static int cleanup_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup_0 (const int i, const double t, Event * _ev) { trace ("cleanup_0", "/home/jiarongw/basilisk/src/layered/hydro.h", 466); 
{
  delete (((scalar []){eta,h,u.x,u.y,{-1}}));
  pfree (tracers,__func__,__FILE__,__LINE__), tracers = NULL;
 end_trace("cleanup_0", "/home/jiarongw/basilisk/src/layered/hydro.h", 470); } return 0; } 
#line 484 "/home/jiarongw/basilisk/src/layered/hydro.h"
double max_slope = 0.577350269189626;
#line 519 "/home/jiarongw/basilisk/src/layered/hydro.h"
void vertical_velocity (scalar w, vector hu, vector hf)
{
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dry = dry;
{ double dry = _dry; NOT_UNUSED(dry);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/hydro.h", .line = 521,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 521 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double dz = _stencil_val(__FILE__,__LINE__,zb,1,0,0) - _stencil_val(__FILE__,__LINE__,zb,-1,0,0);
    double wm = 0.;
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,w,0,0,0) = wm + (_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hu.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry)*
 (dz + _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0))/(2.*Delta);
      IF (point.l > 0)
 {
#line 528

   _stencil_val(__FILE__,__LINE__,w,0,0,0) -= (_stencil_val(__FILE__,__LINE__,hu.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hu.x,1,0,-1))
   /(_stencil_val(__FILE__,__LINE__,hf.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,-1) + dry)*dz/(2.*Delta);
#line 528

   _stencil_val(__FILE__,__LINE__,w,0,0,0) -= (_stencil_val(__FILE__,__LINE__,hu.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hu.y,0,1,-1))
   /(_stencil_val(__FILE__,__LINE__,hf.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,-1) + dry)*dz/(2.*Delta);}
      {
#line 531

 _stencil_val(__FILE__,__LINE__,w,0,0,0) -= (_stencil_val(__FILE__,__LINE__,hu.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,0,0,0))/Delta;
#line 531

 _stencil_val(__FILE__,__LINE__,w,0,0,0) -= (_stencil_val(__FILE__,__LINE__,hu.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,0,0))/Delta;}
      dz += _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0), wm = _stencil_val(__FILE__,__LINE__,w,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/hydro.h", 521, "dry");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 535
foreach(){

#line 521 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double dz = val(zb,1,0,0) - val(zb,-1,0,0);
    double wm = 0.;
     { foreach_block_inner() {
      val(w,0,0,0) = wm + (val(hu.x,0,0,0) + val(hu.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry)*
 (dz + val(h,1,0,0) - val(h,-1,0,0))/(2.*Delta);
      if (point.l > 0)
 {
#line 528

   val(w,0,0,0) -= (val(hu.x,0,0,-1) + val(hu.x,1,0,-1))
   /(val(hf.x,0,0,-1) + val(hf.x,1,0,-1) + dry)*dz/(2.*Delta);
#line 528

   val(w,0,0,0) -= (val(hu.y,0,0,-1) + val(hu.y,0,1,-1))
   /(val(hf.y,0,0,-1) + val(hf.y,0,1,-1) + dry)*dz/(2.*Delta);}
      {
#line 531

 val(w,0,0,0) -= (val(hu.x,1,0,0) - val(hu.x,0,0,0))/Delta;
#line 531

 val(w,0,0,0) -= (val(hu.y,0,1,0) - val(hu.y,0,0,0))/Delta;}
      dz += val(h,1,0,0) - val(h,-1,0,0), wm = val(w,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach(); }
}
#line 546 "/home/jiarongw/basilisk/src/layered/hydro.h"
double _radiation (Point point, double ref, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 547 "/home/jiarongw/basilisk/src/layered/hydro.h"

  double H = 0.;
   { foreach_block()
    H += val(h,0,0,0); end_foreach_block(); }



  return sqrt (G*max(H,0.)) - sqrt(G*max(ref - val(zb,0,0,0), 0.));


#if _call__radiation
}
#define _IN_STENCIL 1

#line 546
static double __radiation (Point point, double ref, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 547 "/home/jiarongw/basilisk/src/layered/hydro.h"

  double H = 0.;
   { foreach_block()
    H += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block(); }



  return sqrt (G*max(H,0.)) - sqrt(G*max(ref - _stencil_val(__FILE__,__LINE__,zb,0,0,0), 0.));


#undef _IN_STENCIL

#endif

#line 556
}
#line 565 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 1 "./elevation.h"
#line 1 "/home/jiarongw/basilisk/src/elevation.h"
#line 133 "/home/jiarongw/basilisk/src/elevation.h"
void conserve_elevation (void) {}
#line 566 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 627 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 1 "./gauges.h"
#line 1 "/home/jiarongw/basilisk/src/gauges.h"
#line 11 "/home/jiarongw/basilisk/src/gauges.h"
typedef struct {
  char * name;
  double x, y;
  char * desc;
  FILE * fp;
} Gauge;

void output_gauges (Gauge * gauges, scalar * list)
{
  int n = 0;
  for (Gauge * g = gauges; g->name; g++, n++);
  coord a[n];
  n = 0;
  for (Gauge * g = gauges; g->name; g++, n++) {
    double xp = g->x, yp = g->y;
    unmap (&xp, &yp);
    a[n].x = xp, a[n].y = yp;
  }
  int len = list_len(list);
  double v[n*len];
  interpolate_array (list, a, n, v, false);

  if (pid() == 0) {
    n = 0;
    for (Gauge * g = gauges; g->name; g++) {
      if (!g->fp) {
 g->fp = fopen (g->name, "w");
 if (g->desc)
   fprintf (g->fp, "%s\n", g->desc);
      }
      if (v[n] != nodata) {
 fprintf (g->fp, "%g", t);
 strongif (list) for (scalar s = *list, *_i75 = list; ((scalar *)&s)->i >= 0; s = *++_i75)
   fprintf (g->fp, " %g", v[n++]);
 fputc ('\n', g->fp);
 fflush (g->fp);
      }
      else
 n += len;
    }
  }
}
#line 628 "/home/jiarongw/basilisk/src/layered/hydro.h"
#line 648 "/home/jiarongw/basilisk/src/layered/hydro.h"
double segment_flux (coord segment[2], double * flux, scalar h, vector u)
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  for (int l = 0; l < nl; l++)
    flux[l] = 0.;
   { 
strongif (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 654
foreach_segment (segment, p){

#line 654 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double dl = 0.;
    {
#line 656
 {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0))/2.;
      dl += sq(dp);
    }
#line 656
 {
      double dp = (p[1].y - p[0].y)*Delta/Delta_y*(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0))/2.;
      dl += sq(dp);
    }}
    dl = sqrt (dl);
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
       { foreach_block_inner()
 flux[_layer] += dl/2.*
 interpolate_linear (point, (struct _interpolate)
       {h, a.x, a.y, 0.})*
 (m.x*interpolate_linear (point, (struct _interpolate)
     {u.x, a.x, a.y, 0.}) +
  m.y*interpolate_linear (point, (struct _interpolate)
     {u.y, a.x, a.y, 0.})); end_foreach_block_inner(); }
    }
  } } end_foreach_segment(); }
strongif (is_constant(fm.x)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 654
foreach_segment (segment, p){

#line 654 "/home/jiarongw/basilisk/src/layered/hydro.h"
 {
    double dl = 0.;
    {
#line 656
 {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0))/2.;
      dl += sq(dp);
    }
#line 656
 {
      double dp = (p[1].y - p[0].y)*Delta/Delta_y*(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0))/2.;
      dl += sq(dp);
    }}
    dl = sqrt (dl);
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
       { foreach_block_inner()
 flux[_layer] += dl/2.*
 interpolate_linear (point, (struct _interpolate)
       {h, a.x, a.y, 0.})*
 (m.x*interpolate_linear (point, (struct _interpolate)
     {u.x, a.x, a.y, 0.}) +
  m.y*interpolate_linear (point, (struct _interpolate)
     {u.y, a.x, a.y, 0.})); end_foreach_block_inner(); }
    }
  } } end_foreach_segment(); } }


  MPI_Allreduce (MPI_IN_PLACE, flux, nl, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double tot = 0.;
  for (int l = 0; l < nl; l++)
    tot += flux[l];
  return tot;
}
#line 691 "/home/jiarongw/basilisk/src/layered/hydro.h"
typedef struct {
  char * name;
  coord s[2];
  char * desc;
  FILE * fp;
} Flux;

void output_fluxes (Flux * fluxes, scalar h, vector u)
{
  for (Flux * f = fluxes; f->name; f++) {
    double flux[nl];
    double tot = segment_flux (f->s, flux, h, u);
    if (pid() == 0) {
      if (!f->fp) {
 f->fp = fopen (f->name, "w");
 if (f->desc)
   fprintf (f->fp, "%s\n", f->desc);
      }
      fprintf (f->fp, "%g %g", t, tot);
      for (int i = 0; i < nl; i++)
 fprintf (f->fp, " %g", flux[i]);
      fputc ('\n', f->fp);
      fflush (f->fp);
    }
  }
}
#line 8 "field_new_PM_restore.c"
#line 1 "layered/nh.h"
#line 1 "/home/jiarongw/basilisk/src/layered/nh.h"
#line 44 "/home/jiarongw/basilisk/src/layered/nh.h"
#line 1 "layered/implicit.h"
#line 1 "/home/jiarongw/basilisk/src/layered/implicit.h"
#line 22 "/home/jiarongw/basilisk/src/layered/implicit.h"
#line 1 "./poisson.h"
#line 1 "/home/jiarongw/basilisk/src/poisson.h"
#line 32 "/home/jiarongw/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 55 "/home/jiarongw/basilisk/src/poisson.h"

 strongif (da) for (scalar s = *da, *_i76 = da; ((scalar *)&s)->i >= 0; s = *++_i76)
    { foreach_blockf (s)
     val(s,0,0,0) = 0.; end_foreach_blockf(); }; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 65 "/home/jiarongw/basilisk/src/poisson.h"

 strongif (da) for (scalar s = *da, *_i77 = da; ((scalar *)&s)->i >= 0; s = *++_i77)
    { foreach_blockf (s)
     val(s,0,0,0) = bilinear (point, s); end_foreach_blockf(); }; } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/poisson.h", .line = 84,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 84 "/home/jiarongw/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i8 = a; scalar * _i9 = da; strongif (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i8, ds = *++_i9)
       { foreach_blockf (s)
 _stencil_val(__FILE__,__LINE__,s,0,0,0) += _stencil_val(__FILE__,__LINE__,ds,0,0,0); end_foreach_blockf(); }
  } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 89
foreach(){

#line 84 "/home/jiarongw/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i8 = a; scalar * _i9 = da; strongif (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i8, ds = *++_i9)
       { foreach_blockf (s)
 val(s,0,0,0) += val(ds,0,0,0); end_foreach_blockf(); }
  } } end_foreach(); }
}
#line 102 "/home/jiarongw/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 125 "/home/jiarongw/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    strongif (da) for (scalar s = *da, *_i78 = da; ((scalar *)&s)->i >= 0; s = *++_i78)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _sum = sum;
{ double sum = _sum; NOT_UNUSED(sum);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/poisson.h", .line = 164,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 164 "/home/jiarongw/basilisk/src/poisson.h"

    strongif (p.b) for (scalar s = *p.b, *_i79 = p.b; ((scalar *)&s)->i >= 0; s = *++_i79)
      sum += _stencil_val(__FILE__,__LINE__,s,0,0,0); } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 166

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:sum)) {

#line 164
foreach (){

#line 164 "/home/jiarongw/basilisk/src/poisson.h"

    strongif (p.b) for (scalar s = *p.b, *_i79 = p.b; ((scalar *)&s)->i >= 0; s = *++_i79)
      sum += val(s,0,0,0); } end_foreach();mpi_all_reduce_array (&sum, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 166
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 199 "/home/jiarongw/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 258 "/home/jiarongw/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 296 "/home/jiarongw/basilisk/src/poisson.h"
  scalar c = a;






   { 
strongif (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/jiarongw/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/jiarongw/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/jiarongw/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 303
foreach_level_or_leaf (l){

#line 303 "/home/jiarongw/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 305
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 305
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }}
#line 319 "/home/jiarongw/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 338 "/home/jiarongw/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;



  double maxres = 0.;
#line 375 "/home/jiarongw/basilisk/src/poisson.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _maxres = maxres;
{ double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/poisson.h", .line = 375,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 375
foreach_stencil(){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((_stencil_val(__FILE__,__LINE__,a,1,0,0) - _stencil_val(__FILE__,__LINE__,a,1 -1,0,0))/Delta))/Delta;
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((_stencil_val(__FILE__,__LINE__,a,0,1,0) - _stencil_val(__FILE__,__LINE__,a,0,1 -1,0))/Delta))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }
strongif (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 375
foreach_stencil(){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((_stencil_val(__FILE__,__LINE__,a,1,0,0) - _stencil_val(__FILE__,__LINE__,a,1 -1,0,0))/Delta))/Delta;
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((_stencil_val(__FILE__,__LINE__,a,0,1,0) - _stencil_val(__FILE__,__LINE__,a,0,1 -1,0))/Delta))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }
strongif (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 375
foreach_stencil(){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((_stencil_val(__FILE__,__LINE__,a,1,0,0) - _stencil_val(__FILE__,__LINE__,a,1 -1,0,0))/Delta))/Delta;
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((_stencil_val(__FILE__,__LINE__,a,0,1,0) - _stencil_val(__FILE__,__LINE__,a,0,1 -1,0))/Delta))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }
strongif (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 375
foreach_stencil(){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,b,0,0,0) - val_lambda(lambda,0,0,0)*_stencil_val(__FILE__,__LINE__,a,0,0,0);
    {
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((_stencil_val(__FILE__,__LINE__,a,1,0,0) - _stencil_val(__FILE__,__LINE__,a,1 -1,0,0))/Delta))/Delta;
#line 377

      _stencil_val(__FILE__,__LINE__,res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,a,0,0,0) - _stencil_val(__FILE__,__LINE__,a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((_stencil_val(__FILE__,__LINE__,a,0,1,0) - _stencil_val(__FILE__,__LINE__,a,0,1 -1,0))/Delta))/Delta;}






    IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
      maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 388

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 375

strongif (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 375
foreach (){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 377

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 377

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
strongif (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 375
foreach (){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 377

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 377

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
strongif (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 375
foreach (){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 377

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 377

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
strongif (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 375
foreach (){

#line 375 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 377

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 377

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > maxres)
      maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 388
 }

  return maxres;
}
#line 402 "/home/jiarongw/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 460 "/home/jiarongw/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/jiarongw/basilisk/src/poisson.h", 470);
  vector uf = q.uf;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/poisson.h", .line = 483,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 483 "/home/jiarongw/basilisk/src/poisson.h"
 {
    _stencil_val(__FILE__,__LINE__,div,0,0,0) = 0.;
    {
#line 485

      _stencil_val(__FILE__,__LINE__,div,0,0,0) += _stencil_val(__FILE__,__LINE__,uf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,uf.x,0,0,0);
#line 485

      _stencil_val(__FILE__,__LINE__,div,0,0,0) += _stencil_val(__FILE__,__LINE__,uf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,uf.y,0,0,0);}
    _stencil_val(__FILE__,__LINE__,div,0,0,0) /= dt*Delta;
  } } end_foreach_stencil(); if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/poisson.h", 483, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 488
foreach(){

#line 483 "/home/jiarongw/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 485

      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
#line 485

      val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 499 "/home/jiarongw/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/poisson.h", .line = 505,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 505
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 506
 }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 505
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    _stencil_val(__FILE__,__LINE__,uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((_stencil_val(__FILE__,__LINE__,p,0,0,0) - _stencil_val(__FILE__,__LINE__,p,0,0 -1,0))/Delta); }  }}  end_foreach_face_stencil()
#line 506
 } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/poisson.h", 505, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 506

strongif (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#line 505
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 506
 end_foreach_face(); }
strongif (is_constant(alpha.x)) {
const struct { double x, y; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#line 505
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 505
{

#line 505 "/home/jiarongw/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  end_foreach_face_generic()
#line 506
 end_foreach_face(); } }

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/jiarongw/basilisk/src/poisson.h", 508);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/jiarongw/basilisk/src/poisson.h", 509); }
#line 23 "/home/jiarongw/basilisk/src/layered/implicit.h"

mgstats mgH;
double theta_H = 0.5;
#line 34 "/home/jiarongw/basilisk/src/layered/implicit.h"
static int defaults0_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults0_0 (const int i, const double t, Event * _ev) { trace ("defaults0_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 34); 
{
  if (CFL_H == 1e40)
    CFL_H = HUGE;
  mgH.nrelax = 4;
 end_trace("defaults0_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 39); } return 0; } 
#line 56 "/home/jiarongw/basilisk/src/layered/implicit.h"

static void relax_hydro (scalar * ql, scalar * rhsl, int lev, void * data)
{ trace ("relax_hydro", "/home/jiarongw/basilisk/src/layered/implicit.h", 58);
  scalar eta = ql[0], rhs_eta = rhsl[0];
  vector alpha = *((vector *)data);
   { 
strongif (!is_constant(cm) && !is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 61
foreach_level_or_leaf (lev){

#line 61 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double d = - val_cm(cm,0,0,0)*Delta;
    double n = d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 65
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 65
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(cm) && !is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 61
foreach_level_or_leaf (lev){

#line 61 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double d = - val_cm(cm,0,0,0)*Delta;
    double n = d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 65
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 65
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (!is_constant(cm) && is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 61
foreach_level_or_leaf (lev){

#line 61 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double d = - val_cm(cm,0,0,0)*Delta;
    double n = d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 65
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 65
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(cm) && is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 61
foreach_level_or_leaf (lev){

#line 61 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double d = - val_cm(cm,0,0,0)*Delta;
    double n = d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 65
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 65
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
 end_trace("relax_hydro", "/home/jiarongw/basilisk/src/layered/implicit.h", 72); }


static double residual_hydro (scalar * ql, scalar * rhsl,
         scalar * resl, void * data)
{ trace ("residual_hydro", "/home/jiarongw/basilisk/src/layered/implicit.h", 77);
  scalar eta = ql[0], rhs_eta = rhsl[0], res_eta = resl[0];
  vector alpha = *((vector *)data);
  double maxres = 0.;
  vector g= new_face_vector("g");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _G = G;
{ double G = _G; NOT_UNUSED(G);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/implicit.h", .line = 82,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 82
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  end_foreach_face_stencil()
#line 83
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 82
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  end_foreach_face_stencil()
#line 83
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 82
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  end_foreach_face_stencil()
#line 83
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 82
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta); }  }}  end_foreach_face_stencil()
#line 83
 } if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 82, "G");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 83

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 82
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.x,0,0,0) = val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.y,0,0,0) = val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta); }  }}  end_foreach_face_generic()
#line 83
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 82
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.x,0,0,0) = val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.y,0,0,0) = val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta); }  }}  end_foreach_face_generic()
#line 83
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 82
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.x,0,0,0) = val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.y,0,0,0) = val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta); }  }}  end_foreach_face_generic()
#line 83
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 82
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.x,0,0,0) = val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 82
{

#line 82 "/home/jiarongw/basilisk/src/layered/implicit.h"

    val(g.y,0,0,0) = val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta); }  }}  end_foreach_face_generic()
#line 83
 end_foreach_face(); } }

   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _maxres = maxres;
{ double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/implicit.h", .line = 85,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 85
foreach_stencil(){

#line 85 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
    {
#line 87

      _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += (_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 87

      _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += (_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
    IF (fabs(_stencil_val(__FILE__,__LINE__,res_eta,0,0,0)) > maxres)
      maxres = fabs(_stencil_val(__FILE__,__LINE__,res_eta,0,0,0));
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 85
foreach_stencil(){

#line 85 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
    {
#line 87

      _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += (_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 87

      _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += (_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
    IF (fabs(_stencil_val(__FILE__,__LINE__,res_eta,0,0,0)) > maxres)
      maxres = fabs(_stencil_val(__FILE__,__LINE__,res_eta,0,0,0));
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 91

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 85

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 85
foreach (){

#line 85 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
    {
#line 87

      val(res_eta,0,0,0) += (val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 87

      val(res_eta,0,0,0) += (val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
    if (fabs(val(res_eta,0,0,0)) > maxres)
      maxres = fabs(val(res_eta,0,0,0));
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 85
foreach (){

#line 85 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
    {
#line 87

      val(res_eta,0,0,0) += (val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 87

      val(res_eta,0,0,0) += (val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
    if (fabs(val(res_eta,0,0,0)) > maxres)
      maxres = fabs(val(res_eta,0,0,0));
  } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 91
 }

  { double _ret =  maxres; delete (((scalar []){g.x,g.y,{-1}}));  end_trace("residual_hydro", "/home/jiarongw/basilisk/src/layered/implicit.h", 93);  return _ret; }
 delete (((scalar []){g.x,g.y,{-1}}));  end_trace("residual_hydro", "/home/jiarongw/basilisk/src/layered/implicit.h", 94); }




scalar res_eta = {-1};

scalar rhs_eta;
vector alpha_eta;
#line 115 "/home/jiarongw/basilisk/src/layered/implicit.h"
static int half_advection_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int half_advection_0 (const int i, const double t, Event * _ev) { trace ("half_advection_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 115);  {
  if (theta_H < 1.)
    advect (tracers, hu, hf, (1. - theta_H)*dt);
 end_trace("half_advection_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 118); } return 0; } 





static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 124); 
{
  vector su= new_face_vector("su");
  alpha_eta = new_face_vector("alpha_eta");
  double C = - sq(theta_H*dt);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _theta_H = theta_H;
 double _G = G;
 double _dry = dry;
 double _dt = dt;
 double _C = C;
{ double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double G = _G; NOT_UNUSED(G);
 double dry = _dry; NOT_UNUSED(dry);
 double dt = _dt; NOT_UNUSED(dt);
 double C = _C; NOT_UNUSED(C);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/implicit.h", .line = 129,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 129
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_stencil()
#line 143
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 129
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_stencil()
#line 143
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 129
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_stencil()
#line 143
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 129
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.x,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.x,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = (1. - theta_H)*(_stencil_val(__FILE__,__LINE__,hu.y,0,0,0) + dt*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax) + theta_H*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += dt*(theta_H*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax);
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) -= _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hu.y,0,0,0);
      _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);
    } end_foreach_block_inner(); }
    _stencil_val(__FILE__,__LINE__,alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_stencil()
#line 143
 } if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 129, "theta_H");
 }
 if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 129, "G");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 129, "dry");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 129, "dt");
 }
 if (_first_call) {
 if (C != _C)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 129, "C");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 143

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 129
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    val(su.x,0,0,0) = val(alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = (1. - theta_H)*(val(hu.x,0,0,0) + dt*val(hf.x,0,0,0)*ax) + theta_H*val(hf.x,0,0,0)*uf;
      val(hu.x,0,0,0) += dt*(theta_H*val(ha.x,0,0,0) - val(hf.x,0,0,0)*ax);
      val(ha.x,0,0,0) -= val(hf.x,0,0,0)*ax;
      val(su.x,0,0,0) += val(hu.x,0,0,0);
      val(alpha_eta.x,0,0,0) += val(hf.x,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    val(su.y,0,0,0) = val(alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = (1. - theta_H)*(val(hu.y,0,0,0) + dt*val(hf.y,0,0,0)*ax) + theta_H*val(hf.y,0,0,0)*uf;
      val(hu.y,0,0,0) += dt*(theta_H*val(ha.y,0,0,0) - val(hf.y,0,0,0)*ax);
      val(ha.y,0,0,0) -= val(hf.y,0,0,0)*ax;
      val(su.y,0,0,0) += val(hu.y,0,0,0);
      val(alpha_eta.y,0,0,0) += val(hf.y,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_generic()
#line 143
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 129
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    val(su.x,0,0,0) = val(alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = (1. - theta_H)*(val(hu.x,0,0,0) + dt*val(hf.x,0,0,0)*ax) + theta_H*val(hf.x,0,0,0)*uf;
      val(hu.x,0,0,0) += dt*(theta_H*val(ha.x,0,0,0) - val(hf.x,0,0,0)*ax);
      val(ha.x,0,0,0) -= val(hf.x,0,0,0)*ax;
      val(su.x,0,0,0) += val(hu.x,0,0,0);
      val(alpha_eta.x,0,0,0) += val(hf.x,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    val(su.y,0,0,0) = val(alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = (1. - theta_H)*(val(hu.y,0,0,0) + dt*val(hf.y,0,0,0)*ax) + theta_H*val(hf.y,0,0,0)*uf;
      val(hu.y,0,0,0) += dt*(theta_H*val(ha.y,0,0,0) - val(hf.y,0,0,0)*ax);
      val(ha.y,0,0,0) -= val(hf.y,0,0,0)*ax;
      val(su.y,0,0,0) += val(hu.y,0,0,0);
      val(alpha_eta.y,0,0,0) += val(hf.y,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_generic()
#line 143
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 129
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    val(su.x,0,0,0) = val(alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = (1. - theta_H)*(val(hu.x,0,0,0) + dt*val(hf.x,0,0,0)*ax) + theta_H*val(hf.x,0,0,0)*uf;
      val(hu.x,0,0,0) += dt*(theta_H*val(ha.x,0,0,0) - val(hf.x,0,0,0)*ax);
      val(ha.x,0,0,0) -= val(hf.x,0,0,0)*ax;
      val(su.x,0,0,0) += val(hu.x,0,0,0);
      val(alpha_eta.x,0,0,0) += val(hf.x,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    val(su.y,0,0,0) = val(alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = (1. - theta_H)*(val(hu.y,0,0,0) + dt*val(hf.y,0,0,0)*ax) + theta_H*val(hf.y,0,0,0)*uf;
      val(hu.y,0,0,0) += dt*(theta_H*val(ha.y,0,0,0) - val(hf.y,0,0,0)*ax);
      val(ha.y,0,0,0) -= val(hf.y,0,0,0)*ax;
      val(su.y,0,0,0) += val(hu.y,0,0,0);
      val(alpha_eta.y,0,0,0) += val(hf.y,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_generic()
#line 143
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 129
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    val(su.x,0,0,0) = val(alpha_eta.x,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = (1. - theta_H)*(val(hu.x,0,0,0) + dt*val(hf.x,0,0,0)*ax) + theta_H*val(hf.x,0,0,0)*uf;
      val(hu.x,0,0,0) += dt*(theta_H*val(ha.x,0,0,0) - val(hf.x,0,0,0)*ax);
      val(ha.x,0,0,0) -= val(hf.x,0,0,0)*ax;
      val(su.x,0,0,0) += val(hu.x,0,0,0);
      val(alpha_eta.x,0,0,0) += val(hf.x,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.x,0,0,0) *= C;
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 129
{

#line 129 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    val(su.y,0,0,0) = val(alpha_eta.y,0,0,0) = 0.;
     { foreach_block_inner() {
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = (1. - theta_H)*(val(hu.y,0,0,0) + dt*val(hf.y,0,0,0)*ax) + theta_H*val(hf.y,0,0,0)*uf;
      val(hu.y,0,0,0) += dt*(theta_H*val(ha.y,0,0,0) - val(hf.y,0,0,0)*ax);
      val(ha.y,0,0,0) -= val(hf.y,0,0,0)*ax;
      val(su.y,0,0,0) += val(hu.y,0,0,0);
      val(alpha_eta.y,0,0,0) += val(hf.y,0,0,0);
    } end_foreach_block_inner(); }
    val(alpha_eta.y,0,0,0) *= C;
  } }  }}  end_foreach_face_generic()
#line 143
 end_foreach_face(); } }
#line 152 "/home/jiarongw/basilisk/src/layered/implicit.h"
  rhs_eta = new_scalar("rhs_eta");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
{ double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/implicit.h", .line = 153,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 153
foreach_stencil(){

#line 153 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,eta,0,0,0);
    {
#line 155

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,su.x,1,0,0) - _stencil_val(__FILE__,__LINE__,su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 155

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,su.y,0,1,0) - _stencil_val(__FILE__,__LINE__,su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 153
foreach_stencil(){

#line 153 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,eta,0,0,0);
    {
#line 155

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,su.x,1,0,0) - _stencil_val(__FILE__,__LINE__,su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 155

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,su.y,0,1,0) - _stencil_val(__FILE__,__LINE__,su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach_stencil(); } if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 153, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 157

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 153
foreach(){

#line 153 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    val(rhs_eta,0,0,0) = val(eta,0,0,0);
    {
#line 155

      val(rhs_eta,0,0,0) -= dt*(val(su.x,1,0,0) - val(su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 155

      val(rhs_eta,0,0,0) -= dt*(val(su.y,0,1,0) - val(su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 153
foreach(){

#line 153 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    val(rhs_eta,0,0,0) = val(eta,0,0,0);
    {
#line 155

      val(rhs_eta,0,0,0) -= dt*(val(su.x,1,0,0) - val(su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 155

      val(rhs_eta,0,0,0) -= dt*(val(su.y,0,1,0) - val(su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach(); } }







  restriction (((scalar []){cm,zb,h,hf.x,hf.y,alpha_eta.x,alpha_eta.y,{-1}}));
#line 176 "/home/jiarongw/basilisk/src/layered/implicit.h"
 delete (((scalar []){su.x,su.y,{-1}}));  end_trace("acceleration_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 176); } return 0; } 







static int pressure_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int pressure_0 (const int i, const double t, Event * _ev) { trace ("pressure_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 184); 
{
  mgH = mg_solve ((struct MGSolve){((scalar []){eta,{-1}}), ((scalar []){rhs_eta,{-1}}), residual_hydro, relax_hydro, &alpha_eta,
    .res = res_eta.i >= 0 ? ((scalar []){res_eta,{-1}}) : NULL,
    .nrelax = 4, .minlevel = 1,
    .tolerance = TOLERANCE});
  delete (((scalar []){rhs_eta,alpha_eta.x,alpha_eta.y,{-1}}));
#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _theta_H = theta_H;
 double _G = G;
 double _dry = dry;
 double _dt = dt;
{ double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double G = _G; NOT_UNUSED(G);
 double dry = _dry; NOT_UNUSED(dry);
 double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/implicit.h", .line = 211,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 211
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_stencil()
#line 220
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 211
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_stencil()
#line 220
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 211
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_stencil()
#line 220
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 211
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,-1,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,-1,0,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.x,-1,0,0) + hr*_stencil_val(__FILE__,__LINE__,u.x,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*ax;
      double hl = _stencil_val(__FILE__,__LINE__,h,0,-1,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,-1,0) : 0.;
      double hr = _stencil_val(__FILE__,__LINE__,h,0,0,0) > dry ? _stencil_val(__FILE__,__LINE__,h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*_stencil_val(__FILE__,__LINE__,u.y,0,-1,0) + hr*_stencil_val(__FILE__,__LINE__,u.y,0,0,0))/(hl + hr) : 0.;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) = theta_H*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*uf + dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0)) - dt*_stencil_val(__FILE__,__LINE__,ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_stencil()
#line 220
 } if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 211, "theta_H");
 }
 if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 211, "G");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 211, "dry");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/implicit.h", 211, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 220

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 211
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.x,0,0,0) += val(hf.x,0,0,0)*ax;
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = theta_H*(val(hf.x,0,0,0)*uf + dt*val(ha.x,0,0,0)) - dt*val(ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.y,0,0,0) += val(hf.y,0,0,0)*ax;
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = theta_H*(val(hf.y,0,0,0)*uf + dt*val(ha.y,0,0,0)) - dt*val(ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_generic()
#line 220
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 211
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.x,0,0,0) += val(hf.x,0,0,0)*ax;
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = theta_H*(val(hf.x,0,0,0)*uf + dt*val(ha.x,0,0,0)) - dt*val(ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.y,0,0,0) += val(hf.y,0,0,0)*ax;
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = theta_H*(val(hf.y,0,0,0)*uf + dt*val(ha.y,0,0,0)) - dt*val(ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_generic()
#line 220
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 211
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.x,0,0,0) += val(hf.x,0,0,0)*ax;
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = theta_H*(val(hf.x,0,0,0)*uf + dt*val(ha.x,0,0,0)) - dt*val(ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.y,0,0,0) += val(hf.y,0,0,0)*ax;
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = theta_H*(val(hf.y,0,0,0)*uf + dt*val(ha.y,0,0,0)) - dt*val(ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_generic()
#line 220
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 211
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.x,0,0,0) += val(hf.x,0,0,0)*ax;
      double hl = val(h,-1,0,0) > dry ? val(h,-1,0,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.x,-1,0,0) + hr*val(u.x,0,0,0))/(hl + hr) : 0.;
      val(hu.x,0,0,0) = theta_H*(val(hf.x,0,0,0)*uf + dt*val(ha.x,0,0,0)) - dt*val(ha.x,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 211
{

#line 211 "/home/jiarongw/basilisk/src/layered/implicit.h"
 {
    double ax = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
     { foreach_block_inner() {
      val(ha.y,0,0,0) += val(hf.y,0,0,0)*ax;
      double hl = val(h,0,-1,0) > dry ? val(h,0,-1,0) : 0.;
      double hr = val(h,0,0,0) > dry ? val(h,0,0,0) : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*val(u.y,0,-1,0) + hr*val(u.y,0,0,0))/(hl + hr) : 0.;
      val(hu.y,0,0,0) = theta_H*(val(hf.y,0,0,0)*uf + dt*val(ha.y,0,0,0)) - dt*val(ha.y,0,0,0);
    } end_foreach_block_inner(); }
  } }  }}  end_foreach_face_generic()
#line 220
 end_foreach_face(); } }
 end_trace("pressure_0", "/home/jiarongw/basilisk/src/layered/implicit.h", 221); } return 0; } 
#line 45 "/home/jiarongw/basilisk/src/layered/nh.h"

scalar w, phi;
mgstats mgp;
double breaking = HUGE;







static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/jiarongw/basilisk/src/layered/nh.h", 56); 
{
  hydrostatic = false;
  mgp.nrelax = 4;

  if (!(nl > 0)) qassert ("/home/jiarongw/basilisk/src/layered/nh.h", 61, "nl > 0");
  w = new_block_scalar("w", "", nl);
  phi = new_block_scalar("phi", "", nl);
  reset (((scalar []){w,phi,{-1}}), 0.);

  if (!linearised)
    tracers = list_append (tracers, w);
 end_trace("defaults_1", "/home/jiarongw/basilisk/src/layered/nh.h", 68); } return 0; } 







static int viscous_term_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term_0 (const int i, const double t, Event * _ev) { trace ("viscous_term_0", "/home/jiarongw/basilisk/src/layered/nh.h", 76); 
{
  if (nu > 0.)
     { 
#define vertical_diffusion _vertical_diffusion
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dt = dt;
 double _nu = nu;
{ double dt = _dt; NOT_UNUSED(dt);
 double nu = _nu; NOT_UNUSED(nu);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 79,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 79 "/home/jiarongw/basilisk/src/layered/nh.h"

      vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.); } end_foreach_stencil(); if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 79, "dt");
 }
 if (_first_call) {
 if (nu != _nu)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 79, "nu");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#undef vertical_diffusion
#line 80
foreach(){

#line 79 "/home/jiarongw/basilisk/src/layered/nh.h"

      vertical_diffusion (point, h, w, dt, nu, 0., 0., 0.); } end_foreach(); }
 end_trace("viscous_term_0", "/home/jiarongw/basilisk/src/layered/nh.h", 81); } return 0; } 
#line 112 "/home/jiarongw/basilisk/src/layered/nh.h"
static void box_matrix (Point point, scalar phi, scalar rhs,
   vector hf, scalar eta,
   double * H, double * d)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 115 "/home/jiarongw/basilisk/src/layered/nh.h"

strongif (!is_constant(cm) && !is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = val(zb,0,0,0) - val(zb,-1,0,0), dzp.x = val(zb,1,0,0) - val(zb,0,0,0);
#line 117

    dz.y = val(zb,0,0,0) - val(zb,0,-1,0), dzp.y = val(zb,0,1,0) - val(zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += val(h,0,0,0) - val(h,-1,0,0), dzp.x += val(h,1,0,0) - val(h,0,0,0);
#line 120

      dz.y += val(h,0,0,0) - val(h,0,-1,0), dzp.y += val(h,0,1,0) - val(h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = val(h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = val(rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) < max_slope ? ((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) : (((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) : (((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) - s)*val(phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) + sp)*val(phi,1,0,m) +
   2.*theta_H*Delta*(val(hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) -
       val(hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) < max_slope ? ((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) : (((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) : (((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) - s)*val(phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) + sp)*val(phi,0,1,m) +
   2.*theta_H*Delta*(val(hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) -
       val(hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) + s)*val(phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) - sp)*val(phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) + s)*val(phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) - sp)*val(phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = val(h,0,0,nl-1-k);
      if (hk > dry) {
 H[l*nl + k] -= 8.*s*val(h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*val(h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= val(h,0,0,m) - val(h,-1,0,m), dzp.x -= val(h,1,0,m) - val(h,0,0,m);
#line 156

      dz.y -= val(h,0,0,m) - val(h,0,-1,m), dzp.y -= val(h,0,1,m) - val(h,0,0,m);}
  }
 }
strongif (is_constant(cm) && !is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = val(zb,0,0,0) - val(zb,-1,0,0), dzp.x = val(zb,1,0,0) - val(zb,0,0,0);
#line 117

    dz.y = val(zb,0,0,0) - val(zb,0,-1,0), dzp.y = val(zb,0,1,0) - val(zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += val(h,0,0,0) - val(h,-1,0,0), dzp.x += val(h,1,0,0) - val(h,0,0,0);
#line 120

      dz.y += val(h,0,0,0) - val(h,0,-1,0), dzp.y += val(h,0,1,0) - val(h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = val(h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = val(rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) < max_slope ? ((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) : (((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) : (((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) - s)*val(phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) + sp)*val(phi,1,0,m) +
   2.*theta_H*Delta*(val(hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) -
       val(hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) < max_slope ? ((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) : (((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) : (((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) - s)*val(phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) + sp)*val(phi,0,1,m) +
   2.*theta_H*Delta*(val(hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) -
       val(hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) + s)*val(phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) - sp)*val(phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) + s)*val(phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) - sp)*val(phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = val(h,0,0,nl-1-k);
      if (hk > dry) {
 H[l*nl + k] -= 8.*s*val(h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*val(h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= val(h,0,0,m) - val(h,-1,0,m), dzp.x -= val(h,1,0,m) - val(h,0,0,m);
#line 156

      dz.y -= val(h,0,0,m) - val(h,0,-1,m), dzp.y -= val(h,0,1,m) - val(h,0,0,m);}
  }
 }
strongif (!is_constant(cm) && is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = val(zb,0,0,0) - val(zb,-1,0,0), dzp.x = val(zb,1,0,0) - val(zb,0,0,0);
#line 117

    dz.y = val(zb,0,0,0) - val(zb,0,-1,0), dzp.y = val(zb,0,1,0) - val(zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += val(h,0,0,0) - val(h,-1,0,0), dzp.x += val(h,1,0,0) - val(h,0,0,0);
#line 120

      dz.y += val(h,0,0,0) - val(h,0,-1,0), dzp.y += val(h,0,1,0) - val(h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = val(h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = val(rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) < max_slope ? ((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) : (((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) : (((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) - s)*val(phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) + sp)*val(phi,1,0,m) +
   2.*theta_H*Delta*(val(hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) -
       val(hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) < max_slope ? ((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) : (((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) : (((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) - s)*val(phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) + sp)*val(phi,0,1,m) +
   2.*theta_H*Delta*(val(hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) -
       val(hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) + s)*val(phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) - sp)*val(phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) + s)*val(phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) - sp)*val(phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = val(h,0,0,nl-1-k);
      if (hk > dry) {
 H[l*nl + k] -= 8.*s*val(h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*val(h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= val(h,0,0,m) - val(h,-1,0,m), dzp.x -= val(h,1,0,m) - val(h,0,0,m);
#line 156

      dz.y -= val(h,0,0,m) - val(h,0,-1,m), dzp.y -= val(h,0,1,m) - val(h,0,0,m);}
  }
 }
strongif (is_constant(cm) && is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = val(zb,0,0,0) - val(zb,-1,0,0), dzp.x = val(zb,1,0,0) - val(zb,0,0,0);
#line 117

    dz.y = val(zb,0,0,0) - val(zb,0,-1,0), dzp.y = val(zb,0,1,0) - val(zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += val(h,0,0,0) - val(h,-1,0,0), dzp.x += val(h,1,0,0) - val(h,0,0,0);
#line 120

      dz.y += val(h,0,0,0) - val(h,0,-1,0), dzp.y += val(h,0,1,0) - val(h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = val(h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = val(rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) < max_slope ? ((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) : (((dz.x - val(h,0,0,m) + val(h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) : (((dzp.x - val(h,1,0,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) - s)*val(phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) + sp)*val(phi,1,0,m) +
   2.*theta_H*Delta*(val(hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) -
       val(hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) < max_slope ? ((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) : (((dz.y - val(h,0,0,m) + val(h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) < max_slope ? ((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) : (((dzp.y - val(h,0,1,m) + val(h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) - s)*val(phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) + sp)*val(phi,0,1,m) +
   2.*theta_H*Delta*(val(hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) -
       val(hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    if (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,-1,0,m) + s)*val(phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,1,0,m) - sp)*val(phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,-1,m) + s)*val(phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,1,m) - sp)*val(phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = val(h,0,0,nl-1-k);
      if (hk > dry) {
 H[l*nl + k] -= 8.*s*val(h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*val(h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= val(h,0,0,m) - val(h,-1,0,m), dzp.x -= val(h,1,0,m) - val(h,0,0,m);
#line 156

      dz.y -= val(h,0,0,m) - val(h,0,-1,m), dzp.y -= val(h,0,1,m) - val(h,0,0,m);}
  }
 }
#if _call_box_matrix
}
#define _IN_STENCIL 1

#line 112
static void _box_matrix (Point point, scalar phi, scalar rhs,
   vector hf, scalar eta,
   double * H, double * d)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 115 "/home/jiarongw/basilisk/src/layered/nh.h"

strongif (!is_constant(cm) && !is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,-1,0,0), dzp.x = _stencil_val(__FILE__,__LINE__,zb,1,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);
#line 117

    dz.y = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,-1,0), dzp.y = _stencil_val(__FILE__,__LINE__,zb,0,1,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0), dzp.x += _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);
#line 120

      dz.y += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,-1,0), dzp.y += _stencil_val(__FILE__,__LINE__,h,0,1,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = _stencil_val(__FILE__,__LINE__,h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = _stencil_val(__FILE__,__LINE__,rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) < max_slope ? ((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) : (((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) - s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,1 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) < max_slope ? ((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) : (((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) - s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,1 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    IF (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) + s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) + s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = _stencil_val(__FILE__,__LINE__,h,0,0,nl-1-k);
      IF (hk > dry) {
 H[l*nl + k] -= 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,-1,0,m), dzp.x -= _stencil_val(__FILE__,__LINE__,h,1,0,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);
#line 156

      dz.y -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,0,-1,m), dzp.y -= _stencil_val(__FILE__,__LINE__,h,0,1,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);}
  }
 }
strongif (is_constant(cm) && !is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,-1,0,0), dzp.x = _stencil_val(__FILE__,__LINE__,zb,1,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);
#line 117

    dz.y = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,-1,0), dzp.y = _stencil_val(__FILE__,__LINE__,zb,0,1,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0), dzp.x += _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);
#line 120

      dz.y += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,-1,0), dzp.y += _stencil_val(__FILE__,__LINE__,h,0,1,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = _stencil_val(__FILE__,__LINE__,h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = _stencil_val(__FILE__,__LINE__,rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) < max_slope ? ((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) : (((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) - s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,1 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) < max_slope ? ((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) : (((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) - s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,1 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    IF (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) + s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) + s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = _stencil_val(__FILE__,__LINE__,h,0,0,nl-1-k);
      IF (hk > dry) {
 H[l*nl + k] -= 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,-1,0,m), dzp.x -= _stencil_val(__FILE__,__LINE__,h,1,0,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);
#line 156

      dz.y -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,0,-1,m), dzp.y -= _stencil_val(__FILE__,__LINE__,h,0,1,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);}
  }
 }
strongif (!is_constant(cm) && is_constant(fm.x)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,-1,0,0), dzp.x = _stencil_val(__FILE__,__LINE__,zb,1,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);
#line 117

    dz.y = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,-1,0), dzp.y = _stencil_val(__FILE__,__LINE__,zb,0,1,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0), dzp.x += _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);
#line 120

      dz.y += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,-1,0), dzp.y += _stencil_val(__FILE__,__LINE__,h,0,1,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = _stencil_val(__FILE__,__LINE__,h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = _stencil_val(__FILE__,__LINE__,rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) < max_slope ? ((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) : (((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) - s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,1 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) < max_slope ? ((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) : (((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) - s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,1 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    IF (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) + s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) + s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = _stencil_val(__FILE__,__LINE__,h,0,0,nl-1-k);
      IF (hk > dry) {
 H[l*nl + k] -= 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,-1,0,m), dzp.x -= _stencil_val(__FILE__,__LINE__,h,1,0,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);
#line 156

      dz.y -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,0,-1,m), dzp.y -= _stencil_val(__FILE__,__LINE__,h,0,1,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);}
  }
 }
strongif (is_constant(cm) && is_constant(fm.x)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#line 115

  coord dz, dzp;
  {
#line 117

    dz.x = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,-1,0,0), dzp.x = _stencil_val(__FILE__,__LINE__,zb,1,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);
#line 117

    dz.y = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,-1,0), dzp.y = _stencil_val(__FILE__,__LINE__,zb,0,1,0) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);}
   { foreach_block()
    {
#line 120

      dz.x += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,-1,0,0), dzp.x += _stencil_val(__FILE__,__LINE__,h,1,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);
#line 120

      dz.y += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,-1,0), dzp.y += _stencil_val(__FILE__,__LINE__,h,0,1,0) - _stencil_val(__FILE__,__LINE__,h,0,0,0);} end_foreach_block(); }
  for (int l = 0, m = nl - 1; l < nl; l++, m--) {
    double a = _stencil_val(__FILE__,__LINE__,h,0,0,m)/(sq(Delta)*val_cm(cm,0,0,0));
    d[l] = _stencil_val(__FILE__,__LINE__,rhs,0,0,m);
    for (int k = 0; k < nl; k++)
      H[l*nl + k] = 0.;
    {
#line 127
 {
      double s = Delta*(fabs((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) < max_slope ? ((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) : (((dz.x - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,-1,0,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.x - _stencil_val(__FILE__,__LINE__,h,1,0,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) - s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.x,0,0,m)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.x,1,0,m)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,1 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,1,0,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }
#line 127
 {
      double s = Delta*(fabs((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) < max_slope ? ((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) : (((dz.y - _stencil_val(__FILE__,__LINE__,h,0,0,m) + _stencil_val(__FILE__,__LINE__,h,0,-1,m))/Delta) > 0. ? max_slope : - max_slope));
      double sp = Delta*(fabs((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) < max_slope ? ((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) : (((dzp.y - _stencil_val(__FILE__,__LINE__,h,0,1,m) + _stencil_val(__FILE__,__LINE__,h,0,0,m))/Delta) > 0. ? max_slope : - max_slope));
      d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) - s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) + sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m) +
   2.*theta_H*Delta*(_stencil_val(__FILE__,__LINE__,hf.y,0,0,m)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta) -
       _stencil_val(__FILE__,__LINE__,hf.y,0,1,m)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,1 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,1,0))/Delta)));
      H[l*nl + l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + s) +
   (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - sp));
    }}
    H[l*nl + l] -= 4.;
    IF (l > 0) {
      H[l*(nl + 1) - 1] = 4.;
      {
#line 140
 {
        double s = Delta*(fabs(dz.x/Delta) < max_slope ? (dz.x/Delta) : ((dz.x/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.x/Delta) < max_slope ? (dzp.x/Delta) : ((dzp.x/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,-1,0,m) + s)*_stencil_val(__FILE__,__LINE__,phi,-1,0,m+1) +
     (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,1,0,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,1,0,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }
#line 140
 {
        double s = Delta*(fabs(dz.y/Delta) < max_slope ? (dz.y/Delta) : ((dz.y/Delta) > 0. ? max_slope : - max_slope));
        double sp = Delta*(fabs(dzp.y/Delta) < max_slope ? (dzp.y/Delta) : ((dzp.y/Delta) > 0. ? max_slope : - max_slope));
 d[l] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,-1,m) + s)*_stencil_val(__FILE__,__LINE__,phi,0,-1,m+1) +
     (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,1,m) - sp)*_stencil_val(__FILE__,__LINE__,phi,0,1,m+1));
 H[l*(nl + 1) - 1] -= a*((2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) - s) +
    (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(_stencil_val(__FILE__,__LINE__,h,0,0,m) + sp));
      }}
    }
    for (int k = l + 1, s = -1; k < nl; k++, s = -s) {
      double hk = _stencil_val(__FILE__,__LINE__,h,0,0,nl-1-k);
      IF (hk > dry) {
 H[l*nl + k] -= 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
 H[l*nl + k - 1] += 8.*s*_stencil_val(__FILE__,__LINE__,h,0,0,m)/hk;
      }
    }
    {
#line 156

      dz.x -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,-1,0,m), dzp.x -= _stencil_val(__FILE__,__LINE__,h,1,0,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);
#line 156

      dz.y -= _stencil_val(__FILE__,__LINE__,h,0,0,m) - _stencil_val(__FILE__,__LINE__,h,0,-1,m), dzp.y -= _stencil_val(__FILE__,__LINE__,h,0,1,m) - _stencil_val(__FILE__,__LINE__,h,0,0,m);}
  }
 }
#undef _IN_STENCIL

#endif

#line 159
}




#line 1 "./hessenberg.h"
#line 1 "/home/jiarongw/basilisk/src/hessenberg.h"
#line 27 "/home/jiarongw/basilisk/src/hessenberg.h"
static inline void givens (double x, double y, double * c, double * s)
{
#line 47 "/home/jiarongw/basilisk/src/hessenberg.h"
  double t = sqrt (sq(x) + sq(y));
  *c = x/t, *s = -y/t;

}

void solve_hessenberg (double * H, double * x, int n)
{
  double v[n], c[n], s[n];
  for (int i = 0; i < n; i++)
    v[i] = H[n*(i + 1) - 1];
  for (int k = n - 1; k >= 1; k--) {
    double a = H[k*n + k - 1];
    givens (v[k], a, &c[k], &s[k]);
    x[k] /= c[k]*v[k] - s[k]*a;
    double ykck = x[k]*c[k], yksk = x[k]*s[k];
    for (int l = 0; l <= k - 2; l++) {
      a = H[l*n + k - 1];
      x[l] -= ykck*v[l] - yksk*a;
      v[l] = c[k]*a + s[k]*v[l];
    }
    a = H[(k - 1)*n + k - 1];
    x[k-1] -= ykck*v[k-1] - yksk*a;
    v[k-1] = c[k]*a + s[k]*v[k-1];
  }
  double tau1 = x[0]/v[0];
  for (int k = 1; k < n; k++) {
    double tau2 = x[k];
    x[k-1] = c[k]*tau1 - s[k]*tau2;
    tau1 = c[k]*tau2 + s[k]*tau1;
  }
  x[n-1] = tau1;
}
#line 165 "/home/jiarongw/basilisk/src/layered/nh.h"

vector hf;


static void relax_nh (scalar * phil, scalar * rhsl, int lev, void * data)
{ trace ("relax_nh", "/home/jiarongw/basilisk/src/layered/nh.h", 170);
  scalar phi = phil[0], rhs = rhsl[0];
  scalar eta = phil[1], rhs_eta = rhsl[1];
  vector alpha = *((vector *)data);
   { 
strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 174
foreach_level_or_leaf (lev){

#line 174 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 184 "/home/jiarongw/basilisk/src/layered/nh.h"
    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, hf, eta, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
     { foreach_block_inner()
      val(phi,0,0,0) = b[l--]; end_foreach_block_inner(); }






    double n = 0.;
    {
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,1,0,0) - val(zb,1 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,1,0,0) + val(h,1 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) + s)*val(phi,1,0,0) - (val(h,1 -1,0,0) - s)*val(phi,1 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) < max_slope ? ((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) : (((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) - s)*val(phi,1,0,1) - (val(h,1 -1,0,0) + s)*val(phi,1 -1,0,1); } pg *= (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*val(hf.x,1,0,0)/(Delta*(val(h,1,0,0) + val(h,1 -1,0,0))); }
 n += pg;
      dz += val(h,1,0,0) - val(h,1 -1,0,0); } end_foreach_block_inner(); } };
    }
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,0,1,0) - val(zb,0,1 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,1,0) + val(h,0,1 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) + s)*val(phi,0,1,0) - (val(h,0,1 -1,0) - s)*val(phi,0,1 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) < max_slope ? ((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) : (((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) - s)*val(phi,0,1,1) - (val(h,0,1 -1,0) + s)*val(phi,0,1 -1,1); } pg *= (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*val(hf.y,0,1,0)/(Delta*(val(h,0,1,0) + val(h,0,1 -1,0))); }
 n += pg;
      dz += val(h,0,1,0) - val(h,0,1 -1,0); } end_foreach_block_inner(); } };
    }}
    n *= theta_H*sq(dt);

    double d = - val_cm(cm,0,0,0)*Delta;
    n += d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 211
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 211
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 174
foreach_level_or_leaf (lev){

#line 174 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 184 "/home/jiarongw/basilisk/src/layered/nh.h"
    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, hf, eta, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
     { foreach_block_inner()
      val(phi,0,0,0) = b[l--]; end_foreach_block_inner(); }






    double n = 0.;
    {
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,1,0,0) - val(zb,1 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,1,0,0) + val(h,1 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) + s)*val(phi,1,0,0) - (val(h,1 -1,0,0) - s)*val(phi,1 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) < max_slope ? ((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) : (((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) - s)*val(phi,1,0,1) - (val(h,1 -1,0,0) + s)*val(phi,1 -1,0,1); } pg *= (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*val(hf.x,1,0,0)/(Delta*(val(h,1,0,0) + val(h,1 -1,0,0))); }
 n += pg;
      dz += val(h,1,0,0) - val(h,1 -1,0,0); } end_foreach_block_inner(); } };
    }
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,0,1,0) - val(zb,0,1 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,1,0) + val(h,0,1 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) + s)*val(phi,0,1,0) - (val(h,0,1 -1,0) - s)*val(phi,0,1 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) < max_slope ? ((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) : (((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) - s)*val(phi,0,1,1) - (val(h,0,1 -1,0) + s)*val(phi,0,1 -1,1); } pg *= (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*val(hf.y,0,1,0)/(Delta*(val(h,0,1,0) + val(h,0,1 -1,0))); }
 n += pg;
      dz += val(h,0,1,0) - val(h,0,1 -1,0); } end_foreach_block_inner(); } };
    }}
    n *= theta_H*sq(dt);

    double d = - val_cm(cm,0,0,0)*Delta;
    n += d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 211
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 211
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 174
foreach_level_or_leaf (lev){

#line 174 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 184 "/home/jiarongw/basilisk/src/layered/nh.h"
    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, hf, eta, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
     { foreach_block_inner()
      val(phi,0,0,0) = b[l--]; end_foreach_block_inner(); }






    double n = 0.;
    {
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,1,0,0) - val(zb,1 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,1,0,0) + val(h,1 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) + s)*val(phi,1,0,0) - (val(h,1 -1,0,0) - s)*val(phi,1 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) < max_slope ? ((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) : (((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) - s)*val(phi,1,0,1) - (val(h,1 -1,0,0) + s)*val(phi,1 -1,0,1); } pg *= (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*val(hf.x,1,0,0)/(Delta*(val(h,1,0,0) + val(h,1 -1,0,0))); }
 n += pg;
      dz += val(h,1,0,0) - val(h,1 -1,0,0); } end_foreach_block_inner(); } };
    }
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,0,1,0) - val(zb,0,1 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,1,0) + val(h,0,1 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) + s)*val(phi,0,1,0) - (val(h,0,1 -1,0) - s)*val(phi,0,1 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) < max_slope ? ((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) : (((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) - s)*val(phi,0,1,1) - (val(h,0,1 -1,0) + s)*val(phi,0,1 -1,1); } pg *= (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*val(hf.y,0,1,0)/(Delta*(val(h,0,1,0) + val(h,0,1 -1,0))); }
 n += pg;
      dz += val(h,0,1,0) - val(h,0,1 -1,0); } end_foreach_block_inner(); } };
    }}
    n *= theta_H*sq(dt);

    double d = - val_cm(cm,0,0,0)*Delta;
    n += d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 211
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 211
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 174
foreach_level_or_leaf (lev){

#line 174 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 184 "/home/jiarongw/basilisk/src/layered/nh.h"
    double H[nl*nl], b[nl];
    box_matrix (point, phi, rhs, hf, eta, H, b);
    solve_hessenberg (H, b, nl);
    int l = nl - 1;
     { foreach_block_inner()
      val(phi,0,0,0) = b[l--]; end_foreach_block_inner(); }






    double n = 0.;
    {
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,1,0,0) - val(zb,1 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,1,0,0) + val(h,1 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) + s)*val(phi,1,0,0) - (val(h,1 -1,0,0) - s)*val(phi,1 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) < max_slope ? ((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) : (((dz + val(h,1,0,0) - val(h,1 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,1,0,0) - s)*val(phi,1,0,1) - (val(h,1 -1,0,0) + s)*val(phi,1 -1,0,1); } pg *= (2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*val(hf.x,1,0,0)/(Delta*(val(h,1,0,0) + val(h,1 -1,0,0))); }
 n += pg;
      dz += val(h,1,0,0) - val(h,1 -1,0,0); } end_foreach_block_inner(); } };
    }
#line 197
 {
      double pg;
      { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
 n -= pg;
      dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
      { double dz = val(zb,0,1,0) - val(zb,0,1 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,1,0) + val(h,0,1 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) + s)*val(phi,0,1,0) - (val(h,0,1 -1,0) - s)*val(phi,0,1 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) < max_slope ? ((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) : (((dz + val(h,0,1,0) - val(h,0,1 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,1,0) - s)*val(phi,0,1,1) - (val(h,0,1 -1,0) + s)*val(phi,0,1 -1,1); } pg *= (2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*val(hf.y,0,1,0)/(Delta*(val(h,0,1,0) + val(h,0,1 -1,0))); }
 n += pg;
      dz += val(h,0,1,0) - val(h,0,1 -1,0); } end_foreach_block_inner(); } };
    }}
    n *= theta_H*sq(dt);

    double d = - val_cm(cm,0,0,0)*Delta;
    n += d*val(rhs_eta,0,0,0);
    val(eta,0,0,0) = 0.;
    {
#line 211
 {
      n += val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val(eta,1 -1,0,0) - val(eta,1,0,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.x,0,0,0)*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val_diagonal(eta,0 -1,0,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.x,1,0,0)*(G*(2.*val_fm_x(fm.x,1,0,0)/(val_cm(cm,1,0,0) + val_cm(cm,1 -1,0,0)))*(val_diagonal(eta,1 -1,0,0) - val_diagonal(eta,1,0,0))/Delta);
    }
#line 211
 {
      n += val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val(eta,0,1 -1,0) - val(eta,0,1,0))/Delta);
      diagonalize (eta)
 d -= val(alpha.y,0,0,0)*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val_diagonal(eta,0,0 -1,0) - val_diagonal(eta,0,0,0))/Delta) - val(alpha.y,0,1,0)*(G*(2.*val_fm_y(fm.y,0,1,0)/(val_cm(cm,0,1,0) + val_cm(cm,0,1 -1,0)))*(val_diagonal(eta,0,1 -1,0) - val_diagonal(eta,0,1,0))/Delta);
    }}
    val(eta,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
 end_trace("relax_nh", "/home/jiarongw/basilisk/src/layered/nh.h", 218); }





static double residual_nh (scalar * phil, scalar * rhsl,
      scalar * resl, void * data)
{ trace ("residual_nh", "/home/jiarongw/basilisk/src/layered/nh.h", 226);
  scalar phi = phil[0], rhs = rhsl[0], res = resl[0];
  scalar eta = phil[1], rhs_eta = rhsl[1], res_eta = resl[1];
  double maxres = 0.;

  vector g = new_block_face_vector("g", nl);
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _theta_H = theta_H;
 double _G = G;
 double _dry = dry;
 double _max_slope = max_slope;
 int _nl = nl;
{ double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double G = _G; NOT_UNUSED(G);
 double dry = _dry; NOT_UNUSED(dry);
 double max_slope = _max_slope; NOT_UNUSED(max_slope);
 int nl = _nl; NOT_UNUSED(nl);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 232,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 232
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); }
      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); }
      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 238
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 232
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); }
      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); }
      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 238
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 232
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); }
      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); }
      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 238
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 232
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(_stencil_val(__FILE__,__LINE__,eta,0 -1,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); }
      _stencil_val(__FILE__,__LINE__,g.x,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.x,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(_stencil_val(__FILE__,__LINE__,eta,0,0 -1,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0))/Delta);
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); }
      _stencil_val(__FILE__,__LINE__,g.y,0,0,0) = - 2.*(pg + _stencil_val(__FILE__,__LINE__,hf.y,0,0,0)*pgh);
    dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 238
 } if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 232, "theta_H");
 }
 if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 232, "G");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 232, "dry");
 }
 if (_first_call) {
 if (max_slope != _max_slope)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 232, "max_slope");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 232, "nl");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 238

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 232
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
      val(g.x,0,0,0) = - 2.*(pg + val(hf.x,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
      val(g.y,0,0,0) = - 2.*(pg + val(hf.y,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 238
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 232
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
      val(g.x,0,0,0) = - 2.*(pg + val(hf.x,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
      val(g.y,0,0,0) = - 2.*(pg + val(hf.y,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 238
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 232
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
      val(g.x,0,0,0) = - 2.*(pg + val(hf.x,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
      val(g.y,0,0,0) = - 2.*(pg + val(hf.y,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 238
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 232
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*(val(eta,0 -1,0,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); }
      val(g.x,0,0,0) = - 2.*(pg + val(hf.x,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 232
{

#line 232 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double pgh = theta_H*(G*(2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*(val(eta,0,0 -1,0) - val(eta,0,0,0))/Delta);
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); }
      val(g.y,0,0,0) = - 2.*(pg + val(hf.y,0,0,0)*pgh);
    dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 238
 end_foreach_face(); } }

   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dry = dry;
 double _max_slope = max_slope;
 int _nl = nl;
 double _theta_H = theta_H;
 double _dt = dt;
 double _maxres = maxres;
{ double dry = _dry; NOT_UNUSED(dry);
 double max_slope = _max_slope; NOT_UNUSED(max_slope);
 int nl = _nl; NOT_UNUSED(nl);
 double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double dt = _dt; NOT_UNUSED(dt);
 double maxres = _maxres; NOT_UNUSED(maxres);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 240,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 240
foreach_stencil(){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs,0,0,0) + 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,0);
      {
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry)*
   (fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.x,1,0,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry)*
   (fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.y,0,1,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      IF (point.l < nl - 1)
        _stencil_val(__FILE__,__LINE__,res,0,0,0) -= 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = _stencil_val(__FILE__,__LINE__,h,0,0,k);
 IF (hk > dry)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) += 8.*s*(_stencil_val(__FILE__,__LINE__,phi,0,0,k) - _stencil_val(__FILE__,__LINE__,phi,0,0,k+1))*_stencil_val(__FILE__,__LINE__,h,0,0,0)/hk;
      }
      IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
      {
#line 276

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 276

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 240
foreach_stencil(){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs,0,0,0) + 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,0);
      {
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry)*
   (fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.x,1,0,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry)*
   (fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.y,0,1,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      IF (point.l < nl - 1)
        _stencil_val(__FILE__,__LINE__,res,0,0,0) -= 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = _stencil_val(__FILE__,__LINE__,h,0,0,k);
 IF (hk > dry)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) += 8.*s*(_stencil_val(__FILE__,__LINE__,phi,0,0,k) - _stencil_val(__FILE__,__LINE__,phi,0,0,k+1))*_stencil_val(__FILE__,__LINE__,h,0,0,0)/hk;
      }
      IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
      {
#line 276

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 276

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 240
foreach_stencil(){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs,0,0,0) + 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,0);
      {
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry)*
   (fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.x,1,0,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry)*
   (fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.y,0,1,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      IF (point.l < nl - 1)
        _stencil_val(__FILE__,__LINE__,res,0,0,0) -= 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = _stencil_val(__FILE__,__LINE__,h,0,0,k);
 IF (hk > dry)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) += 8.*s*(_stencil_val(__FILE__,__LINE__,phi,0,0,k) - _stencil_val(__FILE__,__LINE__,phi,0,0,k+1))*_stencil_val(__FILE__,__LINE__,h,0,0,0)/hk;
      }
      IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
      {
#line 276

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 276

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 240
foreach_stencil(){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,res,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs,0,0,0) + 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,0);
      {
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,0) + _stencil_val(__FILE__,__LINE__,g.x,1,0,0))/(_stencil_val(__FILE__,__LINE__,hf.x,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) + dry)*
   (fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.x,1,0,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.x,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 _stencil_val(__FILE__,__LINE__,res,0,0,0) += _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,0) + _stencil_val(__FILE__,__LINE__,g.y,0,1,0))/(_stencil_val(__FILE__,__LINE__,hf.y,0,0,0) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) + dry)*
   (fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 IF (point.l > 0)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) -= _stencil_val(__FILE__,__LINE__,h,0,0,0)*(_stencil_val(__FILE__,__LINE__,g.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,g.y,0,1,-1))/
     (_stencil_val(__FILE__,__LINE__,hf.y,0,0,-1) + _stencil_val(__FILE__,__LINE__,hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      IF (point.l < nl - 1)
        _stencil_val(__FILE__,__LINE__,res,0,0,0) -= 4.*_stencil_val(__FILE__,__LINE__,phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = _stencil_val(__FILE__,__LINE__,h,0,0,k);
 IF (hk > dry)
   _stencil_val(__FILE__,__LINE__,res,0,0,0) += 8.*s*(_stencil_val(__FILE__,__LINE__,phi,0,0,k) - _stencil_val(__FILE__,__LINE__,phi,0,0,k+1))*_stencil_val(__FILE__,__LINE__,h,0,0,0)/hk;
      }
      IF (fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0)) > maxres)
 maxres = fabs (_stencil_val(__FILE__,__LINE__,res,0,0,0));
      {
#line 276

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 276

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) = _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) - _stencil_val(__FILE__,__LINE__,eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.x,1,0,0) - _stencil_val(__FILE__,__LINE__,g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        _stencil_val(__FILE__,__LINE__,res_eta,0,0,0) += theta_H*sq(dt)/2.*(_stencil_val(__FILE__,__LINE__,g.y,0,1,0) - _stencil_val(__FILE__,__LINE__,g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach_stencil(); } if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 240, "dry");
 }
 if (_first_call) {
 if (max_slope != _max_slope)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 240, "max_slope");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 240, "nl");
 }
 if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 240, "theta_H");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 240, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 293

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(max:maxres)) {

#line 240

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 240
foreach (){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(res,0,0,0) = val(rhs,0,0,0) + 4.*val(phi,0,0,0);
      {
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.x,0,0,0) + val(g.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry)*
   (fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,0,0,-1) + val(g.x,1,0,-1))/
     (val(hf.x,0,0,-1) + val(hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.y,0,0,0) + val(g.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry)*
   (fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,0,-1) + val(g.y,0,1,-1))/
     (val(hf.y,0,0,-1) + val(hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      if (point.l < nl - 1)
        val(res,0,0,0) -= 4.*val(phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = val(h,0,0,k);
 if (hk > dry)
   val(res,0,0,0) += 8.*s*(val(phi,0,0,k) - val(phi,0,0,k+1))*val(h,0,0,0)/hk;
      }
      if (fabs (val(res,0,0,0)) > maxres)
 maxres = fabs (val(res,0,0,0));
      {
#line 276

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 276

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 240
foreach (){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(res,0,0,0) = val(rhs,0,0,0) + 4.*val(phi,0,0,0);
      {
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.x,0,0,0) + val(g.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry)*
   (fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,0,0,-1) + val(g.x,1,0,-1))/
     (val(hf.x,0,0,-1) + val(hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.y,0,0,0) + val(g.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry)*
   (fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,0,-1) + val(g.y,0,1,-1))/
     (val(hf.y,0,0,-1) + val(hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      if (point.l < nl - 1)
        val(res,0,0,0) -= 4.*val(phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = val(h,0,0,k);
 if (hk > dry)
   val(res,0,0,0) += 8.*s*(val(phi,0,0,k) - val(phi,0,0,k+1))*val(h,0,0,0)/hk;
      }
      if (fabs (val(res,0,0,0)) > maxres)
 maxres = fabs (val(res,0,0,0));
      {
#line 276

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 276

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 240
foreach (){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(res,0,0,0) = val(rhs,0,0,0) + 4.*val(phi,0,0,0);
      {
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.x,0,0,0) + val(g.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry)*
   (fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,0,0,-1) + val(g.x,1,0,-1))/
     (val(hf.x,0,0,-1) + val(hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.y,0,0,0) + val(g.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry)*
   (fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,0,-1) + val(g.y,0,1,-1))/
     (val(hf.y,0,0,-1) + val(hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      if (point.l < nl - 1)
        val(res,0,0,0) -= 4.*val(phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = val(h,0,0,k);
 if (hk > dry)
   val(res,0,0,0) += 8.*s*(val(phi,0,0,k) - val(phi,0,0,k+1))*val(h,0,0,0)/hk;
      }
      if (fabs (val(res,0,0,0)) > maxres)
 maxres = fabs (val(res,0,0,0));
      {
#line 276

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 276

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 240
foreach (){

#line 240 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
#line 254 "/home/jiarongw/basilisk/src/layered/nh.h"
    coord dz;
    {
#line 255

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 255

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(res,0,0,0) = val(rhs,0,0,0) + 4.*val(phi,0,0,0);
      {
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.x,0,0,0) + val(g.x,1,0,0))/(val(hf.x,0,0,0) + val(hf.x,1,0,0) + dry)*
   (fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.x,0,0,-1) + val(g.x,1,0,-1))/
     (val(hf.x,0,0,-1) + val(hf.x,1,0,-1) + dry)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }
#line 259
 {
 val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));
 val(res,0,0,0) += val(h,0,0,0)*(val(g.y,0,0,0) + val(g.y,0,1,0))/(val(hf.y,0,0,0) + val(hf.y,0,1,0) + dry)*
   (fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
 if (point.l > 0)
   val(res,0,0,0) -= val(h,0,0,0)*(val(g.y,0,0,-1) + val(g.y,0,1,-1))/
     (val(hf.y,0,0,-1) + val(hf.y,0,1,-1) + dry)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
      }}
      if (point.l < nl - 1)
        val(res,0,0,0) -= 4.*val(phi,0,0,1);
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s) {
 double hk = val(h,0,0,k);
 if (hk > dry)
   val(res,0,0,0) += 8.*s*(val(phi,0,0,k) - val(phi,0,0,k+1))*val(h,0,0,0)/hk;
      }
      if (fabs (val(res,0,0,0)) > maxres)
 maxres = fabs (val(res,0,0,0));
      {
#line 276

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 276

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
    } end_foreach_block_inner(); }
#line 289 "/home/jiarongw/basilisk/src/layered/nh.h"
    val(res_eta,0,0,0) = val(rhs_eta,0,0,0) - val(eta,0,0,0);
     { foreach_block_inner()
      {
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.x,1,0,0) - val(g.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 291

        val(res_eta,0,0,0) += theta_H*sq(dt)/2.*(val(g.y,0,1,0) - val(g.y,0,0,0))/(Delta*val_cm(cm,0,0,0));} end_foreach_block_inner(); }
  } } end_foreach(); }mpi_all_reduce_array (&maxres, double, MPI_MAX, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 293
 }

  delete ((scalar *)((vector []){{g.x,g.y},{{-1},{-1}}}));
  { double _ret =  maxres; end_trace("residual_nh", "/home/jiarongw/basilisk/src/layered/nh.h", 296);  return _ret; }
 end_trace("residual_nh", "/home/jiarongw/basilisk/src/layered/nh.h", 297); }







static int pressure_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int pressure_1 (const int i, const double t, Event * _ev) { trace ("pressure_1", "/home/jiarongw/basilisk/src/layered/nh.h", 305); 
{
#line 350 "/home/jiarongw/basilisk/src/layered/nh.h"
  scalar rhs = new_block_scalar("rhs", "", nl);
  double h1 = 0., v1 = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _max_slope = max_slope;
 double _theta_H = theta_H;
 double _dt = dt;
 double _h1 = h1;
 double _v1 = v1;
{ double max_slope = _max_slope; NOT_UNUSED(max_slope);
 double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double dt = _dt; NOT_UNUSED(dt);
 double h1 = _h1; NOT_UNUSED(h1);
 double v1 = _v1; NOT_UNUSED(v1);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 352,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 352
foreach_stencil(){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) = 2.*_stencil_val(__FILE__,__LINE__,w,0,0,0);
      {
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.x,0,0,0)*(fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.y,0,0,0)*(fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      IF (point.l > 0)
 {
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += 4.*s*_stencil_val(__FILE__,__LINE__,w,0,0,k);
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) *= 2.*_stencil_val(__FILE__,__LINE__,h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 367

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 352
foreach_stencil(){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) = 2.*_stencil_val(__FILE__,__LINE__,w,0,0,0);
      {
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.x,0,0,0)*(fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.y,0,0,0)*(fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      IF (point.l > 0)
 {
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += 4.*s*_stencil_val(__FILE__,__LINE__,w,0,0,k);
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) *= 2.*_stencil_val(__FILE__,__LINE__,h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 367

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 352
foreach_stencil(){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) = 2.*_stencil_val(__FILE__,__LINE__,w,0,0,0);
      {
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.x,0,0,0)*(fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.y,0,0,0)*(fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      IF (point.l > 0)
 {
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += 4.*s*_stencil_val(__FILE__,__LINE__,w,0,0,k);
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) *= 2.*_stencil_val(__FILE__,__LINE__,h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 367

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 352
foreach_stencil(){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(_stencil_val(__FILE__,__LINE__,zb,1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,-1,0,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(_stencil_val(__FILE__,__LINE__,zb,0,1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(_stencil_val(__FILE__,__LINE__,zb,0,-1,0) + _stencil_val(__FILE__,__LINE__,zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) = 2.*_stencil_val(__FILE__,__LINE__,w,0,0,0);
      {
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.x,0,0,0)*(fabs((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += (_stencil_val(__FILE__,__LINE__,hu.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   _stencil_val(__FILE__,__LINE__,u.y,0,0,0)*(fabs((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      IF (point.l > 0)
 {
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += _stencil_val(__FILE__,__LINE__,u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 _stencil_val(__FILE__,__LINE__,rhs,0,0,0) += 4.*s*_stencil_val(__FILE__,__LINE__,w,0,0,k);
      _stencil_val(__FILE__,__LINE__,rhs,0,0,0) *= 2.*_stencil_val(__FILE__,__LINE__,h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += _stencil_val(__FILE__,__LINE__,hf.x,1,0,0) - _stencil_val(__FILE__,__LINE__,hf.x,0,0,0);
#line 367

 dz.y += _stencil_val(__FILE__,__LINE__,hf.y,0,1,0) - _stencil_val(__FILE__,__LINE__,hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*_stencil_val(__FILE__,__LINE__,h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); } if (_first_call) {
 if (max_slope != _max_slope)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 352, "max_slope");
 }
 if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 352, "theta_H");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 352, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 372

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:h1)  reduction(+:v1)) {

#line 352

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 352
foreach (){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(rhs,0,0,0) = 2.*val(w,0,0,0);
      {
#line 358

        val(rhs,0,0,0) += (val(hu.x,1,0,0) - val(hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.x,0,0,0)*(fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        val(rhs,0,0,0) += (val(hu.y,0,1,0) - val(hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.y,0,0,0)*(fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      if (point.l > 0)
 {
#line 362

   val(rhs,0,0,0) += val(u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   val(rhs,0,0,0) += val(u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 val(rhs,0,0,0) += 4.*s*val(w,0,0,k);
      val(rhs,0,0,0) *= 2.*val(h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 367

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*val(h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 352
foreach (){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(rhs,0,0,0) = 2.*val(w,0,0,0);
      {
#line 358

        val(rhs,0,0,0) += (val(hu.x,1,0,0) - val(hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.x,0,0,0)*(fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        val(rhs,0,0,0) += (val(hu.y,0,1,0) - val(hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.y,0,0,0)*(fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      if (point.l > 0)
 {
#line 362

   val(rhs,0,0,0) += val(u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   val(rhs,0,0,0) += val(u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 val(rhs,0,0,0) += 4.*s*val(w,0,0,k);
      val(rhs,0,0,0) *= 2.*val(h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 367

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*val(h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 352
foreach (){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(rhs,0,0,0) = 2.*val(w,0,0,0);
      {
#line 358

        val(rhs,0,0,0) += (val(hu.x,1,0,0) - val(hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.x,0,0,0)*(fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        val(rhs,0,0,0) += (val(hu.y,0,1,0) - val(hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.y,0,0,0)*(fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      if (point.l > 0)
 {
#line 362

   val(rhs,0,0,0) += val(u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   val(rhs,0,0,0) += val(u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 val(rhs,0,0,0) += 4.*s*val(w,0,0,k);
      val(rhs,0,0,0) *= 2.*val(h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 367

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*val(h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 352
foreach (){

#line 352 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    coord dz;
    {
#line 354

      dz.x = (val_fm_x(fm.x,1,0,0)*(val(zb,1,0,0) + val(zb,0,0,0)) - val_fm_x(fm.x,0,0,0)*(val(zb,-1,0,0) + val(zb,0,0,0)))/2.;
#line 354

      dz.y = (val_fm_y(fm.y,0,1,0)*(val(zb,0,1,0) + val(zb,0,0,0)) - val_fm_y(fm.y,0,0,0)*(val(zb,0,-1,0) + val(zb,0,0,0)))/2.;}
     { foreach_block_inner() {
      val(rhs,0,0,0) = 2.*val(w,0,0,0);
      {
#line 358

        val(rhs,0,0,0) += (val(hu.x,1,0,0) - val(hu.x,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.x,0,0,0)*(fabs((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.x + val(hf.x,1,0,0) - val(hf.x,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 358

        val(rhs,0,0,0) += (val(hu.y,0,1,0) - val(hu.y,0,0,0))/(Delta*val_cm(cm,0,0,0)) -
   val(u.y,0,0,0)*(fabs((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) < max_slope ? ((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) : (((dz.y + val(hf.y,0,1,0) - val(hf.y,0,0,0))/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      if (point.l > 0)
 {
#line 362

   val(rhs,0,0,0) += val(u.x,0,0,-1)*(fabs(dz.x/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.x/(Delta*val_cm(cm,0,0,0))) : ((dz.x/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));
#line 362

   val(rhs,0,0,0) += val(u.y,0,0,-1)*(fabs(dz.y/(Delta*val_cm(cm,0,0,0))) < max_slope ? (dz.y/(Delta*val_cm(cm,0,0,0))) : ((dz.y/(Delta*val_cm(cm,0,0,0))) > 0. ? max_slope : - max_slope));}
      for (int k = - 1, s = -1; k >= - point.l; k--, s = -s)
 val(rhs,0,0,0) += 4.*s*val(w,0,0,k);
      val(rhs,0,0,0) *= 2.*val(h,0,0,0)/(theta_H*dt);
      {
#line 367

 dz.x += val(hf.x,1,0,0) - val(hf.x,0,0,0);
#line 367

 dz.y += val(hf.y,0,1,0) - val(hf.y,0,0,0);}
      h1 += (sq(Delta)*val_cm(cm,0,0,0))*val(h,0,0,0);
      v1 += (sq(Delta)*val_cm(cm,0,0,0));
    } end_foreach_block_inner(); }
  } } end_foreach(); }mpi_all_reduce_array (&h1, double, MPI_SUM, 1);
mpi_all_reduce_array (&v1, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 372
 }






  scalar res;
  if (res_eta.i >= 0)
    res = new_block_scalar("res", "", nl);
  mgp = mg_solve ((struct MGSolve){((scalar []){phi,eta,{-1}}), ((scalar []){rhs,rhs_eta,{-1}}), residual_nh, relax_nh, &alpha_eta,
    .res = res_eta.i >= 0 ? ((scalar []){res,res_eta,{-1}}) : NULL,
    .nrelax = 4, .minlevel = 1,
    .tolerance = TOLERANCE*sq(h1/(dt*v1))});
  delete (((scalar []){rhs,{-1}}));
  if (res_eta.i >= 0)
    delete (((scalar []){res,{-1}}));





  vector su= new_face_vector("su");
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dry = dry;
 double _max_slope = max_slope;
 int _nl = nl;
 double _theta_H = theta_H;
 double _dt = dt;
{ double dry = _dry; NOT_UNUSED(dry);
 double max_slope = _max_slope; NOT_UNUSED(max_slope);
 int nl = _nl; NOT_UNUSED(nl);
 double theta_H = _theta_H; NOT_UNUSED(theta_H);
 double dt = _dt; NOT_UNUSED(dt);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 395,
    .each = "foreach_face", .first = _first_call
  };

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 403
 }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 403
 }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 403
 }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 395
foreach_face_stencil() { int ig = -1; VARIABLES;  strongif (is_stencil_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.x,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0 -1,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*_stencil_val(__FILE__,__LINE__,hf.x,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0 -1,0,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.x,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.x,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.x,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_stencil_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    _stencil_val(__FILE__,__LINE__,su.y,0,0,0) = 0.;
    double pg;
    { double dz = _stencil_val(__FILE__,__LINE__,zb,0,0,0) - _stencil_val(__FILE__,__LINE__,zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0,0) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,0); IF (point.l < nl - 1) { double s = Delta*(fabs((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) < max_slope ? ((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) : (((dz + _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (_stencil_val(__FILE__,__LINE__,h,0,0,0) - s)*_stencil_val(__FILE__,__LINE__,phi,0,0,1) - (_stencil_val(__FILE__,__LINE__,h,0,0 -1,0) + s)*_stencil_val(__FILE__,__LINE__,phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*_stencil_val(__FILE__,__LINE__,hf.y,0,0,0)/(Delta*(_stencil_val(__FILE__,__LINE__,h,0,0,0) + _stencil_val(__FILE__,__LINE__,h,0,0 -1,0))); } {
      _stencil_val(__FILE__,__LINE__,ha.y,0,0,0) += pg;
      _stencil_val(__FILE__,__LINE__,su.y,0,0,0) -= pg;
      _stencil_val(__FILE__,__LINE__,hu.y,0,0,0) += theta_H*dt*pg;
    } dz += _stencil_val(__FILE__,__LINE__,h,0,0,0) - _stencil_val(__FILE__,__LINE__,h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_stencil()
#line 403
 } if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 395, "dry");
 }
 if (_first_call) {
 if (max_slope != _max_slope)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 395, "max_slope");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 395, "nl");
 }
 if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 395, "theta_H");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 395, "dt");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 403

strongif (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.x,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); } {
      val(ha.x,0,0,0) += pg;
      val(su.x,0,0,0) -= pg;
      val(hu.x,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.y,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); } {
      val(ha.y,0,0,0) += pg;
      val(su.y,0,0,0) -= pg;
      val(hu.y,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 403
 end_foreach_face(); }
strongif (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.x,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); } {
      val(ha.x,0,0,0) += pg;
      val(su.x,0,0,0) -= pg;
      val(hu.x,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.y,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); } {
      val(ha.y,0,0,0) += pg;
      val(su.y,0,0,0) -= pg;
      val(hu.y,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 403
 end_foreach_face(); }
strongif (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.x,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); } {
      val(ha.x,0,0,0) += pg;
      val(su.x,0,0,0) -= pg;
      val(hu.x,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.y,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); } {
      val(ha.y,0,0,0) += pg;
      val(su.y,0,0,0) -= pg;
      val(hu.y,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 403
 end_foreach_face(); }
strongif (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 395
foreach_face_generic() { int ig = -1; VARIABLES;  strongif (is_face_x()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.x,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0 -1,0,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0 -1,0,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0 -1,0,0) - s)*val(phi,0 -1,0,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0 -1,0,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0 -1,0,0) + s)*val(phi,0 -1,0,1); } pg *= (2.*val_fm_x(fm.x,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0 -1,0,0)))*val(hf.x,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0 -1,0,0))); } {
      val(ha.x,0,0,0) += pg;
      val(su.x,0,0,0) -= pg;
      val(hu.x,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0 -1,0,0); } end_foreach_block_inner(); } };
  } }  }}  { int jg = -1; VARIABLES;  strongif (is_face_y()) {
#line 395
{

#line 395 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    val(su.y,0,0,0) = 0.;
    double pg;
    { double dz = val(zb,0,0,0) - val(zb,0,0 -1,0);  { foreach_block_inner() { pg = 0.; if (val(h,0,0,0) + val(h,0,0 -1,0) > dry) { double s = Delta*(fabs(dz/Delta) < max_slope ? (dz/Delta) : ((dz/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) + s)*val(phi,0,0,0) - (val(h,0,0 -1,0) - s)*val(phi,0,0 -1,0); if (point.l < nl - 1) { double s = Delta*(fabs((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) < max_slope ? ((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) : (((dz + val(h,0,0,0) - val(h,0,0 -1,0))/Delta) > 0. ? max_slope : - max_slope)); pg -= (val(h,0,0,0) - s)*val(phi,0,0,1) - (val(h,0,0 -1,0) + s)*val(phi,0,0 -1,1); } pg *= (2.*val_fm_y(fm.y,0,0,0)/(val_cm(cm,0,0,0) + val_cm(cm,0,0 -1,0)))*val(hf.y,0,0,0)/(Delta*(val(h,0,0,0) + val(h,0,0 -1,0))); } {
      val(ha.y,0,0,0) += pg;
      val(su.y,0,0,0) -= pg;
      val(hu.y,0,0,0) += theta_H*dt*pg;
    } dz += val(h,0,0,0) - val(h,0,0 -1,0); } end_foreach_block_inner(); } };
  } }  }}  end_foreach_face_generic()
#line 403
 end_foreach_face(); } }
#line 423 "/home/jiarongw/basilisk/src/layered/nh.h"
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _breaking = breaking;
 double _G = G;
 double _dry = dry;
 int _nl = nl;
 double _dt = dt;
 double _theta_H = theta_H;
{ double breaking = _breaking; NOT_UNUSED(breaking);
 double G = _G; NOT_UNUSED(G);
 double dry = _dry; NOT_UNUSED(dry);
 int nl = _nl; NOT_UNUSED(nl);
 double dt = _dt; NOT_UNUSED(dt);
 double theta_H = _theta_H; NOT_UNUSED(theta_H);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/nh.h", .line = 423,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 423
foreach_stencil(){

#line 423 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double wmax = HUGE;
    IF (breaking < HUGE) {
      wmax = 0.;
       { foreach_block_inner()
 wmax += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block_inner(); }
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
     { foreach_block_inner()
      IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) > dry) {
 IF (point.l == nl - 1)
   _stencil_val(__FILE__,__LINE__,w,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,phi,0,0,0)/_stencil_val(__FILE__,__LINE__,h,0,0,0);
 
   _stencil_val(__FILE__,__LINE__,w,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,phi,0,0,1) - _stencil_val(__FILE__,__LINE__,phi,0,0,0))/_stencil_val(__FILE__,__LINE__,h,0,0,0);
 IF (fabs(_stencil_val(__FILE__,__LINE__,w,0,0,0)) > wmax)
   _stencil_val(__FILE__,__LINE__,w,0,0,0) = (_stencil_val(__FILE__,__LINE__,w,0,0,0) > 0. ? 1. : -1.)*wmax;
      } end_foreach_block_inner(); }
#line 448 "/home/jiarongw/basilisk/src/layered/nh.h"
    {
#line 448

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) += theta_H*sq(dt)*(_stencil_val(__FILE__,__LINE__,su.x,1,0,0) - _stencil_val(__FILE__,__LINE__,su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 448

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) += theta_H*sq(dt)*(_stencil_val(__FILE__,__LINE__,su.y,0,1,0) - _stencil_val(__FILE__,__LINE__,su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 423
foreach_stencil(){

#line 423 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double wmax = HUGE;
    IF (breaking < HUGE) {
      wmax = 0.;
       { foreach_block_inner()
 wmax += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block_inner(); }
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
     { foreach_block_inner()
      IF (_stencil_val(__FILE__,__LINE__,h,0,0,0) > dry) {
 IF (point.l == nl - 1)
   _stencil_val(__FILE__,__LINE__,w,0,0,0) += dt*_stencil_val(__FILE__,__LINE__,phi,0,0,0)/_stencil_val(__FILE__,__LINE__,h,0,0,0);
 
   _stencil_val(__FILE__,__LINE__,w,0,0,0) -= dt*(_stencil_val(__FILE__,__LINE__,phi,0,0,1) - _stencil_val(__FILE__,__LINE__,phi,0,0,0))/_stencil_val(__FILE__,__LINE__,h,0,0,0);
 IF (fabs(_stencil_val(__FILE__,__LINE__,w,0,0,0)) > wmax)
   _stencil_val(__FILE__,__LINE__,w,0,0,0) = (_stencil_val(__FILE__,__LINE__,w,0,0,0) > 0. ? 1. : -1.)*wmax;
      } end_foreach_block_inner(); }
#line 448 "/home/jiarongw/basilisk/src/layered/nh.h"
    {
#line 448

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) += theta_H*sq(dt)*(_stencil_val(__FILE__,__LINE__,su.x,1,0,0) - _stencil_val(__FILE__,__LINE__,su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 448

      _stencil_val(__FILE__,__LINE__,rhs_eta,0,0,0) += theta_H*sq(dt)*(_stencil_val(__FILE__,__LINE__,su.y,0,1,0) - _stencil_val(__FILE__,__LINE__,su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach_stencil(); } if (_first_call) {
 if (breaking != _breaking)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "breaking");
 }
 if (_first_call) {
 if (G != _G)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "G");
 }
 if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "dry");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "nl");
 }
 if (_first_call) {
 if (dt != _dt)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "dt");
 }
 if (_first_call) {
 if (theta_H != _theta_H)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/nh.h", 423, "theta_H");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 450

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 423
foreach(){

#line 423 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double wmax = HUGE;
    if (breaking < HUGE) {
      wmax = 0.;
       { foreach_block_inner()
 wmax += val(h,0,0,0); end_foreach_block_inner(); }
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
     { foreach_block_inner()
      if (val(h,0,0,0) > dry) {
 if (point.l == nl - 1)
   val(w,0,0,0) += dt*val(phi,0,0,0)/val(h,0,0,0);
 else
   val(w,0,0,0) -= dt*(val(phi,0,0,1) - val(phi,0,0,0))/val(h,0,0,0);
 if (fabs(val(w,0,0,0)) > wmax)
   val(w,0,0,0) = (val(w,0,0,0) > 0. ? 1. : -1.)*wmax;
      } end_foreach_block_inner(); }
#line 448 "/home/jiarongw/basilisk/src/layered/nh.h"
    {
#line 448

      val(rhs_eta,0,0,0) += theta_H*sq(dt)*(val(su.x,1,0,0) - val(su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 448

      val(rhs_eta,0,0,0) += theta_H*sq(dt)*(val(su.y,0,1,0) - val(su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 423
foreach(){

#line 423 "/home/jiarongw/basilisk/src/layered/nh.h"
 {
    double wmax = HUGE;
    if (breaking < HUGE) {
      wmax = 0.;
       { foreach_block_inner()
 wmax += val(h,0,0,0); end_foreach_block_inner(); }
      wmax = wmax > 0. ? breaking*sqrt(G*wmax) : 0.;
    }
     { foreach_block_inner()
      if (val(h,0,0,0) > dry) {
 if (point.l == nl - 1)
   val(w,0,0,0) += dt*val(phi,0,0,0)/val(h,0,0,0);
 else
   val(w,0,0,0) -= dt*(val(phi,0,0,1) - val(phi,0,0,0))/val(h,0,0,0);
 if (fabs(val(w,0,0,0)) > wmax)
   val(w,0,0,0) = (val(w,0,0,0) > 0. ? 1. : -1.)*wmax;
      } end_foreach_block_inner(); }
#line 448 "/home/jiarongw/basilisk/src/layered/nh.h"
    {
#line 448

      val(rhs_eta,0,0,0) += theta_H*sq(dt)*(val(su.x,1,0,0) - val(su.x,0,0,0))/(Delta*val_cm(cm,0,0,0));
#line 448

      val(rhs_eta,0,0,0) += theta_H*sq(dt)*(val(su.y,0,1,0) - val(su.y,0,0,0))/(Delta*val_cm(cm,0,0,0));}
  } } end_foreach(); } }
 delete (((scalar []){su.x,su.y,{-1}}));  end_trace("pressure_1", "/home/jiarongw/basilisk/src/layered/nh.h", 451); } return 0; } 






static int cleanup_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup_1 (const int i, const double t, Event * _ev) { trace ("cleanup_1", "/home/jiarongw/basilisk/src/layered/nh.h", 458);  {
  delete (((scalar []){w,phi,{-1}}));
 end_trace("cleanup_1", "/home/jiarongw/basilisk/src/layered/nh.h", 460); } return 0; } 
#line 9 "field_new_PM_restore.c"
#line 1 "layered/remap.h"
#line 1 "/home/jiarongw/basilisk/src/layered/remap.h"
#line 11 "/home/jiarongw/basilisk/src/layered/remap.h"
#line 1 "./ppr/ppr.h"
#line 1 "/home/jiarongw/basilisk/src/ppr/ppr.h"
#line 12 "/home/jiarongw/basilisk/src/ppr/ppr.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#line 34 "/home/jiarongw/basilisk/src/ppr/ppr.h"
void my_remap (int * npos, int * nnew, int * nvar, int * ndof,
        double * xpos, double * xnew,
        double * fdat, double * fnew,
        int * edge_meth, int * cell_meth, int * cell_lim);
#line 12 "/home/jiarongw/basilisk/src/layered/remap.h"


int edge_meth = 101, cell_meth = 202, cell_lim = 300;
#line 23 "/home/jiarongw/basilisk/src/layered/remap.h"
double * beta = NULL;

static int defaults_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_2 (const int i, const double t, Event * _ev) { trace ("defaults_2", "/home/jiarongw/basilisk/src/layered/remap.h", 25); 
{
  beta = pmalloc (nl*sizeof(double),__func__,__FILE__,__LINE__);
  for (int l = 0; l < nl; l++)
    beta[l] = 1./nl;
 end_trace("defaults_2", "/home/jiarongw/basilisk/src/layered/remap.h", 30); } return 0; } 
#line 41 "/home/jiarongw/basilisk/src/layered/remap.h"
void geometric_beta (double rmin, bool top)
{
  if (rmin <= 0. || rmin >= 1. || nl < 2)
    return;
  double r = 1. + 2.*(1./rmin - 1.)/(nl - 1.);
  double hmin = (r - 1.)/(pow(r, nl) - 1.);
  for (int l = 0; l < nl; l++)
    beta[l] = hmin*pow(r, top ? nl - 1 - l : l);
}







void vertical_remapping (scalar h, scalar * tracers)
{ trace ("vertical_remapping", "/home/jiarongw/basilisk/src/layered/remap.h", 58);
  int nvar = list_len(tracers), ndof = 1, npos = nl + 1;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _dry = dry;
 int _npos = npos;
 int _nvar = nvar;
 int _nl = nl;
 double _beta[10] = {0};
 int _ndof = ndof;
 int _edge_meth = edge_meth;
 int _cell_meth = cell_meth;
 int _cell_lim = cell_lim;
{ double dry = _dry; NOT_UNUSED(dry);
 int npos = _npos; NOT_UNUSED(npos);
 int nvar = _nvar; NOT_UNUSED(nvar);
 int nl = _nl; NOT_UNUSED(nl);
 double * beta = _beta; NOT_UNUSED(beta);
 int ndof = _ndof; NOT_UNUSED(ndof);
 int edge_meth = _edge_meth; NOT_UNUSED(edge_meth);
 int cell_meth = _cell_meth; NOT_UNUSED(cell_meth);
 int cell_lim = _cell_lim; NOT_UNUSED(cell_lim);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "/home/jiarongw/basilisk/src/layered/remap.h", .line = 60,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 60 "/home/jiarongw/basilisk/src/layered/remap.h"
 {
#line 71 "/home/jiarongw/basilisk/src/layered/remap.h"
    double H = 0.;
     { foreach_block_inner()
      H += _stencil_val(__FILE__,__LINE__,h,0,0,0); end_foreach_block_inner(); }


    IF (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
       { foreach_block_inner() {
 zpos[point.l+1] = zpos[point.l] + max(_stencil_val(__FILE__,__LINE__,h,0,0,0),dry);
 int i = nvar*point.l;
 strongif (tracers) for (scalar s = *tracers, *_i80 = tracers; ((scalar *)&s)->i >= 0; s = *++_i80)
   fdat[i++] = _stencil_val(__FILE__,__LINE__,s,0,0,0);






 _stencil_val(__FILE__,__LINE__,h,0,0,0) = H*beta[point.l];

 znew[point.l+1] = znew[point.l] + _stencil_val(__FILE__,__LINE__,h,0,0,0);
      } end_foreach_block_inner(); }

      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
  &edge_meth, &cell_meth, &cell_lim);

       { foreach_block_inner() {
 int i = nvar*point.l;
 strongif (tracers) for (scalar s = *tracers, *_i81 = tracers; ((scalar *)&s)->i >= 0; s = *++_i81)
   _stencil_val(__FILE__,__LINE__,s,0,0,0) = fnew[i++];
      } end_foreach_block_inner(); }
    }
  } } end_foreach_stencil(); if (_first_call) {
 if (dry != _dry)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "dry");
 }
 if (_first_call) {
 if (npos != _npos)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "npos");
 }
 if (_first_call) {
 if (nvar != _nvar)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "nvar");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "nl");
 }
 if (_first_call) {
 for (int i = 0; i < (10*sizeof(double)); i++)
   if (((char *)_beta)[i] != 0) {
     reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "beta");
     break; }
 }
 if (_first_call) {
 if (ndof != _ndof)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "ndof");
 }
 if (_first_call) {
 if (edge_meth != _edge_meth)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "edge_meth");
 }
 if (_first_call) {
 if (cell_meth != _cell_meth)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "cell_meth");
 }
 if (_first_call) {
 if (cell_lim != _cell_lim)
   reduction_warning ("/home/jiarongw/basilisk/src/layered/remap.h", 60, "cell_lim");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 105
foreach(){

#line 60 "/home/jiarongw/basilisk/src/layered/remap.h"
 {
#line 71 "/home/jiarongw/basilisk/src/layered/remap.h"
    double H = 0.;
     { foreach_block_inner()
      H += val(h,0,0,0); end_foreach_block_inner(); }


    if (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
       { foreach_block_inner() {
 zpos[point.l+1] = zpos[point.l] + max(val(h,0,0,0),dry);
 int i = nvar*point.l;
 strongif (tracers) for (scalar s = *tracers, *_i80 = tracers; ((scalar *)&s)->i >= 0; s = *++_i80)
   fdat[i++] = val(s,0,0,0);






 val(h,0,0,0) = H*beta[point.l];

 znew[point.l+1] = znew[point.l] + val(h,0,0,0);
      } end_foreach_block_inner(); }

      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
  &edge_meth, &cell_meth, &cell_lim);

       { foreach_block_inner() {
 int i = nvar*point.l;
 strongif (tracers) for (scalar s = *tracers, *_i81 = tracers; ((scalar *)&s)->i >= 0; s = *++_i81)
   val(s,0,0,0) = fnew[i++];
      } end_foreach_block_inner(); }
    }
  } } end_foreach(); }
 end_trace("vertical_remapping", "/home/jiarongw/basilisk/src/layered/remap.h", 106); }




static int remap_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int remap_0 (const int i, const double t, Event * _ev) { trace ("remap_0", "/home/jiarongw/basilisk/src/layered/remap.h", 111);  {
  if (nl > 1)
    vertical_remapping (h, tracers);
 end_trace("remap_0", "/home/jiarongw/basilisk/src/layered/remap.h", 114); } return 0; } 




static int cleanup_2_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup_2 (const int i, const double t, Event * _ev) { trace ("cleanup_2", "/home/jiarongw/basilisk/src/layered/remap.h", 119); 
{
  pfree (beta,__func__,__FILE__,__LINE__), beta = NULL;
 end_trace("cleanup_2", "/home/jiarongw/basilisk/src/layered/remap.h", 122); } return 0; } 
#line 10 "field_new_PM_restore.c"
#line 1 "layered/perfs.h"
#line 1 "/home/jiarongw/basilisk/src/layered/perfs.h"






static int perfs_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i += 1);   *ip = i; *tp = t;   return ret; } static int perfs (const int i, const double t, Event * _ev) { trace ("perfs", "/home/jiarongw/basilisk/src/layered/perfs.h", 7);  {
  static FILE * fp =NULL; strongif (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("perfs", "w");
  if (i == 0)
    fprintf (fp,
      "t dt mgp.i mgp.nrelax grid->tn perf.t perf.speed npe\n");
  fprintf (fp, "%g %g %d %d %ld %g %g %d\n",
    t, dt,

    mgp.i, mgp.nrelax,





    grid->tn*nl, perf.t, perf.speed*nl, npe());
  fflush (fp);
 end_trace("perfs", "/home/jiarongw/basilisk/src/layered/perfs.h", 23); } return 0; } 






static int perf_plot_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 10);   *ip = i; *tp = t;   return ret; } static int perf_plot (const int i, const double t, Event * _ev) { trace ("perf_plot", "/home/jiarongw/basilisk/src/layered/perfs.h", 30);  {
  if (getenv ("DISPLAY"))
    qpopen ("gnuplot -e 'set term x11 noraise title perfs' "
    "$BASILISK/layered/perfs.plot 2> /dev/null "
    "& read dummy; kill $!", "w");
 end_trace("perf_plot", "/home/jiarongw/basilisk/src/layered/perfs.h", 35); } return 0; } 
#line 11 "field_new_PM_restore.c"

#line 1 "output_mpi.h"
#line 1 "./output_mpi.h"





void output_matrix_mpi (struct OutputMatrix p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  float fn = p.n, Delta = L0/fn;
  float ** field = matrix_new (p.n, p.n, sizeof(float));

  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate ((struct _interpolate){p.f, xp, yp});
      }
      else {
        Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 21 "./output_mpi.h"

        field[i][j] = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif


    fwrite (&fn, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < p.n; i++){
      float xp = Delta*i + X0 + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < p.n; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
#if 1
  else
  MPI_Reduce (field[0], NULL, p.n*p.n, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif

  matrix_free (field);
}

struct OutputMatrix_part {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};

void output_matrix_part_mpi (struct OutputMatrix_part p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  float fn = p.n,
  Delta = (p.box[1][0] - p.box[0][0])/fn;
  float ny = (p.box[1][1] - p.box[0][1])/Delta;
  float nx = (p.box[1][0] - p.box[0][0])/Delta;

  float ** field = matrix_new ((int)nx,(int)ny, sizeof(float));

  for (int i = 0; i < nx; i++) {
    float xp = Delta*i + p.box[0][0] + Delta/2.;
    for (int j = 0; j < ny; j++) {
      float yp = Delta*j + p.box[0][1] + Delta/2.;
      if (p.linear) {
        field[i][j] = interpolate ((struct _interpolate){p.f, xp, yp});
      }
      else {
        Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); POINT_VARIABLES; 
#line 88 "./output_mpi.h"

        field[i][j] = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], nx*ny, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif


    fwrite (&nx, sizeof(float), 1, p.fp);
 fwrite (&ny, sizeof(float), 1, p.fp);
    for (int j = 0; j < ny; j++) {
      float yp = Delta*j + p.box[0][1] + Delta/2.;
      fwrite (&yp, sizeof(float), 1, p.fp);
    }

    for (int i = 0; i < nx; i++){
      float xp = Delta*i + p.box[0][0] + Delta/2.;
      fwrite (&xp, sizeof(float), 1, p.fp);
      for (int j = 0; j < ny; j++) {
        fwrite (&field[i][j], sizeof(float), 1, p.fp);
      }
    }
    fflush (p.fp);
  }
#if 1
  else
  MPI_Reduce (field[0], NULL, nx*ny, MPI_FLOAT, MPI_MIN, 0,MPI_COMM_WORLD);
#endif

  matrix_free (field);
}
#line 13 "field_new_PM_restore.c"






double h_ = 10;
double gpe_base = 0;
double TEND = 50.;
int NLAYER = 10;
int LEVEL_data = 7;

#line 1 "./spectrum.h"
#line 1 "././spectrum.h"






double kp_ = 2.*pi/10.;
double F_kxky_[32*(32 +1)], omega[32*(32 +1)], phase[32*(32 +1)];
double kx_[32], ky_[32 +1];
double dkx_, dky_;
int RANDOM;



double randInRange (int min, int max)
{
  return min + (rand() / (double) (RAND_MAX) * (max - min + 1));
}





void power_input () {
#line 34 "././spectrum.h"
  int length1D, length2D;
  char message[20];
  int i, rank, size;
  MPI_Status status;
  int root = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  if (rank == root) {

    length2D = 32*(32 +1);
    float * a = (float*) pmalloc (sizeof(float)*length2D,__func__,__FILE__,__LINE__);
    char filename[100];
    sprintf (filename, "F_kxky");
    FILE * fp = fopen (filename, "rb");
    fread (a, sizeof(float), length2D, fp);
    for (int i=0;i<length2D;i++) {
      F_kxky_[i] = (double)a[i];
    }
    fclose (fp);


    length1D = 32;
    float * b1 = (float*) pmalloc (sizeof(float)*length1D,__func__,__FILE__,__LINE__);
    sprintf (filename, "kx");
    FILE *fp1 = fopen (filename, "rb");
    fread (b1, sizeof(float), length1D, fp1);
    for (int i=0;i<length1D;i++) {
      kx_[i] = (double)b1[i];
    }
    fclose (fp1);


    float * b2 = (float*) pmalloc (sizeof(float)*(length1D+1),__func__,__FILE__,__LINE__);
    sprintf (filename, "ky");
    FILE *fp2 = fopen (filename, "rb");
    fread (b2, sizeof(float), length1D+1, fp2);
    for (int i=0;i<length1D+1;i++) {
      ky_[i] = (double)b2[i];
    }
    fclose (fp2);


    double kmod = 0;
    int index = 0;
    srand(RANDOM);
    for (int i=0; i<32; i++) {
      for (int j=0; j<32 +1; j++) {
 index = j*32 + i;
 kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
 omega[index] = sqrt(9.8*kmod);
 phase[index] = randInRange (0, 2.*pi);
      }
    }
  }

  MPI_Bcast(&kx_, length1D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&ky_, length1D+1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&F_kxky_, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&omega, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&phase, length2D, MPI_DOUBLE, root, MPI_COMM_WORLD);

  char checkout[100];
  sprintf (checkout, "F-%d", pid());
  FILE * fout = fopen (checkout, "w");
  for (int i=0; i<length2D; i++)
    fprintf (fout, "%g ", F_kxky_[i]);
  fclose (fout);
  sprintf (checkout, "ky-%d", pid());
  fout = fopen (checkout, "w");
  for (int i=0; i<length1D+1; i++)
    fprintf (fout, "%g ", ky_[i]);
  fclose (fout);
#line 157 "././spectrum.h"
}



double wave (double x, double y)
{
  double eta = 0;
  double ampl = 0, a = 0;
  int index = 0;
  for (int i=0; i<32; i++) {
    for (int j=0; j<32 +1; j++) {
      index = j*32 + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      eta += ampl*cos(a);
    }
  }
  return eta;
}
double u_x (double x, double y, double z) {
  int index = 0;
  double u_x = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<32; i++) {
    for (int j=0; j<32 +1; j++) {
      index = j*32 + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);

      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_x += sqrt(9.8*kmod)*ampl*exp(kmod*z_actual)*cos(a)*cos(theta);
    }
  }
  return u_x;
}

double u_y (double x, double y, double z) {
  int index = 0;
  double u_y = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0, theta = 0;
  for (int i=0; i<32; i++) {
    for (int j=0; j<32 +1; j++) {
      index = j*32 + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      theta = atan(ky_[j]/kx_[i]);
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_y += sqrt(9.8*kmod)*ampl*exp(kmod*z_actual)*cos(a)*sin(theta);
    }
  }
  return u_y;
}

double u_z (double x, double y, double z) {
  int index = 0;
  double u_z = 0;
  double ampl = 0, a = 0;
  double z_actual = 0, kmod = 0;
  for (int i=0; i<32; i++) {
    for (int j=0; j<32 +1; j++) {
      index = j*32 + i;
      ampl = sqrt(2.*F_kxky_[index]*dkx_*dky_);
      z_actual = (z < ampl ? (z) : ampl);
      kmod = sqrt(sq(kx_[i]) + sq(ky_[j]));
      a = (kx_[i]*x + ky_[j]*y + phase[index]);
      u_z += sqrt(9.8*kmod)*ampl*exp(kmod*z_actual)*sin(a);
    }
  }
  return u_z;
}
#line 26 "field_new_PM_restore.c"
#line 36 "field_new_PM_restore.c"
int main (int argc, char * argv[])
{ _init_solver();
  if (argc > 1)
    NLAYER = atoi(argv[1]);
  if (argc > 2)
    LEVEL_data = atoi(argv[2]);
  if (argc > 3)
    TEND = atof(argv[3]);
  if (argc > 4)
    nu = atof(argv[4]);
  else
    nu = 0.;
  if (argc > 5)
    RANDOM = atoi(argv[5]);
  if (argc > 6)
    L0 = atof(argv[6]);
  else
    L0 = 50.;
  if (argc > 7)
    kp_ = 2.*pi/atof(argv[7]);
  else
    kp_ = 2.*pi/(L0/5.);
  if (argc > 8)
    theta_H = atof(argv[8]);
  else
    theta_H = 0.5;
  origin ((struct _origin){-L0/2., -L0/2.});
  periodic (right);
  periodic (top);
  N = 1 << LEVEL_data;
  nl = NLAYER;
  G = 9.8;
  h_ = 2.*pi/kp_;
#line 78 "field_new_PM_restore.c"
  gpe_base = -0.5*sq(h_)*sq(L0)*9.8;



  CFL_H = 1;


  run();
 free_solver(); }






int phony = 1;
int j = 0;

static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "field_new_PM_restore.c", 96); 
{
  if (!restore ((struct Dump){"restart"})) {
    power_input();
    dkx_ = kx_[1] - kx_[0];
    dky_ = ky_[1] - ky_[0];
    fprintf (ferr, "dkx = %g, dky = %g\n", dkx_, dky_);
    geometric_beta (1./3., true);
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _h_ = h_;
 int _nl = nl;
{ double h_ = _h_; NOT_UNUSED(h_);
 int nl = _nl; NOT_UNUSED(nl);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "field_new_PM_restore.c", .line = 104,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 104 "field_new_PM_restore.c"
 {
      _stencil_val(__FILE__,__LINE__,zb,0,0,0) = -h_;
      _stencil_val(__FILE__,__LINE__,eta,0,0,0) = wave(x, y);
      double H = wave(x, y) - _stencil_val(__FILE__,__LINE__,zb,0,0,0);
       { foreach_block_inner() {
       _stencil_val(__FILE__,__LINE__,h,0,0,0) = H/nl;
      } end_foreach_block_inner(); }
    } } end_foreach_stencil(); if (_first_call) {
 if (h_ != _h_)
   reduction_warning ("field_new_PM_restore.c", 104, "h_");
 }
 if (_first_call) {
 if (nl != _nl)
   reduction_warning ("field_new_PM_restore.c", 104, "nl");
 }
  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 111
foreach(){

#line 104 "field_new_PM_restore.c"
 {
      val(zb,0,0,0) = -h_;
      val(eta,0,0,0) = wave(x, y);
      double H = wave(x, y) - val(zb,0,0,0);
       { foreach_block_inner() {
       val(h,0,0,0) = H/nl;
      } end_foreach_block_inner(); }
    } } end_foreach(); }

    vertical_remapping (h, tracers);
     { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{ {  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "field_new_PM_restore.c", .line = 114,
    .each = "foreach", .first = _first_call
  };
foreach_stencil(){

#line 114 "field_new_PM_restore.c"
 {
      double z = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
       { foreach_block_inner() {
        z += _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.;
        _stencil_val(__FILE__,__LINE__,u.x,0,0,0) = u_x(x, y, z);
        _stencil_val(__FILE__,__LINE__,u.y,0,0,0) = u_y(x, y, z);
        _stencil_val(__FILE__,__LINE__,w,0,0,0) = u_z(x, y, z);
        z += _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.;
      } end_foreach_block_inner(); }
    } } end_foreach_stencil();  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 123
foreach(){

#line 114 "field_new_PM_restore.c"
 {
      double z = val(zb,0,0,0);
       { foreach_block_inner() {
        z += val(h,0,0,0)/2.;
        val(u.x,0,0,0) = u_x(x, y, z);
        val(u.y,0,0,0) = u_y(x, y, z);
        val(w,0,0,0) = u_z(x, y, z);
        z += val(h,0,0,0)/2.;
      } end_foreach_block_inner(); }
    } } end_foreach(); }
    fprintf (ferr,"Done initialization!\n");
    dump((struct Dump){"initial"});
  }
  else {
    phony = 0;
  }
 end_trace("init_0", "field_new_PM_restore.c", 130); } return 0; } 






static int energy_before_remap_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int energy_before_remap (const int i, const double t, Event * _ev) { trace ("energy_before_remap", "field_new_PM_restore.c", 137); 
{
  if (i==10) {
    fprintf(ferr, "energy output before remap!\n");
    fflush(ferr);
  }
  double ke = 0., gpe = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _ke = ke;
 double _gpe = gpe;
{ double ke = _ke; NOT_UNUSED(ke);
 double gpe = _gpe; NOT_UNUSED(gpe);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "field_new_PM_restore.c", .line = 144,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 144
foreach_stencil(){

#line 144 "field_new_PM_restore.c"
 {
    double zc = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(_stencil_val(__FILE__,__LINE__,w,0,0,0));
      {
#line 148

        norm2 += sq(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));
#line 148

        norm2 += sq(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));}
        ke += norm2*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        gpe += (zc + _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.)*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        zc += _stencil_val(__FILE__,__LINE__,h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 144
foreach_stencil(){

#line 144 "field_new_PM_restore.c"
 {
    double zc = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(_stencil_val(__FILE__,__LINE__,w,0,0,0));
      {
#line 148

        norm2 += sq(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));
#line 148

        norm2 += sq(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));}
        ke += norm2*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        gpe += (zc + _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.)*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        zc += _stencil_val(__FILE__,__LINE__,h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 154

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:ke)  reduction(+:gpe)) {

#line 144

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 144
foreach (){

#line 144 "field_new_PM_restore.c"
 {
    double zc = val(zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(val(w,0,0,0));
      {
#line 148

        norm2 += sq(val(u.x,0,0,0));
#line 148

        norm2 += sq(val(u.y,0,0,0));}
        ke += norm2*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        gpe += (zc + val(h,0,0,0)/2.)*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        zc += val(h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 144
foreach (){

#line 144 "field_new_PM_restore.c"
 {
    double zc = val(zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(val(w,0,0,0));
      {
#line 148

        norm2 += sq(val(u.x,0,0,0));
#line 148

        norm2 += sq(val(u.y,0,0,0));}
        ke += norm2*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        gpe += (zc + val(h,0,0,0)/2.)*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
        zc += val(h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach(); }mpi_all_reduce_array (&ke, double, MPI_SUM, 1);
mpi_all_reduce_array (&gpe, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 154
 }
  static FILE * fp =NULL; strongif (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen("energy_before_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., 9.8*gpe - gpe_base);
  fflush (fp);
 end_trace("energy_before_remap", "field_new_PM_restore.c", 158); } return 0; } 

static int energy_after_remap_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int energy_after_remap (const int i, const double t, Event * _ev) { trace ("energy_after_remap", "field_new_PM_restore.c", 160); 
{
  if (i==10) {
    fprintf(ferr, "energy output after remap!\n");
    fflush(ferr);
  }
  double ke = 0., gpe = 0.;
   { 
disable_fpe (FE_DIVBYZERO|FE_INVALID);
{  double _ke = ke;
 double _gpe = gpe;
{ double ke = _ke; NOT_UNUSED(ke);
 double gpe = _gpe; NOT_UNUSED(gpe);
  static bool _first_call = true;
  ForeachData _foreach_data = {
    .fname = "field_new_PM_restore.c", .line = 167,
    .each = "foreach", .first = _first_call
  };

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) _stencil_val(__FILE__,__LINE__,a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) _stencil_fine(__FILE__,__LINE__,a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) _stencil_coarse(__FILE__,__LINE__,a,i,j,k)
#line 167
foreach_stencil(){

#line 167 "field_new_PM_restore.c"
 {
    double zc = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(_stencil_val(__FILE__,__LINE__,w,0,0,0));
      {
#line 171

 norm2 += sq(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));
#line 171

 norm2 += sq(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));}
      ke += norm2*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      gpe += (zc + _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.)*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      zc += _stencil_val(__FILE__,__LINE__,h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 167
foreach_stencil(){

#line 167 "field_new_PM_restore.c"
 {
    double zc = _stencil_val(__FILE__,__LINE__,zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(_stencil_val(__FILE__,__LINE__,w,0,0,0));
      {
#line 171

 norm2 += sq(_stencil_val(__FILE__,__LINE__,u.x,0,0,0));
#line 171

 norm2 += sq(_stencil_val(__FILE__,__LINE__,u.y,0,0,0));}
      ke += norm2*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      gpe += (zc + _stencil_val(__FILE__,__LINE__,h,0,0,0)/2.)*_stencil_val(__FILE__,__LINE__,h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      zc += _stencil_val(__FILE__,__LINE__,h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach_stencil(); }  _first_call = false;
}}
enable_fpe (FE_DIVBYZERO|FE_INVALID);
#line 177

#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel reduction(+:ke)  reduction(+:gpe)) {

#line 167

strongif (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 167
foreach (){

#line 167 "field_new_PM_restore.c"
 {
    double zc = val(zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(val(w,0,0,0));
      {
#line 171

 norm2 += sq(val(u.x,0,0,0));
#line 171

 norm2 += sq(val(u.y,0,0,0));}
      ke += norm2*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      gpe += (zc + val(h,0,0,0)/2.)*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      zc += val(h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach(); }
strongif (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 167
foreach (){

#line 167 "field_new_PM_restore.c"
 {
    double zc = val(zb,0,0,0);
     { foreach_block_inner() {
      double norm2 = sq(val(w,0,0,0));
      {
#line 171

 norm2 += sq(val(u.x,0,0,0));
#line 171

 norm2 += sq(val(u.y,0,0,0));}
      ke += norm2*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      gpe += (zc + val(h,0,0,0)/2.)*val(h,0,0,0)*(sq(Delta)*val_cm(cm,0,0,0));
      zc += val(h,0,0,0);
    } end_foreach_block_inner(); }
  } } end_foreach(); }mpi_all_reduce_array (&ke, double, MPI_SUM, 1);
mpi_all_reduce_array (&gpe, double, MPI_SUM, 1);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 177
 }
  static FILE * fp =NULL; strongif (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen("energy_after_remap.dat","w");
  fprintf (fp, "%g %g %g\n", t, ke/2., 9.8*gpe - gpe_base);
  fflush (fp);
 end_trace("energy_after_remap", "field_new_PM_restore.c", 181); } return 0; } 





static int movie_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += 0.1);   *ip = i; *tp = t;   return ret; } static int movie_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t <= TEND);   *ip = i; *tp = t;   return ret; } static int movie (const int i, const double t, Event * _ev) { trace ("movie", "field_new_PM_restore.c", 187); 
{
  char s[80];
  view ((struct _view_set){.fov = 20, .quat = {0.475152,0.161235,0.235565,0.832313}, .width = 800, .height = 600});
  sprintf (s, "t = %.2f", t);
  draw_string ((struct _draw_string){s, .size = 30});
  sprintf (s, "u%d.x", nl-1);
  squares ((struct _squares){s, .linear = true, .z = "eta", .min = -1.6*sqrt(1./kp_), .max = 1.6*sqrt(1./kp_)});
  {
  static FILE * fp =NULL; strongif (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("ux" ".ppm", "a");
  save ((struct _save){.fp = fp});
  }
#line 211 "field_new_PM_restore.c"
  char filename1[50], filename2[50], filename3[50];
  sprintf (filename1, "surface/eta_matrix_%g", t);
  sprintf (filename2, "surface/ux_matrix_%g", t);
  sprintf (filename3, "surface/uy_matrix_%g", t);
  FILE * feta = fopen (filename1, "w");

  output_matrix_mpi ((struct OutputMatrix){eta, feta, N, .linear = true});
  fclose (feta);
  sprintf (s, "u%d", nl-1);
  vector u_temp = lookup_vector (s);
  FILE * fux = fopen (filename2, "w");
  output_matrix_mpi ((struct OutputMatrix){u_temp.x, fux, N, .linear = true});
  fclose (fux);
  FILE * fuy = fopen (filename3, "w");
  output_matrix_mpi ((struct OutputMatrix){u_temp.y, fuy, N, .linear = true});
  fclose (fuy);
 end_trace("movie", "field_new_PM_restore.c", 227); } return 0; } 


int writefields (double t, const char *suffix) {
  char s[80];
  char filename1[50], filename2[50], filename3[50], filename4[50];
  vector u_temp;
  scalar w_temp, h_temp;
  for (int j=0; j<nl; ++j) {
    sprintf (filename1, "field/ux_%s_t%g_l%d", suffix, t, j);
    sprintf (filename2, "field/uy_%s_t%g_l%d", suffix, t, j);
    sprintf (filename3, "field/uz_%s_t%g_l%d", suffix, t, j);
    sprintf (filename4, "field/h_%s_t%g_l%d", suffix, t, j);
    if (j==0) {

      sprintf (s, "u");
      u_temp = lookup_vector (s);
      sprintf (s, "w");
      w_temp = lookup_field (s);
      sprintf (s, "h");
      h_temp = lookup_field (s);
    }
    else {
      sprintf (s, "u%d", j);
      u_temp = lookup_vector (s);
      sprintf (s, "w%d", j);
      w_temp = lookup_field (s);
      sprintf (s, "h%d", j);
      h_temp = lookup_field (s);
    }
    FILE * fux = fopen (filename1, "w");
    output_matrix_mpi ((struct OutputMatrix){u_temp.x, fux, N, .linear = true});
    fclose (fux);
    FILE * fuy = fopen (filename2, "w");
    output_matrix_mpi ((struct OutputMatrix){u_temp.y, fuy, N, .linear = true});
    fclose (fuy);
    FILE * fuz = fopen (filename3, "w");
    output_matrix_mpi ((struct OutputMatrix){w_temp, fuz, N, .linear = true});
    fclose (fuz);
    FILE * fh = fopen (filename4, "w");
    output_matrix_mpi ((struct OutputMatrix){h_temp, fh, N, .linear = true});
    fclose (fh);
  }
  return 0;
}

static int output_before_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int output_before (const int i, const double t, Event * _ev) { trace ("output_before", "field_new_PM_restore.c", 273);  {



  if (phony == 0) {
    j = i;
    phony = 1;
    fprintf (ferr, "Index starts at to j = %d!\n", j);
    char *suffix = "matrix_before";
    writefields (t, suffix);
  }
  if (i == j+1) {
    char *suffix = "matrix_after";
    writefields (t, suffix);
  }
 end_trace("output_before", "field_new_PM_restore.c", 288); } return 0; } 
#line 309 "field_new_PM_restore.c"
static int endrun_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = TEND);   *ip = i; *tp = t;   return ret; } static int endrun (const int i, const double t, Event * _ev) { trace ("endrun", "field_new_PM_restore.c", 309);  {
  dump ((struct Dump){0});
 end_trace("endrun", "field_new_PM_restore.c", 311); } return 0; } 
size_t datasize = 1*sizeof (double);
static int defaults0 (const int i, const double t, Event * _ev);
static int defaults0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup (const int i, const double t, Event * _ev);
static int cleanup_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int face_fields (const int i, const double t, Event * _ev);
static int face_fields_expr0 (int * ip, double * tp, Event * _ev);
static int half_advection (const int i, const double t, Event * _ev);
static int half_advection_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int pressure (const int i, const double t, Event * _ev);
static int pressure_expr0 (int * ip, double * tp, Event * _ev);
static int update_eta (const int i, const double t, Event * _ev);
static int update_eta_expr0 (int * ip, double * tp, Event * _ev);
static int remap (const int i, const double t, Event * _ev);
static int remap_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup_0 (const int i, const double t, Event * _ev);
static int cleanup_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults0_0 (const int i, const double t, Event * _ev);
static int defaults0_0_expr0 (int * ip, double * tp, Event * _ev);
static int half_advection_0 (const int i, const double t, Event * _ev);
static int half_advection_0_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int pressure_0 (const int i, const double t, Event * _ev);
static int pressure_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term_0 (const int i, const double t, Event * _ev);
static int viscous_term_0_expr0 (int * ip, double * tp, Event * _ev);
static int pressure_1 (const int i, const double t, Event * _ev);
static int pressure_1_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup_1 (const int i, const double t, Event * _ev);
static int cleanup_1_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_2 (const int i, const double t, Event * _ev);
static int defaults_2_expr0 (int * ip, double * tp, Event * _ev);
static int remap_0 (const int i, const double t, Event * _ev);
static int remap_0_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup_2 (const int i, const double t, Event * _ev);
static int cleanup_2_expr0 (int * ip, double * tp, Event * _ev);
static int perfs (const int i, const double t, Event * _ev);
static int perfs_expr0 (int * ip, double * tp, Event * _ev);
static int perf_plot (const int i, const double t, Event * _ev);
static int perf_plot_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int energy_before_remap (const int i, const double t, Event * _ev);
static int energy_before_remap_expr0 (int * ip, double * tp, Event * _ev);
static int energy_after_remap (const int i, const double t, Event * _ev);
static int energy_after_remap_expr0 (int * ip, double * tp, Event * _ev);
static int movie (const int i, const double t, Event * _ev);
static int movie_expr0 (int * ip, double * tp, Event * _ev);
static int movie_expr1 (int * ip, double * tp, Event * _ev);
static int output_before (const int i, const double t, Event * _ev);
static int output_before_expr0 (int * ip, double * tp, Event * _ev);
static int endrun (const int i, const double t, Event * _ev);
static int endrun_expr0 (int * ip, double * tp, Event * _ev);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults0, {defaults0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 87, "defaults0"});
  event_register ((Event){ 0, 1, defaults0_0, {defaults0_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/implicit.h", 34, "defaults0"});
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/run.h", 42, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 119, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/nh.h", 56, "defaults"});
  event_register ((Event){ 0, 1, defaults_2, {defaults_2_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/remap.h", 25, "defaults"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 159, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 96, "init"});
  event_register ((Event){ 0, 1, perfs, {perfs_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/perfs.h", 7, "perfs"});
  event_register ((Event){ 0, 1, perf_plot, {perf_plot_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/perfs.h", 30, "perf_plot"});
  event_register ((Event){ 0, 2, movie, {movie_expr0, movie_expr1}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 187, "movie"});
  event_register ((Event){ 0, 1, output_before, {output_before_expr0}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 273, "output_before"});
  event_register ((Event){ 0, 1, endrun, {endrun_expr0}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 309, "endrun"});
  event_register ((Event){ 0, 1, cleanup, {cleanup_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/run.h", 50, "cleanup"});
  event_register ((Event){ 0, 1, cleanup_0, {cleanup_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 466, "cleanup"});
  event_register ((Event){ 0, 1, cleanup_1, {cleanup_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/nh.h", 458, "cleanup"});
  event_register ((Event){ 0, 1, cleanup_2, {cleanup_2_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/remap.h", 119, "cleanup"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 173, "set_dtmax"});
  event_register ((Event){ 0, 1, face_fields, {face_fields_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 196, "face_fields"});
  event_register ((Event){ 0, 1, half_advection, {half_advection_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 382, "half_advection"});
  event_register ((Event){ 0, 1, half_advection_0, {half_advection_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/implicit.h", 115, "half_advection"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/diffusion.h", 144, "viscous_term"});
  event_register ((Event){ 0, 1, viscous_term_0, {viscous_term_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/nh.h", 76, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 396, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/implicit.h", 124, "acceleration"});
  event_register ((Event){ 0, 1, pressure, {pressure_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 398, "pressure"});
  event_register ((Event){ 0, 1, pressure_0, {pressure_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/implicit.h", 184, "pressure"});
  event_register ((Event){ 0, 1, pressure_1, {pressure_1_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/nh.h", 305, "pressure"});
  event_register ((Event){ 0, 1, update_eta, {update_eta_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 438, "update_eta"});
  event_register ((Event){ 0, 1, remap, {remap_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/hydro.h", 454, "remap"});
  event_register ((Event){ 0, 1, remap_0, {remap_0_expr0}, ((int *)0), ((double *)0),
    "/home/jiarongw/basilisk/src/layered/remap.h", 111, "remap"});
  event_register ((Event){ 0, 1, energy_before_remap, {energy_before_remap_expr0}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 137, "energy_before_remap"});
  event_register ((Event){ 0, 1, energy_after_remap, {energy_after_remap_expr0}, ((int *)0), ((double *)0),
    "field_new_PM_restore.c", 160, "energy_after_remap"});
  void allocate_globals (int);
  allocate_globals (1);
  set_fpe();
  multigrid_methods();
  init_scalar ((scalar){0}, "zb");
  init_const_scalar ((scalar){_NVARMAX+5}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+4}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+2},{_NVARMAX+3}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1}}, "zerof", (double []) {0.,0.,0.});
}
