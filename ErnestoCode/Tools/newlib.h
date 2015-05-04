#include "mex.h"
#include "math.h"

#ifndef  real 
#define  real   double
#endif

#define NDIMS_Odims  50
#define NW           64

/* funciones ya definidas pero dentro de los codigos */
/* de InterpPointsOn3DGrid */
/* #define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L) */
/* No sacada a libreria pues en InterpOn3DGrid esta definida distinta! */

  
/* functions */
int  newGetInterval( real , real *, int, int);
void prodvm3x3(real *, real *, real *);
void prodm4x4(real *, real *, real *);
void prodm   (real *, real *, real *, int, int, int);
void prodm_c (real *, real *, real *, int, int, int);
void prodmpv (real *, real *, real *, real *, int, int, int);
void interpv3(real *, real *, real *, real);
void interpv2(real *, real *, real *, real);
void print_matrix (char *, real *, int, int);
void print_matrixt(char *, real *, int, int);
void diff_matrix  (real *, real *, real *, int, int);
real paco(void);
real det2x2   (real *);
real det3x3   (real *);
real det4x4mh (real *);
void invt2x2m (real *, real *);
void invt3x3m (real *, real *);
void inv3x3m  (real *, real *);
void invt4x4mh(real *, real *);
void checkm4x4ish(char *, real *, int);
real bspline3(real *, real *, real *, real *, real *);
int checkIsSortedNonrep_int(int *, int);
int checkIsSortedNr(real *, int);
