#include "newlib.h"

real paco(void){

  return 3.9;
}

int newGetInterval( real x , real *G , int I , int lid ) /* gav */
{
  int  il, im, ih;
  real xm;
  
  if( x<G[0]   ){ return( -2 ); }
  if( x==G[0]  ){ return(  0 ); }
  if( x>G[I-1] ){ return(-101); }
  if( x==G[I-1]){ return( I-2); }
  
  il=0;
  ih=I;
  
  while( ih-il > 1 ){
    
    im = (il+ih) >> 1;
    xm = G[im];
    
    if(x < xm)
      ih = im;
    else if (x > xm)
      il = im;
    else
      return(im);
  }
  return il;
}

/* c=a*b (a, b, c = 4x4) */
void prodm4x4(real *c, real *a, real *b) /* gav */
{
  int i, j, k;
  
  for(i=0; i<4; i++)
    for(k=0; k<4; k++)
      for(j=0; j<4; j++)
	*(c+i*4+j) += *(a+i*4+k) * *(b+k*4+j);
}

/* c=a*b (a=mxn; b=nxp; c=mxp) */
void prodm(real *c, real *a, real *b, int m, int n, int p) /* gav */
{
  int i, j, k;

  for(i=0; i<m; i++)
    for(k=0; k<n; k++)
      for(j=0; j<p; j++)
	*(c+i*p+j) += *(a+i*n+k) * *(b+k*p+j);
}
/* c=a*b (a=mxn; b=nxp; c=mxp) */
void prodm_c(real *c, real *a, real *b, int m, int n, int p) /* gav */
{
  int i, j, k;
  
  /* inicializa */
  for(i=0; i<m; i++)
    for(j=0; j<p; j++)
      *(c+i*p+j) = 0;

  /* multiplica */
  for(i=0; i<m; i++)
    for(k=0; k<n; k++)
      for(j=0; j<p; j++)
	*(c+i*p+j) += *(a+i*n+k) * *(b+k*p+j);
}

void prodvm3x3(real *c, real *v, real *a)
{
  c[0] = v[0] * a[0] + v[1] * a[3] + v[2] * a[6] ;
  c[1] = v[0] * a[1] + v[1] * a[4] + v[2] * a[7] ;
  c[2] = v[0] * a[2] + v[1] * a[5] + v[2] * a[8] ;
}

/* c=a*b+v (a=mxn; b=nxp; v=mxp; c=mxp) */
void prodmpv(real *c, real *a, real *b, real *v, int m, int n, int p) /* gav */
{
  int i, j, k;

  for(i=0; i<m; i++)
    for(k=0; k<n; k++)
      for(j=0; j<p; j++)
	*(c+i*p+j) += *(a+i*n+k) * *(b+k*p+j);
  
  for(i=0; i<m*p; i++) 
    *(c+i) += *(v+i);
  
}

/* asi la usaremos tambien para el evol4D llamandola dos veces consecutivas */
void interpv3(real *vel, real *CVs, real *uu, real norm)
{
  int i, j;
  int dim_CVs, dimvel;
  
  dim_CVs = 8;
  dimvel = 3;

  for(i=0; i<dimvel; i++){
    for(j=0; j<dim_CVs; j++){
      vel[i] += CVs[i*dim_CVs+j] * uu[j];
    }
  }

  for(i=0; i<dimvel; i++){
    vel[i] *= norm;
  }

}

/* con esta nos evitamos tener que multiplicar y sumar 12 ceros 
pero a cambio tenermos que ir a buscar el indice
*/
void interpv2(real *vel, real *CVs, real *uu, real norm)
{
  int i, j;
  int dim_CVs, dim_CVsmed, dimvel;
  
  dim_CVs = 8;   dim_CVsmed = 4;
  dimvel = 2;

  for(i=0; i<dimvel; i++){
    for(j=0; j<dim_CVsmed; j++){
      vel[i] += CVs[i*dim_CVs+j] * uu[j];
    }
  }

  for(i=0; i<dimvel; i++){
    vel[i] *= norm;
  }
  vel[2] = 0.;
  
}


void print_matrix(char *mname, real *mx, int m, int n) /* gav */
{
  int i, j, index;

  /* printing matrices */
  /* mexPrintf("\nMatrix %s =[\n", mname); */
  mexPrintf("\n%s =[\n", mname);
  for(i=0; i<m; i++) 
  {
      for(j=0; j<n; j++) 
      {
	index = i*n + j;
	mexPrintf("%+-lf ", mx[index]);
      }
      if(i!=m-1)
	mexPrintf(";\n");
      else
	mexPrintf("]\n");
  }
}
void print_matrixt(char *mname, real *mx, int m, int n) /* gav */
{
  int i, j, index;

  /* mexPrintf("\nMatrixt %s =[\n", mname); */
  mexPrintf("\n%s =[\n", mname);
  for(i=0; i<m; i++) 
  {
      for(j=0; j<n; j++) 
      {
	index = i + j*m;
	mexPrintf("%+-lf ", mx[index]);
      }
      if(i!=m-1)
	mexPrintf(";\n");
      else
	mexPrintf("];\n");
  }
}

void diff_matrix(real *diffmx, real *mx1, real *mx2, int m, int n) /* gav */
{
  int i;

  for(i=0; i<m*n; i++){
    *(diffmx + i) = *(mx1 + i) - *(mx2 + i);    
  }

}

real det2x2(real *m) /* gav */
{

  real det;
  
  det = m[0]*m[3] - m[1]*m[2];

  return det;
}


real det3x3(real *m) /* gav */
{

  real det;
  det = 0;
  
  det += m[0] * (m[4]*m[8] - m[5]*m[7]);
  det -= m[1] * (m[3]*m[8] - m[5]*m[6]);
  det += m[2] * (m[3]*m[7] - m[4]*m[6]);

  /* mexPrintf("det3x3=%lf\n", det); */

  return det;
}

/* matriz homogenea (last row: 0 0 0 1) */
real det4x4mh(real *m) /* gav */
{
  real det;
  
  det = 0;
  
  det -= m[0] * (m[5]*m[10] - m[6]*m[9]);
  det += m[1] * (m[4]*m[10] - m[6]*m[8]);
  det -= m[2] * (m[4]*m[ 9] - m[5]*m[8]);

  return det;
}

void inv3x3m(real *invm, real *m)
{
  real det, Inv_det;

  det = det3x3(m);
  Inv_det = 1./det;

  invm[0] = ( m[4]*m[8] - m[5]*m[7]) * Inv_det ;
  invm[1] = (-m[1]*m[8] + m[2]*m[7]) * Inv_det ;
  invm[2] = ( m[1]*m[5] - m[2]*m[4]) * Inv_det ;
  invm[3] = (-m[3]*m[8] + m[5]*m[6]) * Inv_det ;
  invm[4] = ( m[0]*m[8] - m[2]*m[6]) * Inv_det ;
  invm[5] = (-m[0]*m[5] + m[2]*m[3]) * Inv_det ;
  invm[6] = ( m[3]*m[7] - m[4]*m[6]) * Inv_det ;
  invm[7] = (-m[0]*m[7] + m[1]*m[6]) * Inv_det ;
  invm[8] = ( m[0]*m[4] - m[1]*m[3]) * Inv_det ;
}

void invt2x2m(real *invm, real *m)
{
  real det, Inv_det;

  det = det2x2(m);
  Inv_det = 1./det;

  invm[0] = (-m[3]) * Inv_det ;
  invm[2] = ( m[2]) * Inv_det ;
  invm[1] = ( m[1]) * Inv_det ;
  invm[3] = (-m[0]) * Inv_det ;
}

void invt3x3m(real *invm, real *m)
{
  real det, Inv_det;

  det = det3x3(m);
  Inv_det = 1./det;

  invm[0] = (-m[4]*m[8] + m[5]*m[7]) * Inv_det ;
  invm[1] = ( m[1]*m[8] - m[2]*m[7]) * Inv_det ;
  invm[2] = (-m[1]*m[5] + m[2]*m[4]) * Inv_det ;
  invm[3] = ( m[3]*m[8] - m[5]*m[6]) * Inv_det ;
  invm[4] = (-m[0]*m[8] + m[2]*m[6]) * Inv_det ;
  invm[5] = ( m[0]*m[5] - m[2]*m[3]) * Inv_det ;
  invm[6] = (-m[3]*m[7] + m[4]*m[6]) * Inv_det ;
  invm[7] = ( m[0]*m[7] - m[1]*m[6]) * Inv_det ;
  invm[8] = (-m[0]*m[4] + m[1]*m[3]) * Inv_det ;
}

/* inversa de la traspuesta de matriz homogenea 4x4, (A')^{-1}=(adjA')^{T}/det(A) */
void invt4x4mh(real *invm, real *m) /* gav */
{
  real det, Inv_det;

  det = det4x4mh(m);
  Inv_det=1./det;

  invm[ 0] = (-m[ 5]*m[10] + m[ 6]*m[ 9]) * Inv_det ;
  invm[ 1] = ( m[ 1]*m[10] - m[ 2]*m[ 9]) * Inv_det ;
  invm[ 2] = (-m[ 1]*m[ 6] + m[ 2]*m[ 5]) * Inv_det ;
  invm[ 3] = 0.0;
          
  invm[ 4] = ( m[ 4]*m[10] - m[ 6]*m[ 8]) * Inv_det ;
  invm[ 5] = (-m[ 0]*m[10] + m[ 2]*m[ 8]) * Inv_det ;
  invm[ 6] = ( m[ 0]*m[ 6] - m[ 2]*m[ 4]) * Inv_det ;
  invm[ 7] = 0.0;
          
  invm[ 8] = (-m[ 4]*m[ 9] + m[ 5]*m[ 8]) * Inv_det ;
  invm[ 9] = ( m[ 0]*m[ 9] - m[ 1]*m[ 8]) * Inv_det ;
  invm[10] = (-m[ 0]*m[ 5] + m[ 1]*m[ 4]) * Inv_det ;
  invm[11] = 0.0;
  
  invm[12] = ( 
           m[ 4]*(m[ 9]*m[14] - m[10]*m[13]) -
           m[ 8]*(m[ 5]*m[14] - m[ 6]*m[13]) + 
           m[12]*(m[ 5]*m[10] - m[ 6]*m[ 9]) ) * Inv_det ;
  invm[13] = (- 
           m[ 0]*(m[ 9]*m[14] - m[10]*m[13]) +
           m[ 8]*(m[ 1]*m[14] - m[ 2]*m[13]) - 
           m[12]*(m[ 1]*m[10] - m[ 2]*m[ 9]) ) * Inv_det ;
  invm[14] = ( 
           m[ 0]*(m[ 5]*m[14] - m[ 6]*m[13]) -
           m[ 4]*(m[ 1]*m[14] - m[ 2]*m[13]) + 
           m[12]*(m[ 1]*m[ 6] - m[ 2]*m[ 5]) ) * Inv_det ;
  invm[15] = 1.0;

}

/* *************************************** */


void checkm4x4ish(char *mname, real *mx, int nsd)  /* nsd: number spatial dimensions */ /* gav */
{
  int dim=4;
  
  /* mexPrintf("checking if %s is a homogeneous matix: ", mname); */
  
  if( nsd == 3 ){
    if( mx[3] != 0 || mx[7] != 0 || mx[11] != 0 || mx[15] != 1 ){
      
      mexPrintf("\nInvalid homogeneous matrix %dD.\n", nsd);
      print_matrix(mname, mx, dim, dim);
      mexErrMsgTxt("Exiting ...\n\n");
    }
  } 
  else if( nsd == 2 ){
    if( mx[ 2] != 0 || mx[ 3] != 0 || mx[ 6] != 0 || mx[ 7] != 0 || mx[ 8] != 0 || 
	mx[ 9] != 0 || mx[10] != 1 || mx[11] != 0 || mx[14] != 0 || mx[15] != 1 ){

      mexPrintf("\nInvalid homogeneous matrix %dD.\n", nsd);
      /* print_matrix(mname, mx, dim, dim); */
      mexErrMsgTxt("Exiting ...\n\n");
    }
  }
  else if( nsd == 1 ){
    if( mx[ 1] != 0 || mx[ 2] != 0 || mx[ 3] != 0 || mx[ 4] != 0 || mx[ 5] != 1 || mx[ 6] != 0 || mx[ 7] != 0 || 
	mx[ 8] != 0 || mx[ 9] != 0 || mx[10] != 1 || mx[11] != 0 || mx[13] != 0 || mx[14] != 0 || mx[15] != 1 ){

      mexPrintf("\nInvalid homogeneous matrix %dD.\n", nsd);
      /* print_matrix(mname, mx, dim, dim); */
      mexErrMsgTxt("Exiting ...\n\n");
    }
  }
  else{
    mexPrintf("\nError: nsd=%d not valid (permited values 2 and 3 only).\n", nsd);
    mexErrMsgTxt("Exiting ...\n\n");
  }

  /* mexPrintf("ok\n"); */
}

int checkIsSortedNonrep_int(int *G, int I)
{
  int    i;
  
  if( I <= 1 ){ return(1); }
  for( i = 1 ; i < I ; i++ ){
    if( G[i-1] >= G[i] ){
      return(0);
    }
  }
  return(1);
}
int checkIsSortedNr( real *G, int I ){
  int    i;
  
  if( I <= 1 ){ return(1); }
  for( i = 1 ; i < I ; i++ ){
    if( G[i-1] >= G[i] ){
      return(0);
    }
  }
  return(1);
}



/*
 *x   :  8 elementos
 *y   :  8 elementos
 *z   :  8 elementos
 *IM  : 64 intensidades (4x4x4)
 *XYZ : pto donde interpolamos

     x x x x     
     x x x x 
x x  o o o o  x x
x x  o o o o  x x
x x  o o o o  x x
x x  o o o o  x x
     x x x x 
     x x x x     
*/


double bspline3(real *x, real *y, real *z, real *IM, real *XYZ)
{
  int i;
  real W[NW];
  real IMinter;
  real sumW;
  
  real X, Y, Z;
  real XX, YY, ZZ;
  real XXX, YYY, ZZZ;
  real dX01, dX02, dX03, dX04, dX05, dX06, dX07, dX12, dX13, dX14, dX15, dX16, dX17, dX23, dX24, dX25, dX26, dX27, dX34, dX35, dX36, dX37, dX45, dX46, dX47, dX56, dX57, dX67; 
  real dY01, dY02, dY03, dY04, dY05, dY06, dY07, dY12, dY13, dY14, dY15, dY16, dY17, dY23, dY24, dY25, dY26, dY27, dY34, dY35, dY36, dY37, dY45, dY46, dY47, dY56, dY57, dY67; 
  real dZ01, dZ02, dZ03, dZ04, dZ05, dZ06, dZ07, dZ12, dZ13, dZ14, dZ15, dZ16, dZ17, dZ23, dZ24, dZ25, dZ26, dZ27, dZ34, dZ35, dZ36, dZ37, dZ45, dZ46, dZ47, dZ56, dZ57, dZ67; 
  real dX3, dX4; 
  real dY3, dY4; 
  real dZ3, dZ4; 

  /* // calcula weights */
  X=XYZ[0];
  Y=XYZ[1];
  Z=XYZ[2];
  
  dX3 = X - x[3]; dX4 = x[4] - X;
  dX01 = x[1]-x[0]; dX02 = x[2]-x[0]; dX03 = x[3]-x[0]; dX04 = x[4]-x[0]; dX05 = x[5]-x[0]; dX06 = x[6]-x[0]; dX07 = x[7]-x[0];
  dX12 = x[2]-x[1]; dX13 = x[3]-x[1]; dX14 = x[4]-x[1]; dX15 = x[5]-x[1]; dX16 = x[6]-x[1]; dX17 = x[7]-x[1];
  dX23 = x[3]-x[2]; dX24 = x[4]-x[2]; dX25 = x[5]-x[2]; dX26 = x[6]-x[2]; dX27 = x[7]-x[2];
  dX34 = x[4]-x[3]; dX35 = x[5]-x[3]; dX36 = x[6]-x[3]; dX37 = x[7]-x[3];
  dX45 = x[5]-x[4]; dX46 = x[6]-x[4]; dX47 = x[7]-x[4];
  dX56 = x[6]-x[5]; dX57 = x[7]-x[5];
  dX67 = x[7]-x[6];

  dY3 = Y - y[3]; dY4 = y[4] - Y;
  dY01 = y[1]-y[0]; dY02 = y[2]-y[0]; dY03 = y[3]-y[0]; dY04 = y[4]-y[0]; dY05 = y[5]-y[0]; dY06 = y[6]-y[0]; dY07 = y[7]-y[0];
  dY12 = y[2]-y[1]; dY13 = y[3]-y[1]; dY14 = y[4]-y[1]; dY15 = y[5]-y[1]; dY16 = y[6]-y[1]; dY17 = y[7]-y[1];
  dY23 = y[3]-y[2]; dY24 = y[4]-y[2]; dY25 = y[5]-y[2]; dY26 = y[6]-y[2]; dY27 = y[7]-y[2];
  dY34 = y[4]-y[3]; dY35 = y[5]-y[3]; dY36 = y[6]-y[3]; dY37 = y[7]-y[3];
  dY45 = y[5]-y[4]; dY46 = y[6]-y[4]; dY47 = y[7]-y[4];
  dY56 = y[6]-y[5]; dY57 = y[7]-y[5];
  dY67 = y[7]-y[6];

  dZ3 = Z - z[3]; dZ4 = z[4] - Z;
  dZ01 = z[1]-z[0]; dZ02 = z[2]-z[0]; dZ03 = z[3]-z[0]; dZ04 = z[4]-z[0]; dZ05 = z[5]-z[0]; dZ06 = z[6]-z[0]; dZ07 = z[7]-z[0];
  dZ12 = z[2]-z[1]; dZ13 = z[3]-z[1]; dZ14 = z[4]-z[1]; dZ15 = z[5]-z[1]; dZ16 = z[6]-z[1]; dZ17 = z[7]-z[1];
  dZ23 = z[3]-z[2]; dZ24 = z[4]-z[2]; dZ25 = z[5]-z[2]; dZ26 = z[6]-z[2]; dZ27 = z[7]-z[2];
  dZ34 = z[4]-z[3]; dZ35 = z[5]-z[3]; dZ36 = z[6]-z[3]; dZ37 = z[7]-z[3];
  dZ45 = z[5]-z[4]; dZ46 = z[6]-z[4]; dZ47 = z[7]-z[4];
  dZ56 = z[6]-z[5]; dZ57 = z[7]-z[5];
  dZ67 = z[7]-z[6];

  XX = X*X;    XXX = XX*X;
  YY = Y*Y;    YYY = YY*Y;
  ZZ = Z*Z;    ZZZ = ZZ*Z;

  /* // weights */
  /* // ---------------------------------------------------------- */
  W[0]=(8*dX13*pow(dX4,3)*dY13*pow(dY4,3)*dZ13*pow(dZ4,3))/(dX04*dX14*dX24*dX34*dY04*dY14*dY24*dY34*dZ04*dZ14*dZ24*dZ34);
    
    W[1]=(8*dY13*pow(dY4,3)*dZ13*pow(dZ4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX))/(dX14*dX15*dX25*dX34*dX35*dY04*dY14*dY24*dY34*dZ04*dZ14*dZ24*dZ34);

    W[2]=(-8*dY13*pow(dY4,3)*dZ13*pow(dZ4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX))/(dX24*dX25*dX26*dX34*dX36*dY04*dY14*dY24*dY34*dZ04*dZ14*dZ24*dZ34);

    W[3]=(8*pow(dX3,3)*dX46*dY13*pow(dY4,3)*dZ13*pow(dZ4,3))/(dX34*dX35*dX36*dX37*dY04*dY14*dY24*dY34*dZ04*dZ14*dZ24*dZ34);

    W[4]=(8*dX13*pow(dX4,3)*dZ13*pow(dZ4,3)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX04*dX14*dX24*dX34*dY14*dY15*dY25*dY34*dY35*dZ04*dZ14*dZ24*dZ34);

    W[5]=(8*dZ13*pow(dZ4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX14*dX15*dX25*dX34*dX35*dY14*dY15*dY25*dY34*dY35*dZ04*dZ14*dZ24*dZ34);

    W[6]=(-8*dZ13*pow(dZ4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX24*dX25*dX26*dX34*dX36*dY14*dY15*dY25*dY34*dY35*dZ04*dZ14*dZ24*dZ34);

    W[7]=(8*pow(dX3,3)*dX46*dZ13*pow(dZ4,3)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX34*dX35*dX36*dX37*dY14*dY15*dY25*dY34*dY35*dZ04*dZ14*dZ24*dZ34);

    W[8]=(-8*dX13*pow(dX4,3)*dZ13*pow(dZ4,3)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX04*dX14*dX24*dX34*dY24*dY25*dY26*dY34*dY36*dZ04*dZ14*dZ24*dZ34);

    W[9]=(-8*dZ13*pow(dZ4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX14*dX15*dX25*dX34*dX35*dY24*dY25*dY26*dY34*dY36*dZ04*dZ14*dZ24*dZ34);

    W[10]=(8*dZ13*pow(dZ4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX24*dX25*dX26*dX34*dX36*dY24*dY25*dY26*dY34*dY36*dZ04*dZ14*dZ24*dZ34);

    W[11]=(-8*pow(dX3,3)*dX46*dZ13*pow(dZ4,3)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX34*dX35*dX36*dX37*dY24*dY25*dY26*dY34*dY36*dZ04*dZ14*dZ24*dZ34);

    W[12]=(8*dX13*pow(dX4,3)*pow(dY3,3)*dY46*dZ13*pow(dZ4,3))/(dX04*dX14*dX24*dX34*dY34*dY35*dY36*dY37*dZ04*dZ14*dZ24*dZ34);

    W[13]=(8*pow(dY3,3)*dY46*dZ13*pow(dZ4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX))/(dX14*dX15*dX25*dX34*dX35*dY34*dY35*dY36*dY37*dZ04*dZ14*dZ24*dZ34);

    W[14]=(-8*pow(dY3,3)*dY46*dZ13*pow(dZ4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX))/(dX24*dX25*dX26*dX34*dX36*dY34*dY35*dY36*dY37*dZ04*dZ14*dZ24*dZ34);

    W[15]=(8*pow(dX3,3)*dX46*pow(dY3,3)*dY46*dZ13*pow(dZ4,3))/(dX34*dX35*dX36*dX37*dY34*dY35*dY36*dY37*dZ04*dZ14*dZ24*dZ34);

    W[16]=(8*dX13*pow(dX4,3)*dY13*pow(dY4,3)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX04*dX14*dX24*dX34*dY04*dY14*dY24*dY34*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[17]=(8*dY13*pow(dY4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY04*dY14*dY24*dY34*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[18]=(-8*dY13*pow(dY4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY04*dY14*dY24*dY34*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[19]=(8*pow(dX3,3)*dX46*dY13*pow(dY4,3)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX34*dX35*dX36*dX37*dY04*dY14*dY24*dY34*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[20]=(8*dX13*pow(dX4,3)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX04*dX14*dX24*dX34*dY14*dY15*dY25*dY34*dY35*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[21]=(8*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY14*dY15*dY25*dY34*dY35*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[22]=(-8*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY14*dY15*dY25*dY34*dY35*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[23]=(8*pow(dX3,3)*dX46*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX34*dX35*dX36*dX37*dY14*dY15*dY25*dY34*dY35*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[24]=(-8*dX13*pow(dX4,3)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX04*dX14*dX24*dX34*dY24*dY25*dY26*dY34*dY36*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[25]=(-8*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY24*dY25*dY26*dY34*dY36*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[26]=(8*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY24*dY25*dY26*dY34*dY36*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[27]=(-8*pow(dX3,3)*dX46*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX34*dX35*dX36*dX37*dY24*dY25*dY26*dY34*dY36*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[28]=(8*dX13*pow(dX4,3)*pow(dY3,3)*dY46*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX04*dX14*dX24*dX34*dY34*dY35*dY36*dY37*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[29]=(8*pow(dY3,3)*dY46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY34*dY35*dY36*dY37*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[30]=(-8*pow(dY3,3)*dY46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY34*dY35*dY36*dY37*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[31]=(8*pow(dX3,3)*dX46*pow(dY3,3)*dY46*(z[2]*z[3]*pow(z[4],2)*z[5] + z[2]*z[3]*z[4]*pow(z[5],2) - z[2]*pow(z[4],2)*pow(z[5],2) - z[3]*pow(z[4],2)*pow(z[5],2) + 3*Z*(z[4]*z[5]*(-(z[2]*z[3]) + z[4]*z[5]) + z[1]*(-(z[3]*z[4]*z[5]) + z[2]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])))) - 3*dZ14*z[4]*z[5]*ZZ + 3*z[2]*z[4]*z[5]*ZZ + 3*z[3]*z[4]*z[5]*ZZ - 3*z[4]*pow(z[5],2)*ZZ + z[1]*(z[4]*z[5]*(-(z[4]*z[5]) + z[3]*(z[4] + z[5])) - z[2]*(-(z[4]*z[5]*(z[4] + z[5])) + z[3]*(pow(z[4],2) + z[4]*z[5] + pow(z[5],2) + 3*ZZ))) + dZ14*dZ24*ZZZ - dZ14*z[3]*ZZZ + z[2]*z[3]*ZZZ + dZ14*z[5]*ZZZ - z[2]*z[5]*ZZZ - z[3]*z[5]*ZZZ + pow(z[5],2)*ZZZ))/(dX34*dX35*dX36*dX37*dY34*dY35*dY36*dY37*dZ14*dZ15*dZ25*dZ34*dZ35);

    W[32]=(-8*dX13*pow(dX4,3)*dY13*pow(dY4,3)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX04*dX14*dX24*dX34*dY04*dY14*dY24*dY34*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[33]=(-8*dY13*pow(dY4,3)*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY04*dY14*dY24*dY34*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[34]=(8*dY13*pow(dY4,3)*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY04*dY14*dY24*dY34*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[35]=(-8*pow(dX3,3)*dX46*dY13*pow(dY4,3)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX34*dX35*dX36*dX37*dY04*dY14*dY24*dY34*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[36]=(-8*dX13*pow(dX4,3)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX04*dX14*dX24*dX34*dY14*dY15*dY25*dY34*dY35*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[37]=(-8*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY14*dY15*dY25*dY34*dY35*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[38]=(8*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY14*dY15*dY25*dY34*dY35*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[39]=(-8*pow(dX3,3)*dX46*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX34*dX35*dX36*dX37*dY14*dY15*dY25*dY34*dY35*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[40]=(8*dX13*pow(dX4,3)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX04*dX14*dX24*dX34*dY24*dY25*dY26*dY34*dY36*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[41]=(8*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY24*dY25*dY26*dY34*dY36*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[42]=(-8*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY24*dY25*dY26*dY34*dY36*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[43]=(8*pow(dX3,3)*dX46*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX34*dX35*dX36*dX37*dY24*dY25*dY26*dY34*dY36*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[44]=(-8*dX13*pow(dX4,3)*pow(dY3,3)*dY46*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX04*dX14*dX24*dX34*dY34*dY35*dY36*dY37*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[45]=(-8*pow(dY3,3)*dY46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX14*dX15*dX25*dX34*dX35*dY34*dY35*dY36*dY37*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[46]=(8*pow(dY3,3)*dY46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX24*dX25*dX26*dX34*dX36*dY34*dY35*dY36*dY37*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[47]=(-8*pow(dX3,3)*dX46*pow(dY3,3)*dY46*(-(pow(z[3],2)*z[4]*z[5]*z[6]) - pow(z[2],2)*(z[4]*z[5]*z[6] + pow(z[3],2)*(z[4] + z[5] + z[6]) - z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]))) + 3*Z*(pow(z[2],2)*pow(z[3],2) + z[3]*z[4]*z[5]*z[6] - z[2]*(-(z[4]*z[5]*z[6]) + z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6])))) - 3*z[4]*z[5]*z[6]*ZZ + z[2]*z[3]*(z[3]*(z[5]*z[6] + z[4]*(z[5] + z[6]) - 3*ZZ) + 3*(dZ25 + z[6])*ZZ + z[4]*(-(z[5]*z[6]) + 3*ZZ)) + dZ25*dZ26*ZZZ - dZ25*z[3]*ZZZ + pow(z[3],2)*ZZZ + dZ25*z[4]*ZZZ - z[3]*z[4]*ZZZ - z[3]*z[6]*ZZZ + z[4]*z[6]*ZZZ))/(dX34*dX35*dX36*dX37*dY34*dY35*dY36*dY37*dZ24*dZ25*dZ26*dZ34*dZ36);

    W[48]=(8*dX13*pow(dX4,3)*dY13*pow(dY4,3)*pow(dZ3,3)*dZ46)/(dX04*dX14*dX24*dX34*dY04*dY14*dY24*dY34*dZ34*dZ35*dZ36*dZ37);

    W[49]=(8*dY13*pow(dY4,3)*pow(dZ3,3)*dZ46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX))/(dX14*dX15*dX25*dX34*dX35*dY04*dY14*dY24*dY34*dZ34*dZ35*dZ36*dZ37);

    W[50]=(-8*dY13*pow(dY4,3)*pow(dZ3,3)*dZ46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX))/(dX24*dX25*dX26*dX34*dX36*dY04*dY14*dY24*dY34*dZ34*dZ35*dZ36*dZ37);

    W[51]=(8*pow(dX3,3)*dX46*dY13*pow(dY4,3)*pow(dZ3,3)*dZ46)/(dX34*dX35*dX36*dX37*dY04*dY14*dY24*dY34*dZ34*dZ35*dZ36*dZ37);

    W[52]=(8*dX13*pow(dX4,3)*pow(dZ3,3)*dZ46*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX04*dX14*dX24*dX34*dY14*dY15*dY25*dY34*dY35*dZ34*dZ35*dZ36*dZ37);

    W[53]=(8*pow(dZ3,3)*dZ46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX14*dX15*dX25*dX34*dX35*dY14*dY15*dY25*dY34*dY35*dZ34*dZ35*dZ36*dZ37);

    W[54]=(-8*pow(dZ3,3)*dZ46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX24*dX25*dX26*dX34*dX36*dY14*dY15*dY25*dY34*dY35*dZ34*dZ35*dZ36*dZ37);

    W[55]=(8*pow(dX3,3)*dX46*pow(dZ3,3)*dZ46*(y[2]*y[3]*pow(y[4],2)*y[5] + y[2]*y[3]*y[4]*pow(y[5],2) - y[2]*pow(y[4],2)*pow(y[5],2) - y[3]*pow(y[4],2)*pow(y[5],2) + 3*Y*(y[4]*y[5]*(-(y[2]*y[3]) + y[4]*y[5]) + y[1]*(-(y[3]*y[4]*y[5]) + y[2]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])))) - 3*dY14*y[4]*y[5]*YY + 3*y[2]*y[4]*y[5]*YY + 3*y[3]*y[4]*y[5]*YY - 3*y[4]*pow(y[5],2)*YY + y[1]*(y[4]*y[5]*(-(y[4]*y[5]) + y[3]*(y[4] + y[5])) - y[2]*(-(y[4]*y[5]*(y[4] + y[5])) + y[3]*(pow(y[4],2) + y[4]*y[5] + pow(y[5],2) + 3*YY))) + dY14*dY24*YYY - dY14*y[3]*YYY + y[2]*y[3]*YYY + dY14*y[5]*YYY - y[2]*y[5]*YYY - y[3]*y[5]*YYY + pow(y[5],2)*YYY))/(dX34*dX35*dX36*dX37*dY14*dY15*dY25*dY34*dY35*dZ34*dZ35*dZ36*dZ37);

    W[56]=(-8*dX13*pow(dX4,3)*pow(dZ3,3)*dZ46*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX04*dX14*dX24*dX34*dY24*dY25*dY26*dY34*dY36*dZ34*dZ35*dZ36*dZ37);

    W[57]=(-8*pow(dZ3,3)*dZ46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX14*dX15*dX25*dX34*dX35*dY24*dY25*dY26*dY34*dY36*dZ34*dZ35*dZ36*dZ37);

    W[58]=(8*pow(dZ3,3)*dZ46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX)*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX24*dX25*dX26*dX34*dX36*dY24*dY25*dY26*dY34*dY36*dZ34*dZ35*dZ36*dZ37);

    W[59]=(-8*pow(dX3,3)*dX46*pow(dZ3,3)*dZ46*(-(pow(y[3],2)*y[4]*y[5]*y[6]) - pow(y[2],2)*(y[4]*y[5]*y[6] + pow(y[3],2)*(y[4] + y[5] + y[6]) - y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]))) + 3*Y*(pow(y[2],2)*pow(y[3],2) + y[3]*y[4]*y[5]*y[6] - y[2]*(-(y[4]*y[5]*y[6]) + y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6])))) - 3*y[4]*y[5]*y[6]*YY + y[2]*y[3]*(y[3]*(y[5]*y[6] + y[4]*(y[5] + y[6]) - 3*YY) + 3*(dY25 + y[6])*YY + y[4]*(-(y[5]*y[6]) + 3*YY)) + dY25*dY26*YYY - dY25*y[3]*YYY + pow(y[3],2)*YYY + dY25*y[4]*YYY - y[3]*y[4]*YYY - y[3]*y[6]*YYY + y[4]*y[6]*YYY))/(dX34*dX35*dX36*dX37*dY24*dY25*dY26*dY34*dY36*dZ34*dZ35*dZ36*dZ37);

    W[60]=(8*dX13*pow(dX4,3)*pow(dY3,3)*dY46*pow(dZ3,3)*dZ46)/(dX04*dX14*dX24*dX34*dY34*dY35*dY36*dY37*dZ34*dZ35*dZ36*dZ37);

    W[61]=(8*pow(dY3,3)*dY46*pow(dZ3,3)*dZ46*(x[2]*x[3]*pow(x[4],2)*x[5] + x[2]*x[3]*x[4]*pow(x[5],2) - x[2]*pow(x[4],2)*pow(x[5],2) - x[3]*pow(x[4],2)*pow(x[5],2) + 3*X*(x[4]*x[5]*(-(x[2]*x[3]) + x[4]*x[5]) + x[1]*(-(x[3]*x[4]*x[5]) + x[2]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])))) - 3*dX14*x[4]*x[5]*XX + 3*x[2]*x[4]*x[5]*XX + 3*x[3]*x[4]*x[5]*XX - 3*x[4]*pow(x[5],2)*XX + x[1]*(x[4]*x[5]*(-(x[4]*x[5]) + x[3]*(x[4] + x[5])) - x[2]*(-(x[4]*x[5]*(x[4] + x[5])) + x[3]*(pow(x[4],2) + x[4]*x[5] + pow(x[5],2) + 3*XX))) + dX14*dX24*XXX - dX14*x[3]*XXX + x[2]*x[3]*XXX + dX14*x[5]*XXX - x[2]*x[5]*XXX - x[3]*x[5]*XXX + pow(x[5],2)*XXX))/(dX14*dX15*dX25*dX34*dX35*dY34*dY35*dY36*dY37*dZ34*dZ35*dZ36*dZ37);

    W[62]=(-8*pow(dY3,3)*dY46*pow(dZ3,3)*dZ46*(-(pow(x[3],2)*x[4]*x[5]*x[6]) - pow(x[2],2)*(x[4]*x[5]*x[6] + pow(x[3],2)*(x[4] + x[5] + x[6]) - x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]))) + 3*X*(pow(x[2],2)*pow(x[3],2) + x[3]*x[4]*x[5]*x[6] - x[2]*(-(x[4]*x[5]*x[6]) + x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6])))) - 3*x[4]*x[5]*x[6]*XX + x[2]*x[3]*(x[3]*(x[5]*x[6] + x[4]*(x[5] + x[6]) - 3*XX) + 3*(dX25 + x[6])*XX + x[4]*(-(x[5]*x[6]) + 3*XX)) + dX25*dX26*XXX - dX25*x[3]*XXX + pow(x[3],2)*XXX + dX25*x[4]*XXX - x[3]*x[4]*XXX - x[3]*x[6]*XXX + x[4]*x[6]*XXX))/(dX24*dX25*dX26*dX34*dX36*dY34*dY35*dY36*dY37*dZ34*dZ35*dZ36*dZ37);

    W[63]=(8*pow(dX3,3)*dX46*pow(dY3,3)*dY46*pow(dZ3,3)*dZ46)/(dX34*dX35*dX36*dX37*dY34*dY35*dY36*dY37*dZ34*dZ35*dZ36*dZ37);
  /* // ---------------------------------------------------------- */

  sumW = 0.;
  for(i=0; i<NW; i++){
    sumW += W[i];
  }
  mexPrintf("sumW = %lf  (%d values)\n", sumW, NW);

  
  /* // interpola */
  /* mexPrintf("# i       W[i]        IM[i]        W[i]*IM[i]    IMinter\n"); */
  IMinter = 0.;
  for(i=0; i<NW; i++){
    IMinter += W[i] * IM[i];
    /* mexPrintf(" %02d  %+-12.6lf %+-12.6lf (%+-12.6lf): %+-12.6lf\n", i, W[i], IM[i], (W[i]*IM[i]),  IMinter); */
  }
  /* mexPrintf("Done: IMinter=%lf\n", IMinter); */

  return( IMinter );
}


/* *************************************** */

/* funciones ya definidas pero dentro de los codigos */

/* de InterpPointsOn3DGrid */
/* #define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L) */
/* No sacada a libreria pues en InterpOn3DGrid esta definida distinta! */
