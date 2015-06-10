#ifndef MZ_INV4X4_H
#define MZ_INV4X4_H


#if !defined( real )
  #define   real       double
#endif


/*a =
     0     2
     1     3*/
real det2x2(real *m) /* gav */
{

  real det;

  det = m[0]*m[3] - m[1]*m[2];
  /*mexPrintf("det2x2=%lf\n", det);*/

  return det;
}
/*a =
     0     3     6
     1     4     7
     2     5     8*/
real det3x3(real *m) /* gav */
{

  real det;
  det = 0;

  det += m[0] * (m[4]*m[8] - m[5]*m[7]);
  det -= m[1] * (m[3]*m[8] - m[5]*m[6]);
  det += m[2] * (m[3]*m[7] - m[4]*m[6]);

  /*mexPrintf("det3x3=%lf\n", det);*/

  return det;
}
/*a =
     0     4     8    12
     1     5     9    13
     2     6    10    14
     3     7    11    15*/
real det4x4m(real *m)
{
  real det;

  det = 0;

  det += m[ 0] * ( m[ 5]*m[10]*m[15] + m[ 7]*m[ 9]*m[14] + m[ 6]*m[11]*m[13] -
	  	   	   	   m[ 7]*m[10]*m[13] - m[ 6]*m[ 9]*m[15] - m[ 5]*m[11]*m[14] );
  det -= m[ 4] * ( m[ 1]*m[10]*m[15] + m[ 3]*m[ 9]*m[14] + m[ 2]*m[11]*m[13] -
 	   	   	   	   m[ 3]*m[10]*m[13] - m[ 2]*m[ 9]*m[15] - m[ 1]*m[11]*m[14] );
  det += m[ 8] * ( m[ 1]*m[ 6]*m[15] + m[ 3]*m[ 5]*m[14] + m[ 2]*m[ 7]*m[ 13] -
 	   	   	   	   m[ 3]*m[ 6]*m[13] - m[ 2]*m[ 5]*m[15] - m[ 1]*m[ 7]*m[14] );
  det -= m[12] * ( m[ 1]*m[ 6]*m[11] + m[ 3]*m[ 5]*m[10] + m[ 2]*m[ 7]*m[ 9] -
 	   	   	   	   m[ 3]*m[ 6]*m[ 9] - m[ 2]*m[ 5]*m[11] - m[ 1]*m[ 7]*m[10] );

  /*mexPrintf("det4x4=%lf\n", det);*/

  return det;
}
real det4x4mh(real *m)
{
  real det;

  det = 0;

  det -= m[0] * (m[5]*m[10] - m[6]*m[9]);
  det += m[1] * (m[4]*m[10] - m[6]*m[8]);
  det -= m[2] * (m[4]*m[ 9] - m[5]*m[8]);

  return det;
}
/*a =
     0     2
     1     3*/
/*
 * A^1 = (A^*)^t / |A|
 */
void invt2x2m(real *invm, real *m)
{
  real det, Inv_det;

  det = det2x2(m);
  Inv_det = 1./det;

  invm[0] = (m[3]) * Inv_det ;
  invm[1] = (-m[1]) * Inv_det ;
  invm[2] = (-m[2]) * Inv_det ;
  invm[3] = (m[0]) * Inv_det ;
}

void invt3x3m(real *invm, real *m)
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
void invt4x4m(real *invm, real *m)
{
  real det, Inv_det;

  det = det4x4m(m);
  Inv_det=1./det;

  invm[ 0] = ( m[ 5]*m[10]*m[15] + m[ 9]*m[ 7]*m[14] + m[ 6]*m[11]*m[13] -
		  	   m[ 7]*m[10]*m[13] - m[ 6]*m[ 9]*m[15] - m[ 5]*m[11]*m[14]) * Inv_det ;
  invm[ 1] =-( m[ 1]*m[10]*m[15] + m[ 3]*m[ 9]*m[14] + m[ 2]*m[11]*m[13] -
	  	   	   m[ 3]*m[10]*m[13] - m[ 2]*m[ 9]*m[15] - m[ 1]*m[11]*m[14]) * Inv_det ;
  invm[ 2] = ( m[ 1]*m[ 6]*m[15] + m[ 3]*m[ 5]*m[14] + m[ 2]*m[ 7]*m[13] -
	  	   	   m[ 3]*m[ 6]*m[13] - m[ 2]*m[ 5]*m[15] - m[ 1]*m[ 7]*m[14]) * Inv_det ;
  invm[ 3] =-( m[ 1]*m[ 6]*m[11] + m[ 3]*m[ 5]*m[10] + m[ 2]*m[ 7]*m[ 9] -
	  	   	   m[ 3]*m[ 6]*m[ 9] - m[ 2]*m[ 5]*m[11] - m[ 1]*m[ 7]*m[10]) * Inv_det ;

  invm[ 4] =-( m[ 4]*m[10]*m[15] + m[ 7]*m[ 8]*m[14] + m[ 6]*m[11]*m[12] -
		  	   m[ 7]*m[10]*m[12] - m[ 6]*m[ 8]*m[15] - m[ 4]*m[11]*m[14]) * Inv_det ;
  invm[ 5] = ( m[ 0]*m[10]*m[15] + m[ 3]*m[ 8]*m[14] + m[ 2]*m[11]*m[12] -
		  	   m[ 3]*m[10]*m[12] - m[ 2]*m[ 8]*m[15] - m[ 0]*m[11]*m[14]) * Inv_det ;
  invm[ 6] =-( m[ 0]*m[ 6]*m[15] + m[ 3]*m[ 4]*m[14] + m[ 2]*m[ 7]*m[12] -
		  	   m[ 3]*m[ 6]*m[12] - m[ 2]*m[ 4]*m[15] - m[ 0]*m[ 7]*m[14]) * Inv_det ;
  invm[ 7] = ( m[ 0]*m[ 6]*m[11] + m[ 3]*m[ 4]*m[10] + m[ 2]*m[ 7]*m[ 8] -
		  	   m[ 3]*m[ 6]*m[ 8] - m[ 2]*m[ 4]*m[11] - m[ 0]*m[ 7]*m[10]) * Inv_det ;

  invm[ 8] = ( m[ 4]*m[ 9]*m[15] + m[ 7]*m[ 8]*m[13] + m[ 5]*m[11]*m[12] -
		  	   m[ 7]*m[ 9]*m[12] - m[ 5]*m[ 8]*m[15] - m[ 4]*m[11]*m[13]) * Inv_det ;
  invm[ 9] =-( m[ 0]*m[ 9]*m[15] + m[ 3]*m[ 8]*m[13] + m[ 1]*m[11]*m[12] -
		  	   m[ 3]*m[ 9]*m[12] - m[ 1]*m[ 8]*m[15] - m[ 0]*m[11]*m[13]) * Inv_det ;
  invm[10] = ( m[ 0]*m[ 5]*m[15] + m[ 3]*m[ 4]*m[13] + m[ 1]*m[ 7]*m[12] -
		  	   m[ 3]*m[ 5]*m[12] - m[ 1]*m[ 4]*m[15] - m[ 0]*m[ 7]*m[13]) * Inv_det ;
  invm[11] =-( m[ 0]*m[ 5]*m[11] + m[ 3]*m[ 4]*m[ 9] + m[ 1]*m[ 7]*m[ 8] -
		  	   m[ 3]*m[ 5]*m[ 8] - m[ 1]*m[ 4]*m[11] - m[ 0]*m[ 7]*m[ 9]) * Inv_det ;

  invm[12] =-( m[ 4]*m[ 9]*m[14] + m[ 6]*m[ 8]*m[13] + m[ 5]*m[10]*m[12] -
		  	   m[ 6]*m[ 9]*m[12] - m[ 5]*m[ 8]*m[14] - m[ 4]*m[10]*m[13]) * Inv_det ;
  invm[13] = ( m[ 0]*m[ 9]*m[14] + m[ 2]*m[ 8]*m[13] + m[ 1]*m[10]*m[12] -
		  	   m[ 2]*m[ 9]*m[12] - m[ 1]*m[ 8]*m[14] - m[ 0]*m[10]*m[13]) * Inv_det ;
  invm[14] =-( m[ 0]*m[ 5]*m[14] + m[ 2]*m[ 4]*m[13] + m[ 1]*m[ 6]*m[12] -
		  	   m[ 2]*m[ 5]*m[12] - m[ 1]*m[ 4]*m[14] - m[ 0]*m[ 6]*m[13]) * Inv_det ;
  invm[15] = ( m[ 0]*m[ 5]*m[10] + m[ 2]*m[ 4]*m[ 9] + m[ 1]*m[ 6]*m[ 8] -
		  	   m[ 2]*m[ 5]*m[ 8] - m[ 1]*m[ 4]*m[10] - m[ 0]*m[ 6]*m[ 9]) * Inv_det ;

}

void invt4x4mh(real *invm, real *m)
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

#ifdef __cplusplus
//Para cualquier tipo de dato
template <typename T>
bool isMH2x2_T(T* m)
{
	if( (m[1] == 0.0) && (m[2] == 0.0) && (m[3] == 1.0) )
		return true;
	else
		return false;
}
template <typename T>
bool isMH3x3_T(T* m)
{
	if( (m[2] == 0.0) && (m[5] == 0.0) && (m[6] == 0.0) && (m[7] == 0.0) && (m[8] == 1.0) )
		return true;
	else
		return false;
}
template <typename T>
bool isMH4x4_T(T* m)
{
	if( (m[ 3] == 0.0) && (m[ 7] == 0.0) && (m[11] == 0.0) && (m[12] == 0.0) && (m[13] == 0.0) && (m[14] == 0.0) && (m[15] == 1.0) )
		return true;
	else
		return false;
}

template<typename T>
T det4x4mh_T(T *m)
{
  T det;

  det = 0;

  det += m[0] * (m[5]*m[10] - m[6]*m[9]);
  det -= m[1] * (m[4]*m[10] - m[6]*m[8]);
  det += m[2] * (m[4]*m[ 9] - m[5]*m[8]);

  return det;
}
template <typename T>
void invt4x4mh_T(T *invm, T *m)
{
	T det, Inv_det;

  det = det4x4mh_T<T>(m);
  Inv_det=1./det;

  invm[ 0] = -(-m[ 5]*m[10] + m[ 6]*m[ 9]) * Inv_det ;
  invm[ 1] = -( m[ 1]*m[10] - m[ 2]*m[ 9]) * Inv_det ;
  invm[ 2] = -(-m[ 1]*m[ 6] + m[ 2]*m[ 5]) * Inv_det ;
  invm[ 3] = 0.0;

  invm[ 4] = -( m[ 4]*m[10] - m[ 6]*m[ 8]) * Inv_det ;
  invm[ 5] = -(-m[ 0]*m[10] + m[ 2]*m[ 8]) * Inv_det ;
  invm[ 6] = -( m[ 0]*m[ 6] - m[ 2]*m[ 4]) * Inv_det ;
  invm[ 7] = 0.0;

  invm[ 8] = -(-m[ 4]*m[ 9] + m[ 5]*m[ 8]) * Inv_det ;
  invm[ 9] = -( m[ 0]*m[ 9] - m[ 1]*m[ 8]) * Inv_det ;
  invm[10] = -(-m[ 0]*m[ 5] + m[ 1]*m[ 4]) * Inv_det ;
  invm[11] = 0.0;

  invm[12] = -(
           m[ 4]*(m[ 9]*m[14] - m[10]*m[13]) -
           m[ 8]*(m[ 5]*m[14] - m[ 6]*m[13]) +
           m[12]*(m[ 5]*m[10] - m[ 6]*m[ 9]) ) * Inv_det ;
  invm[13] = -(-
           m[ 0]*(m[ 9]*m[14] - m[10]*m[13]) +
           m[ 8]*(m[ 1]*m[14] - m[ 2]*m[13]) -
           m[12]*(m[ 1]*m[10] - m[ 2]*m[ 9]) ) * Inv_det ;
  invm[14] = -(
           m[ 0]*(m[ 5]*m[14] - m[ 6]*m[13]) -
           m[ 4]*(m[ 1]*m[14] - m[ 2]*m[13]) +
           m[12]*(m[ 1]*m[ 6] - m[ 2]*m[ 5]) ) * Inv_det ;
  invm[15] = 1.0;

}

template<typename T>
T det4x4m_T(T *m)
{
  T det;

  det = 0;

  det += m[ 0] * ( m[ 5]*m[10]*m[15] + m[ 7]*m[ 9]*m[14] + m[ 6]*m[11]*m[13] -
	  	   	   	   m[ 7]*m[10]*m[13] - m[ 6]*m[ 9]*m[15] - m[ 5]*m[11]*m[14] );
  det -= m[ 4] * ( m[ 1]*m[10]*m[15] + m[ 3]*m[ 9]*m[14] + m[ 2]*m[11]*m[13] -
 	   	   	   	   m[ 3]*m[10]*m[13] - m[ 2]*m[ 9]*m[15] - m[ 1]*m[11]*m[14] );
  det += m[ 8] * ( m[ 1]*m[ 6]*m[15] + m[ 3]*m[ 5]*m[14] + m[ 2]*m[ 7]*m[ 13] -
 	   	   	   	   m[ 3]*m[ 6]*m[13] - m[ 2]*m[ 5]*m[15] - m[ 1]*m[ 7]*m[14] );
  det -= m[12] * ( m[ 1]*m[ 6]*m[11] + m[ 3]*m[ 5]*m[10] + m[ 2]*m[ 7]*m[ 9] -
 	   	   	   	   m[ 3]*m[ 6]*m[ 9] - m[ 2]*m[ 5]*m[11] - m[ 1]*m[ 7]*m[10] );

  /*mexPrintf("det4x4=%lf\n", det);*/

  return det;
}
template<typename T>
void invt4x4m_T(T *invm, T *m)
{
  T det, Inv_det;

  det = det4x4m_T<T>(m);
  Inv_det=1./det;

  invm[ 0] = ( m[ 5]*m[10]*m[15] + m[ 9]*m[ 7]*m[14] + m[ 6]*m[11]*m[13] -
		  	   m[ 7]*m[10]*m[13] - m[ 6]*m[ 9]*m[15] - m[ 5]*m[11]*m[14]) * Inv_det ;
  invm[ 1] =-( m[ 1]*m[10]*m[15] + m[ 3]*m[ 9]*m[14] + m[ 2]*m[11]*m[13] -
	  	   	   m[ 3]*m[10]*m[13] - m[ 2]*m[ 9]*m[15] - m[ 1]*m[11]*m[14]) * Inv_det ;
  invm[ 2] = ( m[ 1]*m[ 6]*m[15] + m[ 3]*m[ 5]*m[14] + m[ 2]*m[ 7]*m[13] -
	  	   	   m[ 3]*m[ 6]*m[13] - m[ 2]*m[ 5]*m[15] - m[ 1]*m[ 7]*m[14]) * Inv_det ;
  invm[ 3] =-( m[ 1]*m[ 6]*m[11] + m[ 3]*m[ 5]*m[10] + m[ 2]*m[ 7]*m[ 9] -
	  	   	   m[ 3]*m[ 6]*m[ 9] - m[ 2]*m[ 5]*m[11] - m[ 1]*m[ 7]*m[10]) * Inv_det ;

  invm[ 4] =-( m[ 4]*m[10]*m[15] + m[ 7]*m[ 8]*m[14] + m[ 6]*m[11]*m[12] -
		  	   m[ 7]*m[10]*m[12] - m[ 6]*m[ 8]*m[15] - m[ 4]*m[11]*m[14]) * Inv_det ;
  invm[ 5] = ( m[ 0]*m[10]*m[15] + m[ 3]*m[ 8]*m[14] + m[ 2]*m[11]*m[12] -
		  	   m[ 3]*m[10]*m[12] - m[ 2]*m[ 8]*m[15] - m[ 0]*m[11]*m[14]) * Inv_det ;
  invm[ 6] =-( m[ 0]*m[ 6]*m[15] + m[ 3]*m[ 4]*m[14] + m[ 2]*m[ 7]*m[12] -
		  	   m[ 3]*m[ 6]*m[12] - m[ 2]*m[ 4]*m[15] - m[ 0]*m[ 7]*m[14]) * Inv_det ;
  invm[ 7] = ( m[ 0]*m[ 6]*m[11] + m[ 3]*m[ 4]*m[10] + m[ 2]*m[ 7]*m[ 8] -
		  	   m[ 3]*m[ 6]*m[ 8] - m[ 2]*m[ 4]*m[11] - m[ 0]*m[ 7]*m[10]) * Inv_det ;

  invm[ 8] = ( m[ 4]*m[ 9]*m[15] + m[ 7]*m[ 8]*m[13] + m[ 5]*m[11]*m[12] -
		  	   m[ 7]*m[ 9]*m[12] - m[ 5]*m[ 8]*m[15] - m[ 4]*m[11]*m[13]) * Inv_det ;
  invm[ 9] =-( m[ 0]*m[ 9]*m[15] + m[ 3]*m[ 8]*m[13] + m[ 1]*m[11]*m[12] -
		  	   m[ 3]*m[ 9]*m[12] - m[ 1]*m[ 8]*m[15] - m[ 0]*m[11]*m[13]) * Inv_det ;
  invm[10] = ( m[ 0]*m[ 5]*m[15] + m[ 3]*m[ 4]*m[13] + m[ 1]*m[ 7]*m[12] -
		  	   m[ 3]*m[ 5]*m[12] - m[ 1]*m[ 4]*m[15] - m[ 0]*m[ 7]*m[13]) * Inv_det ;
  invm[11] =-( m[ 0]*m[ 5]*m[11] + m[ 3]*m[ 4]*m[ 9] + m[ 1]*m[ 7]*m[ 8] -
		  	   m[ 3]*m[ 5]*m[ 8] - m[ 1]*m[ 4]*m[11] - m[ 0]*m[ 7]*m[ 9]) * Inv_det ;

  invm[12] =-( m[ 4]*m[ 9]*m[14] + m[ 6]*m[ 8]*m[13] + m[ 5]*m[10]*m[12] -
		  	   m[ 6]*m[ 9]*m[12] - m[ 5]*m[ 8]*m[14] - m[ 4]*m[10]*m[13]) * Inv_det ;
  invm[13] = ( m[ 0]*m[ 9]*m[14] + m[ 2]*m[ 8]*m[13] + m[ 1]*m[10]*m[12] -
		  	   m[ 2]*m[ 9]*m[12] - m[ 1]*m[ 8]*m[14] - m[ 0]*m[10]*m[13]) * Inv_det ;
  invm[14] =-( m[ 0]*m[ 5]*m[14] + m[ 2]*m[ 4]*m[13] + m[ 1]*m[ 6]*m[12] -
		  	   m[ 2]*m[ 5]*m[12] - m[ 1]*m[ 4]*m[14] - m[ 0]*m[ 6]*m[13]) * Inv_det ;
  invm[15] = ( m[ 0]*m[ 5]*m[10] + m[ 2]*m[ 4]*m[ 9] + m[ 1]*m[ 6]*m[ 8] -
		  	   m[ 2]*m[ 5]*m[ 8] - m[ 1]*m[ 4]*m[10] - m[ 0]*m[ 6]*m[ 9]) * Inv_det ;

}

template <typename T>
void prodm4x4(T *c, T *a, T *b)
{
  int i, j, k;

  for(i=0; i<4; i++)
    for(k=0; k<4; k++)
      for(j=0; j<4; j++)
	*(c+i*4+j) += *(a+i*4+k) * *(b+k*4+j);
}
/* c=a*b (a=mxn; b=nxp; c=mxp) */
template <typename T>
void prodm(T *c, T *a, T *b, int m, int n, int p)
{
  int i, j, k;

  for(i=0; i<m; i++)
    for(k=0; k<n; k++)
      for(j=0; j<p; j++)
	*(c+i*p+j) += *(a+i*n+k) * *(b+k*p+j);
}
/* c=a*b (a=mxn; b=nxp; c=mxp) */
template <typename T>
void prodm_matlab(T *c, T *a, T *b, int m, int n, int p)
{
  int i, j, k;

  for(i=0; i<m; i++)
    for(k=0; k<n; k++)
      for(j=0; j<p; j++)
	*(c+i+j*n) += *(a+i+k*m) * *(b+k+j*n);
}
#endif

#endif
