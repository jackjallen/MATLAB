/*
  o = Solve3x3SymSys( M, b )

*/

#include "mex.h"

#define O(t,c)            O[ (t) + (c)*IJK]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double  *IM, v1, v2, v3, *V1, *V2, *V3,  *O, det, *A, *B, *C, *E, *F, *L, a, b, c, e, f, l;
  int     I,J,K;
  int     Odims[4];
  long    IJ,IJK,t;

  IM   = (double *)mxGetData(prhs[0]);
  V1   = (double *)mxGetData(prhs[1]);
  I    = (int) *( mxGetDimensions(prhs[0]) + 0 );
  J    = (int) *( mxGetDimensions(prhs[0]) + 1 );
  K    = (int) *( mxGetDimensions(prhs[0]) + 2 );
/*  ndims = mxGetNumberOfDimensions(prhs[0]);
  if( ndims > 2){
    K= (int) *( mxGetDimensions(prhs[0]) + 2 );
  } else {
    K=1;
  }*/
  IJK= I*J*K;

  A = IM + IJK*0;  
  B = IM + IJK*1;  
  C = IM + IJK*2;
  E = IM + IJK*3;
  F = IM + IJK*4;
  L = IM + IJK*5;
  
  V2 = V1 + IJK;  V3 = V2 + IJK;
  
  Odims[0]= I;
  Odims[1]= J;
  Odims[2]= K;
  Odims[3]= 3;
  

  plhs[0] = mxCreateNumericArray( 4 , Odims , mxDOUBLE_CLASS , mxREAL );


  O = (double *)mxGetData( plhs[0] );

  for( t =0 ; t <IJK ; t++  ){
      a = A[t]; b = B[t]; c = C[t]; e = E[t]; f = F[t];  l = L[t]; 
      v1 = V1[t];  v2 = V2[t];  v3 = V3[t];
      
            
      det = -c*e*c + b*f*c + c*b*f - a*f*f - b*b*l + a*e*l;
      
      det = 1/det;
      
      O(t,0) = ( (e*l - f*f)*v1 + (c*f - b*l)*v2 + (b*f - c*e)*v3 )*det;
      O(t,1) = ( (f*c - b*l)*v1 + (a*l - c*c)*v2 + (c*b - a*f)*v3 )*det;
      O(t,2) = ( (b*f - e*c)*v1 + (b*c - a*f)*v2 + (a*e - b*b)*v3 )*det;

  }
  
}
