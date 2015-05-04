/*
  o = Solve3x3Sys( M, b )

*/

#include "mex.h"

#define O(t,c)            O[ (t) + (c)*IJK]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double  *IM, v1, v2, v3, *V1, *V2, *V3,  *O, det, *A, *B, *C, *D, *E, *F, *G, *H, *L, a, b, c, d, e, f, g, h, l;
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

  A = IM + IJK*0;  B = IM + IJK*3;  C = IM + IJK*6;
  D = IM + IJK*1;  E = IM + IJK*4;  F = IM + IJK*7;
  G = IM + IJK*2;  H = IM + IJK*5;  L = IM + IJK*8;
  
  V2 = V1 + IJK;  V3 = V2 + IJK;
  
  Odims[0]= I;
  Odims[1]= J;
  Odims[2]= K;
  Odims[3]= 3;
  

  plhs[0] = mxCreateNumericArray( 4 , Odims , mxDOUBLE_CLASS , mxREAL );


  O = (double *)mxGetData( plhs[0] );

  for( t =0 ; t <IJK ; t++  ){
      a = A[t]; b = B[t]; c = C[t]; d = D[t]; e = E[t]; f = F[t];  g = G[t];  h = H[t];  l = L[t]; 
      v1 = V1[t];  v2 = V2[t];  v3 = V3[t];
      
            
      det = -c*e*g + b*f*g + c*d*h - a*f*h - b*d*l + a*e*l;
      
      det = 1/det;
      
      O(t,0) = ( (e*l - f*h)*v1 + (c*h - b*l)*v2 + (b*f - c*e)*v3 )*det;
      O(t,1) = ( (f*g - d*l)*v1 + (a*l - c*g)*v2 + (c*d - a*f)*v3 )*det;
      O(t,2) = ( (d*h - e*g)*v1 + (b*g - a*h)*v2 + (a*e - b*d)*v3 )*det;

  }
  
}
