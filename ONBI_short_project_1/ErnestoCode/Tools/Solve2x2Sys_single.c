/*
  o = Solve2x2Sys( M, b )

*/

#include "mex.h"

#define O(t,c)            O[ (t) + (c)*IJ]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  float  *IM, v1, v2, *V1, *V2, *O, det, *A, *B, *C, *D, a, b, c, d;
  int     I,J;
  int     Odims[3];
  long    IJ, t;

  IM   = (float *)mxGetData(prhs[0]);
  V1   = (float *)mxGetData(prhs[1]);
  I    = (int) *( mxGetDimensions(prhs[0]) + 0 );
  J    = (int) *( mxGetDimensions(prhs[0]) + 1 );

  IJ= I*J;

  A = IM + IJ*0;  B = IM + IJ*2;  
  C = IM + IJ*1;  D = IM + IJ*3;
  
  V2 = V1 + IJ;
  
  Odims[0]= I;
  Odims[1]= J;
  Odims[2]= 2;

  plhs[0] = mxCreateNumericArray( 3 , Odims , mxSINGLE_CLASS , mxREAL );


  O = (float *)mxGetData( plhs[0] );

  for( t =0 ; t <IJ ; t++  ){
      a = A[t]; b = B[t]; c = C[t]; d = D[t]; 
      v1 = V1[t];  v2 = V2[t];
      
            
      det = a*d-b*c;
      
      det = 1/det;
      
      O(t,0) = (  d*v1 - b*v2  )*det;
      O(t,1) = ( -c*v1 + a*v2  )*det;

  }
  
}
