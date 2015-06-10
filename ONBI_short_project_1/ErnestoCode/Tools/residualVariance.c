/*

  residualVariance( X , M , C ) == var( X - M*C , 1 )

*/

#include "myMEX.h"

#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  ALLOCATES();
  real    *M;
  real    *C;
  real    *X;
  
  int     i,j,c, nC, N, nV;
  int     Odims[3];
  
  real    *VAR;
  real    r, sumx2, sumx;
  
  if( mySize(prhs[0],0) != mySize(prhs[1],0) ){ myErrMsgTxt("size(X,1) have to be equal to size(M,1)."); }
  if( mySize(prhs[1],1) != mySize(prhs[2],0) ){ myErrMsgTxt("size(M,2) have to be equal to size(C,1)."); }

  X = myGetPr( prhs[0] );
  M = myGetPr( prhs[1] );
  C = myGetPr( prhs[2] );

  nV = mySize( prhs[0] , 1 );
  N  = mySize( prhs[0] , 0 );
  nC = mySize( prhs[2] , 0 );

  Odims[0]= 1;
  Odims[1]= nV;
  plhs[0] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
  VAR = (real *) mxGetData( plhs[0] );

  for( j = 0 ; j < nV ; j++ ){
    
    sumx2 = 0;
    sumx  = 0;
    for( i = 0 ; i < N ; i++ ){
    
      r = X[i];
      for( c = 0 ; c < nC ; c++ ){
        r -= M[ i + c*N ]*C[c];
      }
      
      sumx  += r;
      sumx2 += r*r;
      
    }
    X = X + N;
    C = C + nC;

    VAR[j] = ( sumx2 - sumx*sumx/N )/( N-1 );
    
  }

  
  EXIT: myFreeALLOCATES();
}

