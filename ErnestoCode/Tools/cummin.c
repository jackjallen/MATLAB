
#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  ALLOCATES(); 
  real    *X, *out;
  int     n,N;
  
  X = myGetPr( prhs[0] ); 
  N = mxGetNumberOfElements( prhs[0] );
  

  
  plhs[0]  = mxCreateNumericArray( mxGetNumberOfDimensions(prhs[0]) , mxGetDimensions(prhs[0]) , mxREAL_CLASS , mxREAL );
  out = (real *) mxGetData( plhs[0] );

  out[0] = X[0];
  
    for( n = 1 ; n<N ; n++ ){
      
       out[n] = out[n-1]<X[n]?out[n-1]:X[n];
    }
    
  EXIT:
    myFreeALLOCATES();
}
