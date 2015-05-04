#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  
  real *cumprod, *size;
  
  real *O, *idxs;
  real IND, r;
  
  int N, n, S, s;
  
  size = myGetPr( prhs[0] );
  S    = myNumel( prhs[0] );

  cumprod = mxMalloc( S * sizeof( real ) );
  cumprod[0] = 1;
  for( s = 0 ; s < S-1 ; s++ ){
    cumprod[s+1] = cumprod[ s ] * size[ s ];
  }
  
  
  N    = mxGetNumberOfElements( prhs[1] );
  idxs = myGetPr( prhs[1] );
  
  plhs[0] = myCreateDoubleMatrix_E( N , S , mxREAL );
  O = mxGetPr( plhs[0] );
  

  for( n = 0 ; n < N ; n++ ){
    IND = idxs[n] - 1;
    
    for( s = S-1 ; s > 0 ; s-- ){
      r = (int) ( IND / cumprod[s] );
      O[ s*N + n ] = ( (double) r ) + 1;
      IND = IND - cumprod[s]*r;
    }
    O[ n ] = IND + 1;
  }
  
  mxFree( cumprod );
  EXIT: myFreeALLOCATES();
}
