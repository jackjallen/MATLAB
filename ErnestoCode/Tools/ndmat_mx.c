#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real    *O, *INPUTS[50];
  int     n, j, c, i, e, r, R, S[50], E[50];

  E[0] = 1;
  for( c = 0 ; c < nrhs ; c++ ){
    INPUTS[c] = myGetPr( prhs[c] );
    S[c] = myNumel( prhs[c] );
    E[c+1] = E[c] * S[c];
  }
  R = E[c];
  
  
  plhs[0] = myCreateDoubleMatrix_E( R , nrhs , mxREAL );
  O = (real *) mxGetData( plhs[0] );


  for( c = 0 ; c < nrhs ; c++ ){
    for( i = 0 ; i < S[c] ; i++ ){

      n = E[c];
      e = E[c+1] - n;
      
      for( r = n*i ; r < R ; r += e ){
        for( j = 0 ; j < n ; j++ ){
          O[r++] = INPUTS[c][i];
        }
      }
    }
    
    O += R;
  }
    
  
  EXIT: myFreeALLOCATES();
}
