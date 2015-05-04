#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real    *G, *DG;
  int     I;
  int     Odims[50], ndims, d, vector;


  ndims = myGetSizes( prhs[0] , Odims );
  
  vector = 0;
  for( d = 0 ; d < ndims ; d++ ){
    if( Odims[d] > 1 ){
      if( vector ){
        myErrMsgTxt("G is not a singleton.");
      } else {
        vector = 1;
        Odims[d]++;
      }
    }
  }
  if( !vector ){
    Odims[0]++;
  }
    
  G = myGetPr( prhs[0] ); 
  I = myNumel( prhs[0] );

  if( ! checkIsSorted(G,I) ){
    myErrMsgTxt("G is not sorted.");
  }

  plhs[0] = mxCreateNumericArray( 0 , 0 , mxDOUBLE_CLASS , mxREAL );
  mxSetDimensions( plhs[0] , Odims , ndims );


  if( I == 1 ){
    mxSetData( plhs[0] , mxMalloc( 16 ) );
    
    DG = mxGetPr( plhs[0] );
    DG[0] = G[0] - 0.5;
    DG[1] = G[0] + 0.5;

  } else {

    DG = DualGrid(G,I);
    mxSetPr( plhs[0] , DG );

  }
  
  EXIT: myFreeALLOCATES();
}
