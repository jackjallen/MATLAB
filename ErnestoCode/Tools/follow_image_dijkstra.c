#include "mex.h"
#include "stdlib.h"
#include "math.h"

// #define RRAND(a,b)      ( rand()*1.0/RAND_MAX*((b)-(a)) + (a) )
#define mxFlush()           mexEvalString("drawnow expose;");

#define XY(i,j)       XY[ ((j)-1)*MN + ((i)-1) ]
#define  C(i,j)        C[ ((j)-1)*mC + ((i)-1) ]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  long      i, nj, j, l, MN, mC;
  double    *F, *XY, *C;
  double    _NAN_;


  XY = mxGetPr( prhs[1] );

  MN =  mxGetNumberOfElements( prhs[0] );
  
  plhs[0] = mxCreateDoubleMatrix( 2 , 5*MN  , mxREAL );
  C  = mxGetPr( plhs[0] );
  mC = mxGetM( plhs[0] );

  if( nrhs == 2 ) {
    _NAN_ = log(-1);
    F  = mxGetPr( mxDuplicateArray( prhs[0] ));
    l = 1;
    for( i = 1 ; i <= MN ; i++ ){
      if( !F[i-1] ){ continue; }
      j = i;
      C(1,l) = XY(j,2); C(2,l)= XY(j,1); l++;
      while( F[j-1] ){
        nj = F[j-1]; F[j-1] = 0;
        j = nj;
        C(1,l) = XY(j,2); C(2,l)= XY(j,1); l++;
      }
      C(1,l) = _NAN_; C(2,l)= _NAN_; l++;
    }
  } else {
    F  = mxGetPr( prhs[0] );
    j = (long) *mxGetPr( prhs[2] );
    l = 1;
    
    C(1,l) = XY(j,2); C(2,l)= XY(j,1); l++;
    while( F[j-1] ){
      j = F[j-1]; 
      C(1,l) = XY(j,2); C(2,l)= XY(j,1); l++;
    }
    l++;
  }
  mxSetN( plhs[0] , l-2 );

}

