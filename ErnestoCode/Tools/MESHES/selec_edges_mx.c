#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int       M, e, i;
  double    *col1, *col2;
  double    node1, node2;
  
  
  plhs[0] = mxDuplicateArray( prhs[0] );

  M = mxGetM( plhs[0] );
  
  col1 = mxGetPr( plhs[0] );
  col2 = col1 + M;
  
  for( e = 0 ; e < M ; e++ ){
    node1 = col1[ e ];  if( !node1 ){ continue; }
    node2 = col2[ e ];  if( !node2 ){ continue; }
    
    for( i = e+1 ; i < M ; i ++ ){
      if(  col1[i] == node1  ||  col1[i] == node2  ){
        col1[i] = 0;
      }

      if(  col2[i] == node1  ||  col2[i] == node2  ){
        col2[i] = 0;
      }
    }
  }

}

