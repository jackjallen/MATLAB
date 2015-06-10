#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include "math.h"

/* //nlhs number of output args
//plhs pointer to the output args
//nrhs number of input args
//prhs pointer to the input args
//M is the number of rows
//N is the number of columns */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  double *points;
  double *o;
  double *W;
  double *newpoints;
  double k;
  
  char basis[50];
  int nsd;
  int d, dd;
  
  long n_points;
  long n_o, n_w;
  long i;
  long p;
  
  
  mxGetString( prhs[3], basis, ( mxGetM(prhs[3]) * mxGetN(prhs[3]) ) + 1 );
  
  nsd     = mxGetN(  prhs[0] );
  n_points= mxGetM(  prhs[0] );
  n_o     = mxGetM(  prhs[1] );
  n_w     = mxGetM(  prhs[2] );
  points  = mxGetPr( prhs[0] );
  o       = mxGetPr( prhs[1] );
  W       = mxGetPr( prhs[2] );
  
  plhs[0]= mxCreateDoubleMatrix( n_points , nsd , mxREAL );
  newpoints= mxGetPr( plhs[0] );
  
  for( p=0; p<n_points; p++ ){
    for( i=0; i<n_o; i++ ){
      k= 0;
      for( d=0; d<nsd; d++ ){
        k += (points[p + d*n_points] - o[i + d*n_o ])*
             (points[p + d*n_points] - o[i + d*n_o ]);
      }
      if( !strcmp(basis,"r") ){
        k=sqrt(k);
      } else {
        if( k < 1e-15 ) { k= 0;        } 
                   else { k= k*log(k); }
      }
      for( d=0; d<nsd ; d++){
        newpoints[ p+d*n_points ] += W[ i+d*n_w ]*k;
      }
    }
    for( d=0 ; d<nsd ; d++ ){
      for( dd=0 ; dd<nsd ; dd++ ){
        newpoints[ p + d*n_points ] += points[ p + dd*n_points ]*W[ n_o + dd + d*n_w ];
      }
      newpoints[ p + d*n_points ] += W[ n_o + dd + d*n_w ];
    }
  }
 
}


