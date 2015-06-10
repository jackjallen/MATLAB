#include "mex.h"
#include <stdlib.h>
#include "string.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  long          n_xyz, n_tri, n_points;
  long          n, p, t;
  long          nn_xyz, nn_tri;
  long          *new_pids, *new_tids;
  double        *tri;
  const char    *field_name;  
  int           f, n_fields;
  int           c, n_cols;
  mxArray       *DATA;
  const char    *names[] = {""};
  const int     dims[1] = {1};
  mxArray       *INPUT;
  double        *data;

  if(nrhs == 0)
  {
	  mexPrintf("Error calling function.\n");
	  return;
  }

  INPUT= mxDuplicateArray( prhs[0] );
  
  n_xyz   =  mxGetM( mxGetField( INPUT, 0, "xyz") );
  n_tri   =  mxGetM( mxGetField( INPUT, 0, "tri") );
  tri     = mxGetPr( mxGetField( INPUT, 0, "tri") );
  n_points=  mxGetM( prhs[1] )*mxGetN( prhs[1] );

  new_pids= mxMalloc( (n_xyz+1)*sizeof( long ) );
  for(n=1 ; n <= n_xyz ; n++ ){
    new_pids[n]= 1;
  }
  
  for( p=0 ; p< n_points ; p++ ){
    new_pids[ (long)*( mxGetPr( prhs[1] ) + p ) ]= 0;
  }
  nn_xyz=0;
  for(n=1 ; n <= n_xyz ; n++ ){
    if ( new_pids[n] ) {
      nn_xyz++;
      new_pids[n]= nn_xyz;
    }
  }

  new_tids= mxMalloc( (n_tri+1)*sizeof( long ) );
  nn_tri=0;
  for( t=1 ; t <= n_tri ; t++ ){
    if ( new_pids[ (long)tri[t-1] ] &&
         new_pids[ (long)tri[t-1+n_tri] ] &&
         new_pids[ (long)tri[t-1+2*n_tri] ] ) {
      nn_tri++;
      new_tids[t]= nn_tri;
      tri[t-1        ]= new_pids[ (long)tri[t-1        ] ];
      tri[t-1+  n_tri]= new_pids[ (long)tri[t-1+  n_tri] ];
      tri[t-1+2*n_tri]= new_pids[ (long)tri[t-1+2*n_tri] ];
    } else {
      new_tids[t]=0;
    }
  }
  
  plhs[0] = mxCreateStructArray(1, dims, 0, names);
  n_fields= mxGetNumberOfFields( INPUT );
  for( f=0 ; f < n_fields ; f++ ){
    field_name= mxGetFieldNameByNumber( INPUT , f );
    n_cols= mxGetN( mxGetField( INPUT, 0, field_name ) );
    if        ( !strncmp( field_name,"xyz",3 ) || !strcmp( field_name, "uv") ){
      DATA= mxCreateDoubleMatrix( nn_xyz , n_cols , mxREAL );
      data= mxGetPr( DATA );
      for( n=1 ; n <= n_xyz ; n++ ){
        if( new_pids[n] ){
          for( c=0 ; c<n_cols ; c++){
            data[ new_pids[n]-1+c*nn_xyz ]=
              *(mxGetPr( mxGetField( INPUT, 0, field_name ) ) + n-1+c*n_xyz );
          } 
        }
      }
    } else {
      if  ( !strncmp( field_name,"tri",3 ) ) {
        DATA= mxCreateDoubleMatrix( nn_tri , n_cols , mxREAL );
        data= mxGetPr( DATA );
        for( t=1 ; t <= n_tri ; t++ ){
          if( new_tids[t] ){
            for( c=0 ; c<n_cols ; c++){
              data[ new_tids[t]-1+c*nn_tri ]= 
                *(mxGetPr( mxGetField( INPUT, 0, field_name ) ) + t-1+c*n_tri);
            } 
          }
        }
      } else {
         DATA= mxDuplicateArray( mxGetField( INPUT, 0, field_name ));
      }
    }
    mxAddField( plhs[0], field_name );
    mxSetField( plhs[0], 0, field_name, DATA );
  }      
  
  mxDestroyArray( INPUT );
  
  mxFree(new_pids);
  mxFree(new_tids);
}
