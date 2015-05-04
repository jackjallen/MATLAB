#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include "math.h"
#include "SCAMCOORD.h"
#define myFlush()                   mexEvalString("drawnow expose;")

#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  int           n_fields, f;
  const char    *field_name;
  int           n_dims= 0, c;
  long          n_tri, t, n_xyz, n, p;
  double        *xyz_e;
  int           *ijk;
  long          *new_ids;
  double        threshold;
  long          last_id;
  double        *DATA;
  double        *newDATA;
  int           n_cols;
  double        distance;
  double        d;
  long          nn_xyz;
  mxArray       *mxNewDATA;
  int           N[3];

  threshold= *( mxGetPr( prhs[1] ) );
  threshold= threshold*threshold;

  n_xyz =  mxGetM( mxGetField( prhs[0], 0, "xyz") );
/* //   mexPrintf("%d\n",n_xyz); myFlush(); */

  ijk= mxMalloc( 3 * n_xyz * sizeof(int) );


  SCAMCOORD( ijk , prhs[0], mxGetPr( mxGetField( prhs[0], 0, "xyz") ) , n_xyz, N);
/* //   for( n=0; n<n_xyz; n++){
//     mexPrintf( "%d  %d  %d\n",
//         ijk[n],ijk[n+n_xyz],ijk[n+2*n_xyz]);
//   } */
  

  xyz_e = mxMalloc( n_xyz*sizeof( double ) );
  n_fields= mxGetNumberOfFields( prhs[0] );
  n_dims= 0;
  for( f=0 ; f < n_fields ; f++ ){
/* //     mexPrintf("%d\n",f); myFlush(); */
    field_name= mxGetFieldNameByNumber( prhs[0] , f );
    if( !strncmp( field_name,"xyz",3 ) || !strcmp( field_name, "uv") ){
      DATA  = mxGetPr( mxGetField( prhs[0], 0, field_name ) );
      n_cols=  mxGetN( mxGetField( prhs[0], 0, field_name ) );
      xyz_e= mxRealloc( xyz_e, (n_dims+n_cols)*n_xyz*sizeof( double ) );
      for( c=0 ; c < n_cols ; c++ ){
        for( n= 0; n< n_xyz ; n++ ){
          xyz_e[ n_dims*n_xyz + n ] = DATA[ c*n_xyz + n ];
        }
        n_dims++;
      }
    }
  }
/* //   mexPrintf("\n\nPorAca1\n\n"); */

  new_ids= mxMalloc( n_xyz*sizeof(long) );
  for( n=0 ; n < n_xyz ; n++ ){
    new_ids[n]= -1;
  }
/* //   mexPrintf("\n\nPorAca2\n\n"); */

  last_id= 0;
  for( n=0; n < n_xyz ; n++ ){
 /* //     mexPrintf( "n: %d of %d\n" , n , n_xyz ); myFlush(); */
    
    if( new_ids[n] != -1 ){
      continue;
    }
    new_ids[n] = last_id;
    for( p=n+1 ; p<n_xyz ; p++ ){
/* //       mexPrintf("n: %d  p: %d   de  %d", n,p,n_xyz);  myFlush(); */
 
/* // mexPrintf("%d %d  %f\n", ijk[p+2*n_xyz],ijk[n+2*n_xyz],fabs(1)); */
      if( new_ids[p] != -1 )                          { continue; }
      if( fabs(ijk[p        ] - ijk[n        ]) > 1 ) { continue; }
      if( fabs(ijk[p+  n_xyz] - ijk[n+  n_xyz]) > 1 ) { continue; }
      if( fabs(ijk[p+2*n_xyz] - ijk[n+2*n_xyz]) > 1 ) { continue; }
  
/* //       mexPrintf(" distance0 - " );  myFlush(); */
      distance=0;
      for( c=0 ; c<n_dims & distance<=threshold ; c++ ){
          d = ( xyz_e[ c*n_xyz + p] - xyz_e[ c*n_xyz + n] );
          distance += d*d;
      }
/* //       mexPrintf(" distance: %f ", distance);  myFlush(); */
      
      
/* //       mexPrintf(" new_ids: %d ", new_ids[ p ]);  myFlush(); */
      if( distance <= threshold ){
        new_ids[ p ]= last_id;
/* //         mexPrintf("last_id : %d\n",last_id); myFlush(); */
      }
/* //       mexPrintf("\n"); myFlush(); */
    }
    last_id++;
  }
  nn_xyz= last_id;
  
  plhs[0]= mxDuplicateArray( prhs[0] );
  for( f=0 ; f<n_fields ; f++ ){
    field_name= mxGetFieldNameByNumber( plhs[0] , f );
/* // mexPrintf("%s  ",field_name); */

    if( !strncmp( field_name,"xyz",3 ) || !strcmp( field_name, "uv") ){
      DATA    = mxGetPr( mxGetField( plhs[0], 0, field_name ) );
      n_cols  = mxGetN(  mxGetField( plhs[0], 0, field_name ) );    
      mxNewDATA= mxCreateDoubleMatrix( nn_xyz , n_cols , mxREAL );
      newDATA  = mxGetPr( mxNewDATA );
      last_id= 0;
      for( n=0 ; n<n_xyz ; n++ ){
        if( new_ids[n] != last_id ) { continue; }
        for( c=0 ; c<n_cols ; c++ ){
          newDATA[ c*nn_xyz + last_id ]= DATA[ c*n_xyz + n ];
        }
        last_id++;
      }
      mxSetField( plhs[0], 0, field_name, mxNewDATA  );
    }

    if( !strcmp( field_name,"tri" ) ){
      DATA   = mxGetPr( mxGetField( plhs[0], 0, field_name ) );
      n_tri  =  mxGetM( mxGetField( plhs[0], 0, field_name ) );    
      for( t=0 ; t<n_tri ; t++ ){
        DATA[t        ]=new_ids[(long)(DATA[t        ]-1)]+1;
        DATA[t+  n_tri]=new_ids[(long)(DATA[t+  n_tri]-1)]+1;
        DATA[t+2*n_tri]=new_ids[(long)(DATA[t+2*n_tri]-1)]+1;
      }
    }
    if( !strcmp( field_name,"PLOC" ) ){
      mxRemoveField( plhs[0],f);
      n_fields--; f--;
    }
    if( !strcmp( field_name,"PSUP" ) ){
      mxRemoveField( plhs[0],f);
      n_fields--; f--;
    }
    if( !strcmp( field_name,"ESUP" ) ){
      mxRemoveField( plhs[0],f);
      n_fields--; f--;
    }
/* // mexPrintf("ok\n");   */
  }
  
  mxFree(ijk);
  mxFree(new_ids);
  mxFree(xyz_e);
}
