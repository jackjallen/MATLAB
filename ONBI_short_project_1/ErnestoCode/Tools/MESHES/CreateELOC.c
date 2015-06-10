#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include "SCAMCOORD.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  long      n_xyz;
  long      n_tri;
  double    *tri;
  long      id;
  int       i, j, k;
  int       mini, minj, mink;
  int       maxi, maxj, maxk;
  long      t;
  mxArray   *ELOC;
  mxArray   **lists;
  long      *counter;
  int       N[3];
  int       *ijk;
  long      *sizes;

  n_xyz =  mxGetM( mxGetField( prhs[0], 0, "xyz") );
  n_tri =  mxGetM( mxGetField( prhs[0], 0, "tri") );
  tri   = mxGetPr( mxGetField( prhs[0], 0, "tri") );

  /* //SCAM coordinates of the nodes */
  ijk   = mxMalloc( 3*n_xyz*sizeof( int ) );
  SCAMCOORD( ijk, prhs[0] , mxGetPr( mxGetField( prhs[0], 0, "xyz") ), n_xyz , N);

/*   //Create the output */
  plhs[0]= mxDuplicateArray( prhs[0] );
  if( mxGetField( plhs[0], 0, "ELOC") == NULL ){
    mxAddField( plhs[0],"ELOC" );
  }
  ELOC= mxCreateCellArray( 3 , N );   /* //Create the N[0]xN[1]xN[2] cell array in ELOC */

/*   //lists point to each list of bucket */
  lists= (mxArray **)mxCalloc( N[0]*N[1]*N[2] , sizeof( mxArray * ));
/*   //number of triangles in each bucket */
  sizes = mxMalloc( N[0]*N[1]*N[2]*sizeof( long ) );
/*  //an auxiliar counter */
  counter   = mxMalloc( N[0]*N[1]*N[2]*sizeof( long ) );
/*   //setting to zero sizes and counter */
  for( i= 0 ; i < N[0]*N[1]*N[2] ; i++ ){
    sizes[i]= 0;
    counter[i]  = 0;
  }

/*   //run over the triangles to count the sizes of the lists */
  for( t=0 ; t<n_tri; t++ ){
    mini= MMIN( ijk[(int)tri[t]-1] , ijk[(int)tri[t+n_tri]-1] , ijk[(int)tri[t+2*n_tri]-1] );
    maxi= MMAX( ijk[(int)tri[t]-1] , ijk[(int)tri[t+n_tri]-1] , ijk[(int)tri[t+2*n_tri]-1] );
    minj= MMIN( ijk[(int)tri[t]-1+n_xyz] , ijk[(int)tri[t+n_tri]-1+n_xyz] , ijk[(int)tri[t+2*n_tri]-1+n_xyz] );
    maxj= MMAX( ijk[(int)tri[t]-1+n_xyz] , ijk[(int)tri[t+n_tri]-1+n_xyz] , ijk[(int)tri[t+2*n_tri]-1+n_xyz] );
    mink= MMIN( ijk[(int)tri[t]-1+2*n_xyz] , ijk[(int)tri[t+n_tri]-1+2*n_xyz] , ijk[(int)tri[t+2*n_tri]-1+2*n_xyz] );
    maxk= MMAX( ijk[(int)tri[t]-1+2*n_xyz] , ijk[(int)tri[t+n_tri]-1+2*n_xyz] , ijk[(int)tri[t+2*n_tri]-1+2*n_xyz] );
    
    for( i=mini-1; i<maxi; i++ ) {
      for( j=minj-1; j<maxj; j++ ) {
        for( k=mink-1; k<maxk; k++ ) {
          id= k*N[0]*N[1] + j*N[0] + i;
          sizes[ id ]++;
        }
      }
    }
  }
  
/*   //create the output lists in the ELOC cell array */
  for( i= 0 ; i < N[0]*N[1]*N[2] ; i++ ){
    lists[i] = mxCreateDoubleMatrix( 1, sizes[i] , mxREAL );
    mxSetCell( ELOC , i , lists[i]);
  }

/*   //filling the list */
  for( t=0 ; t<n_tri; t++ ){
    mini= MMIN( ijk[(int)tri[t]-1] , ijk[(int)tri[t+n_tri]-1] , ijk[(int)tri[t+2*n_tri]-1] );
    maxi= MMAX( ijk[(int)tri[t]-1] , ijk[(int)tri[t+n_tri]-1] , ijk[(int)tri[t+2*n_tri]-1] );
    minj= MMIN( ijk[(int)tri[t]-1+n_xyz] , ijk[(int)tri[t+n_tri]-1+n_xyz] , ijk[(int)tri[t+2*n_tri]-1+n_xyz] );
    maxj= MMAX( ijk[(int)tri[t]-1+n_xyz] , ijk[(int)tri[t+n_tri]-1+n_xyz] , ijk[(int)tri[t+2*n_tri]-1+n_xyz] );
    mink= MMIN( ijk[(int)tri[t]-1+2*n_xyz] , ijk[(int)tri[t+n_tri]-1+2*n_xyz] , ijk[(int)tri[t+2*n_tri]-1+2*n_xyz] );
    maxk= MMAX( ijk[(int)tri[t]-1+2*n_xyz] , ijk[(int)tri[t+n_tri]-1+2*n_xyz] , ijk[(int)tri[t+2*n_tri]-1+2*n_xyz] );
    
    for( i=mini-1; i<maxi; i++ ) {
      for( j=minj-1; j<maxj; j++ ) {
        for( k=mink-1; k<maxk; k++ ) {
          id= k*N[0]*N[1] + j*N[0] + i;
          *( mxGetPr( lists[id] ) + (counter[id]++) )= t+1;
        }
      }
    }
  }

/*   //fixing the output */
  mxSetField( plhs[0], 0, "ELOC", ELOC );

  mxFree(ijk);
/* //   for( i=0 ; i < N[0]*N[1]*N[2] ; i++ ){ 
//     mxDestroyArray(lists[i]);
//   } */
  mxFree(lists);
  mxFree(counter);
  mxFree(sizes);
}
