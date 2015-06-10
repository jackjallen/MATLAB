#include "mex.h"
#include <stdlib.h>
#include <string.h>
#include "SCAMCOORD.h"
/*
//nlhs number of output args
//plhs pointer to the output args
//nrhs number of input args
//prhs pointer to the input args
//M is the number of rows
//N is the number of columns
*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
  long      n_xyz;
  long      id;
  int       i, j, k;
  long      n;
  mxArray   *PLOC;
  mxArray   **lists;
  long      *counter;
  int       N[3];
  int       *ijk;
  long      *sizes;

  n_xyz = mxGetM( mxGetField( prhs[0], 0, "xyz") );

  /* //SCAM coordinates of the nodes */
  ijk   = mxMalloc( 3*n_xyz*sizeof( int ) );
  SCAMCOORD( ijk, prhs[0] , mxGetPr( mxGetField( prhs[0], 0, "xyz") ), n_xyz , N);

/*   //Create the output */
  plhs[0]= mxDuplicateArray( prhs[0] );
  if( mxGetField( plhs[0], 0, "PLOC") == NULL ){
    mxAddField( plhs[0],"PLOC" );
  }
  PLOC= mxCreateCellArray( 3 , N );   /* //Create the N[0]xN[1]xN[2] cell array in PLOC */

/*   //lists point to each list of nodes */
  lists= (mxArray **)mxCalloc( N[0]*N[1]*N[2] , sizeof( mxArray * ));
/*   //number of nodes in each bucket */
  sizes = mxMalloc( N[0]*N[1]*N[2]*sizeof( long ) );
/*  //an auxiliar counter */
  counter   = mxMalloc( N[0]*N[1]*N[2]*sizeof( long ) );
/*   //setting to zero sizes and counter */
  for( i= 0 ; i < N[0]*N[1]*N[2] ; i++ ){
    sizes[i]= 0;
    counter[i]  = 0;
  }

/*   //run over the nodes to count the sizes of the lists */
  for( n=0 ; n<n_xyz; n++ ){
    i = ijk[ n           ]-1;
    j = ijk[ n + n_xyz   ]-1;
    k = ijk[ n + 2*n_xyz ]-1;
    id= k*N[0]*N[1] + j*N[0] + i;  /*     //like sub2ind index */
    sizes[ id ]++;
  }

/*   //create the output lists in the PLOC cell array */
  for( i= 0 ; i < N[0]*N[1]*N[2] ; i++ ){
    lists[i] = mxCreateDoubleMatrix( 1, sizes[i] , mxREAL );
    mxSetCell( PLOC , i , lists[i]);
  }

/*   //filling the list */
  for( n=0 ; n<n_xyz; n++ ){
    i = ijk[ n           ]-1;
    j = ijk[ n + n_xyz   ]-1;
    k = ijk[ n + 2*n_xyz ]-1;
    id= k*N[0]*N[1] + j*N[0] + i;
    *( mxGetPr( lists[id] ) + (counter[id]++) )= n+1;
  }

 /*   //fixing the output */
  mxSetField( plhs[0], 0, "PLOC", PLOC );

  mxFree(ijk);
/* //   for( i=0 ; i < N[0]*N[1]*N[2] ; i++ ){ 
//     mxDestroyArray(lists[i]);
//   } */
  mxFree(lists);
  mxFree(counter);
  mxFree(sizes);
}
