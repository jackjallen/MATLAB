/*
 [ labels , num_labels ] = toboggan_image( I );
*/

#include "mex.h"
#include "stdlib.h"

#define L(i,j)        L[ (j)*M + (i) ]
#define I(i,j)        I[ (j)*M + (i) ]
#define Z(i,j)        Z[ (j)*M + (i) ]

int       M,N;
double    *I, *L, *Z;
double    lastL;
FILE      *f;

void toboggan(int i, int j){
  double v, vmin;
  int    ni,nj,ii,jj;
  v = I(i,j);
  vmin = v + 1;

  L(i,j) = -1;
  
  for( ii=i-1 ; ii<=i+1 ; ii++ ){
    for( jj=j-1 ; jj<=j+1 ; jj++ ){
      if( !(ii==i && jj==j) && ii>=0 && jj>=0 && ii<M && jj<N && L(ii,jj) >= 0 && I(ii,jj) < vmin ) {
        ni = ii;
        nj = jj;
        vmin = I(ii,jj);
      }
    }
  }
  
  if( vmin <= v ){
//     if( vmin==v ){ Z(ni,nj)=0; }
    if( L(ni,nj) <= 0 ){ 
      toboggan(ni,nj);
    }
    L(i,j) = L(ni,nj);
  } else {
    lastL = lastL + 1;
//     mexPrintf("nuevo L en %d,%d  (%f)\n",i,j,lastL);
    L(i,j) = lastL;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int       i, j, ii, jj;
  double    v;

  lastL = 0;
  
  M =  mxGetM( prhs[0] );
  N =  mxGetN( prhs[0] );
  I = mxGetPr( prhs[0] );
  
  plhs[0] = mxCreateDoubleMatrix( M , N , mxREAL );
  plhs[1] = mxCreateDoubleMatrix( 1 , 1 , mxREAL );
  plhs[2] = mxCreateDoubleMatrix( M , N , mxREAL );
 
  L = mxGetPr( plhs[0] );
  Z = mxGetPr( plhs[2] );
  
  for(j=0;j<N;j++){
    for(i=0;i<M;i++){
      L(i,j)= 0;
      Z(i,j)= 1;
    }
  }

  for(j=0;j<N;j++){
    for(i=0;i<M;i++){
      v = I(i,j);
      for( ii=i-1 ; ii<=i+1 ; ii++ ){
        for( jj=j-1 ; jj<=j+1 ; jj++ ){
          if( !(ii==i && jj==j) && ii>=0 && jj>=0 && ii<M && jj<N && I(ii,jj) != v ) {
            Z(i,j)=0;
          }
        }
      }
    }
  }
  
  for(j=0;j<N;j++){
    for(i=0;i<M;i++){
      if( L(i,j) <= 0 ){
        toboggan(i,j);
      }
    }
  }

  *mxGetPr( plhs[1] ) = lastL;

}
