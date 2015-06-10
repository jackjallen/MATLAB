/*

  [id,distance] = ClosestPoint( X , N , maxd )

 
  X = randn( 100000  ,2);
  N = randn(   5000  ,2);
  tic
  id = ClosestPoint( X , N );
  toc
 
  figure
  line(X(:,1),X(:,2),'marker','.','linestyle','none');
  line(N(:,1),N(:,2),'marker','+','linestyle','none','color',[1 0 0]);
  for i = 1:numel(id)
    line( [ N(i,1) X(id(i),1) ],[ N(i,2) X(id(i),2) ],'color',[1 0 1] );
  end
  axis equal


*/

#define   real               double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  ALLOCATES();
  
  int     nX, nN;
  
  int     nsd, c, n, x;
  double  *X, *N, *O, *D, *P;
  double  minO, minD, dist, xx;
  double  maximalDIST;
  
  P = NULL;
  
  nsd = mySize( prhs[0] , 1 );
  if( mySize( prhs[1] , 1 ) != nsd ){ myErrMsgTxt("number of spatial coordinates (number of column) have to be the same."); }
  
  
  nX = mySize( prhs[0] , 0 );
  X = myGetPr( prhs[0] );
  
  nN = mySize( prhs[1] , 0 );
  N = myGetPr( prhs[1] );
  
  
  if( nrhs > 2 ){
    maximalDIST = *( myGetPr( prhs[2] ) );
    maximalDIST = maximalDIST*maximalDIST;
//     mexPrintf("minimal %g\n\n",maximalDIST);
  } else {
    maximalDIST = INF;
  }
  
  plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
  mxSetM( plhs[0] , nN ); mxSetN( plhs[0] , 1 );  mxSetData( plhs[0] , mxMalloc( nN << 3 ) );
  O = mxGetPr( plhs[0] );
  
  if( nlhs > 1 ){
    plhs[1] = mxCreateDoubleMatrix(0,0,mxREAL);
    mxSetM( plhs[1] , nN ); mxSetN( plhs[1] , 1 );  mxSetData( plhs[1] , mxMalloc( nN << 3 ) );
    D = mxGetPr( plhs[1] );
  }
  
  P = mxMalloc( nsd << 3 );
  
  for( n = 0 ; n < nN ; n++ ){

    for( c = 0 ; c < nsd ; c++ ){ P[c] = N[ c*nN + n ]; }

    dist = 0;
    for( c = 0 ; ( c < nsd ) && ( dist < maximalDIST ); c++ ){
      xx = X[ c*nX + 0 ] - P[c];
      xx = xx*xx;
      dist += xx;
    }
    
    minO = 0;
    minD = dist;
//     for( x = 1 ; x < nX  && ( minD > maximalDIST ); x++ ){
    for( x = 1 ; x < nX  ; x++ ){
      dist = 0;
      for( c = 0 ; c < nsd   &&  ( dist < minD ) && ( dist < maximalDIST ); c++ ){
        xx = X[ c*nX + x ] - P[c];
//         if( fabs(xx) > minD ){ dist = minD + 1e10; break; }
        xx = xx*xx;
        dist += xx;
      }
      
      if( dist < minD ){
        minD = dist;
        minO = x;
      }
    }
    
    O[n] = minO + 1;
    if( nlhs > 1 ){
      D[n] = sqrt( minD );
    }
  }
  

  EXIT:
    if( P != NULL ){ mxFree(P); }
    myFreeALLOCATES();
}
