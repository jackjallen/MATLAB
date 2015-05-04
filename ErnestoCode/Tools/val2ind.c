/*

 G = sort( (1:20) + ( rand(1,20)-0.5 )*5 );
 X = randn(10000,1)*20;
 figure; plot( X , val2ind( G , X ) , '.' ); set(gca,'xtick',sort(G),'xgrid','on')
 figure; plot( X , getInterval( X , G ) , '.' ); set(gca,'xtick',sort(G),'xgrid','on')
 
 

 if G is sorted 
        1            2         k               n-1              n
 ...............)[........)[.........       ..........)[................... 
            *        *           *              *               *
           G1       G2          Gk             Gn-1            Gn
 
 
 else
 
 
        k1           k2         kk              kn-1            kn
 ...............)[........)[.........       ..........)[................... 
            *        *           *              *               *
           Gk1      Gk2         Gkk            Gkn-1           Gkn
        

  
 G = sort( (1:2000000) + ( rand(1,2000000)-0.5 )*5 );
 X = sort( randn(1000,1)*20 );

 Gns = G(randperm(numel(G)));
 tic;
 idxs1 = val2ind( Gns , X );
 toc;
 
 tic;
 idxs1 = val2ind( G , X );
 toc;

 tic;
 idxs2 = val2ind( G , X , 'sorted' );
 toc;
 maxnorm( idxs1 - idxs2 )

 perm = randperm( numel(X) );
 Xns = X(perm);
 tic;
 idxs3 = val2ind( G , Xns , 'sorted' );
 toc;
 maxnorm( idxs3 - idxs2(perm) )

 *
 *
 *
 
 with 'inside' option
   -Inf       1      2         k               n-1         n        +Inf
 ..........)[...)[........)[.........       ..........)[........](......... 
            *        *           *              *               *
           G1       G2          Gk             Gn-1            Gn

  
*/

#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real    *G, *X, *out, x, distance;
  int     i,I,n,N,closest;
  int     issorted = 0, InfOutside = 0;
  char    STR[100];
  
  G = myGetPr( prhs[0] ); 
  I = myNumel( prhs[0] );
  
  if( I == 0 ){
    plhs[0]  = mxCreateDoubleMatrix( 0, 0, mxREAL );
    goto EXIT;
  }
  
  X = myGetPr( prhs[1] ); 
  N = myNumel( prhs[1] );
  
  
  if( nrhs > 2 ){
    if( !mxIsChar(prhs[2]) ){ myErrMsgTxt("String expected\n"); }
    mxGetString( prhs[2], STR, 100 );
    if( ! myStrcmpi( STR , "sorted" ) ){ issorted = 1;   }
    if( ! myStrcmpi( STR , "inside" ) ){ InfOutside = 1; }
  }
  if( nrhs > 3 ){
    if( !mxIsChar(prhs[3]) ){ myErrMsgTxt("String expected\n"); }
    mxGetString( prhs[3], STR, 100 );
    if( ! myStrcmpi( STR , "sorted" ) ){ issorted = 1;   }
    if( ! myStrcmpi( STR , "inside" ) ){ InfOutside = 1; }
  }

  
  plhs[0]  = mxCreateNumericArray( mxGetNumberOfDimensions(prhs[1]) , mxGetDimensions(prhs[1]) , mxREAL_CLASS , mxREAL );
  out = (real *) mxGetData( plhs[0] );

  if( !issorted ){
    issorted = checkIsSorted(G,I);
  }
  
  if( issorted ) {
    
    i = -1;
    for( n = 0 ; n<N ; n++ ){
      x = X[n];

      if( InfOutside ){
        if( x < G[ 0 ] ){ out[n] = INFn; continue; }
        if( x > G[I-1] ){ out[n] = INFp; continue; }
      }
      
/*****************************************************************************

I = numel(G)

... -2   ) [.  0   ) [.  1   ) [.  2   ) ... [. I-3  ) [. I-2  ] (   -101  ...
          *         *         *             *         *         *
         G0        G1        G2            GI-3      GI-2      GI-1

*****************************************************************************/
      i = GetInterval( x , G , I , i );
      if( i == -2 ){
        i = 0;
      } else if( i == -101 ){
        i = I-1;
      } else if( fabs( x - G[ i ] ) > fabs( x - G[i+1] ) ){
        i++;
      }        
      out[n] = i+1;
    }
    
  } else {
    
    if( InfOutside ){ myErrMsgTxt("No se puede usar la opcion 'inside' si la grilla no esta sorted\n"); }
    
    for( n = 0 ; n<N ; n++ ){
      x = X[n];
      closest  = 0;
      distance = fabs( x - G[0] );
      for( i = I-1 ; i ; i-- ){
//       for( i = 0 ; i < I ; i++ ){
        if( fabs( x - G[i] ) < distance ){
          closest  = i;
          distance = fabs( x - G[i] );
        }
        if( distance == 0 ){ break; }
      }
      out[n] = closest+1;
    }
    
  }

  
  EXIT:
    myFreeALLOCATES();
    
}
