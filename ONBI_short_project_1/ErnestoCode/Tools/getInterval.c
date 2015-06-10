/*****************************************************************************

I = numel(G)
 
...  0   ) [.  1   ) [.  2   ) [.  3   ) ... [. I-2  ) [. I-1  .] (   I  ....
          *         *         *             *         *          *
         G1        G2        G3            GI-2      GI-1       GI

tic
G = 0:20;
X = [ randn( 1 , 100000 ) linspace(-1,20,10000) -10:.1:30 ];
plot( X , getInterval(X,G) , '.' )
toc

 *
 *
 *
 *
hay que arreglar    getInterval( 0 , [-10 10])  , a veces da 1 y a veces 3
 *
 *
*****************************************************************************/

#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real  *G, *X, *O;
  int     I, n, N;
  real  boundary_size;
  
  char boundary_str[100];
  
  X = myGetPr( prhs[0] );
  N = myNumel( prhs[0] );
  
  G = myGetPr( prhs[1] );
  I = myNumel( prhs[1] );


  if( !checkIsSorted(G,I) ){ myErrMsgTxt("G  Coordinates are not sorted."); }

  plhs[0]  = myDuplicateArrayWithClass( prhs[0] , mxREAL_CLASS , mxREAL );
  O = (real *)mxGetData( plhs[0] );

  for(n=0;n<N;n++){
    O[n] = GetInterval( X[n] , G , I , -1 )+1;
    if( O[n] ==  -1  ){ 
      O[n] = 0; 
    } else if( O[n] == -100 ){ 
      O[n] = I; 
    }
//     O[n] = GetInterval( X[n] , G , I , -1 );
  }
  
  EXIT: myFreeALLOCATES();
}
