/*

v = MixtureOfGaussians( C , w , h , Y , cutoff )

 C: centers of gaussians
 w: weight of each gaussian ( each column, for each dim )
    if size(w,1) == 1 , all gaussians with the same weight
 h: sigma of each gaussian
    if numel(h) == 1 , all gaussians with the same sigma
 Y: points where eval
 
 v_i = sum_{c=1}^{#C}  w_c exp( -| y_i - C_c |^2 / h_c^2 )

 if  h  isscalar  ->  h_c = h  \forall c
 
 if  size(w) == [  nC    vec ]
 
 if  size( w , 1 ) == 1  ->  w_c = w  \forall c
 
 
 size( v ) = [ size(Y,1)   size(q,2) ]

*/

#include "myMEX.h"

#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real    *C, *Y, *W , *H, *V;
  int     nC, nY, cc , yy;
  int     uniqueW, uniqueH;
  int     Odims[2];
  int     nsd, vec, d; 
  
  real    cutoff, this_cutoff, dist2 , this_dist , this_w, this_h, expo;
  
  
  if( nrhs < 4 || nrhs > 5 ){
    myErrMsgTxt( "4 or 5 arguments expected." );
  }
  
  
  nsd = mySize( prhs[0] , 1 );  //number of spatial dimensions
  nC  = mySize( prhs[0] , 0 );  //number of center points
  
  vec = mySize( prhs[1] , 1 );  //vectorial dimension
  if( mySize( prhs[1] , 0 ) == 1 ){
    uniqueW = 1;
  } else {
    uniqueW = 0;
    if( mySize( prhs[1],0 ) != nC ){
      myErrMsgTxt("numel(W,1) do not coincide with size(C,1).");
    }
  }
    
  if( mySize( prhs[2] , 1 ) > 1 ){
    myErrMsgTxt("size(H,2) has to be 1.");
  }
  if( mySize( prhs[2] , 0 ) == 1 ){
    uniqueH = 1;
  } else {
    uniqueH = 0;
    if( mySize( prhs[2],0 ) != nC ){
      myErrMsgTxt("numel(H,1) do not coincide with size(C,1).");
    }
  }
  
  if( mySize( prhs[3] , 1 ) != nsd ){
    myErrMsgTxt("numel(Y,2) do not coincide with size(C,2).");
  }
  nY = mySize( prhs[3] , 0 );    //number of points where eval
//   DISP( nY );
  
  
  C = myGetPr( prhs[0] );
  W = myGetPr( prhs[1] );
  H = myGetPr( prhs[2] );
  Y = myGetPr( prhs[3] );
  
  
  if( nrhs == 5 ){
    cutoff = myGetValue( prhs[4] );
    cutoff = MIN( cutoff , 7 );
    cutoff = cutoff * cutoff ;
  } else {
    cutoff = 49;
  }
  

  /*Creating output*/
  Odims[0]= nY;
  Odims[1]= vec;
  plhs[0] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
  V = (real*) mxGetData( plhs[0] );
  /*END Creating output*/

  this_h = H[0];
  this_h = this_h * this_h;
  this_cutoff = this_h * cutoff;
  this_w = W[0];  

//   for( yy = 0 ; yy < nY ; yy++ ){
//     for( cc = 0; cc < nC ; cc++ ){
//       dist = 0;
//       for( d = 0 ; d < nsd ; d++ ){
//         this_dist = C[ cc + d*nC ] - Y[ yy + d*nY ];
//         dist +=  this_dist * this_dist;
//       }
//       V[ yy ] += this_w * exp( - dist / this_h );
//     }
//   }
  
  for( cc = 0; cc < nC ; cc++ ){
    
    if( !uniqueH ){
      this_h      = H[cc];
      this_h      = this_h * this_h;
      this_cutoff = this_h * cutoff;
    };
    if( !uniqueW ){
      this_w      = W[cc];
    }
    
    for( yy = 0; yy < nY ; yy++ ){
      dist2 = 0;
      for( d = 0 ; d < nsd   &&  dist2 < this_cutoff  ; d++ ){
        this_dist = C[ cc + d*nC ] - Y[ yy + d*nY ];
        dist2 +=  this_dist * this_dist;
      }

      if( dist2 < this_cutoff ){
        expo = exp( - dist2 / this_h );
        
        if( vec == 1 ){
          V[ yy ] += this_w * expo;
        } else {
          for( d = 0 ; d < vec ; d++ ){
            this_w      = uniqueW ? W[d] : W[cc + d*nC ];
            V[ yy + d*nY ] += this_w * expo;
          }
        }
      }
      
    }
    
  }
  

  EXIT: myFreeALLOCATES();
}

