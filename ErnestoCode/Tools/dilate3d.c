#include "myMEX.h"

#define O(i,j,k)  O[ (i) + (j)*I + (k)*IJ ]
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  int             n, i, j ,k, I, J, K, IJ, I1, J1, K1;
  int             i0,i1,j0,j1,k0,k1,slice;
  mwSize          sizes[3];
  mxLogical       *V, *O;
  int             connect;

  if( nrhs < 1 ){
    mxErrMsgTxt("sintax error\n");
  }
  
  if( !mxIsLogical( prhs[0] ) ){
    mxErrMsgTxt("array have to be logical.\n");
  }
  if( myNDims( prhs[0] ) > 3 ){
    mxErrMsgTxt("bigger than 3d arrays is not allowed.\n");
  }
  
  if( nrhs > 1 ){
    connect = myGetValue( prhs[1] );
    if( connect != 26 && connect != 18 && connect != 6 ){
      mxErrMsgTxt("valid connects are 6, 18, 26.\n");
    }
  } else {
    connect = 26;
  }
  

  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );
  K  = mySize( prhs[0] , 2 );
  IJ = I*J;
  I1 = I-1;
  J1 = J-1;
  K1 = K-1;
  
  
  sizes[0] = I;
  sizes[1] = J;
  sizes[2] = K;
  
  V = (mxLogical *) mxGetData( prhs[0] );

  plhs[0] = mxCreateLogicalArray( 3 , sizes );
  O = (mxLogical *) mxGetData( plhs[0] );
  

  if( connect == 26 ){

    n = 0;
    for( k = 0 ; k < K ; k++ ){
      k0 = k > 0  ? k - 1 : k;
      k1 = k < K1 ? k + 1 : k;

      for( j = 0 ; j < J ; j++ ){
        j0 = j > 0  ? j - 1 : j;
        j1 = j < J1 ? j + 1 : j;

        for( i = 0 ; i < I ; i++ ){

          if( V[n] ){
            i0 = i > 0  ? i - 1 : i;
            i1 = i < I1 ? i + 1 : i;

            O( i0 , j0 , k0 ) = 1;
            O( i  , j0 , k0 ) = 1;
            O( i1 , j0 , k0 ) = 1;
            O( i0 , j  , k0 ) = 1;
            O( i  , j  , k0 ) = 1;
            O( i1 , j  , k0 ) = 1;
            O( i0 , j1 , k0 ) = 1;
            O( i  , j1 , k0 ) = 1;
            O( i1 , j1 , k0 ) = 1;

            O( i0 , j0 , k  ) = 1;
            O( i  , j0 , k  ) = 1;
            O( i1 , j0 , k  ) = 1;
            O( i0 , j  , k  ) = 1;
            O[ n            ] = 1;
            O( i1 , j  , k  ) = 1;
            O( i0 , j1 , k  ) = 1;
            O( i  , j1 , k  ) = 1;
            O( i1 , j1 , k  ) = 1;

            O( i0 , j0 , k1 ) = 1;
            O( i  , j0 , k1 ) = 1;
            O( i1 , j0 , k1 ) = 1;
            O( i0 , j  , k1 ) = 1;
            O( i  , j  , k1 ) = 1;
            O( i1 , j  , k1 ) = 1;
            O( i0 , j1 , k1 ) = 1;
            O( i  , j1 , k1 ) = 1;
            O( i1 , j1 , k1 ) = 1;

          }
          n++;
        }
      }
    }

  } else if( connect == 18 ){

    n = 0;
    for( k = 0 ; k < K ; k++ ){
      k0 = k > 0  ? k - 1 : k;
      k1 = k < K1 ? k + 1 : k;

      for( j = 0 ; j < J ; j++ ){
        j0 = j > 0  ? j - 1 : j;
        j1 = j < J1 ? j + 1 : j;

        for( i = 0 ; i < I ; i++ ){

          if( V[n] ){
            i0 = i > 0  ? i - 1 : i;
            i1 = i < I1 ? i + 1 : i;

            O( i  , j0 , k0 ) = 1;
            O( i0 , j  , k0 ) = 1;
            O( i  , j  , k0 ) = 1;
            O( i1 , j  , k0 ) = 1;
            O( i  , j1 , k0 ) = 1;

            O( i0 , j0 , k  ) = 1;
            O( i  , j0 , k  ) = 1;
            O( i1 , j0 , k  ) = 1;
            O( i0 , j  , k  ) = 1;
            O[ n            ] = 1;
            O( i1 , j  , k  ) = 1;
            O( i0 , j1 , k  ) = 1;
            O( i  , j1 , k  ) = 1;
            O( i1 , j1 , k  ) = 1;

            O( i  , j0 , k1 ) = 1;
            O( i0 , j  , k1 ) = 1;
            O( i  , j  , k1 ) = 1;
            O( i1 , j  , k1 ) = 1;
            O( i  , j1 , k1 ) = 1;

          }
          n++;
        }
      }
    }

  } else if( connect == 6 ){

    n = 0;
    for( k = 0 ; k < K ; k++ ){
      k0 = k > 0   ? k - 1 : k;
      k1 = k < K-1 ? k + 1 : k;

      for( j = 0 ; j < J ; j++ ){
        j0 = j > 0   ? j - 1 : j;
        j1 = j < J-1 ? j + 1 : j;

        for( i = 0 ; i < I ; i++ ){

          if( V[n] ){
            i0 = i > 0   ? i - 1 : i;
            i1 = i < I-1 ? i + 1 : i;

            O( i0 , j  , k  ) = 1;
            O( i  , j  , k  ) = 1;
            O( i1 , j  , k  ) = 1;

            O( i  , j0 , k  ) = 1;
            O[ n            ] = 1;
            O( i  , j1 , k  ) = 1;

            O( i  , j  , k0 ) = 1;
            O( i  , j  , k  ) = 1;
            O( i  , j  , k1 ) = 1;

          }
          n++;
        }
      }
    }

  }

  
}
