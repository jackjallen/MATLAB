/*

M = I3D('X',-20:20,'Y',-20:20,'Z',-20:20);
M.data(:) = sqrt( sum( M.XYZ.^2 , 2 ) );

I = single( randn( [220 200 180 ] ) );
 
tic;
[m,E] = max_grayscale_erode_3d( I , crop( M.data <= 4.8 ) );
toc

tic;
mm = max_grayscale_erode_3d( I , crop( M.data <= 4.8 ) );
toc

[ mm - m , mm - max(E(:)) ]

tic;
EE = imerode( I , M.data <= 4.8 );
toc

maxnorm( E - EE )
  
*/


#include "myMEX.h"


#define  INFp   1000000
#define  INFn  -1000000

#define V(i,j,k)  V[ (i) + (j)*I  + (k)*IJ  ]
#define E(i,j,k)  E[ (i) + (j)*I  + (k)*IJ  ]
#define M(i,j,k)  M[ (i) + (j)*MI + (k)*MIJ ]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ ALLOCATES();

  int             I, J, K, IJ, MI, MJ, MK, MIJ, RI, RJ, RK, LI, LJ, LK, FI, FJ, FK;
  int             ci, cj, ck, ii, jj, kk;
  mwSize          sizes[3];
  mxLogical       *M;
  double          *E, *V;
  float           *Vs;
  double          thisMIN, thisV, maximum;
  float           thisMINs, thisVs, maximums;
  int             *MASK_indexs, N_MASK, m, *V_indexs, N_V, v;
  unsigned char   skip;
  

  if( nlhs > 2 ){
    myErrMsgTxt("too much outputs\m");
  }

  if( nrhs < 2 ){
    myErrMsgTxt("sintax error. max_grayscale_erode_3d( IMAGE , MASK )\m");
  }

  if( !mxIsLogical( prhs[1] ) ){
    myErrMsgTxt("MASK have to be logical.\m");
  }

  if( myNDims( prhs[0] ) > 3  || myNDims( prhs[1] ) > 3 ){
    myErrMsgTxt("bigger than 3d arrays is not allowed.\m");
  }

  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );
  K  = mySize( prhs[0] , 2 );
  IJ = I*J;


  MI  = mySize( prhs[1] , 0 );  if( MI % 2 == 0 ){ myErrMsgTxt("MASK has to have odd size along 1st dim.\m"); }
  MJ  = mySize( prhs[1] , 1 );  if( MJ % 2 == 0 ){ myErrMsgTxt("MASK has to have odd size along 2nd dim.\m"); }
  MK  = mySize( prhs[1] , 2 );  if( MK % 2 == 0 ){ myErrMsgTxt("MASK has to have odd size along 3rd dim.\m"); }
  MIJ = MI*MJ;

  
  RI = (MI-1)/2;  LI = I - RI; FI = MIN(RI,I); LI = MAX(LI,0);
  RJ = (MJ-1)/2;  LJ = J - RJ; FJ = MIN(RJ,J); LJ = MAX(LJ,0);
  RK = (MK-1)/2;  LK = K - RK; FK = MIN(RK,K); LK = MAX(LK,0);



  
  M = (mxLogical *) mxGetData( prhs[1] );
  MASK_indexs = (int *) mxMalloc( MIJ * MK *sizeof( int ) );
  m = 0;
  for( kk = 0 ; kk < MK ; kk++ ){  for( jj = 0 ; jj < MJ ; jj++ ){  for( ii = 0 ; ii < MI ; ii++ ){
    if( M(ii,jj,kk) ){
      MASK_indexs[ m ] = ( ii - RI ) + ( jj - RJ )*I + ( kk - RK )*IJ;
      m++;
    }
  }}}
  N_MASK = m;


  if( nlhs < 2 ) {

    V_indexs = (int *) mxMalloc( ( LI - FI + 1 )*( LJ - FJ + 1 )*( LK - FK + 1 )*sizeof( int ) );
    v = 0;
    for( ck = FK ; ck < LK ; ck++ ){ for( cj = FJ ; cj < LJ ; cj++ ){ for( ci = FI ; ci < LI ; ci++ ){
      V_indexs[ v ] = ci + cj*I  + ck*IJ;
      v++;
    }}}
    N_V = v;
    

    if( mxGetClassID( prhs[0] ) == mxSINGLE_CLASS ){
      
      Vs = (float *) mxGetData( prhs[0] );

      maximums = INFn;
      
      if( N_V ){ 
        v = N_V/2;
        thisMINs = INFp; skip = 0;
        for( m = 0 ; m < N_MASK ; m++ ){
          thisVs = Vs[ V_indexs[v] + MASK_indexs[m] ];
          if( thisVs != thisVs   ){ skip = 1; break; }
          if( thisVs <  thisMINs ){ thisMINs = thisVs; }
        }
        if( !skip ){ maximums = thisMINs; }
      }
      

      for( v = 0 ; v < N_V ; v++ ){
        thisMINs = INFp;  skip    = 0;

        for( m = 0 ; m < N_MASK ; m++ ){
          thisVs = Vs[ V_indexs[v] + MASK_indexs[m] ];
          if( thisVs != thisVs   ||  
              thisVs <  maximums   ){ skip = 1; break; }
          if( thisVs <  thisMINs   ){ thisMINs = thisVs; }
        }

        if( !skip ){ maximums = thisMINs; }
      }
      
      
    } else {
    
      V = myGetPr( prhs[0] );

      maximum = INFn;
      
      if( N_V ){ 
        v = N_V/2;
        thisMIN = INFp; skip = 0;
        for( m = 0 ; m < N_MASK ; m++ ){
          thisV = V[ V_indexs[v] + MASK_indexs[m] ];
          if( thisV != thisV   ){ skip = 1; break; }
          if( thisV <  thisMIN ){ thisMIN = thisV; }
        }
        if( !skip ){ maximum = thisMIN; }
      }
      

      for( v = 0 ; v < N_V ; v++ ){
        thisMIN = INFp;  skip    = 0;

        for( m = 0 ; m < N_MASK ; m++ ){
          thisV = V[ V_indexs[v] + MASK_indexs[m] ];
          if( thisV != thisV   ||  
              thisV <  maximum   ){ skip = 1; break; }
          if( thisV <  thisMIN   ){ thisMIN = thisV; }
        }

        if( !skip ){ maximum = thisMIN; }
      }
      
    }
    
    mxFree( V_indexs );
    
    
  } else {

    V = myGetPr( prhs[0] );
    maximum = INFn;
    
    sizes[0] = I;
    sizes[1] = J;
    sizes[2] = K;

    plhs[1] = mxCreateNumericArray( 3 , sizes , mxREAL_CLASS , mxREAL );
    E = (real *) mxGetPr( plhs[1] );
    
    
    for( ck = 0 ; ck < FK ; ck++ ){ for( cj = 0 ; cj < J ; cj++ ){ for( ci = 0 ; ci < I ; ci++ ){       
          E( ci , cj , ck ) = NAN;
    }}}
    for( ck = LK ; ck < K ; ck++ ){ for( cj = 0 ; cj < J ; cj++ ){ for( ci = 0 ; ci < I ; ci++ ){       
          E( ci , cj , ck ) = NAN;
    }}}
    for( ck = 0 ; ck < K ; ck++ ){
      for( cj = 0 ; cj < FJ ; cj++ ){
        for( ci = 0 ; ci < I ; ci++ ){       
          E( ci , cj , ck ) = NAN;
        }
      }
      for( cj = LJ ; cj < J ; cj++ ){
        for( ci = 0 ; ci < I ; ci++ ){       
          E( ci , cj , ck ) = NAN;
        }
      }
    }
    for( ck = 0 ; ck < K ; ck++ ){
      for( cj = 0 ; cj < J ; cj++ ){
        for( ci = 0  ; ci < FI ; ci++ ){ E( ci , cj , ck ) = NAN; }
        for( ci = LI ; ci < I  ; ci++ ){ E( ci , cj , ck ) = NAN; }
      }
    }
    
    
    for( ck = FK ; ck < LK ; ck++ ){ for( cj = FJ ; cj < LJ ; cj++ ){ for( ci = FI ; ci < LI ; ci++ ){
      thisMIN = INFp; skip = 0;

      for( m = 0 ; m < N_MASK ; m++ ){
        thisV = V( ci + MASK_indexs[m] , cj , ck );
        if( thisV != thisV   ){ skip = 1; break; }
        if( thisV  < thisMIN ){ thisMIN = thisV; }
      }

      if( skip ){
        E( ci , cj , ck ) = NAN;
      } else {
        if( thisMIN > maximum ){ maximum = thisMIN; }
        E( ci , cj , ck ) = thisMIN;
      }
      
    }}}  

    maximums = (float) maximum;
    
  }

  mxFree( MASK_indexs );
  
  if( mxGetClassID( prhs[0] ) == mxSINGLE_CLASS ){

    plhs[0] = mxCreateNumericMatrix( 1 , 1 , mxSINGLE_CLASS , mxREAL );
    *( (float *) mxGetData( plhs[0] ) ) = maximums;

  } else {

    plhs[0] = mxCreateNumericMatrix( 1 , 1 , mxDOUBLE_CLASS , mxREAL );
    *( (double *) mxGetData( plhs[0] ) ) = maximum;

  }
  
  EXIT: myFreeALLOCATES();
}
