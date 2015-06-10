/*

dist_int = @(x,x0) 1/(1-tanh(abs(x-x0)/(2*(22.4./(abs(x0)+ 4)))).^3);
dist_mat = reshape(sqrt(sum(ndmat(-5:5,-5:5,-5:5).^2,2)),[11 11 11]);
maximos = max_min_multiples_erodes( T , @(ijk) Dist_Fun_in_c( a.data , dist_mat , dist_int , ijk ) , 1:1000 )
 
*/
#include "myMEX.h"



#define  EVERY  1000

struct triplet {
  unsigned short i;
  unsigned short j;
  unsigned short k;
  char isnan;
}; typedef struct triplet triplet;


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ ALLOCATES();

  CreateTicTacToc( CallMatlab );
  CreateTicTacToc( callSort   );
  int               I, J, K, ii, jj, kk;
  int               IJ, IJK, IJK_1;
  int               DI, DJ, DK, DIJK, DIJK_1;
  int               CDI, CDJ, CDK;
  int               result, fevals = 0;
  int               NVOLS, NVOLS_1, n, s, s_start, s_end, v_init;
  real              *volumes, *V, x, y, *DIST, *order, last_distance;
  int               *VV=NULL, nV, v, vv;
  char              skip;
  triplet           *TS=NULL, *DTS=NULL, T;
  mxArray           *INPUT[2]={NULL,NULL}, *OUTPUT[3]={NULL,NULL,NULL};
  double            *MAXs, LAST_MAX;
  double            thisMINx, thisMINy;
  double            *idxs;
  double            *vols;
  double            *ijk;
  char              callSort;
  mwSize            toVec[2]={1,1};
  char              VERBOSE = 0;
  char              STR[1024];
  
  
  if( nlhs > 1 ){
    mxErrMsgTxt("too much outputs");
  }


  if( mxIsChar( prhs[nrhs-1] ) ){
    mxGetString( prhs[nrhs-1], STR, 100 );
    if( ! myStrcmpi(STR,"verbose") ){ 
      VERBOSE = 1;
    } else {
      mxErrMsgTxt("only 'verbose' option allowed.");
    }
    nrhs = nrhs-1;
  }
  
  
  if( nrhs != 3 ){
    mxErrMsgTxt("sintax error. max_min_multiples_erodes( V , F , volumes )");
  }

  if( mxGetClassID( prhs[1] ) != mxFUNCTION_CLASS ){
    mxErrMsgTxt("F have to be a function_handle.");
  }

  if( myNDims( prhs[0] ) > 3  ){
    mxErrMsgTxt("bigger than 3d arrays is not allowed.");
  }
  
  NVOLS   = myNumel( prhs[2] );
  NVOLS_1 = NVOLS - 1;
  volumes = myGetPr( prhs[2] );
  
  
  I     = mySize( prhs[0] , 0 );
  J     = mySize( prhs[0] , 1 );
  K     = mySize( prhs[0] , 2 );
  IJ    = I*J;
  IJK   = IJ*K;
  

  VV = (int     *) mxMalloc( IJK*sizeof( int     ) );
  TS     = (triplet *) mxMalloc( IJK*sizeof( triplet ) );

  V = myGetPr( prhs[0] );

  v  = 0; 
  nV = 0;
  for( kk = 0 ; kk < K ; kk++ ){ for( jj = 0 ; jj < J ; jj++ ){ for( ii = 0 ; ii < I ; ii++ ){
    x = V[ v ];
    if( x == x ){
      VV[ nV ] = v;
      TS[ v ].isnan = 0;
      TS[ v ].i     = ii;
      TS[ v ].j     = jj;
      TS[ v ].k     = kk;
      nV++;
    } else {
      TS[ v ].isnan = 1;
    }
    v++;
  }}}


  INPUT[0] = prhs[1];
  INPUT[1] = mxCreateNumericMatrix( 1 , 3 , mxDOUBLE_CLASS , mxREAL );
  ijk = (double *) mxGetData( INPUT[1] );
  
  
  ijk[0] = TS[ VV[ nV/2] ].i + 1;
  ijk[1] = TS[ VV[ nV/2] ].j + 1;
  ijk[2] = TS[ VV[ nV/2] ].k + 1;
  


  OUTPUT[2] = mexCallMATLABWithTrap( 2 , OUTPUT , 2 , INPUT , "feval" );
  
  if( OUTPUT[2] == NULL ){
    
    callSort = 0;
    if( mxGetClassID( OUTPUT[0] ) != mxDOUBLE_CLASS ){
      if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1] );  INPUT[1]=NULL;  }
      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
      mxErrMsgTxt("F debe retornar un double en el primer output.");
    }
    if( mxGetClassID( OUTPUT[1] ) != mxDOUBLE_CLASS ){
      if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1] );  INPUT[1]=NULL;  }
      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
      mxErrMsgTxt("F debe retornar un double en el segundo output.");
    }
    
  } else {

    callSort = 1;
    if( VERBOSE ){
      mexPrintf("sort has to be called\n");
    }
    
    mxDestroyArray( OUTPUT[2] ); OUTPUT[2] = NULL;

    result = mexCallMATLAB( 1 , OUTPUT , 2 , INPUT , "feval" );
    if( result ){ mxErrMsgTxt("error computing la funcion."); }

    if( mxGetClassID( OUTPUT[0] ) != mxDOUBLE_CLASS ){
      if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1] );  INPUT[1]=NULL;  }
      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
      mxErrMsgTxt("F debe retornar un double en el primer output.");
    }

  }

  DI   = mySize( OUTPUT[0] , 0 );
  DJ   = mySize( OUTPUT[0] , 1 );
  DK   = mySize( OUTPUT[0] , 2 );
  
  DTS  = (triplet *) mxMalloc( 2*DI*DJ*DK*sizeof( triplet ) );
  

  plhs[0] = mxCreateNumericMatrix( NVOLS , 1 , mxREAL_CLASS , mxREAL );
  MAXs    = (real *) mxGetData( plhs[0] );
  for( n = 0 ; n < NVOLS ; n++ ){
    MAXs[n] = -10000;
  }

  
  LAST_MAX = MAXs[ NVOLS_1 ];
  for( v_init = 0 ; v_init < EVERY ; v_init++ ){
    if( utIsInterruptPending() ){
      if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1] );  INPUT[1]=NULL;  }
      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
      mexPrintf("USER INTERRUP!!!\n");
      mxErrMsgTxt("USER INTERRUP!!!");
    }
    if( VERBOSE ){
      mexPrintf("v_init:  %d  (%g)  of  %d\n", v_init , LAST_MAX , EVERY );
    }

    for( v = v_init ; v < nV ; v += EVERY ){
      vv = VV[ v ];
      
      thisMINx =   V[ vv ];
      thisMINy =  -thisMINx;
      if( ( thisMINx < LAST_MAX )  && ( thisMINy < LAST_MAX ) ){
        continue;
      }

      T = TS[ vv ];
      
      ijk[0] = T.i + 1;
      ijk[1] = T.j + 1;
      ijk[2] = T.k + 1;
      

      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }

      if( !callSort ){
        tic( CallMatlab );
        result = mexCallMATLAB( 2 , OUTPUT , 2 , INPUT , "feval" ); fevals++;
        tac( CallMatlab );
      } else {
        tic( CallMatlab );
        result = mexCallMATLAB( 1 , OUTPUT , 2 , INPUT , "feval" ); fevals++;
        tac( CallMatlab );
      }
      
      DI   = mySize( OUTPUT[0] , 0 );
      DJ   = mySize( OUTPUT[0] , 1 );
      DK   = mySize( OUTPUT[0] , 2 );

      DIJK    = DI*DJ*DK;

      if( volumes[ NVOLS_1 ] > DIJK ){
      if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1] );  INPUT[1]=NULL;  }
      if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
      if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
      if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
        mxErrMsgTxt("el maximo volumen debe ser menor que numel(DIST)");
      }

      DIJK_1  = DIJK - 1;

      DIST  = (double  *) mxGetData( OUTPUT[0] );
      
      DTS  = (triplet *) mxRealloc( DTS , DIJK*sizeof( triplet ) );
      s = 0;
      for( kk = 0 ; kk < DK ; kk++ ){ for( jj = 0 ; jj < DJ ; jj++ ){ for( ii = 0 ; ii < DI ; ii++ ){
        DTS[ s ].i = ii;
        DTS[ s ].j = jj;
        DTS[ s ].k = kk;
        s++;
      }}}


      if( !callSort ){
        order = (double  *) mxGetData( OUTPUT[1] );
      } else {
        toVec[0] = mxGetNumberOfElements( OUTPUT[0] );
        mxSetDimensions( OUTPUT[0] ,  toVec  , 2 );
        
        tic( callSort );
        result = mexCallMATLAB( 2 , OUTPUT+1 , 1 , OUTPUT , "sort" );
        tac( callSort );
      
        order = (double  *) mxGetData( OUTPUT[2] );
      }
      
      CDI = DTS[ (int) ( order[0] - 1 ) ].i;
      CDJ = DTS[ (int) ( order[0] - 1 ) ].j;
      CDK = DTS[ (int) ( order[0] - 1 ) ].k;
      
      
      skip   = 0;

      s = 0;
      for( n = 0 ; n < NVOLS ; n++ ){
        s_end = (int) ( volumes[n] - 1 );
        last_distance = DIST[ (int) order[ s_end ] - 1 ];
        
        while( s_end < DIJK_1 && DIST[ (int) ( order[ s_end + 1 ] - 1 ) ] == last_distance ){
          s_end++;
        }
        s_end++;
        
        for( ; s < s_end ; s++ ){
          vv = (int) ( order[ s ] - 1 );
          
          ii = T.i + DTS[ vv ].i - CDI;  if( ii < 0 || ii > I ){ skip = 1; break; }
          jj = T.j + DTS[ vv ].j - CDJ;  if( jj < 0 || jj > J ){ skip = 1; break; }
          kk = T.k + DTS[ vv ].k - CDK;  if( kk < 0 || kk > K ){ skip = 1; break; }
            
          vv = ii + jj*I + kk*IJ;
          if( TS[ vv ].isnan ){  skip = 1; break; }
          x =  V[ vv ];  if( x < thisMINx ){ thisMINx = x; }
          y = -x;        if( y < thisMINy ){ thisMINy = y; }
          if( ( thisMINx < LAST_MAX )  && ( thisMINy < LAST_MAX ) ){
            skip = 1; break; 
          }
        }
        if( skip ){  break; }
        if( thisMINx > MAXs[n] ){ MAXs[n] = thisMINx; }
        if( thisMINy > MAXs[n] ){ MAXs[n] = thisMINy; }
      }
      LAST_MAX = MAXs[ NVOLS_1 ];

    }
    
  }
  if( INPUT[1]  != NULL ){ mxDestroyArray( INPUT[1]  ); INPUT[1] =NULL; }
  if( OUTPUT[0] != NULL ){ mxDestroyArray( OUTPUT[0] ); OUTPUT[0]=NULL; }
  if( OUTPUT[1] != NULL ){ mxDestroyArray( OUTPUT[1] ); OUTPUT[1]=NULL; }
  if( OUTPUT[2] != NULL ){ mxDestroyArray( OUTPUT[2] ); OUTPUT[2]=NULL; }
  
  if( VERBOSE ){
    mexPrintf( "\nfevals: %d  en  tiempo:   CallMatlab: %20.30g    sorting: %20.30g\n" , fevals , toc( CallMatlab ) , toc( callSort ) );
  }
  
  
  if(  VV != NULL ){  mxFree(  VV ); }
  if(  TS != NULL ){  mxFree(  TS ); }
  if( DTS != NULL ){  mxFree( DTS ); }

  myFreeALLOCATES();

}
