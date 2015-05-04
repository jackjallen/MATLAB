/*

addpath g:\MoreWork\MatlabTools\gdft\

x = rand(6,2);
frt_mx(x,'dct' ) ./ dtt(x,'dct2e',1)   ,  frt_mx(x,'dct'  ) ./ dct(x)
frt_mx(x,'idct') ./ dtt(x,'dct3e',1)   ,  frt_mx(x,'idct' ) ./ idct(x)

frt_mx(x,'dct-i'  ) ./ dtt(x,'dct1e',1)
frt_mx(x,'dct-iv' ) ./ dtt(x,'dct4e',1)


frt_mx(x,'dst-i'  ) ./ dtt(x,'dst1e',1)
frt_mx(x,'dst-ii' ) ./ dtt(x,'dst2e',1)
frt_mx(x,'dst-iii') ./ dtt(x,'dst3e',1)
frt_mx(x,'dst-iv' ) ./ dtt(x,'dst4e',1)

frt_mx(x,'hartley') ./ ( gallery('orthog', size(x,1) ,5)*x )
 
 
[
 frt_mx( x(:,1) , 'dct1no' ) ./ dctn( x(:,1) , [] , '1no' ) ...
 frt_mx( x(:,2) , 'dct1no' ) ./ dctn( x(:,2) , [] , '1no' ) ...
 frt_mx( x(:,3) , 'dct1no' ) ./ dctn( x(:,3) , [] , '1no' ) ...
 frt_mx( x(:,4) , 'dct1no' ) ./ dctn( x(:,4) , [] , '1no' ) ...
] 
*/
 
#include "mex.h"
#include "fftw3.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  fftw_plan       plan;
  fftw_r2r_kind   TYPE;
  int             M , N, i, j, numel;
  double          f, f2, *O;
  char            TYPE_STR[100];

  M     = mxGetM( prhs[0] );
  N     = mxGetN( prhs[0] );
  
  mxGetString( prhs[1], TYPE_STR, 100 );
       if( !strcmp( TYPE_STR , "dct-ii"  ) || !strcmp( TYPE_STR , "dct"   ) ){   TYPE = FFTW_REDFT10; }
  else if( !strcmp( TYPE_STR , "dct-iii" ) || !strcmp( TYPE_STR , "idct"  ) ){   TYPE = FFTW_REDFT01; }
  else if( !strcmp( TYPE_STR , "dct-i"   ) ){   TYPE = FFTW_REDFT00; }
  else if( !strcmp( TYPE_STR , "dct-iv"  ) ){   TYPE = FFTW_REDFT11; }
  
  else if( !strcmp( TYPE_STR , "dst-i"   ) ){   TYPE = FFTW_RODFT00; }
  else if( !strcmp( TYPE_STR , "dst-ii"  ) ){   TYPE = FFTW_RODFT10; }
  else if( !strcmp( TYPE_STR , "dst-iii" ) ){   TYPE = FFTW_RODFT01; }
  else if( !strcmp( TYPE_STR , "dst-iv"  ) ){   TYPE = FFTW_RODFT11; }
  
  else if( !strcmp( TYPE_STR , "hartley" ) ){   TYPE = FFTW_DHT;     }
  
  else if( !strcmp( TYPE_STR , "dct1no"  ) ){   TYPE = 100;     }

  
  switch( TYPE ){

    case 100:      // dct-i no ortonormal
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );

      TYPE = FFTW_REDFT00;
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      break;

    
    case FFTW_REDFT00:   // frt_mx(x,'dct-i') ./ dtt(x,'dct1e',1)
      plhs[0] = mxDuplicateArray( prhs[0] );
      
      O = mxGetPr( plhs[0] ) - M;
      for( j = 0 ; j < N ; j++ ){
        O = O + M;
        O[ 0 ] *= 1.4142135623730951;
        O[M-1] *= 1.4142135623730951;
      }
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( plhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2 * (M-1) );
      
      O = mxGetPr( plhs[0] ) - M;
      for( j = 0 ; j < N ; j++ ){
        O = O + M;

        O[ 0 ] *= 0.7071067811865475;
        O[M-1] *= 0.7071067811865475;
        for ( i = 0 ; i < M  ; i++ ){ 
          O[i] *= f; 
        }
      }
      break;


    case FFTW_REDFT10:   // frt_mx(x,'dct-ii') ./ dtt(x,'dct2e',1)
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2.0 * M );
      f2 = f*0.7071067811865475;
      
      O = mxGetPr( plhs[0] ) - M;
      for( j = 0 ; j < N ; j++ ){
        O = O + M;

        O[0] *= f2;
        for ( i = 1 ; i < M  ; i++ ){ 
          O[i] *= f; 
        }
      }
      break;


    case FFTW_REDFT01:   // frt_mx(x,'dct-iii') ./ dtt(x,'dct3e',1)
      plhs[0] = mxDuplicateArray( prhs[0] );
      
      O = mxGetPr( plhs[0] ) - M;
      for( j = 0 ; j < N ; j++ ){
        O = O + M;
        O[0] = O[0]*1.4142135623730951;
      }
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( plhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2.0 * M );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;
      

    case FFTW_REDFT11:   // frt_mx(x,'dct-iv') ./ dtt(x,'dct4e',1)
                         // gallery('orthog', size(x,1) ,6)*x ./ frt_mx( x , 'dct-iv' )
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2.0 * M );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;


  
      
      
    case FFTW_RODFT00:   // frt_mx(x,'dst-i') ./ dtt(x,'dst1e',1)
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      

      f = 1/sqrt( 2*(M+1) );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;

      
    case FFTW_RODFT10:   // frt_mx(x,'dst-ii') ./ dtt(x,'dst2e',1)
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2.0 * M );
      
      O = mxGetPr( plhs[0] ) - M;
      for( j = 0 ; j < N ; j++ ){
        O = O + M;

        O[M-1] *= 0.7071067811865475;
        for ( i = 0 ; i < M  ; i++ ){ 
          O[i] *= f; 
        }
      }
      break;
  
      
    case FFTW_RODFT01:   // frt_mx(x,'dst-iii') ./ dtt(x,'dst3e',1)
      plhs[0] = mxDuplicateArray( prhs[0] );
      
      O = mxGetPr( plhs[0] );
      for( j = 0 ; j < N ; j++ ){
        O = O + M;
        O[-1] *= 1.4142135623730951;
      }
      
//       O = mxGetPr( plhs[0] ) - M;
//       for( j = 0 ; j < N ; j++ ){
//         O = O + M;
//         O[M-1] *= 1.4142135623730951;
//       }

      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( plhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      
      f  = 1.0/sqrt( 2.0 * M );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;
      

    case FFTW_RODFT11:   // frt_mx(x,'dst-iv') ./ dtt(x,'dst4e',1)
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      

      f = 1/sqrt( 2*M );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;
      
      
    case FFTW_DHT:   // gallery('orthog', size(x,1) ,5)*x ./ frt_mx( x , 'hartley' )
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetN( plhs[0] , N ); mxSetM( plhs[0] , M );
      mxSetData( plhs[0] , mxMalloc( mxGetNumberOfElements(prhs[0]) << 3 ) );
      
      plan = fftw_plan_many_r2r( 1 , &M , N , mxGetPr( prhs[0] ) , NULL , 1 , M , mxGetPr( plhs[0] ) , NULL , 1 , M , &TYPE , FFTW_ESTIMATE );
      fftw_execute(plan); fftw_destroy_plan(plan);
      

      f = 1/sqrt( M );
      
      O = mxGetPr( plhs[0] );
      numel = mxGetNumberOfElements( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f; 
      }
      break;
      
  }

  
}




// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
// 
//   fftw_plan       plan;
//   fftw_r2r_kind   types[10];
//   double          f2, f, *O;
//   int             N, M, j, i;
//   
//   
// 
//   f  = 1.0/sqrt( 2 * M );
//   f2 = f*0.7071067811865475;
//   for( j = 0 ; j < N ; j++ ){
//     O = mxGetPr( plhs[0] ) + j*M;
//     
//     O[0] *= f2;
//     for ( i = 1 ; i < M  ; i++ ){ O[i] *= f; }
//   }
//   
//   
//   
//   fftw_destroy_plan(plan);
//   




