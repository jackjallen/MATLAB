/*

addpath g:\MoreWork\MatlabTools\gdft\

x = rand(3,5,2,2,2,7);
clc
maxnorm(unique(vec(   frtn_mx(x,'dct'    ) ./ dttn(x,'dct2e',1)   )).'-1)  %,  frtn_mx(x,'dct'  ) ./ dctn(x)
maxnorm(unique(vec(   frtn_mx(x,'idct'   ) ./ dttn(x,'dct3e',1)   )).'-1)  %,  frtn_mx(x,'idct' ) ./ idctn(x)
maxnorm(unique(vec(   frtn_mx(x,'dct-i'  ) ./ dttn(x,'dct1e',1)   )).'-1)
maxnorm(unique(vec(   frtn_mx(x,'dct-iv' ) ./ dttn(x,'dct4e',1)   )).'-1)

maxnorm(unique(vec(   frtn_mx(x,'dst-i'  ) ./ dttn(x,'dst1e',1)   )).'-1)
maxnorm(unique(vec(   frtn_mx(x,'dst-ii' ) ./ dttn(x,'dst2e',1)   )).'-1)
maxnorm(unique(vec(   frtn_mx(x,'dst-iii') ./ dttn(x,'dst3e',1)   )).'-1)
maxnorm(unique(vec(   frtn_mx(x,'dst-iv' ) ./ dttn(x,'dst4e',1)   )).'-1)

maxnorm(unique(vec(   frtn_mx(x,'hartley') ./ dttn(x,'hartley')   )).'-1)
 
*/


#include "mex.h"
#include "fftw3.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  fftw_plan       plan;
  fftw_r2r_kind   types[64], TYPE;
  fftw_iodim      iodims[64];
  int             numel, i, c, j, id;
  int             cumprod[64];
  double          f, *O;
  int             ndims;
  const int       *dims;
  char            TYPE_STR[100];

  
  numel = mxGetNumberOfElements( prhs[0] );
  dims  = mxGetDimensions( prhs[0] );
  ndims = mxGetNumberOfDimensions( prhs[0] );
  

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

  
  
  f = 1;
  cumprod[0] = 1;
  for( c = 0; c < ndims ; c++ ){ 
    types[c] = TYPE; 

    iodims[ c ].n  = dims[ c ];
    iodims[ c ].is = iodims[ c ].os = cumprod[c];

    cumprod[c+1]  = cumprod[c] * dims[c];
    
    f  =  f * 2 * dims[c];
  }
  
  
  switch( TYPE ){
    
    case 100:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      for( c = 0; c < ndims ; c++ ){ 
        types[c] = FFTW_REDFT00; 
      }      
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);
      
      break;
      
      
    
    
    case FFTW_REDFT00:
      plhs[0] = mxDuplicateArray( prhs[0] );
      O = mxGetPr( plhs[0] );

      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( !id                 ){  O[i] *= 1.4142135623730951;  }
          if( id == ( dims[c]-1 ) ){  O[i] *= 1.4142135623730951;  }
          j -= id*cumprod[c];
        }
      }

      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( plhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1;
      for( c = 0 ; c < ndims ; c++ ){
        f = f * 2 * ( dims[c] - 1 );
      }
      
      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( !id                 ){  O[i] *= 0.7071067811865475;  }
          if( id == ( dims[c]-1 ) ){  O[i] *= 0.7071067811865475;  }
          j -= id*cumprod[c];
        }
        O[i] *= f;
      }
      break;
  
    case FFTW_REDFT10:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1./sqrt(f);

      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( !id ){
            O[i] *= 0.7071067811865475;
          }
          j -= id*cumprod[c];
        }
        O[i] *= f;
      }
      break;

      
    case FFTW_REDFT01:
      plhs[0] = mxDuplicateArray( prhs[0] );
      O = mxGetPr( plhs[0] );

      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( !id ){
            O[i] *= 1.4142135623730951;
          }
          j -= id*cumprod[c];
        }
      }
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( plhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;
      

    case FFTW_REDFT11:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;

 
    case FFTW_RODFT00:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1;
      for( c = 0 ; c < ndims ; c++ ){
        f = f * 2 * ( dims[c] + 1 );
      }
      
      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;
  
      
    case FFTW_RODFT10:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( id == ( dims[c]-1 ) ){  O[i] *= 0.7071067811865475;  }
          j -= id*cumprod[c];
        }
        O[i] *= f;
      }
      break;
      
      
    case FFTW_RODFT01:
      plhs[0] = mxDuplicateArray( prhs[0] );
      O = mxGetPr( plhs[0] );

      for( i = 0 ; i < numel ; i++ ){
        j = i;
        for( c = ndims-1 ; c >= 0 ; c-- ){
          id = (int) ( j / cumprod[c] );
          if( id == ( dims[c]-1 ) ){  O[i] *= 1.4142135623730951;  }
          j -= id*cumprod[c];
        }
      }

      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( plhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;
  

    case FFTW_RODFT11:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);


      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;
  
      
    case FFTW_DHT:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
      mxSetDimensions( plhs[0] , dims , ndims );
      mxSetData( plhs[0] , mxMalloc( numel << 3 ) );
      
      plan = fftw_plan_guru_r2r( ndims , iodims , 0 , NULL , mxGetPr( prhs[0] ) ,  mxGetPr( plhs[0] ) , types , FFTW_ESTIMATE );
      fftw_execute(plan);   fftw_destroy_plan(plan);

      f = 1;
      for( c = 0 ; c < ndims ; c++ ){
        f = f * dims[c];
      }

      f = 1./sqrt(f);
      O = mxGetPr( plhs[0] );
      for( i = 0 ; i < numel ; i++ ){
        O[i] *= f;
      }
      break;
      
  }
  
}
