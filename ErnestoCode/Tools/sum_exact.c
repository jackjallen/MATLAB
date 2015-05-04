/*

rand('state',0);
z = nonans(noinfs( [ 10.^([-300:0.01:300]) -rand(1)*10.^([-300:0.01:300]) 1 drand(1,10000)/100 ] ));
z = z( NTHrandperm(0,numel(z) ) );
[ss,rr] = sum_exact_m( z );
 

tic; [r,r] = sum_exact( z , 0 ); toc; uneval(r,r)

 *
 *
 * 
 
 rand('state',0)
 z = noinfs( nonans( [ 10.^(-300:.1:300) drand(1,100) ] ));
 z0 = [ z 1 -0.5 -0.5 -z 2e-4*ones(1,1e4) 1e14 1];

 clc 
 for p = 1:100
   z = z0( NTHrandperm( p , numel(z0) ) );

   for K = [1:50 0]
    [r,r] = sum_exact(z,K); 
    disp( uneval( K , r , r ) );
   end
   disp('');
   disp('');
 end
 
 
 *
 *
 z=zeros(1,100);r=[1e307 10.^(-323:0.020:304)  ]; b = 1e10;
 z = [b z r z r z z b z r z r z];
 for K = 0:50
 [r,r] = sum_exact( z , K ); disp( uneval(K,r,r) );
 end 
 
*/

#include "mex.h"

double TwoSum( double a , double b , double *y );

#define isPOS( x )    ( mxIsInf(x)  &&  (x) > 0  )
#define isNEG( x )    ( mxIsInf(x)  &&  (x) < 0  )
#define isNAN( x )    ( !( x == x ) )
#define SWAP( z , w )    temp = X[z]; X[z] = X[w]; X[w] = temp;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  mxClassID         CLASS;
  int               N, K, K1, K2, i, k;
  double            *Pd;
  float             *Ps;
  double            S, r, a;
  double            *X = NULL;
  double            temp;
  int               IT = 0;
  int               first_no_zero, last_zero, last_no_zero;
  
  
  
  CLASS = mxGetClassID( prhs[0] );
  switch( CLASS ){
    case mxDOUBLE_CLASS:
      Pd = (double *) mxGetData( prhs[0] );
      break;

    case mxSINGLE_CLASS:
      Ps = (float  *) mxGetData( prhs[0] );
      break;

    default:
      mxErrMsgTxt("first arguments must be a float");
  }
  
  
  N = mxGetNumberOfElements( prhs[0] );
  
  switch( N ){
    case 0:
      S = 0;
      r   = 0;
      break;
      
    case 1:
      if( CLASS == mxDOUBLE_CLASS ){
        S = (double) Pd[0];
      } else {
        S = (double) Ps[0];
      }
      r   = 0;
      break;
  
    case 2:
      if( CLASS == mxDOUBLE_CLASS ){
        S = TwoSum( (double) Pd[0] , (double) Pd[1] , &r );
      } else {
        S = TwoSum( (double) Ps[0] , (double) Ps[1] , &r );
      }
      break;


    default:
      if( nrhs > 1 ){
        K = *( mxGetPr( prhs[1] ) );
      } else {
        K = 0;
      }
      
      if( K < 0 ){ mxErrMsgTxt("K >=0 expected"); }
      
      if( K == 0 ){
        
        X = mxMalloc( sizeof( double )* N );
        
        if( CLASS == mxDOUBLE_CLASS ){
          for( i = 0 ; i < N ; i++ ){ X[i] = (double) Pd[i]; }
        } else {
          for( i = 0 ; i < N ; i++ ){ X[i] = (double) Ps[i]; }
        }
        
        
        first_no_zero = 0;
        for( IT = 0 ; IT < 45 ; IT++ ){
          while( X[first_no_zero] == 0 && first_no_zero < N-1 ){ first_no_zero++; }
          
//           mexPrintf( "IT: %5d    -  first_no_zero = %5d\n", IT , first_no_zero);
          
          for( k = first_no_zero + 1 ; k < N ; k++ ){
            X[k] = TwoSum( X[k-1] , X[k] , X + k - 1 );
          }
          

          if( isPOS( X[N-1] ) || isNEG( X[N-1] ) || isNAN( X[N-1] ) ){
            first_no_zero = N-1;
            break;
          }
          


          
          
          last_zero     = N-1;
          while( last_zero > first_no_zero  &&  X[last_zero] != 0 ){
            last_zero--;
          }
          
          last_no_zero = last_zero-1;
          
          while( last_zero > first_no_zero  &&  last_no_zero > first_no_zero ){
            
            while( last_no_zero >  first_no_zero  &&  X[last_no_zero] == 0  ){
              last_no_zero--;
            }

            SWAP( last_zero , last_no_zero );
            last_zero--; last_no_zero--;
            
            while( last_zero >  first_no_zero  &&  X[last_zero] != 0  ){
              last_zero--;
            }
            
          }
          
          
          
          
        }
        
        
        
        
        while( X[first_no_zero] == 0 && first_no_zero < N-1 ){ first_no_zero++; }
        S = X[N-1];

        
      } else {

        if( N < K ){ 
          K = N; 
        }

        K1 = K-1;
        K2 = K-2;

        X = mxMalloc( sizeof( double )* K );
        for( i = 0 ; i < K ; i++ ){
          X[i] = 0;
        }

        if( CLASS == mxDOUBLE_CLASS ){

          for( i = 0 ; i < K1 ; i++ ){
            r = (double) Pd[i];
            for( k = 0 ; k < i-1 ; k++ ){
              X[k] = TwoSum( X[k] , r , &r ); 
            }
            X[i] = r;
          }


          for( i = K1 ; i < N ; i++ ){
            a = (double) Pd[i];
            for( k = 0 ; k < K1 ; k++ ){
              X[k] = TwoSum( X[k] , a , &a );
            }
            r = r + a;
          }

        } else {


          for( i = 0 ; i < K1 ; i++ ){
            r = (double) Ps[i];
            for( k = 0 ; k < i-1 ; k++ ){
              X[k] = TwoSum( X[k] , r , &r ); 
            }
            X[i] = r;
          }


          for( i = K1 ; i < N ; i++ ){
            a = (double) Ps[i];
            for( k = 0 ; k < K1 ; k++ ){
              X[k] = TwoSum( X[k] , a , &a );
            }
            r = r + a;
          }

        }


        for( i = 0 ; i < K2 ; i++ ){
          a = X[i];
          for( k = i+1 ; k < K1 ; k++ ){
            X[k] = TwoSum( X[k] , a , &a );
          }
          r = r + a;
        }

        S = TwoSum( r , X[K2] , &r );
        
      }

  }
  
    
  plhs[0] = mxCreateDoubleScalar( S );

  if( nlhs > 1 ){
    if( K > 0 ){
      plhs[1] = mxCreateDoubleScalar( r );
    } else {
      N = N - first_no_zero;
      
      if( N > 1 ){
        plhs[1] = mxCreateDoubleMatrix( 1 , N - 1 , mxREAL );

        Pd = mxGetPr( plhs[1] );
        for( k = 0; k < N - 1 ; k++ ){
          Pd[ N - 2 - k ] = X[ k + first_no_zero ];
        }
      } else {
        plhs[1] = mxCreateDoubleScalar( 0 );
      }
    }
  }
  
  
  if(  X != NULL ){  mxFree( X ); }
  
}


double TwoSum( double a , double b , double *y ){
  double x;
  double z;
  
  x  = a + b;
  z  = x - a;
  *y = ( a - (x-z) ) + ( b - z );
  
  return( x );
}

