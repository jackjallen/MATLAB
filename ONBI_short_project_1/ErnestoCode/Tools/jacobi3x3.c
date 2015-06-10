#include "mex.h"
#include "math.h"

#define N  3
#define NN 9

void sorteigs(double *d, double *v) {
  int      k,j,i;
  double   p;

  for( i=0 ; i<N ; i++ ){
    k = i;
    p = d[ i+N*i ];
    for( j=i+1 ; j<N ; j++ ){
      if( d[ j+N*j ] < p ){
        k = j;
        p = d[ j+N*j ];
      }
    }
    if( k != i ){
      d[k+N*k] = d[i+N*i];
      d[i+N*i] = p;
      for( j=0 ; j<N ; j++ ) {
        p       = v[j+N*i];
        v[j+N*i]= v[j+N*k];
        v[j+N*k]= p;
      }
    }
  }
}

void diagonalize( double *a , double *v, double *b, double *z) {
	int    k,j,i,iter;
	double theta,tau,t,s,h,g,c;

	for( i=0 ; i<NN ; i++ ){
    v[i] = 0.0;
  }
	for( i=0 ; i<N ; i++ ){
    v[ i+N*i ] = 1.0;
		b[i]       = a[ i+N*i ];
	}

  iter= 1;
  while( iter ) {
    iter = 0;
		for( i=0 ; i < N-1 ; i++ ) {
			for( j=i+1 ; j < N ; j++ ) {
        if( fabs( a[i+N*j] ) > 1e-20 ) {
          iter= 1;
          h     = a[j+N*j]-a[i+N*i];
          theta = 0.5*h/a[i+N*j];
          if( theta >= 0.0 ){
            t =  1.0/(fabs(theta)+sqrt(1.0+theta*theta));
          } else {
            t = -1.0/(fabs(theta)+sqrt(1.0+theta*theta));
          }
          c = 1.0/sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*a[i+N*j];
          z[i] -= h;
          z[j] += h;
          a[i+N*i] -= h;
          a[j+N*j] += h;
          a[i+N*j] =  0.0;
          for( k=0    ; k < i ; k++ ) {
            g = a[ k+N*i ];
            h = a[ k+N*j ];
            a[ k+N*i ] = g-s*(h+g*tau);
            a[ k+N*j ] = h+s*(g-h*tau);
          }
          for( k = i+1  ; k<j ; k++ ) {
            g = a[ i+N*k ];
            h = a[ k+N*j ];
            a[ i+N*k ] = g-s*(h+g*tau);
            a[ k+N*j ] = h+s*(g-h*tau);
          }
          for( k=j+1  ; k<N   ; k++ ) {			
            g = a[ i+N*k ];
            h = a[ j+N*k ];
            a[ i+N*k ] = g-s*(h+g*tau);
            a[ j+N*k ] = h+s*(g-h*tau);
          }
          for( k=0    ; k<N   ; k++ ) {		
            g = v[ k+N*i ];
            h = v[ k+N*j ];
            v[ k+N*i ] = g-s*(h+g*tau);
            v[ k+N*j ] = h+s*(g-h*tau);
          }
        }
			}
		}
		for( i=0 ; i<N ; i++ ){
			z[i]  = 0.0;
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double    *a, *v;
  double    z[N], b[N];
  int       i,k;
  long      numel;
  
  numel= mxGetNumberOfElements( prhs[0] );
  numel= numel/NN;

  for( i=0 ; i<N ; i++ ){  z[i]= 0.0;  }

  plhs[0]= mxDuplicateArray( prhs[0] );
  plhs[1]= mxDuplicateArray( prhs[0] );

  v = mxGetPr( plhs[0] );
  a = mxGetPr( plhs[1] );
  for( k=1 ; k<=numel ; k++ ){
//     mexPrintf("%d\n" , k );

    diagonalize( a , v , b , z );
    sorteigs( a , v );

    a[1]= 0.0;
    a[2]= 0.0;
    a[3]= 0.0;
    a[5]= 0.0;
    a[6]= 0.0;
    a[7]= 0.0;

    v += NN;
    a += NN;
  }
}
