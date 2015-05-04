/*


 Las sumas se hacen sin tener en cuenta los Infs y NaNs !!!
 
 [R,r] = dot_exact( a[N] , b[N] )          \sum_i^N   a_i * b_i
 [R,r] = dot_exact( a[N] , a[N] )          \sum_i^N   a_i ^ 2
 [R,r] = dot_exact( a[N] , 'sq' )          \sum_i^N   a_i ^ 2
 [R,r] = dot_exact( a[N] , s[1] )          s * \sum_i^N a_i
 [R,r] = dot_exact( a[N] )                 \sum_i^N a_i     ==  dot_exact( a[N] , 1 )
 
*/

#include "myMEX.h"


#define isNAN( x )    ( !( x == x ) )
#define SWAP( z , w )    temp = S->v[z]; S->v[z] = S->v[w]; S->v[w] = temp;



#define  SIZE_ACCUMULATOR  200
struct accumulator {
    double  v[ SIZE_ACCUMULATOR ];
    int     first_nz;
};
typedef struct accumulator accumulator;



void ADD( accumulator *S , double t );
void TIMES( accumulator *S , double s );
void SUM( accumulator *S );
void FULLSUM( accumulator *S );
double TwoSum( double a , double b , double *y );
double TwoSumWithCheck( double a , double b , double *y , char *updated );
double TwoProduct( double a , double b , double *y );
double Square( double a , double *y );




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum               OPERATORS { NONE , SQUARE };

	double             *A, *B, *OUTPUT;
  double             AB , rAB, Ai, Bi;
  accumulator        S;
  int                N, i;
  char               STR[100];
  
  enum    OPERATORS OPERATOR;



  
  S.first_nz = SIZE_ACCUMULATOR;
  for( i = 0 ; i < SIZE_ACCUMULATOR ; i++ ){ S.v[i] = 0; }
  
  
  N = myNumel( prhs[0] );
  
  
  OPERATOR = NONE;
  if( nrhs == 2 &&   mxIsChar( prhs[1] ) ){
    mxGetString( prhs[1], STR, 100 );
    if( ! myStrcmpi(STR,"sq")     || ! myStrcmpi(STR,"square")   ){ 
      
      OPERATOR = SQUARE;

    } else {
     
      myErrMsgTxt("invalid operator.");
      
    }
    
    
  }
  
  if(   (  nrhs == 2    &&   OPERATOR == SQUARE  )  ||
        (  nrhs == 2    &&   OPERATOR == NONE   &&   mxGetPr( prhs[0] ) == mxGetPr( prhs[1] )  )
    ){
    
    A = myGetPr( prhs[0] );

    for( i = 0 ; i < N  ; i++ ){
      Ai = A[i];
      if( Ai == 0 || mxIsInf( Ai ) || isNAN( Ai ) ){ continue; }
      AB = Square( Ai , &rAB );
      ADD( &S ,  AB );
      ADD( &S , rAB );
    }
    
  } else if( nrhs == 2    &&   OPERATOR == NONE   &&    myNumel( prhs[1] ) != 1 ){
  
    if( N  !=  myNumel( prhs[1] ) ){ myErrMsgTxt("Incorrect sizes."); }
   
    A = myGetPr( prhs[0] );
    B = myGetPr( prhs[1] );

    for( i = 0 ; i < N  ; i++ ){
      Ai = A[i]; Bi = B[i];

      if( Ai == 0       || Bi == 0       ||
          mxIsInf( Ai ) || mxIsInf( Bi ) || 
          isNAN( Ai )   || isNAN( Bi )   ){
            continue;
      }

      AB = TwoProduct( Ai , Bi , &rAB );

      ADD( &S ,  AB );
      ADD( &S , rAB );
    }
    
  }  else if( nrhs == 2    &&   OPERATOR == NONE   &&    myNumel( prhs[1] ) == 1 ){
    
    A = myGetPr( prhs[0] );

    for( i = 0 ; i < N  ; i++ ){
      Ai = A[i];

      if( Ai == 0 || mxIsInf( Ai ) || isNAN( Ai ) ){ continue; }

      ADD( &S ,  Ai );
    }
    
    Bi  = *myGetPr( prhs[1] );
    if( Bi != 1 ){
      TIMES( &S , Bi );
    }
    
  } else if( nrhs == 1 ){
    
    A = myGetPr( prhs[0] );

    for( i = 0 ; i < N  ; i++ ){
      Ai = A[i];
      if( Ai == 0 || mxIsInf( Ai ) || isNAN( Ai ) ){ continue; }
      ADD( &S ,  Ai );
    }
    
    
  } else {
    myErrMsgTxt("Incorrect syntax.");
  }

  FULLSUM( &S );
  
  
  
  
  
  
  
  plhs[0] = mxCreateDoubleScalar( S.v[ SIZE_ACCUMULATOR - 1 ] );

  if( nlhs > 1 ){

    N = SIZE_ACCUMULATOR - 1 - S.first_nz;
    

    if( N > 0 ){
    
      plhs[1] = mxCreateDoubleMatrix( 1 , N , mxREAL );

      OUTPUT = mxGetPr( plhs[1] );
      for( i = 0 ; i < N ; i++ ){
        OUTPUT[ i ] = S.v[ SIZE_ACCUMULATOR - 2 - i ];
      }
      
    } else {
      
      plhs[1] = mxCreateDoubleScalar( 0 );
    
    }

  }
  
    
  EXIT:
    myFreeALLOCATES();

}



double TwoSum( double a , double b , double *y ){
  double x;
  double z;
  
  x  = a + b;
  z  = x - a;
  *y = ( a - (x-z) ) + ( b - z );

  
  return( x );
}


double TwoSumWithCheck( double a , double b , double *y , char *updated ){
  double x;
  double z;
  
  x  = a + b;
  z  = x - a;
  *y = ( a - (x-z) ) + ( b - z );
  
  if( !*updated && x != b ){ 
    *updated = 1; 
  }

  return( x );
}


double TwoProduct( double a , double b , double *y ){
  double  x;
  double  a1, a2, b1, b2, c;
  
  
  c  = 134217729.0 * a;
  a1 = c - ( c - a );
  a2 = a - a1;
  
  
  c  = 134217729.0 * b;
  b1 = c - ( c - b );
  b2 = b - b1;
  
  x = a*b;
  
  *y = ( a2 * b2 ) - ( ( ( x - a1*b1) - a2*b1 ) - a1*b2 );
  
  return( x );
}


double Square( double a , double *y ){
  double  x;
  double  a1, a2, c;
  
  
  c  = 134217729.0 * a;
  a1 = c - ( c - a );
  a2 = a - a1;

  c  = a1*a2;
  
  x = a*a;
  
  *y = ( a2 * a2 ) - ( ( ( x - a1*a1) - c ) - c );
  
  return( x );
}


void ADD( accumulator *S , double t ){
  
  S->v[ S->first_nz - 1 ] = t;
  S->first_nz--;

  if( S->first_nz < 2 ){
    SUM( S );
  }
  
}


void TIMES( accumulator *S , double s ){
  int     k, first_nz;
  double  r;
  
  
  FULLSUM( S );

  first_nz = S->first_nz;
  for( k = SIZE_ACCUMULATOR - 1 ; k >= first_nz ; k-- ){
    
    S->v[k] = TwoProduct( S->v[k] , s , &r );
    
    S->v[ S->first_nz - 1 ] = r;
    S->first_nz--;
    
  }
  
  
}

void SUM( accumulator *S ){
  int  k, last_z, last_nz;
  double temp;
  
  while( S->first_nz < SIZE_ACCUMULATOR - 50 ){
    while( S->v[ S->first_nz ] == 0  &&  S->first_nz < SIZE_ACCUMULATOR - 1 ){ S->first_nz++; }

    for( k = S->first_nz + 1 ; k < SIZE_ACCUMULATOR ; k++ ){
      S->v[k] = TwoSum( S->v[k-1] , S->v[k] ,  ( S->v ) + k - 1 );
    }

  
  
//           last_z     = SIZE_ACCUMULATOR-1;
//           while( last_z > S->first_nz  &&  S->v[last_z] != 0 ){
//             last_z--;
//           }
//           
//           last_nz = last_z-1;
//           
//           while( last_z > S->first_nz  &&  last_nz > S->first_nz ){
//             
//             while( last_nz >  S->first_nz  &&  S->v[last_nz] == 0  ){
//               last_nz--;
//             }
// 
//             SWAP( last_z , last_nz );
//             last_z--; last_nz--;
//             
//             while( last_z >  S->first_nz  &&  S->v[last_z] != 0  ){
//               last_z--;
//             }
//             
//           }
  
  
  
  
  
  
  
  
  }
  while( S->v[ S->first_nz ] == 0  &&  S->first_nz < SIZE_ACCUMULATOR - 1 ){ S->first_nz++; }

}




void FULLSUM( accumulator *S ){
  int  k, IT;
  char updated;

  for( IT = 0 ; IT < 1500 ; IT++ ){
    while( S->v[ S->first_nz ] == 0  &&  S->first_nz < SIZE_ACCUMULATOR - 1 ){ S->first_nz++; }

    updated = 0;
    for( k = S->first_nz + 1 ; k < SIZE_ACCUMULATOR ; k++ ){
      S->v[k] = TwoSumWithCheck( S->v[k-1] , S->v[k] ,  ( S->v ) + k - 1 , &updated );
    }
    
    if( updated == 0 ){ 
      break; 
    }
  }

  while( S->v[ S->first_nz ] == 0  &&  S->first_nz < SIZE_ACCUMULATOR - 1 ){ S->first_nz++; }

}

