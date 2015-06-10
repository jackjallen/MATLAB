
#define MAX(a,b) (a)>(b)?(a):(b)
#define MIN(a,b) (a)<(b)?(a):(b)
#define MMAX(a,b,c) (a)>(b)?((a)>(c)?(a):(c)):((b)>(c)?(b):(c))
#define MMIN(a,b,c) (a)<(b)?((a)<(c)?(a):(c)):((b)<(c)?(b):(c))


void SCAMCOORD(int *ijk, const mxArray *m, double *points, const long n , int *N){
  long      p;
  int       d;
  mxArray   *SCAM;
  double    O[3], D;

  SCAM= mxGetField( m, 0, "SCAM");
  
/* //   mexPrintf("EN SCAMCOORD\n\n"); mexEvalString("drawnow expose;"); */
  
  if( SCAM != NULL ) {
/* //     mexPrintf("No es null\n\n"); mexEvalString("drawnow expose;"); */
    
    D= *mxGetPr( mxGetField( SCAM, 0 , "D"));
    for( d=0 ; d<3 ; d++ ){
      O[d]= *(mxGetPr( mxGetField( SCAM, 0 , "O"))+d);
      N[d]= (int)*(mxGetPr( mxGetField( SCAM, 0 , "N"))+d);
    }

/* //     mexPrintf("empieza pintos\n\n"); mexEvalString("drawnow expose;");*/ 
    
    for( d=0 ; d<3 ; d++ ){
      for( p=0 ; p<n ; p++){
        ijk[ p+d*n ]= (int)( (points[ p+d*n ]-O[d])/D)+1;
        ijk[ p+d*n ]= (int)MAX( 1   , ijk[ p+d*n ] );
        ijk[ p+d*n ]= (int)MIN( N[d], ijk[ p+d*n ] );  
      }
    }
  } else {
    for( d=0 ; d<3 ; d++ ){
      N[d]= 1;
      for( p=0 ; p<n ; p++){
        ijk[ p+d*n ]= 1;
        ijk[ p+d*n ]= 1;
        ijk[ p+d*n ]= 1;  
      }
    }
  }
/* //   mexPrintf("Sale de EN SCAMCOORD\n\n"); mexEvalString("drawnow expose;"); */

}
