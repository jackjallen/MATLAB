#ifndef __NEWMYMEX__
#define __NEWMYMEX__

#include "mex.h"
#include "math.h"

#if !defined( real )
  #define   real       double
#endif
#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

/* gav typedef real double gav */

char *uneval( double x ){

  mxArray   *X[1], *STR[1];
  char      *string;

  X[0] = mxCreateDoubleScalar( x );
  mexCallMATLAB( 1 , STR , 1 , X , "uneval" );


  string = mxArrayToString( STR[0] );

  mxDestroyArray(   X[0] );
  mxDestroyArray( STR[0] );


  return( string );
}

#ifdef mxSetLogical
mxArray *mexCallMATLABWithTrap(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], const char *functionName){
  mxArray   *ME = NULL;
  int       result;

  mexSetTrapFlag(1);
  result = mexCallMATLAB(nlhs,plhs,nrhs,prhs,functionName);

  if( result ) {
    mexCallMATLAB( 1 , &ME , 0 , NULL , "lasterror" );
  }

  return( ME );
}
#endif


#if defined (_WIN32)
    #include <windows.h>
#elif defined (__linux__)
    #include <unistd.h>
#endif

#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif



mxArray *myCreateDoubleMatrix_E( mwSize m, mwSize n, mxComplexity ComplexFlag ){
  mxArray *M;

  M = mxCreateDoubleMatrix( 0, 0, ComplexFlag);
  mxSetM( M , m );
  mxSetN( M , n );

  mxSetData( M , mxRealloc( mxGetData( M ) , (m*n) << 3 ) );
  if( ComplexFlag == mxCOMPLEX ){
    mxSetImagData( M , mxRealloc( mxGetImagData( M ) , (m*n) << 3 ) );
  }

  return( M );
}


mxArray *myCreateNumericArray_E(mwSize ndim, const mwSize *dims, mxClassID classid, mxComplexity ComplexFlag){
  mxArray   *M;
  int       numel;

  M = mxCreateNumericArray( 0 , 0 , classid , ComplexFlag );
  mxSetDimensions( M , dims , ndim );

/*
//   numel = mxGetElementSize(M);
//   for( d = 0 ; d < ndim ; d++ ){
//     numel *= dims[d];
//   }
*/

  numel = mxGetElementSize(M)*mxGetNumberOfElements(M);
  mxSetData( M , mxRealloc( mxGetData( M ) , numel ) );

  if( ComplexFlag == mxCOMPLEX ){
    mxSetImagData( M , mxRealloc( mxGetImagData( M ) , numel ) );
  }

  return( M );
}


mxArray *myDuplicateArray_E(const mxArray *A){
  mxArray *M;
  int     numel;

  M = mxCreateNumericArray( 0 , 0 , mxGetClassID(A) , ( mxIsComplex(A) ) ? mxCOMPLEX : mxREAL  );
  mxSetDimensions( M , mxGetDimensions(A) , mxGetNumberOfDimensions(A) );

  numel = mxGetElementSize(A)*mxGetNumberOfElements(A);

  mxSetData( M , mxMalloc( numel ) );

  if( mxIsComplex(A) ){
    mxSetImagData( M , mxMalloc( numel ) );
  }

  return( M );
}


#define ALLOCATES()                 myALLOCS ALLOCS = {0}
#define myGetPr( mxA )              myGetPr_( mxA , &ALLOCS )
#define myFreeALLOCATES()           myFreeALLOCATES_( &ALLOCS );                        \
/*                                       mexPrintf("deleting ALLOCS before ERROR ... ");   \ */

/* #define myErrMsgTxt( err )          sprintf(myERROR, #err ); goto EXIT  */
#define myErrMsgTxt( ... )          mexErrMsgTxt( __VA_ARGS__ )

/* +0/0       -1.#IND  */
const union { int bits[2]; double value; } myNANp = {0,   -524288};
#if !defined( NAN )

/*
//   #define NAN                         myNANp.value
*/

#define NAN                         mxGetNaN()
#endif
#define toNAN(x)                    x = 0/(x*0)
/*
// #define myISNAN(x)                  ~(x==x)
*/
#define myISNAN(x)                  mxIsNaN( x )

/* -0/0        1.#QNAN */
const union { int bits[2]; double value; } myNANn = {0,2146959360};
#if !defined( NANn )
  #define NANn                         myNANn.value
#endif

/* +1/0        1.#INF */
const union { int bits[2]; double value; } myINFp = {0,2146435072};
#if !defined( INFp )
/*
//   #define INFp                         myINFp.value
*/
  #define INFp                         mxGetInf()
#endif
#if !defined( INF )
/*
//   #define INF                          myINFp.value
*/
  #define INF                          mxGetInf()
#endif

/* -1/0       -1.#INF */
const union { int bits[2]; double value; } myINFn = {0,  -1048576};
#if !defined( INFn )
  #define INFn                         myINFn.value
#endif

/*
// #define myISINF( x )                ( ( x ) == myINFp.value || ( x ) == myINFn.value )
// #define myISINFp( x )               ( ( x ) == myINFp.value )
// #define myISINFn( x )               ( ( x ) == myINFn.value )
*/
#define myISINF(  x )                 mxIsInf( x )
#define myISINFp( x )               ( mxIsInf( x ) && ( x > 0 ) )
#define myISINFn( x )               ( mxIsInf( x ) && ( x < 0 ) )


real myEPS( real x ){
  real y;

  union { unsigned int   bits[2]; double value; } ud;
  union { unsigned short bits[2]; float  value; } uf;

  if( sizeof(real) == 8 ) {

    ud.value = x;
    ud.bits[1] &= 0x7fffffff;
    y = ud.value;

    if( ud.bits[0] == 0xffffffff ){
      ud.bits[1] += 1;
      ud.bits[0] = 0x00000000;
    } else {
      ud.bits[0] += 1;
    }

    return( ud.value - y );

  } else {

    uf.value = x;
    uf.bits[1] &= 0x7fff;
    y = uf.value;

    if( uf.bits[0] == 0xffff ){
      uf.bits[1] += 1;
      uf.bits[0] = 0x0000;
    } else {
      uf.bits[0] += 1;
    }

    return( uf.value - y );

  }

}


#define myFlush()                   mexEvalString("drawnow expose;")
/*// #define myFlush()                   mexCallMATLAB(0,NULL,0,NULL,"drawnow expose;")*/
#define DISP(x)                     mexPrintf(". %s : %g\n" , #x , (double) x ); mexCallMATLAB(0,NULL,0,NULL,"drawnow('expose');")
#define myIsEmpty(mxA)              ( mxGetNumberOfElements( (mxA) ) == 0 )
#define myIsParameter( mxA )        ( mxGetNumberOfElements((mxA)) < 2 ) && ( mxIsNumeric((mxA)) || mxIsLogical((mxA)) )
#define myNumel( mxA )              ( mxGetNumberOfElements((mxA)) )
#define myNDims( mxA )              ( mxGetNumberOfDimensions((mxA)) )
#define mySize( mxA , i )           ( mxGetNumberOfDimensions( (mxA) ) > (i) ? (int) *( mxGetDimensions( (mxA) ) + i ) : 1 )



#ifdef _WIN32

  #include <windows.h>
  #define CreateTicToc()              LARGE_INTEGER __tic, __toc, __frec
  #define TIC()                       QueryPerformanceCounter(   &__tic  )
  #define TOC()                       QueryPerformanceCounter(   &__toc  );   \
                                      QueryPerformanceFrequency( &__frec );   \
                                      mexPrintf(".Elapsed time is %.8g seconds.\n", ( (double) (__toc.QuadPart - __tic.QuadPart) )/( (double) __frec.QuadPart )  )


  #define CreateTicTacToc(i)          LARGE_INTEGER __tic##_##i, __toc##_##i, __frec##_##i ; double __tac##_##i=0
  #define tic(i)                      QueryPerformanceCounter(  &__tic##_##i )
  #define tac(i)                      QueryPerformanceCounter(  &__toc##_##i );   \
                                      QueryPerformanceFrequency(&__frec##_##i);   \
                                      __tac##_##i = __tac##_##i + ( (double) __toc##_##i.QuadPart - (double) __tic##_##i.QuadPart )
  #define toc(i)                      __tac##_##i/( (double) __frec##_##i.QuadPart )

#else

  #include <time.h>
  #define CreateTicToc()              clock_t     __tic, __toc
  #define TIC()                       __tic = clock()
  #define TOC()                       __toc = clock();          \
                                      mexPrintf(".Elapsed time is %.8g seconds.\n", ( (double) (__toc - __tic) )/( (double) CLOCKS_PER_SEC )  )


  #define CreateTicTacToc(i)          clock_t __tic##_##i, __toc##_##i; double __tac##_##i=0
  #define tic(i)                      __tic##_##i = clock()
  #define tac(i)                      __tac##_##i = __tac##_##i + ( (double) clock() - (double) __tic##_##i )
  #define toc(i)                      __tac##_##i/( (double) CLOCKS_PER_SEC )


#endif













#define ABS(x)                      fabs( (x) )
#define SWAP(a, b)                  a ^= b; b ^= a; a ^= b



#if 0
mxArray * myCreateRealArray( int kk , ... ){
  mwSize  dims[100];
  mxArray *out;
  int     n;

  va_list sizes;
  va_start( sizes , 0);

  n = -1;
  do{
    dims[++n] = va_arg( sizes , int);
  }while( dims[n] > 0 ){
  va_end( sizes );

  out = mxCreateArray( dims , n , mxREAL_CLASS , mxREAL );
  return( out );
}
#endif


int myGetSizes( mxArray *M , int *D ){
  int n, N;
  N = mxGetNumberOfDimensions( M );
  for( n = 0 ; n < N ; n++ ){
    D[n] = (int) *( mxGetDimensions( M ) + n );
  }
  for( ; n<50 ; n++ ){ D[n] = -1; }
  return( N );
}

struct myALLOCS {
    int     n;
    real    *pointers[100];
}; typedef struct myALLOCS myALLOCS;

int GetInterval( real , real * , int , int );
real *DualGrid( real *, int );
int checkEqualSpaced( real *, int );
real myGetValue( const mxArray *);
real *myGetPr_( const mxArray *, myALLOCS *);
void myFreeALLOCATES_( myALLOCS *);
real myGetValueIth(const mxArray * , int );


#define IS_LOWERCASE(x)  (( (x)>='a') && ( (x) <='z') )
#define TO_UPPERCASE(x)  (IS_LOWERCASE (x)?(x)-'a'+'A':(x))
#define IS_UPPERCASE(x)  (( (x)>='A') && ( (x) <='Z') )
#define TO_LOWERCASE(x)  (IS_UPPERCASE (x) ? (x)-'A' + 'a': (x))
#define IS_CHAR(x)  ( ( (x)>='a') && ( (x) <='z') || ( (x)>='A') && ( (x) <='Z') ) ? true : false
int myStrcmpi( char *str1, char *str2 )
{
  int nn;
  for( nn = 0 ; str1[nn] && str2[nn] ; nn++){
    if( TO_LOWERCASE( str1[nn] ) != TO_LOWERCASE( str2[nn] ) ){
      return(1);
    }
  }
  if( str1[nn] || str2[nn] ){
    return(1);
  }
  return( 0 );
}

/*mmonla*/
void myStrTOUppercase( char **str, int iSizeSrc )
{
  /*int nn;
  for( nn = 0 ; nn < iSizeSrc ; nn++)
  {
	  (*str)[nn] = TO_UPPERCASE( (*str)[nn] );
  }
  return;*/
  int nn;
  nn = 0;
  while( ((*str)[nn] != '\0') && (nn < iSizeSrc) )
  {
	  (*str)[nn] = TO_UPPERCASE( (*str)[nn] );
	  nn++;
  }
  (*str)[nn] = '\0';
  return;
}
void myStrTOLowercase( char **str, int iSizeSrc )
{
  /*int nn;
  for( nn = 0 ; nn < iSizeSrc ; nn++)
  {
	  (*str)[nn] = TO_UPPERCASE( (*str)[nn] );
  }
  return;*/
  int nn;
  nn = 0;
  while( ((*str)[nn] != '\0') && (nn < iSizeSrc) )
  {
	  (*str)[nn] = TO_LOWERCASE( (*str)[nn] );
	  nn++;
  }
  (*str)[nn] = '\0';
  return;
}
int myStrGetUppercase( char **srcDst, char *srcStr, int iSizeSrc)
{
  int nn, iChFounds;
  iChFounds = 0;
  nn = 0;
  while( (srcStr[nn] != '\0') && (nn < iSizeSrc) )
  {
	  if( IS_UPPERCASE( srcStr[nn] ) )
	  {
		  (*srcDst)[iChFounds] = srcStr[nn];
		  iChFounds++;
	  }
	  nn++;
  }
  (*srcDst)[iChFounds] = '\0';

  return iChFounds;
}
/*mmonla*/

/* gav se puede usar directamente strcasecmp(str1, str2) */


real *myGetPr_(const mxArray *p, myALLOCS *ALLOCATES){
/*   CreateTicToc(); */
  long              n,N;
  mxClassID         CLASS;
  real              *P;

  double            *dataD;
  float             *dataF;
  unsigned char     *dataUC;
  char              *dataC;
  unsigned short    *dataUS;
  short             *dataS;
  unsigned int      *dataUI;
  int               *dataI;


  CLASS = mxGetClassID(p);
  if( CLASS == mxDOUBLE_CLASS && sizeof(real) == sizeof(double) ){
    return( (real *) mxGetData( p ) );
  }
  if( CLASS == mxSINGLE_CLASS && sizeof(real) == sizeof(float ) ){
    return( (real  *) mxGetData( p ) );
  }

  N = mxGetNumberOfElements( p );
  P = (real *) mxMalloc( sizeof(real)* N );
  if( P == NULL ){
    mexPrintf("Insuficient Memory to copy array\n");
  }

  switch( CLASS ){
    case mxDOUBLE_CLASS:
      dataD = (double *) mxGetData( p );
/*       TIC(); */
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataD[n]; }
/*       TOC(); */
/*       mexPrintf("Correctamente allocados %d elementos de DOUBLE a real(%d bytes)\n", N,sizeof(real)); myFlush(); */
      break;

    case mxSINGLE_CLASS:
      dataF = (float *) mxGetData( p );
/*       TIC(); */
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataF[n]; }
/*       TOC();  */
/*       mexPrintf("Correctamente allocados %d elementos de SINGLE a real(%d bytes)\n", N,sizeof(real)); myFlush(); */
      break;


    case mxUINT8_CLASS:
      dataUC = (unsigned char *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataUC[n]; }
      break;

    case mxINT8_CLASS:
      dataC = (char *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataC[n]; }
      break;

    case mxUINT16_CLASS:
      dataUS = (unsigned short *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataUS[n]; }
      break;

    case mxINT16_CLASS:
      dataS = (short *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataS[n]; }
      break;

    case mxUINT32_CLASS:
      dataUI = (unsigned int *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataUI[n]; }
      break;

    case mxINT32_CLASS:
      dataI = (int *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataI[n]; }
      break;

    case mxLOGICAL_CLASS:
      dataUC = (unsigned char *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataUC[n]; }
      break;

    case mxCHAR_CLASS:
      dataC = (char *) mxGetData( p );
      for( n=0 ; n<N ; n++ ){ P[n] = (real) dataC[n]; }
      break;
  }

  ALLOCATES->pointers[ALLOCATES->n] = P;
  ALLOCATES->n++;
  return( P );
}

void myFreeALLOCATES_(myALLOCS *ALLOCATES){
  int n;
  for( n = 0 ; n < ALLOCATES->n ; n++ ){
    mxFree( ALLOCATES->pointers[n] );
  }
}

real myGetValueIth(const mxArray *p , int offset ){
  switch( mxGetClassID(p) ){
    case mxDOUBLE_CLASS:    return( (real) *( offset + (double *)          mxGetData(p) ));
    case mxSINGLE_CLASS:    return( (real) *( offset + (float  *)          mxGetData(p) ));
    case mxLOGICAL_CLASS:   return( (real) *( offset + (unsigned char *)   mxGetData(p) ));
    case mxCHAR_CLASS:      return( (real) *( offset + (char *)            mxGetData(p) ));
    case mxUINT8_CLASS:     return( (real) *( offset + (unsigned char *)   mxGetData(p) ));
    case mxINT8_CLASS:      return( (real) *( offset + (char *)            mxGetData(p) ));
    case mxUINT16_CLASS:    return( (real) *( offset + (unsigned short *)  mxGetData(p) ));
    case mxINT16_CLASS:     return( (real) *( offset + (short *)           mxGetData(p) ));
    case mxUINT32_CLASS:    return( (real) *( offset + (unsigned int *)    mxGetData(p) ));
    case mxINT32_CLASS:     return( (real) *( offset + (int *)             mxGetData(p) ));
    case mxUINT64_CLASS:    return( (real) *( offset + (unsigned long *)   mxGetData(p) ));
    case mxINT64_CLASS:     return( (real) *( offset + (long *)            mxGetData(p) ));
  }
  return(0);
}


real myGetValue(const mxArray *p ){
  switch( mxGetClassID(p) ){
    case mxDOUBLE_CLASS:    return( (real) *(double *)          mxGetData(p) );
    case mxSINGLE_CLASS:    return( (real) *(float  *)          mxGetData(p) );
    case mxLOGICAL_CLASS:   return( (real) *(unsigned char *)   mxGetData(p) );
    case mxCHAR_CLASS:      return( (real) *(char *)            mxGetData(p) );
    case mxUINT8_CLASS:     return( (real) *(unsigned char *)   mxGetData(p) );
    case mxINT8_CLASS:      return( (real) *(char *)            mxGetData(p) );
    case mxUINT16_CLASS:    return( (real) *(unsigned short *)  mxGetData(p) );
    case mxINT16_CLASS:     return( (real) *(short *)           mxGetData(p) );
    case mxUINT32_CLASS:    return( (real) *(unsigned int *)    mxGetData(p) );
    case mxINT32_CLASS:     return( (real) *(int *)             mxGetData(p) );
    case mxUINT64_CLASS:    return( (real) *(unsigned long *)   mxGetData(p) );
    case mxINT64_CLASS:     return( (real) *(long *)            mxGetData(p) );
  }
  return(0);
}


mxArray * myDuplicateArrayWithClass( const mxArray *A , mxClassID classid , mxComplexity ComplexFlag ){
  mxArray *M;
  int     numel;
  int     sz;

   M = mxCreateNumericArray( mxGetNumberOfDimensions(A) , mxGetDimensions(A) , classid , ComplexFlag  );

/*
//   M = mxCreateNumericArray( 0 , 0 , classid , ComplexFlag  );
//
//   numel = mxGetNumberOfElements(A);
//   sz    = mxGetElementSize(A);
//
//   mxSetData( M , mxCalloc( numel , sz ) );
//
//   if( ComplexFlag ){
//     mxSetImagData( M , mxCalloc( numel , sz ) );
//   }
//
//   mxSetDimensions( M , mxGetDimensions(A) , mxGetNumberOfDimensions(A) );
*/

  return( M );
}



/*****************************************************************************

I = numel(G)

... -2   ) [.  0   ) [.  1   ) [.  2   ) ... [. I-3  ) [. I-2  ] (   -101  ...
          *         *         *             *         *         *
         G0        G1        G2            GI-3      GI-2      GI-1

*****************************************************************************/
int GetInterval( real x , real *G , int I , int lid ){
  int      id, iid, im, s;
  real     xm;

  if( x<G[0]   ){ return( -2 ); }
  if( x==G[0]  ){ return(  0 ); }
  I--;
  if( x>G[I]   ){ return(-101); }
  if( x==G[I]  ){ return( I-1); }

  if( lid < 0   ){ lid = (I+1) >> 1; }

  s = 1;
  id = lid;
  while( G[id] > x ){
    id -= s;
    if( id <= 0 ){
      id = 0;
      break;
    }
    s <<= 1;
  }
  if( G[id] == x ){ return(id); }

  s = 1;
  iid = lid + 1;
  while( G[iid] < x ){
    iid += s;
    if( iid >= I ){
      iid = I;
      break;
    }
    s <<= 1;
  }
  if( G[iid] == x ){ return(iid); }

  while( iid-id > 1 ){
    im = (id+iid) >> 1;
    xm = G[im];

    if( x < xm ){
      iid = im;
    }
    else if( x > xm ){
      id = im;
    }
    else {
      return(im);
    }
  }
  return(id);
}



real *DualGrid( real *X , int I){
  real    *DX;
  int     i;

  DX = (real *) mxMalloc( (I+1)*sizeof(real) );


  if( I == 1 ){
    DX[0] = X[0] - 0.5*0.1;
    DX[1] = X[0] + 0.5*0.1;
  } else {
    DX[0] = X[0] - (X[1]-X[0])/2;
    for( i = 1 ; i < I ; i++ ){
      DX[i] = (X[i-1]+X[i])/2;
    }
    DX[I] = X[I-1] + ( X[I-1]-X[I-2] )/2;
  }

  return( DX );
}

int checkEqualSpaced( real *G, int I ){
  int    i;
  real   d, eps;

  if( I <= 1 ){
    return(1);
  }
  d = G[1]-G[0];
  eps = d/1000.0;
  I -= 1;
  for( i = 1 ; i < I ; i++ ){
    if( fabs( G[i+1]-G[i] - d ) > eps ){
      return(0);
    }
  }
  return(1);
}


real MAX( real x, real y ){ if(x>y){ return(x); } return(y); }
real MIN( real x, real y ){ if(y>x){ return(x); } return(y); }

real MAX3( real x, real y , real z ){ return( MAX( MAX(x,y) , z ) ); }
real MIN3( real x, real y , real z ){ return( MIN( MIN(x,y) , z ) ); }


int checkIsSorted( real *G, int I ){
  int    i;

  if( I <= 1 ){ return(1); }
  for( i = 1 ; i < I ; i++ ){
    if( G[i-1] > G[i] ){
      return(0);
    }
  }
  return(1);
}

































/*
void myPrintMatrix(double* M, int rows, int cols) {
	mxArray* matrix;

  matrix = mxCreateDoubleMatrix( rows , cols, mxREAL );

  memcpy( mxGetPr( matrix ) , M , rows*cols*sizeof(double) );
	mexCallMATLAB(0, NULL, 1, &matrix, "disp");

  mxFree( matrix );
}
*/












#define isinteger(n) ( (int)(n) == (n) )

#define myIs2DMatrix( m )   ( mxGetNumberOfDimensions(m) == 2 )
#define myIs2DArray( m )    ( mxGetNumberOfDimensions(m) == 2 )
#define myIs3DArray( m )    ( mxGetNumberOfDimensions(m) == 3 )
#define myIsSquare( m )   ( ( mxGetNumberOfDimensions(m) == 2 )  && ( mxGetM(m) == mxGetN(m) ) )
#define myIsVector( m )   ( ( mxGetNumberOfDimensions(m) <= 2 )  && ( mxGetM(m) == 1 || mxGetN(m) == 1 ) && ( mxGetNumberOfElements(m) != 0 ) )
#define myIsScalar( m )   ( ( ~mxIsComplex(m) ) && ( mxGetNumberOfElements(m) == 1 ) )



/*
mwSize myNNZ(const mxArray* arr) {
	if (mxIsCell(arr)) {
		mwSize n = mxGetNumberOfElements(arr);
		mwSize c = 0;
		mwIndex ix;

		for (ix = 0; ix < n; ix++) {
			const mxArray* cell = mxGetCell(arr, ix);
			if (cell != NULL) {
				c += nonzerocount(cell);
			}
		}
		return c;
	} else if (mxIsStruct(arr)) {
		mexErrMsgIdAndTxt(__ARGTYPEMISMATCH__, "Structure arrays are not supported.");
	} else if (mxIsSparse(arr)) {
		mwIndex* jc = mxGetJc(arr);
		mwSize n = mxGetN(arr);
		return jc[n];
	} else {
		mwSize n = mxGetNumberOfElements(arr);
		mwSize c = 0;
		const double* pr = mxGetPr(arr);
		const double* pi = mxGetPi(arr);
		mwSize i;

		if (pr == NULL) {
			mexErrMsgIdAndTxt(__ARGTYPEMISMATCH__, "Operation supported only on numeric arrays.");
		}

		if (pi != NULL) {
			for (i = 0; i < n; i++) {
				if (isnonzero(*(pr++)) || isnonzero(*(pi++))) {
					c++;
				}
			}
		} else {
			for (i = 0; i < n; i++) {
				if (isnonzero(*(pr++))) {
					c++;
				}
			}
		}
		return c;
	}
}
*/



real *newmyGetData_v0(const mxArray *p){

  switch( mxGetClassID(p) ){
    case mxDOUBLE_CLASS:    return(   (double *)          mxGetData(p) );
    /* en estos de abajo se queja, claro, pues la salida es real*/
    /* case mxSINGLE_CLASS:    return(   (float  *)          mxGetData(p) ); */
    /* case mxLOGICAL_CLASS:   return(   (unsigned char *)   mxGetData(p) ); */
    /* case mxCHAR_CLASS:      return(   (char *)            mxGetData(p) ); */
    /* case mxUINT8_CLASS:     return(   (unsigned char *)   mxGetData(p) ); */
    /* case mxINT8_CLASS:      return(   (char *)            mxGetData(p) ); */
    /* case mxUINT16_CLASS:    return(   (unsigned short *)  mxGetData(p) ); */
    /* case mxINT16_CLASS:     return(   (short *)           mxGetData(p) ); */
    /* case mxUINT32_CLASS:    return(   (unsigned int *)    mxGetData(p) ); */
    /* case mxINT32_CLASS:     return(   (int *)             mxGetData(p) ); */
    /* case mxUINT64_CLASS:    return(   (unsigned long *)   mxGetData(p) ); */
    /* case mxINT64_CLASS:     return(   (long *)            mxGetData(p) ); */
  }
  /* return(0); */
}

#endif
