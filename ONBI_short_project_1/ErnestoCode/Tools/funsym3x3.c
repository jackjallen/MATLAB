/*

 xx = rand( [20 5 3 2 40 3 2 2] );
 [yy1,yy2] = funsym3x3( { xx , 3,6} , 'xtx' , 'exp' , 'sqrtm' );
 
*/

#include "myMEX.h"
#include "math.h" 

#if !defined( memcpy )
  #include "string.h"
#endif

#define N  3
#define NN 9


#define MAX9(a0,a1,a2,a3,a4,a5,a6,a7,a8)  MAX(MAX(MAX((a0),(a1)),MAX((a2),(a3))),MAX(MAX((a4),(a5)),MAX((a6),MAX((a7),(a8)))))

#define TWO_THREE             1.2599210498948732
#define TWO_SIX               1.5874010519681996
#define SQRT3                 1.7320508075688772
#define PI                    3.1415926535897931
#define EPS                   1.0e-25



#define Mtimes3x3( o , a , b)     o[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];  \
                                  o[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];  \
                                  o[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];  \
                                  o[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];  \
                                  o[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];  \
                                  o[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];  \
                                  o[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];  \
                                  o[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];  \
                                  o[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8]

#define FillOk    Ok[0] = V[0]*V[0]*fD[0] + V[3]*V[3]*fD[1] + V[6]*V[6]*fD[2];   \
                  Ok[1] = V[1]*V[0]*fD[0] + V[4]*V[3]*fD[1] + V[7]*V[6]*fD[2];   \
                  Ok[2] = V[2]*V[0]*fD[0] + V[5]*V[3]*fD[1] + V[8]*V[6]*fD[2];   \
                  Ok[3] = V[0]*V[1]*fD[0] + V[3]*V[4]*fD[1] + V[6]*V[7]*fD[2];   \
                  Ok[4] = V[1]*V[1]*fD[0] + V[4]*V[4]*fD[1] + V[7]*V[7]*fD[2];   \
                  Ok[5] = V[2]*V[1]*fD[0] + V[5]*V[4]*fD[1] + V[8]*V[7]*fD[2];   \
                  Ok[6] = V[0]*V[2]*fD[0] + V[3]*V[5]*fD[1] + V[6]*V[8]*fD[2];   \
                  Ok[7] = V[1]*V[2]*fD[0] + V[4]*V[5]*fD[1] + V[7]*V[8]*fD[2];   \
                  Ok[8] = V[2]*V[2]*fD[0] + V[5]*V[5]*fD[1] + V[8]*V[8]*fD[2];


#define ToOutput3x3     if( K1 == 1 && K2 == 1 ){                                                   \
                          memcpy( O + k3*NN , Ok , NN*sizeof( real ) );                             \
                        } else {                                                                    \
                          for( jj = 0 ; jj < N ; jj++ ){                                            \
                            for( ii = 0 ; ii < N ; ii++ ){                                          \
                              O[  k1 + K1*( ii + N*( k2 + K2*( jj + N*k3 ))) ] = Ok[ ii + N*jj ];   \
                            }                                                                       \
                          }                                                                         \
                        }                     

#define ToOutput1x1     if( K1 == 1 && K2 == 1 ){                                                   \
                          O[k3] = Ok[0];                                                            \
                        } else {                                                                    \
                          O[  k1 + K1*( k2 + K2* k3 ) ] = Ok[0];                                    \
                        }

#define ToOutput3x1     if( K1 == 1 && K2 == 1 ){                                                   \
                          memcpy( O + k3*N , Ok , N*sizeof( real ) );                               \
                        } else {                                                                    \
                          for( ii = 0 ; ii < N ; ii++ ){                                            \
                            O[  k1 + K1*( ii + N*( k2 + K2* k3 )) ] = Ok[ ii ];                     \
                          }                                                                         \
                        }

          
real eval[8];
real evec[9];

void EigenValuesSym3x3( real *M , real *D );
void EigenVectorsSym3x3( real *M , real *D , real *V );
void GetOrder(real *D, int *P);
void SortEigenValues(real *D, int *P , real *DS );
void SortEigenVectors(real *V, int *P , real *VS );
void JacobiDecomposition( real *M , real *D, real *V );
// void EigenSystemSym3x3( real *M , real *D , real *V );
void MatrixExp3x3( real *A , real *E );


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  real     *X, XX[NN], Xk[NN], XXr[NN], D[N], DS[N], fD[N], V[NN], VS[NN], *O, Ok[NN];
  long     k1, K1, k2, K2, k3, K3, ii, jj;
  int      ORDER[N];
  char     STR[100];

  enum     symmetrize_ { NONE , X_Xt , Xt_X , PLUS };
  enum     symmetrize_ symmetrize;  
  
  enum     outs_types { ID , INVERSE , EIGVALS , EIGVECS , EIGVALSDIAG , MAXDIFF , SQRT , LOG , EXP , evalFUN , DET , TRACE };
  enum     outs_types outs[50];
  int      ndims, Odims[50];
  real     det;
  
  int      oid;
  int      D_computed, V_computed, D_sorted, V_sorted, ORDER_computed;
  int      dim1 , dim2;
  mxArray  *mxX;
  
  if( nlhs == 0 ){ nlhs = 1; }
  if( nlhs > nrhs-2 ){
    myErrMsgTxt("There are no enough inputs.\n");
  }

  if( mxIsCell( prhs[0] ) ){
    mxX = mxGetCell( prhs[0] , 0 );
    dim1 = myGetValue( mxGetCell( prhs[0] , 1 ) );
    dim2 = myGetValue( mxGetCell( prhs[0] , 2 ) );
  } else {
    mxX = prhs[0];
    dim1 = 1;
    dim2 = 2;
  }
  
  if( ( mySize( mxX , dim1-1 ) != N ) || ( mySize( mxX , dim1-1 ) != N ) ){
    myErrMsgTxt("size( X , %d ) and size( X , %d ) have to be equal to three.\n",dim1,dim2);
  }
  

  if( myIsEmpty( prhs[1] ) ){
    symmetrize = NONE;
  } else if( mxIsChar(prhs[1]) ){
    mxGetString( prhs[1], STR, 100 );
    
    if(        !myStrcmpi(STR,"+") ) {
      symmetrize = PLUS;
    } else if( !myStrcmpi(STR,"xt*x") || !myStrcmpi(STR,"xtx") ) {
      symmetrize = Xt_X;
    } else if( !myStrcmpi(STR,"x*xt") || !myStrcmpi(STR,"xxt") ) {
      symmetrize = X_Xt;
    } else {
      myErrMsgTxt("Second argument expected: 'xt*x' , 'x*xt' , '+' , or [].\n");
    }

  } else {
    myErrMsgTxt("Second argument expected: 'xt*x' , 'x*xt' , '+' , or [].\n");
  }
  

  for( oid = 0 ; oid < nlhs ; oid++ ){
    if( ! mxIsChar(prhs[oid+2]) ){
      myErrMsgTxt("Valids arguments are: 'id' 'eigval' 'eigvec' 'log' 'sqrt' 'exp' 'error'.\n");
    }
    mxGetString( prhs[oid+2], STR, 100 );
    
    //{ ID , EIGVALS , EIGVECS , EIGVALSDIAG , MAXDIFF , SQRT , LOG , EXP , evalFUN };
    if(        !myStrcmpi(STR,"id") ) {
      outs[oid] = ID;
    } else if( !myStrcmpi(STR,"inv") || !myStrcmpi(STR,"inverse") ) {
      outs[oid] = INVERSE;
    } else if( !myStrcmpi(STR,"det") ) {
      outs[oid] = DET;
    } else if( !myStrcmpi(STR,"trace") ) {
      outs[oid] = TRACE;
    } else if( !myStrcmpi(STR,"evec") || !myStrcmpi(STR,"v") || !myStrcmpi(STR,"eigenvec") || !myStrcmpi(STR,"eigvec") || !myStrcmpi(STR,"eigenvectors") ) {
      outs[oid] = EIGVECS;
    } else if( !myStrcmpi(STR,"eval") || !myStrcmpi(STR,"d") || !myStrcmpi(STR,"eigenval") || !myStrcmpi(STR,"eigval") || !myStrcmpi(STR,"eigenvalues") ) {
      outs[oid] = EIGVALS;
    } else if( !myStrcmpi(STR,"diag") || !myStrcmpi(STR,"dd") ) {
      outs[oid] = EIGVALSDIAG;
    } else if( !myStrcmpi(STR,"error") ) {
      outs[oid] = MAXDIFF;
    } else if( !myStrcmpi(STR,"sqrt") || !myStrcmpi(STR,"sqrtm") ) {
      outs[oid] = SQRT;
    } else if( !myStrcmpi(STR,"log") || !myStrcmpi(STR,"logm") ) {
      outs[oid] = LOG;
    } else if( !myStrcmpi(STR,"exp") || !myStrcmpi(STR,"expm") ) {
      outs[oid] = EXP;
    } else {
      myErrMsgTxt("Valids arguments are: 'id' 'inv' 'det' 'trace' 'eigval' 'eigvec' 'log' 'sqrt' 'exp' 'error'.\n");
    }
  }
  
  
  X = mxGetPr( mxX );
  ndims = myGetSizes( mxX , Odims );
  
  for( oid = 0 ; oid < nlhs ; oid++ ){
    switch( outs[oid] ){
      case ID:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
        break;
      case INVERSE:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
        break;
      case DET:
        Odims[dim1-1] = 1;
        Odims[dim2-1] = 1;
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        Odims[dim1-1] = N;
        Odims[dim2-1] = N;
        break;
      case TRACE:
        Odims[dim1-1] = 1;
        Odims[dim2-1] = 1;
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        Odims[dim1-1] = N;
        Odims[dim2-1] = N;
        break;
      case EIGVALS:
        Odims[dim2-1] = 1;
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
        Odims[dim2-1] = N;
        break;
      case EIGVECS:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
      case EIGVALSDIAG:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
      case MAXDIFF:
        Odims[dim1-1] = 1;
        Odims[dim2-1] = 1;
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        Odims[dim1-1] = N;
        Odims[dim2-1] = N;
        break;
      case EXP:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
      case LOG:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
      case SQRT:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
      case evalFUN:
        plhs[oid] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );        
        break;
    }
  }
  
  K1 = 1;
  for( k1 = 0    ; k1 < dim1-1 ; k1++ ){  K1 *= Odims[k1]; }
  K2 = 1;
  for( k2 = dim1 ; k2 < dim2-1 ; k2++ ){  K2 *= Odims[k2]; }
  K3 = 1;
  for( k3 = dim2 ; k3 < ndims  ; k3++ ){  K3 *= Odims[k3]; }


  for( k3 = 0 ; k3 < K3 ; k3++ ){
  for( k2 = 0 ; k2 < K2 ; k2++ ){
  for( k1 = 0 ; k1 < K1 ; k1++ ){
    if( K1 == 1 && K2 == 1 ){
      memcpy( Xk , X + k3*NN , NN*sizeof( real ) );
    } else {
      for( jj = 0 ; jj < N ; jj++ ){ for( ii = 0 ; ii < N ; ii++ ){
        Xk[ ii + N*jj ] = X[  k1 + K1*( ii + N*( k2 + K2*( jj + N*k3 ))) ];
      } }
    }
    
    
    switch( symmetrize ){
      case NONE:
        memcpy( XX , Xk , NN*sizeof( real ) );
        break;

      case X_Xt:
        XX[0] = Xk[0]*Xk[0] + Xk[3]*Xk[3] + Xk[6]*Xk[6];
        XX[1] = Xk[1]*Xk[0] + Xk[4]*Xk[3] + Xk[7]*Xk[6];
        XX[2] = Xk[2]*Xk[0] + Xk[5]*Xk[3] + Xk[8]*Xk[6];

        XX[3] = Xk[0]*Xk[1] + Xk[3]*Xk[4] + Xk[6]*Xk[7];
        XX[4] = Xk[1]*Xk[1] + Xk[4]*Xk[4] + Xk[7]*Xk[7];
        XX[5] = Xk[2]*Xk[1] + Xk[5]*Xk[4] + Xk[8]*Xk[7];

        XX[6] = Xk[0]*Xk[2] + Xk[3]*Xk[5] + Xk[6]*Xk[8];
        XX[7] = Xk[1]*Xk[2] + Xk[4]*Xk[5] + Xk[7]*Xk[8];
        XX[8] = Xk[2]*Xk[2] + Xk[5]*Xk[5] + Xk[8]*Xk[8];
        break;

      case Xt_X:
        XX[0] = Xk[0]*Xk[0] + Xk[1]*Xk[1] + Xk[2]*Xk[2];
        XX[1] = Xk[3]*Xk[0] + Xk[4]*Xk[1] + Xk[5]*Xk[2];
        XX[2] = Xk[6]*Xk[0] + Xk[7]*Xk[1] + Xk[8]*Xk[2];

        XX[3] = Xk[0]*Xk[3] + Xk[1]*Xk[4] + Xk[2]*Xk[5];
        XX[4] = Xk[3]*Xk[3] + Xk[4]*Xk[4] + Xk[5]*Xk[5];
        XX[5] = Xk[6]*Xk[3] + Xk[7]*Xk[4] + Xk[8]*Xk[5];

        XX[6] = Xk[0]*Xk[6] + Xk[1]*Xk[7] + Xk[2]*Xk[8];
        XX[7] = Xk[3]*Xk[6] + Xk[4]*Xk[7] + Xk[5]*Xk[8];
        XX[8] = Xk[6]*Xk[6] + Xk[7]*Xk[7] + Xk[8]*Xk[8];
        break;

      case PLUS:
        XX[0]         = Xk[0];
        XX[1] = XX[3] = ( Xk[1] + Xk[3] )/2.0;
        XX[2] = XX[6] = ( Xk[2] + Xk[6] )/2.0;
        XX[4]         = Xk[4];
        XX[5] = XX[7] = ( Xk[5] + Xk[7] )/2.0;
        XX[8]         = Xk[8];
        break;
    }
    
    D_computed      = 0;
    V_computed      = 0;
    ORDER_computed  = 0;
    D_sorted        = 0;
    V_sorted        = 0;
    
    for( oid = 0 ; oid < nlhs ; oid++ ){
      switch( outs[oid] ){
        case ID:
          memcpy( Ok , XX , NN*sizeof( real ) );
          
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;
          
          
        case INVERSE:
          det =  XX[0]*( XX[4]*XX[8] - XX[5]*XX[7] )
               - XX[1]*( XX[3]*XX[8] - XX[5]*XX[6] )
               + XX[2]*( XX[3]*XX[7] - XX[4]*XX[6] );
          det = 1.0/det;
          
          Ok[0] = ( XX[4]*XX[8] - XX[5]*XX[7] )*det;
          Ok[1] = ( XX[2]*XX[7] - XX[1]*XX[8] )*det;
          Ok[2] = ( XX[1]*XX[5] - XX[2]*XX[4] )*det;
          Ok[3] = ( XX[5]*XX[6] - XX[3]*XX[8] )*det;
          Ok[4] = ( XX[0]*XX[8] - XX[2]*XX[6] )*det;
          Ok[5] = ( XX[2]*XX[3] - XX[0]*XX[5] )*det;
          Ok[6] = ( XX[3]*XX[7] - XX[4]*XX[6] )*det;
          Ok[7] = ( XX[1]*XX[6] - XX[0]*XX[7] )*det;
          Ok[8] = ( XX[0]*XX[4] - XX[1]*XX[3] )*det;
          
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;
          

        case DET:
          Ok[0] =   XX[0]*( XX[4]*XX[8] - XX[5]*XX[7] )
                  - XX[1]*( XX[3]*XX[8] - XX[5]*XX[6] )
                  + XX[2]*( XX[3]*XX[7] - XX[4]*XX[6] );
          
          O = mxGetPr( plhs[oid] );
          ToOutput1x1;
          break;

          
        case TRACE:
          Ok[0] =  XX[0] + XX[4] + XX[8];
          
          O = mxGetPr( plhs[oid] );
          ToOutput1x1;
          break;

          
        case EIGVALS:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !D_sorted       ){ SortEigenValues( D, ORDER , DS );  D_sorted = 1;        }
          
          memcpy( Ok , DS , N*sizeof( real ) );
          
          O = mxGetPr( plhs[oid] );
          ToOutput3x1;
          break;


        case EIGVALSDIAG:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !D_sorted       ){ SortEigenValues( D, ORDER , DS );  D_sorted = 1;        }

          Ok[0] = DS[0];
          Ok[4] = DS[1];
          Ok[8] = DS[2];
          Ok[1] = Ok[2] = Ok[3] = Ok[5] = Ok[6] = Ok[7] = 0;

          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;

          
        case EIGVECS:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !V_sorted       ){ SortEigenVectors( V, ORDER , VS);  V_sorted = 1;        }

          memcpy( Ok , VS , NN*sizeof( real ) );
          
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;

          
        case MAXDIFF:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }

          Ok[0] = V[0]*V[0]*D[0] + V[3]*V[3]*D[1] + V[6]*V[6]*D[2];
          Ok[1] = V[1]*V[0]*D[0] + V[4]*V[3]*D[1] + V[7]*V[6]*D[2];
          Ok[2] = V[2]*V[0]*D[0] + V[5]*V[3]*D[1] + V[8]*V[6]*D[2];

          Ok[3] = V[0]*V[1]*D[0] + V[3]*V[4]*D[1] + V[6]*V[7]*D[2];
          Ok[4] = V[1]*V[1]*D[0] + V[4]*V[4]*D[1] + V[7]*V[7]*D[2];
          Ok[5] = V[2]*V[1]*D[0] + V[5]*V[4]*D[1] + V[8]*V[7]*D[2];

          Ok[6] = V[0]*V[2]*D[0] + V[3]*V[5]*D[1] + V[6]*V[8]*D[2];
          Ok[7] = V[1]*V[2]*D[0] + V[4]*V[5]*D[1] + V[7]*V[8]*D[2];
          Ok[8] = V[2]*V[2]*D[0] + V[5]*V[5]*D[1] + V[8]*V[8]*D[2];

          Ok[0] =  MAX9( fabs( Ok[0] - XX[0] ) ,
                         fabs( Ok[1] - XX[1] ) ,
                         fabs( Ok[2] - XX[2] ) ,
                         fabs( Ok[3] - XX[3] ) ,
                         fabs( Ok[4] - XX[4] ) ,
                         fabs( Ok[5] - XX[5] ) ,
                         fabs( Ok[6] - XX[6] ) ,
                         fabs( Ok[7] - XX[7] ) ,
                         fabs( Ok[8] - XX[8] ) );

          O = mxGetPr( plhs[oid] );
          ToOutput1x1;
          break;

          
        case SQRT:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }

          fD[0] = sqrt( D[0] );
          fD[1] = sqrt( D[1] );
          fD[2] = sqrt( D[2] );
          
          FillOk;
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;

          
        case LOG:
          if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
          if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }

          fD[0] = log( D[0] );
          fD[1] = log( D[1] );
          fD[2] = log( D[2] );
          
          FillOk;
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;


        case EXP:
          
          MatrixExp3x3( XX , Ok );
          
//           if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
//           if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
// 
//           fD[0] = exp( D[0] );
//           fD[1] = exp( D[1] );
//           fD[2] = exp( D[2] );
//           
//           FillOk;
          O = mxGetPr( plhs[oid] );
          ToOutput3x3;
          break;

          
          
//         case evalFUN:
//           O = mxGetPr( plhs[oid] ) + k*NN;
//           if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
//           if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
// 
//           fD[0] = exp( D[0] );
//           fD[1] = exp( D[1] );
//           fD[2] = exp( D[2] );
// 
//           mexCallMATLAB(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[], const char *name)
//           
//           FillO
//           break;
      
      }
    }
    
  }}}
  
  EXIT: myFreeALLOCATES();
}


#define v0      eval[0]
#define v1      eval[1]
#define v2      eval[2]
#define tmp     eval[3]
#define theta   eval[4]
#define rho     eval[5]
#define Re      eval[6]
#define Im      eval[7]
void EigenValuesSym3x3( real *M , real *D ) {
  v0 = M[0]+M[4]+M[8];
  v1 = -v0*v0 - 3.0*(M[1]*M[1] + M[2]*M[2] - M[0]*M[4] + M[5]*M[5] - M[0]*M[8] - M[4]*M[8]);
  v2 = 9*(        -M[0]*M[1]*M[1] - M[0]*M[2]*M[2] + M[0]*M[0]*M[4] - 
                   M[1]*M[1]*M[4] + M[0]*M[4]*M[4] - M[4]*M[5]*M[5] + 
                   M[0]*M[0]*M[8] - M[2]*M[2]*M[8] + M[4]*M[4]*M[8] - 
                   M[5]*M[5]*M[8] + M[0]*M[8]*M[8] + M[4]*M[8]*M[8] + 
            2.0* ( M[2]*M[2]*M[4] + M[0]*M[5]*M[5] + M[1]*M[1]*M[8] - 
                                                 3.0*M[1]*M[2]*M[5] ) ) - 2.0*v0*v0*v0;

  tmp = 4.0*v1*v1*v1 + v2*v2;
  if( tmp < 0 ){
    tmp   = sqrt( - tmp );
    theta = atan( tmp / v2 );
    if( theta < 0 ){
      theta =  1.0/3.0  * theta + PI;
    } else {
      theta =  1.0/3.0  * theta;
    }
    rho   = pow( sqrt(tmp*tmp + v2*v2) ,  1.0/3.0  );

    Re = rho * cos(theta);
    Im = rho * sin(theta);		
  } else {
    Re = pow( ( v2 + sqrt(tmp) ),   1.0/3.0  );
    Im = 0;
  }

  tmp  =  1.0/3.0  * ( TWO_SIX - 2.0*TWO_THREE * v1 /(Re*Re +Im*Im) );

  D[0] =  1.0/3.0 *v0 - 0.50*tmp *  Re;
  D[1] =  1.0/3.0 *v0 + 0.25*tmp * ( Re - SQRT3 * Im );
  D[2] =  1.0/3.0 *v0 + 0.25*tmp * ( Re + SQRT3 * Im );
  
  return;
}  

#define ev2x    evec[0]
#define ev1x    evec[1]
#define ev0x    evec[2]
#define ev2y    evec[3]
#define ev1y    evec[4]
#define ev0y    evec[5]
#define n2      evec[6]
#define n1      evec[7]
#define n0      evec[8]
void EigenVectorsSym3x3( real *M , real *D , real *V ) {
  
  if( fabs( D[0] - D[1] ) < 1e-5 || fabs( D[0] - D[2] ) < 1e-5 || fabs( D[1] - D[2] ) < 1e-5 ){
//     mexPrintf("calling jacobi\n\n");
    JacobiDecomposition( M , D , V );
    return;
  }
  
  if( fabs(M[1]) < EPS && fabs(M[2]) < EPS && fabs(M[5]) < EPS) {

    if( M[0] == MAX3( M[0], M[4], M[8] ) ) {
      if(M[4] > M[8]){
        V[0] = 0.0;		V[3] = 0.0;		V[6] = 1.0;
        V[1] = 0.0;		V[4] = 1.0;		V[7] = 0.0;
        V[2] = 1.0;		V[5] = 0.0;		V[8] = 0.0;
//         mexPrintf("case 1\n");
        return;
      }

        V[0] = 0.0;		V[3] = 0.0;		V[6] = 1.0;
        V[1] = 1.0;		V[4] = 0.0;		V[7] = 0.0;
        V[2] = 0.0;		V[5] = 1.0;		V[8] = 0.0;
//         mexPrintf("case 2\n");
        return;
    } 
    
    if( M[4] == MAX3( M[0], M[4], M[8] ) ) {
      if(M[0] > M[8]) {
        V[0] = 0.0;		V[3] = 1.0;		V[6] = 0.0;
        V[1] = 0.0;		V[4] = 0.0;		V[7] = 1.0;
        V[2] = 1.0;		V[5] = 0.0;		V[8] = 0.0;
//         mexPrintf("case 3\n");
        return;
      }
      
        V[0] = 1.0;		V[3] = 0.0;		V[6] = 0.0;
        V[1] = 0.0;		V[4] = 0.0;		V[7] = 1.0;
        V[2] = 0.0;		V[5] = 1.0;		V[8] = 0.0;
//         mexPrintf("case 4\n");
        return;
    } 
    
    if(M[0] > M[4]){
      V[0] = 0.0;		V[3] = 1.0;		V[6] = 0.0;
      V[1] = 1.0;		V[4] = 0.0;		V[7] = 0.0;
      V[2] = 0.0;		V[5] = 0.0;		V[8] = 1.0;
//       mexPrintf("case 6\n");
      return;
    }

      V[0] = 1.0;		V[3] = 0.0;		V[6] = 0.0;
      V[1] = 0.0;		V[4] = 1.0;		V[7] = 0.0;
      V[2] = 0.0;		V[5] = 0.0;		V[8] = 1.0;
//       mexPrintf("case 7\n");
      return;
  } 
  
  if( fabs(M[1]) > EPS && fabs(M[2]) < EPS && fabs(M[5]) < EPS ){
    ev2x = (D[2]-M[4])/M[1];
    ev1x = (D[1]-M[4])/M[1];
    ev0x = (D[0]-M[4])/M[1];

    n2 = 1.0/sqrt(ev2x*ev2x + 1.0);
    n1 = 1.0/sqrt(ev1x*ev1x + 1.0);
    n0 = 1.0/sqrt(ev0x*ev0x + 1.0);

    if( fabs(M[8]-D[2]) < 1e-10 ){
      V[0] = ev0x*n0;		V[3] = ev1x*n1; 	V[6] = 0.0;
      V[1] = n0;        V[4] = n1;        V[7] = 0.0;
      V[2] = 0.0;       V[5] = 0.0;       V[8] = 1.0;
//       mexPrintf("case 8\n");
      return;
    } 

    if( fabs(M[8]-D[1]) < EPS && M[8] > D[0]){
      V[0] = ev0x*n0;		V[3] = 0.0;       V[6] = ev2x*n2;
      V[1] = n0;        V[4] = 0.0;       V[7] = n2;
      V[2] = 0.0;       V[5] = 1.0;       V[8] = 0.0;
//       mexPrintf("case 9\n");
      return;
    }

      V[0] = 0.0;       V[3] = ev1x*n1; 	V[6] = ev2x*n2;
      V[1] = 0.0;       V[4] = n1;        V[7] = n2;
      V[2] = 1.0;       V[5] = 0.0;       V[8] = 0.0;
//       mexPrintf("case 10\n");
      return;
  }
  
  if( fabs(M[2]) > EPS && fabs(M[1]) < EPS && fabs(M[5]) < EPS ){
    ev2x = (D[2]-M[8])/M[2];
    ev1x = (D[1]-M[8])/M[2];
    ev0x = (D[0]-M[8])/M[2];

    n2 = 1.0/sqrt(ev2x*ev2x + 1.0);
    n1 = 1.0/sqrt(ev1x*ev1x + 1.0);
    n0 = 1.0/sqrt(ev0x*ev0x + 1.0);

    if( fabs(M[4]-D[2]) < EPS ){
      V[0] = ev0x*n0;		V[3] = ev1x*n1; 	V[6] = 0.0;
      V[1] = 0.0;       V[4] = 0.0;       V[7] = 1.0;
      V[2] = n0;        V[5] = n1;        V[8] = 0.0;
//       mexPrintf("case 11\n");
      return;
    } 
    
    if( fabs(M[4]-D[1]) < EPS && M[4] > D[0]){
      V[0] = ev0x*n0;		V[3] = 0.0;       V[6] = ev2x*n2;
      V[1] = 0.0;       V[4] = 1.0;       V[7] = 0.0;
      V[2] = n0;        V[5] = 0.0;       V[8] = n2;
//       mexPrintf("case 12\n");
      return;
    } 
    
      V[0] = 0.0;       V[3] = ev1x*n1;		V[6] = ev2x*n2;
      V[1] = 1.0;       V[4] = 0.0;       V[7] = 0.0;
      V[2] = 0.0;       V[5] = n1;        V[8] = n2;
//       mexPrintf("case 13\n");
      return;
  }
  
  if( fabs(M[5]) > EPS && fabs(M[1]) < EPS && fabs(M[2]) < EPS ){
    ev2y = (D[2]-M[8])/M[5];
    ev1y = (D[1]-M[8])/M[5];
    ev0y = (D[0]-M[8])/M[5];

    n2 = 1.0/sqrt(ev2y*ev2y + 1.0);
    n1 = 1.0/sqrt(ev1y*ev1y + 1.0);
    n0 = 1.0/sqrt(ev0y*ev0y + 1.0);

    if( fabs(M[0]-D[2]) < EPS ){
      V[0] = 0.0;       V[3] = 0.0;       V[6] = 1.0;
      V[1] = ev0y*n0;		V[4] = ev1y*n1;		V[7] = 0.0;
      V[2] = n0;        V[5] = n1;        V[8] = 0.0;
//       mexPrintf("case 14\n");
      return;
    } 
    
    if( fabs(M[0]-D[1]) < EPS && M[0] > D[0]){
      V[0] = 0.0;       V[3] = 1.0;       V[6] = 0.0;
      V[1] = ev0y*n0;		V[4] = 0.0;       V[7] = ev2y*n2;
      V[2] = n0;        V[5] = 0.0;       V[8] = n2;
//       mexPrintf("case 15\n");
      return;
    } 
    
      V[0] = 1.0;       V[3] = 0.0;       V[6] = 0.0;
      V[1] = 0.0;       V[4] = ev1y*n1;		V[7] = ev2y*n2;
      V[2] = 0.0;       V[5] = n1;        V[8] = n2;
//       mexPrintf("case 16\n");
      return;
  }
  
  n2 = 1.0/ (M[1]*M[2] + (D[2]-M[0])*M[5]);
  n1 = 1.0/ (M[1]*M[2] + (D[1]-M[0])*M[5]);
  n0 = 1.0/ (M[1]*M[2] + (D[0]-M[0])*M[5]);

  ev2x = (D[2]*M[1] + M[2]*M[5] - M[1]*M[8]) * n2;
  ev2y = (D[2]*D[2] - M[2]*M[2] + M[0]*M[8] - D[2]*(M[0] + M[8])) * n2;

  ev1x = (D[1]*M[1] + M[2]*M[5] - M[1]*M[8]) * n1;
  ev1y = (D[1]*D[1] - M[2]*M[2] + M[0]*M[8] - D[1]*(M[0] + M[8])) * n1;

  ev0x = (D[0]*M[1] + M[2]*M[5] - M[1]*M[8]) * n0;
  ev0y = (D[0]*D[0] - M[2]*M[2] + M[0]*M[8] - D[0]*(M[0] + M[8])) * n0;

  n2 = 1.0/sqrt(ev2x*ev2x + ev2y*ev2y + 1.0);
  n1 = 1.0/sqrt(ev1x*ev1x + ev1y*ev1y + 1.0);
  n0 = 1.0/sqrt(ev0x*ev0x + ev0y*ev0y + 1.0);

  V[0] = ev0x*n0;		V[3] = ev1x*n1;		V[6] = ev2x*n2;
  V[1] = ev0y*n0;		V[4] = ev1y*n1;		V[7] = ev2y*n2;
  V[2] = n0;        V[5] = n1;        V[8] = n2;
//   mexPrintf("case 17\n");
  return;
}


void GetOrder(real *D, int *P) {
  int      smaller_id, i, k;
  real     smaller;

  for( i = 0 ; i < N ; i++ ){ P[i] = i; }
  
  for( i = 0 ; i < N ; i++ ){
    if( i == N-1 ){ continue; }
    smaller = D[ P[i] ];
    smaller_id = i;
    for( k = i + 1 ; k < N ; k ++ ){
      if( D[ P[k] ] < smaller ){
        smaller    = D[ P[k] ];
        smaller_id = k;
      }
    }
    k = P[i];
    P[i] = P[smaller_id];
    P[smaller_id] = k;
  }
}

void SortEigenValues( real *D, int *P , real *DS){
//   int i;
//   for( i = 0 ; i < N ; i++ ){
//     DS[i] = D[ P[i] ];
//   }
  
  DS[0] = D[ P[0] ];
  DS[1] = D[ P[1] ];
  DS[2] = D[ P[2] ];
  
}

void SortEigenVectors(real *V, int *P , real *VS){
//   int i;
//   for( i = 0 ; i < N ; i++ ){
//     memcpy( VS + i*N , V + P[i]*N , N*sizeof( real ) );
//   }
  
  memcpy( VS     , V + P[0]*3 , N*sizeof( real ) );  
  memcpy( VS + 3 , V + P[1]*3 , N*sizeof( real ) );  
  memcpy( VS + 6 , V + P[2]*3 , N*sizeof( real ) );  
  
}

void JacobiDecomposition( real *M , real *D, real *V ) {
  real   z[N], b[N], A[NN], Dnew[N];
	int    k,j,i,iter, ORDERorig[N], ORDERnew[N];
	real   theta,tau,t,s,h,g,c;
  real   x;

  
//   mexPrintf("%g   %g   %g \n", D[0] , D[1] , D[2] );
  
  memcpy( A , M , NN*sizeof( real ) );
  
  for( i=0 ; i<NN ; i++ ){
    V[i] = 0.0;
  }
	for( i=0 ; i<N ; i++ ){
    z[i]       = 0;
    V[ i+N*i ] = 1.0;
		b[i]       = A[ i+N*i ];
	}

//   for( i = 0 ; i < NN ; i++ ){ mexPrintf("A:  %g\n",A[i]); }
//   for( i = 0 ; i < NN ; i++ ){ mexPrintf("V:  %g\n",V[i]); }
//   for( i = 0 ; i < N  ; i++ ){ mexPrintf("b:  %g\n",b[i]); }
//   for( i = 0 ; i < N  ; i++ ){ mexPrintf("z:  %g\n",z[i]); }
  
  iter= 1;
  while( iter ) {
    iter = 0;
		for( i=0 ; i < N-1 ; i++ ) {
			for( j=i+1 ; j < N ; j++ ) {
        if( fabs( A[i+N*j] ) > EPS ) {
          iter= 1;
          h     = A[j+N*j]-A[i+N*i];
          theta = 0.5*h/A[i+N*j];
          if( theta >= 0.0 ){
            t =  1.0/(fabs(theta)+sqrt(1.0+theta*theta));
          } else {
            t = -1.0/(fabs(theta)+sqrt(1.0+theta*theta));
          }
          c = 1.0/sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*A[i+N*j];
          z[i] -= h;
          z[j] += h;
          A[i+N*i] -= h;
          A[j+N*j] += h;
          A[i+N*j] =  0.0;
          for( k=0    ; k < i ; k++ ) {
            g = A[ k+N*i ];
            h = A[ k+N*j ];
            A[ k+N*i ] = g-s*(h+g*tau);
            A[ k+N*j ] = h+s*(g-h*tau);
          }
          for( k = i+1  ; k<j ; k++ ) {
            g = A[ i+N*k ];
            h = A[ k+N*j ];
            A[ i+N*k ] = g-s*(h+g*tau);
            A[ k+N*j ] = h+s*(g-h*tau);
          }
          for( k=j+1  ; k<N   ; k++ ) {			
            g = A[ i+N*k ];
            h = A[ j+N*k ];
            A[ i+N*k ] = g-s*(h+g*tau);
            A[ j+N*k ] = h+s*(g-h*tau);
          }
          for( k=0    ; k<N   ; k++ ) {		
            g = V[ k+N*i ];
            h = V[ k+N*j ];
            V[ k+N*i ] = g-s*(h+g*tau);
            V[ k+N*j ] = h+s*(g-h*tau);
          }
        }
			}
		}
		for( i=0 ; i<N ; i++ ){
			z[i]  = 0.0;
		}
	}


  Dnew[0] = A[0];
  Dnew[1] = A[4];
  Dnew[2] = A[8];
  

  GetOrder( Dnew , ORDERnew  );
  GetOrder( D    , ORDERorig );
  
//   mexPrintf("%g   %g   %g \n", D[0] , D[1] , D[2] );
//   mexPrintf("ORDERorig: %d   %d   %d \n\n", ORDERorig[0] , ORDERorig[1] , ORDERorig[2] );

  
//   mexPrintf("%g   %g   %g \n", Dnew[0] , Dnew[1] , Dnew[2] );
//   mexPrintf("ORDERnew : %d   %d   %d \n\n", ORDERnew[0] , ORDERnew[1] , ORDERnew[2] );
  
  
  for( i = 0 ; i < N-1 ; i++ ){
//     mexPrintf( "ORDERorig: %d  ... ORDERnew:  %d   ", ORDERorig[i] , ORDERnew[i] );
    
    for( j = 0 ; j < N ; j++ ){
      if( ORDERnew[j] == ORDERorig[i] ){
        break;
      }
    }

//     mexPrintf("i: %d   -> j: %d\n",i,j );

    if( i == j ){ continue; }
    
//     mexPrintf("changing cols: %d  -> %d\n", ORDERnew[i] , ORDERnew[j] );
    x = Dnew[ ORDERnew[i] ];
    Dnew[ ORDERnew[i] ] = Dnew[ ORDERnew[j] ];
    Dnew[ ORDERnew[j] ] = x;

    
    for( k = 0 ; k < N ; k++ ){
      x = V[ k + N*ORDERnew[i] ];
      V[ k + N*ORDERnew[i] ] = V[ k + N*ORDERnew[j] ];
      V[ k + N*ORDERnew[j] ] = x;
    }
    
    
    k = ORDERnew[i];
    ORDERnew[i] = ORDERnew[j];
    ORDERnew[j] = k;
  }
  
  
  D[0] = Dnew[0];
  D[1] = Dnew[1];
  D[2] = Dnew[2];
  
//   mexPrintf("\n%g   %g   %g \n", D[0] , D[1] , D[2] );
  
}



void MatrixExp3x3( real *A , real *E ){
  real    normA;
  real    A1[9] = { 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 };
  real    A2[9];
  real    A4[9];
  real    A6[9];
  real    A8[9];
  real    U[9] ;
  real    V[9] ;
  real    aux[9];
  int     i, s;
  real    scale = 2.0;
  real    det;
  
  normA = MAX3(  ABS( A[0] ) + ABS( A[1] ) + ABS( A[2] ) ,
                 ABS( A[3] ) + ABS( A[4] ) + ABS( A[5] ) ,
                 ABS( A[6] ) + ABS( A[7] ) + ABS( A[8] ) );

  
  if( normA <= 1.495585217958292e-2 ){

    Mtimes3x3( A2 , A , A );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A2[i] +  60.0*A1[i]; 
      V[i] = 12.0*A2[i] + 120.0*A1[i];
    }

  } else if( normA <= 0.2539398330063230 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A4[i] +  420.0*A2[i] + 15120.0*A1[i]; 
      V[i] = 30.0*A4[i] + 3360.0*A2[i] + 30240.0*A1[i];
    }

  } else if( normA <= 0.9504178996162932 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A6[i] +  1512.0*A4[i] +  277200.0*A2[i] +  8648640.0*A1[i]; 
      V[i] = 56.0*A6[i] + 25200.0*A4[i] + 1995840.0*A2[i] + 17297280.0*A1[i];
    }

  } else if( normA <= 2.097847961257068 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    Mtimes3x3( A8 , A6 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A8[i] +   3960.0*A6[i] +  2162160.0*A4[i] +  302702400.0*A2[i] +  8821612800.0*A1[i]; 
      V[i] = 90.0*A8[i] + 110880.0*A6[i] + 30270240.0*A4[i] + 2075673600.0*A2[i] + 17643225600.0*A1[i];
    }
    
  } else if( normA <= 5.371920351148152 ){
    
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    
    for(i=0;i<9;i++){ 
      U[i] =   1.0*A6[i] +  16380.0*A4[i] +   40840800.0*A2[i];
      V[i] = 182.0*A6[i] + 960960.0*A4[i] + 1323241920.0*A2[i];
    }
    
    Mtimes3x3( aux , A6 , U ); memcpy( U , aux , 9*sizeof( real ) );
    Mtimes3x3( aux , A6 , V ); memcpy( V , aux , 9*sizeof( real ) );
        
    for(i=0;i<9;i++){ 
      U[i] = U[i] +  33522128640.0*A6[i] +  10559470521600.0*A4[i] + 1187353796428800.0*A2[i] + 32382376266240000.0*A1[i];
      V[i] = V[i] + 670442572800.0*A6[i] + 129060195264000.0*A4[i] + 7771770303897600.0*A2[i] + 64764752532480000.0*A1[i];
    }
    
  } else {

    #define log2   0.69314718055994529
    s = (int) ( log( normA/5.371920351148152 )/log2  +  1 );
    for( i = 1; i<s ; i++ ){
      scale = scale*2.0;
    }
    scale = 1.0/scale;
    
    aux[0] = A[0]*scale;
    aux[1] = A[1]*scale;
    aux[2] = A[2]*scale;
    aux[3] = A[3]*scale;
    aux[4] = A[4]*scale;
    aux[5] = A[5]*scale;
    aux[6] = A[6]*scale;
    aux[7] = A[7]*scale;
    aux[8] = A[8]*scale;
    
    MatrixExp3x3( aux , E );
    
    for( i=0 ; i<s ; i++ ){
      Mtimes3x3( aux , E , E );
      memcpy( E , aux , 9*sizeof( real ) );
    }
    
    return;

  }
  
  Mtimes3x3( aux , A , U ); memcpy( U , aux , 9*sizeof( real ) );
  for(i=0;i<9;i++){ 
    A2[i] = V[i]-U[i];
    A4[i] = V[i]+U[i];
  }

  det =   A2[0]*( A2[4]*A2[8] - A2[5]*A2[7] )
        - A2[1]*( A2[3]*A2[8] - A2[5]*A2[6] )
        + A2[2]*( A2[3]*A2[7] - A2[4]*A2[6] );
  det = 1.0/det;
  
  aux[0] = ( A2[4]*A2[8] - A2[5]*A2[7] )*det;
  aux[1] = ( A2[2]*A2[7] - A2[1]*A2[8] )*det;
  aux[2] = ( A2[1]*A2[5] - A2[2]*A2[4] )*det;
  aux[3] = ( A2[5]*A2[6] - A2[3]*A2[8] )*det;
  aux[4] = ( A2[0]*A2[8] - A2[2]*A2[6] )*det;
  aux[5] = ( A2[2]*A2[3] - A2[0]*A2[5] )*det;
  aux[6] = ( A2[3]*A2[7] - A2[4]*A2[6] )*det;
  aux[7] = ( A2[1]*A2[6] - A2[0]*A2[7] )*det;
  aux[8] = ( A2[0]*A2[4] - A2[1]*A2[3] )*det;
  
  Mtimes3x3( E , aux , A4 );
  
  return;
}



























// void EigenSystemSym3x3( real *M , real *D , real *V ) {
//   EigenValuesSym3x3( M , D );
//   EigenVectorsSym3x3( M , D , V );
//   return;
// }


// void SortEigenValues(real *D, real *V) {
//   int      col_id, smaller_id, v_id, i;
//   real     x;
// 
//   for( col_id = 0 ; col_id < N-1 ; col_id++ ){
//     smaller_id  = col_id;
//     for( v_id = smaller_id + 1 ; v_id < N ; v_id++ ){
//       if( D[ v_id ] < D[ smaller_id ] ){
//         smaller_id = v_id;
//       }
//       if( smaller_id != col_id ){
//         x                     = D[ N*col_id     + i ];
//         D[ N*col_id     + i ] = D[ N*smaller_id + i ];
//         D[ N*smaller_id + i ] = x;
//         if( V != NULL ){
//           for( i = 0 ; i < N ; i++ ){
//             x                     = V[ N*col_id     + i ];
//             V[ N*col_id     + i ] = V[ N*smaller_id + i ];
//             V[ N*smaller_id + i ] = x;
//           }
//         }
//       }
//       if( V != NULL ){
//         if( V[ N*col_id ] < 0 ){
//           for( i = 0 ; i < N ; i++ ){
//             V[ N*col_id + i ] *= -1;
//           }      
//         }
//       }
//     }
//   }
//   if( V != NULL ){
//     col_id ++;
//     if( V[ N*col_id ] < 0 ){
//       for( i = 0 ; i < N ; i++ ){
//         V[ N*col_id + i ] *= -1;
//       }      
//     }
//   }
// }
