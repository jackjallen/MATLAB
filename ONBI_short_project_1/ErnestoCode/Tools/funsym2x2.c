#include "myMEX.h"
#include "math.h" 

#if !defined( memcpy )
  #include "string.h"
#endif

#define N  2
#define NN 4

#define MAX3(a,b,c)           MAX( (a) , MAX( (b) , (c) ) )
#define MAX4(a0,a1,a2,a3)     MAX( MAX( (a0),(a1) ) , MAX( (a2),(a3) ) )


#define FillOk    Ok[0] = V[0]*V[0]*fD[0] + V[3]*V[3]*fD[1] + V[6]*V[6]*fD[2];   \
                  Ok[1] = V[1]*V[0]*fD[0] + V[4]*V[3]*fD[1] + V[7]*V[6]*fD[2];   \
                  Ok[2] = V[2]*V[0]*fD[0] + V[5]*V[3]*fD[1] + V[8]*V[6]*fD[2];   \
                  Ok[3] = V[0]*V[1]*fD[0] + V[3]*V[4]*fD[1] + V[6]*V[7]*fD[2];   \
                  Ok[4] = V[1]*V[1]*fD[0] + V[4]*V[4]*fD[1] + V[7]*V[7]*fD[2];   \
                  Ok[5] = V[2]*V[1]*fD[0] + V[5]*V[4]*fD[1] + V[8]*V[7]*fD[2];   \
                  Ok[6] = V[0]*V[2]*fD[0] + V[3]*V[5]*fD[1] + V[6]*V[8]*fD[2];   \
                  Ok[7] = V[1]*V[2]*fD[0] + V[4]*V[5]*fD[1] + V[7]*V[8]*fD[2];   \
                  Ok[8] = V[2]*V[2]*fD[0] + V[5]*V[5]*fD[1] + V[8]*V[8]*fD[2];


#define ToOutput2x2     if( K1 == 1 && K2 == 1 ){                                                   \
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

#define ToOutput2x1     if( K1 == 1 && K2 == 1 ){                                                   \
                          memcpy( O + k3*N , Ok , N*sizeof( real ) );                               \
                        } else {                                                                    \
                          for( ii = 0 ; ii < N ; ii++ ){                                            \
                            O[  k1 + K1*( ii + N*( k2 + K2* k3 )) ] = Ok[ ii ];                     \
                          }                                                                         \
                        }

          
real eval[4];
real evec[4];

void EigenValues2x2( real *M , real *D );
void EigenVectors2x2( real *M , real *D , real *V );
void GetOrder(real *D, int *P);
void SortEigenValues(real *D, int *P , real *DS );
void SortEigenVectors(real *V, int *P , real *VS );
void JacobiDecomposition( real *M , real *D, real *V );


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
        XX[0] = Xk[0]*Xk[0] + Xk[2]*Xk[2];
        XX[1] = Xk[1]*Xk[0] + Xk[3]*Xk[2];

        XX[2] = Xk[0]*Xk[1] + Xk[2]*Xk[3];
        XX[3] = Xk[1]*Xk[1] + Xk[3]*Xk[3];
        break;

      case Xt_X:
        XX[0] = Xk[0]*Xk[0] + Xk[1]*Xk[1];
        XX[1] = Xk[2]*Xk[0] + Xk[3]*Xk[1];

        XX[2] = Xk[0]*Xk[2] + Xk[1]*Xk[3];
        XX[3] = Xk[2]*Xk[2] + Xk[3]*Xk[3];
        break;

      case PLUS:
        XX[0]         = Xk[0];
        XX[1] = XX[2] = ( Xk[1] + Xk[2] )/2.0;
        XX[3]         = Xk[3];
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
          ToOutput2x2;
          break;
          
          
        case INVERSE:
          det = XX[0]*XX[3] - XX[1]*XX[2];
          det = 1.0/det;
          
          Ok[0] =  XX[3] *det;
          Ok[1] = -XX[1] *det;
          Ok[2] = -XX[2] *det;
          Ok[3] =  XX[0] *det;
          
          O = mxGetPr( plhs[oid] );
          ToOutput2x2;
          break;
          

        case DET:
          Ok[0] = XX[0]*XX[3] - XX[1]*XX[2];
          
          O = mxGetPr( plhs[oid] );
          ToOutput1x1;
          break;

          
        case TRACE:
          Ok[0] =  XX[0] + XX[3];
          
          O = mxGetPr( plhs[oid] );
          ToOutput1x1;
          break;

          
        case EIGVALS:
          if( !D_computed     ){ EigenValues2x2( XX , D );          D_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !D_sorted       ){ SortEigenValues( D, ORDER , DS );  D_sorted = 1;        }
          
          memcpy( Ok , DS , N*sizeof( real ) );
          
          O = mxGetPr( plhs[oid] );
          ToOutput2x1;
          break;


        case EIGVALSDIAG:
          if( !D_computed     ){ EigenValues2x2( XX , D );          D_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !D_sorted       ){ SortEigenValues( D, ORDER , DS );  D_sorted = 1;        }

          Ok[0] = DS[0];
          Ok[3] = DS[1];
          Ok[1] = Ok[2] = 0;

          O = mxGetPr( plhs[oid] );
          ToOutput2x2;
          break;

          
        case EIGVECS:
          if( !D_computed     ){ EigenValues2x2( XX , D );          D_computed = 1;      }
          if( !V_computed     ){ EigenVectors2x2( XX , D , V);      V_computed = 1;      }
          if( !ORDER_computed ){ GetOrder( D, ORDER );              ORDER_computed = 1;  }
          if( !V_sorted       ){ SortEigenVectors( V, ORDER , VS);  V_sorted = 1;        }

          memcpy( Ok , VS , NN*sizeof( real ) );
          
          O = mxGetPr( plhs[oid] );
          ToOutput2x2;
          break;


//         case MAXDIFF:
//           if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
//           if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
// 
//           Ok[0] = V[0]*V[0]*D[0] + V[3]*V[3]*D[1] + V[6]*V[6]*D[2];
//           Ok[1] = V[1]*V[0]*D[0] + V[4]*V[3]*D[1] + V[7]*V[6]*D[2];
//           Ok[2] = V[2]*V[0]*D[0] + V[5]*V[3]*D[1] + V[8]*V[6]*D[2];
// 
//           Ok[3] = V[0]*V[1]*D[0] + V[3]*V[4]*D[1] + V[6]*V[7]*D[2];
//           Ok[4] = V[1]*V[1]*D[0] + V[4]*V[4]*D[1] + V[7]*V[7]*D[2];
//           Ok[5] = V[2]*V[1]*D[0] + V[5]*V[4]*D[1] + V[8]*V[7]*D[2];
// 
//           Ok[6] = V[0]*V[2]*D[0] + V[3]*V[5]*D[1] + V[6]*V[8]*D[2];
//           Ok[7] = V[1]*V[2]*D[0] + V[4]*V[5]*D[1] + V[7]*V[8]*D[2];
//           Ok[8] = V[2]*V[2]*D[0] + V[5]*V[5]*D[1] + V[8]*V[8]*D[2];
// 
//           Ok[0] =  MAX9( fabs( Ok[0] - XX[0] ) ,
//                          fabs( Ok[1] - XX[1] ) ,
//                          fabs( Ok[2] - XX[2] ) ,
//                          fabs( Ok[3] - XX[3] ) ,
//                          fabs( Ok[4] - XX[4] ) ,
//                          fabs( Ok[5] - XX[5] ) ,
//                          fabs( Ok[6] - XX[6] ) ,
//                          fabs( Ok[7] - XX[7] ) ,
//                          fabs( Ok[8] - XX[8] ) );
// 
//           O = mxGetPr( plhs[oid] );
//           ToOutput1x1;
//           break;
// 
//           
//         case SQRT:
//           if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
//           if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
// 
//           fD[0] = sqrt( D[0] );
//           fD[1] = sqrt( D[1] );
//           fD[2] = sqrt( D[2] );
//           
//           FillOk;
//           O = mxGetPr( plhs[oid] );
//           ToOutput3x3;
//           break;
// 
//           
//         case LOG:
//           if( !D_computed     ){ EigenValuesSym3x3( XX , D );       D_computed = 1;      }
//           if( !V_computed     ){ EigenVectorsSym3x3( XX , D , V);   V_computed = 1;      }
// 
//           fD[0] = log( D[0] );
//           fD[1] = log( D[1] );
//           fD[2] = log( D[2] );
//           
//           FillOk;
//           O = mxGetPr( plhs[oid] );
//           ToOutput3x3;
//           break;

        #define r       eval[0]
        #define e       eval[1]
        #define er      eval[2]
        case EXP:
          r = XX[0] - XX[3];
          r = 4*XX[1]*XX[2] + r*r;
          
          if( r >= 0 ){
            r = sqrt(r);
            
            
            e = exp( ( XX[0] + XX[3] - r )/2 )/r/2.0;
            er = exp(r);
            
            Ok[0] = e * ( XX[3] - XX[0] + r + er*(XX[0]-XX[3]+r) );
            
            Ok[1] = e * XX[1] * ( er - 1 ) * 2;
            
            Ok[2] = e * XX[2] * ( er - 1 ) * 2;
            
            Ok[3] = e * ( XX[0] - XX[3] + r + er*(XX[3]-XX[0]+r) );

          } else {
            
            r = sqrt( -r );
            
            e = exp( ( XX[0] + XX[3] )/2 )/r;
            er = sin(r/2);
            
            Ok[0] = e * ( r*cos(r/2) + ( XX[0]-XX[3] )*er );
            
            Ok[1] = 2 * XX[1] * e * er;
            
            Ok[2] = 2 * XX[2] * e * er;            
            
            Ok[3] = e * ( r*cos(r/2) + ( XX[3]-XX[0] )*er );
            
          }
          
          
          O = mxGetPr( plhs[oid] );
          ToOutput2x2;
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


#define r       eval[0]
#define a0pa3   eval[1]
void EigenValues2x2( real *M , real *D ) {
  
  r = M[0] - M[3];
  r = sqrt( 4*M[1]*M[2] + r*r );
  
  a0pa3 = M[0] + M[3];
  
  D[0] = ( a0pa3 - r )/2;
  D[1] = ( a0pa3 + r )/2;
  
  return;
}

#define n       evec[0]
#define a0ma3   evec[1]
void EigenVectors2x2( real *M , real *D , real *V ) {

  r = M[0] - M[3];
  r = sqrt( 4*M[1]*M[2] + r*r );
  
  a0ma3 = M[0] - M[3];
  
  V[0] = a0ma3 - r;
  V[1] = 2*M[1];

  n = 1/sqrt( V[0]*V[0] + V[1]*V[1] );
  V[0] *= n;
  V[1] *= n;
  

  V[2] = a0ma3 + r;
  V[3] = 2*M[1];

  n = 1/sqrt( V[2]*V[2] + V[3]*V[3] );
  V[2] *= n;
  V[3] *= n;

  return;
}


void GetOrder(real *D, int *P) {
  if( D[0] <= D[1] ){
    P[0]=0; P[1]=1;
  } else {
    P[0]=1; P[1]=0;
  }
  
  return;
}

void SortEigenValues( real *D, int *P , real *DS){
  DS[0] = D[ P[0] ];
  DS[1] = D[ P[1] ];
  
  return;
}

void SortEigenVectors(real *V, int *P , real *VS){
  
  if( P[0] == 0 ){
    VS[0] = V[0];
    VS[1] = V[1];
    VS[2] = V[2];
    VS[3] = V[3];
  } else {
    VS[0] = V[2];
    VS[1] = V[3];
    VS[2] = V[0];
    VS[3] = V[1];
  }

  return;
}

