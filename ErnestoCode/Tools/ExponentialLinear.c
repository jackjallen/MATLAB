/*
  phi = ExponentialLinear( V , gX , POINTS , 
                           BOUNDARY_MODE , [ BOUNDARY_SIZE ]
                         );

      BOUNDARY_MODE:
           'value'                          stop the evolution outside the grid
           'closest'
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY)/2

 
  [ phi , jac ] = ExponentialLinear( V , gX , POINTS , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );

      size(jac) = [ size( points ) ]

 */

#include "myMEX.h"

#if !defined( real )
  #define   real       double
#endif
#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

#define V(i)       V[ (i) ]

#define P(p)       P[ (p) ]

#define O(p)       O[ (p) ]


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID , CLOSEST };
  real    *X, dx;                           /*/ grillas */
  int     I,ii;        /*/ sizes, indices   */
  real    X0, X1, V0, V1;
  
  real    *V, A, B, Vx;       /*/input velocities, and v at current point (homogeneous coordinates)  */
  
  real    *P;                   /*/ todos los puntos  */
  int     nP, p, K, ITER;                /*/ numero total puntos, id punto actual  */
  real    x;       /*/ punto actual; (homogeneous coordinates)  */
  
  real    T;
  real    T_evol;
  real    log_arg, J;
  
  real    boundary_size;
  real    *O, *OJ;                   /*/ output   */
  real    OX, LX;
  
  enum boundary_modes boundary_mode;
  
  char STR[100];
  int  argN, computeJAC;

  V  = myGetPr(prhs[0]);
  I  = myNumel(prhs[0]);  
  
  if( myNumel(prhs[1]) != I ){ 
    myErrMsgTxt("numel(oGx) Coordinates do not coincide with numel(V)."); 
  }
  X  = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("gX  Coordinates are not sorted."); }
  
  
  nP = myNumel( prhs[2] );
  P  = myGetPr( prhs[2] );
  

  /*Parsing arguments*/
  /*Defaults*/
  boundary_mode = VALUE;
  boundary_size = X[I-1]-X[0];
  
  argN = 3;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );
      if( ! myStrcmpi(STR,"circular") || ! myStrcmpi(STR,"periodic") ||
          ! myStrcmpi(STR,"circ")     || ! myStrcmpi(STR,"per")    ) { boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"zero")   ) { boundary_mode = VALUE;    argN++; continue; }
      if( ! myStrcmpi(STR,"closest")                               ) { boundary_mode = CLOSEST;  argN++; continue; }

      if( ! myStrcmpi( STR,"decay_to_zero") || 
          ! myStrcmpi( STR,"tozero")  || ! myStrcmpi( STR,"decay") ){
        boundary_mode = DECAY_TO_ZERO; argN++;
        if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
          if( !myIsEmpty( prhs[argN] ) ){ boundary_size = myGetValue( prhs[argN] ); }
          argN++; continue;
        }
        continue;
      }

      mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
    }
  }
  /*END Parsing arguments*/
  
  
  if( boundary_mode == CIRCULAR ){
    OX = X[ 0 ] - ( X[ 1 ] - X[ 0 ] )/2;
    LX = X[I-1] + ( X[I-1] - X[I-2] )/2 - OX;
    
//     mexPrintf("OX: %g   LX: %g\n",OX,LX);
  }
  
  
  /*Creating output*/
  plhs[0]  = myDuplicateArrayWithClass( prhs[2] , mxREAL_CLASS , mxREAL );
  O = (real *)mxGetData( plhs[0] );
  
  if( nlhs > 1 ){
    plhs[1]  = myDuplicateArrayWithClass( prhs[2] , mxREAL_CLASS , mxREAL );
    OJ = (real *) mxGetData( plhs[1] );
    computeJAC = 1;
  } else {
    computeJAC = 0;
  }
  /*END Creating output*/
  
  ii = -1;
  for( p=0 ; p < nP ; p++ ){
    x = P(p);  if( computeJAC ){ J = 1; }
    K = 0;
//     mexPrintf("x: %.8g\n",x); myFlush();
    
    ii = GetInterval( x , X , I , ii );
//     mexPrintf("ii: %d\n",ii); myFlush();

    if( boundary_mode == DECAY_TO_ZERO && ( x < X[ 0 ]-boundary_size || x > X[I-1]+boundary_size ) ){
      O(p) = x;
      if( computeJAC ) { OJ[p] = J; }
      continue;
    }
    if( boundary_mode == VALUE         && ( x < X[ 0 ]               || x > X[I-1]               ) ){
      O(p) = x;
      if( computeJAC ) { OJ[p] = J; }
      continue;
    }
    
    
    T = 0;   ITER = 0;
    
    while( T < 1 ){
      switch( boundary_mode ){
        
/*****************************************************************************

... -2   ) [.  0   ) [.  1   ) [.  2   ) ... [. I-3  ) [. I-2  ] (   -101  ...
          *         *         *             *         *         *
         G0        G1        G2            GI-3      GI-2      GI-1

*****************************************************************************/
        
        case DECAY_TO_ZERO:
          if( ii == -2 || ii == -1 || ( x == X[0]  &&  V[0] < 0 ) ){
            ii =  -1;
            X0 = X1 - boundary_size;        V1 = V[0];
            X1 = X[0];                      V0 = 0;
          } else if( ii == -101 || ii >= I-1 || ( x == X[I-1] && V[I-1] > 0 ) ){
            ii =  I-1;
            X0 =  X[I-1];                   V0 = V[I-1];
            X1 =  X0 + boundary_size;       V1 = 0;
          } else if( x == X[ii] && V[ii] < 0 ){
            ii = ii-1;
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          } else {
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          }
          break;
        case VALUE:
          if( ii == -2 || ii == -1 || ( x == X[0]  &&  V[0] < 0 ) ){
            ii = -1;
            X0 = X1 - 1000000000;           V1 = 0;
            X1 = X[0];                      V0 = 0;
          } else if( ii == -101 || ii >= I-1 || ( x == X[I-1] && V[I-1] > 0 ) ){
            ii =  I-1;
            X0 =  X[I-1];                   V0 = 0;
            X1 =  X0 + 1000000000;          V1 = 0;
          } else if( x == X[ii] && V[ii] < 0 ){
            ii = ii-1;
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          } else {
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          }
          break;
        case CLOSEST:
          if( ii == -2 || ii == -1 || ( x == X[0]  &&  V[0] < 0 ) ){
            ii = -1;
            X0 = X1 - 1000000000;           V1 = V[0];
            X1 = X[0];                      V0 = V[0];
          } else if( ii == -101 || ii >= I-1 || ( x == X[I-1] && V[I-1] > 0 ) ){
            ii =  I-1;
            X0 = X[I-1];                    V0 = V[I-1];
            X1 = X0 + 1000000000;           V1 = V[I-1];
          } else if( x == X[ii] && V[ii] < 0 ){
            ii = ii-1;
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          } else {
            X0 = X[ ii   ];                 V0 = V[ ii   ];
            X1 = X[ ii+1 ];                 V1 = V[ ii+1 ];
          }
          break;
        case CIRCULAR:
          if( ii < 0 || ii > I-2 ){
            K  = floor( ( x - OX )/LX );
            ii = GetInterval( x - K*LX , X , I , ii );
          }
            
          if( ii == -1 || ii == -2 || ( ( x - K*LX ) == X[0]   &&  V[0] < 0 ) ){
            ii = -1;
            X0 = X[I-1] +  K   *LX;         V0 = V[I-1];
            X1 = X[ 0 ] + (K-1)*LX;         V1 = V[ 0 ];
          } else if( ii == -101 || ( ( x - K*LX ) == X[I-1]  &&  V[I-1] > 0 ) ){
            ii = I-1;
            X0 = X[I-1] + (K-1)*LX;         V0 = V[I-1];
            X1 = X[ 0 ] +  K   *LX;         V1 = V[ 0 ];
          } else if( x - K*LX == X[ii] && V[ii] < 0 ){
            ii = ii-1;
            X0 = X[ ii   ]+K*LX;            V0 = V[ ii   ];
            X1 = X[ ii+1 ]+K*LX;            V1 = V[ ii+1 ];
          } else {
            X0 = X[ ii   ]+K*LX;            V0 = V[ ii   ];
            X1 = X[ ii+1 ]+K*LX;            V1 = V[ ii+1 ];
          }
          break;

      }

      
      A = ( V1 - V0 )      /( X1 - X0 );
      B = ( V0*X1 - V1*X0 )/( X1 - X0 );
      
//       if( ii >= 0 && (x - K*LX) == X[ii] ){
//         Vx = V[ ii ];
//       } else if( ii < I-1 && (x - K*LX) == X[ii+1] ) {
//         Vx = V[ii+1];
//       } else {
        Vx = x*A + B;
//       }
      
      if( Vx == 0 ){
        break;
      } else if( Vx > 0 ){

        log_arg = ( X1*A + B )/Vx;
        if( log_arg < 1  &&  A > 0 ){
          T_evol = 2;
        } else {
          T_evol = log( log_arg )/A;
        }
        
        if( T_evol == 0 ){
          ITER = ITER + 1;
          if( ITER > 2 ){ break; }
        } else {
          ITER = 0;
        }
          
        
        if(  ( 1 - T ) > T_evol ){
          x = X1;
          T += T_evol;
          ii++;
          if( computeJAC ) { J = J * exp( A * T_evol ); }
        } else {
          x = ( exp( A*( 1 - T ) )*Vx - B )/A;
          if( computeJAC ) { J = J * exp( A * (1-T) );  }
          break;
        }
        
      } else {
        
        log_arg = ( X0*A + B )/Vx;
        if( log_arg < 1  &&  A > 0 ){
          T_evol = 2;
        } else {
          T_evol = log( log_arg )/A;
        }

        if( T_evol == 0 ){
          ITER = ITER + 1;
          if( ITER > 2 ){ break; }
        } else {
          ITER = 0;
        }

        if(  ( 1 - T ) > T_evol ){
          x = X0;
          T += T_evol;
          ii--;
          if( computeJAC ) { J = J * exp( A * T_evol ); }
        } else {
          x = ( exp( A*( 1 - T ) )*Vx - B )/A;
          if( computeJAC ) { J = J * exp( A * (1-T) );  }
          break;
        }

      }
      
    }
      
    O(p) = x;
    if( computeJAC ) { OJ[p] = J; }
  }

  EXIT: myFreeALLOCATES();
}  
