/*
   phi = EulerBilinear( V , Gx , Gy , points , N , 
                           BOUNDARY_MODE , [ BOUNDARY_SIZE ]
                       );

      if N == []  set N = 500
 
      BOUNDARY_MODE:
           'value'                          stop the evolution outside the grid
           'closest'
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY)/2


   [ phi , jac ] = EulerBilinear( V , Gx , Gy , points , N , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );
 
      size( jac ) = [ size(points)  ,  2 ]

  
 
   [ phi , jac , DETjac ] = EulerBilinear( V , Gx , Gy , points , N , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );
 
      size( DETjac ) = [ numel(points)/2  ,  1 ]

 
*/

#include "myMEX.h"

#if !defined( real )
  #define   real       double
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

#define   N_def         500

#define VX(i,j)    ( ((i)!=-1 && (j)!=-1 ) ? VX[ (j)*I + (i) ] : 0 )
#define VY(i,j)    ( ((i)!=-1 && (j)!=-1 ) ? VY[ (j)*I + (i) ] : 0 )

#define Px(p)       P[ (p) ]
#define Py(p)       P[ (p) + nP ]

#define Ox(p)       O[ (p) ]
#define Oy(p)       O[ (p) + nP ]

#define OJx_x(p)      OJ[ (p)        ]
#define OJy_x(p)      OJ[ (p) +   nP ]
#define OJx_y(p)      OJ[ (p) + 2*nP ]
#define OJy_y(p)      OJ[ (p) + 3*nP ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST };
  
  real  *X, *Y, OX, LX, OY, LY, Jx_x, Jx_y, Jy_x, Jy_y;
  int     I,J,xi,yj, ii0, ii1, jj0, jj1;

  real  *VX, *VY, v[2];
  real  Vx_x, Vy_x, Vx_y, Vy_y ,  deltaJx_x, deltaJy_x, deltaJx_y, deltaJy_y;
  
  real  *P;
  int     nP, p;
  real  XY[2], xyl[2];

  real  dt;
  int     m, N;
  real  D, U0, U1, V0, V1;

  real  X0,X1,Y0,Y1;
  real  VX00, VX10, VX01, VX11;
  real  VY00, VY10, VY01, VY11;
  real  boundary_size;
  enum    boundary_modes boundary_mode;
  int     isOutside;
  
  
  real  *O, *OJ, *DetJ;
  char    STR[100];
  int     argN;
  int     Odims[50], ndims, d, computeJAC;
  
  

  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );

  if( myNumel( prhs[0] ) != I*J*2 ){
    myErrMsgTxt("size(V) has to be [numel(Gx) numel(Gy) 2].");
  }
  
  
  VX = myGetPr( prhs[0] );
  VY = VX + I*J;

  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(Gy) Coordinates do not coincide with size(V,2)."); }

  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("oGx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("oGy  Coordinates are not sorted."); }
  
  nP = myNumel( prhs[3] );
  nP = nP/2;
  P  = myGetPr( prhs[3] );


  if( nrhs > 4 && myIsParameter( prhs[4] ) && !myIsEmpty( prhs[4] ) ){ 
    N = (int) myGetValue( prhs[4] ); 
  } else {
    N = N_def;
  }
  dt = 1.0/N;

  
  /*Parsing arguments*/
  /*Defaults*/
  boundary_mode = VALUE;
  boundary_size = (X[I-1]-X[0]+Y[J-1]-Y[0])/2;
  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );

      if( ! myStrcmpi(STR,"circular") || ! myStrcmpi(STR,"periodic") || 
          ! myStrcmpi(STR,"circ")     || ! myStrcmpi(STR,"per")      ) { boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"zero")     ) { boundary_mode = VALUE;    argN++; continue; }

      if( ! myStrcmpi( STR,"decay_to_zero") || 
          ! myStrcmpi( STR,"tozero")  || ! myStrcmpi( STR,"decay") ){
        boundary_mode = DECAY_TO_ZERO; argN++;
        if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
          if( !myIsEmpty( prhs[argN] ) ){ boundary_size = myGetValue(prhs[argN]); }
          argN++; continue;
        }
        continue;
      }

      mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
    }
  }
  /*END Parsing arguments*/
  
  
  /*Creating output*/
  plhs[0]  = myDuplicateArrayWithClass( prhs[3] , mxREAL_CLASS , mxREAL );
  O = (real *)mxGetData( plhs[0] );
  
  if( nlhs > 1 ){
    ndims = myNDims( prhs[3] );
    for( d=0 ; d<ndims ; d++ ){
      Odims[d] = mySize( prhs[3] , d );
    }
    Odims[d++] = 2;
    plhs[1] = mxCreateNumericArray( d , Odims , mxREAL_CLASS , mxREAL );
    OJ = (real *) mxGetData( plhs[1] );
    computeJAC = 1;
  } else {
    computeJAC = 0;
  }
  if( nlhs > 2 ){
    Odims[0] = nP;
    Odims[1] = 1;
    plhs[2] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
    DetJ = (real *) mxGetData( plhs[2] );
  }
  /*END Creating output*/

  switch( boundary_mode ){
    case VALUE:
      OX = X[0];  LX = X[I-1];
      OY = Y[0];  LY = Y[J-1];
      break;
    case DECAY_TO_ZERO:
      OX = X[0]-boundary_size;  LX = X[I-1]+boundary_size;
      OY = Y[0]-boundary_size;  LY = Y[J-1]+boundary_size;
      break;
    case CIRCULAR:
      OX = X[0] - (X[1]-X[0])/2;  LX = X[I-1] + ( X[I-1]-X[I-2] )/2 - OX;
      OY = Y[0] - (Y[1]-Y[0])/2;  LY = Y[J-1] + ( Y[J-1]-Y[J-2] )/2 - OY;
      break;
  }
  
  xi = -1;  yj = -1;
  for( p=0 ; p<nP ; p++ ){
    XY[0] = Px(p);  XY[1] = Py(p);
    
    if( computeJAC ){
      Jx_x = 1;  Jy_x = 0;
      Jx_y = 0;  Jy_y = 1;
    }

    m = 0;
    while( m < N ){
      switch( boundary_mode ){
        case VALUE:
          if( XY[0] < OX || XY[0] > LX ){ isOutside= 1; break; } 
          if( XY[1] < OY || XY[1] > LY ){ isOutside= 1; break; } 
          isOutside = 0;

          xyl[0] = XY[0];
          xi = GetInterval( xyl[0] , X , I , xi );

          ii0 = xi;             ii1 = xi+1;
          X0  = X[ xi ];         X1 = X[xi+1];
          
          xyl[1] = XY[1];
          yj = GetInterval( xyl[1] , Y , J , yj );
          jj0 = yj;             jj1 = yj+1;
          Y0  = Y[ yj ];         Y1 = Y[yj+1];
          break;

        case CLOSEST:
          xyl[0] = XY[0];
          if( xyl[0] < X[0] ){ xyl[0] = X[0]; } else if( xyl[0] > X[I-1] ){ xyl[0] = X[I-1]; }

          xi = GetInterval( xyl[0] , X , I , xi );
          ii0 = xi;             ii1 = xi+1;
          X0  = X[ xi ];         X1 = X[xi+1];
          
          xyl[1] = XY[1];
          if( xyl[1] < Y[0] ){ xyl[1] = Y[0]; } else if( xyl[1] > Y[J-1] ){ xyl[1] = Y[J-1]; }
            
          yj = GetInterval( xyl[1] , Y , J , yj );
          jj0 = yj;             jj1 = yj+1;
          Y0  = Y[ yj ];         Y1 = Y[yj+1];
          break;

        case CIRCULAR:
          isOutside = 0;
          xyl[0] = PutInside(XY[0],OX,LX);
          xyl[1] = PutInside(XY[1],OY,LY);


          xi = GetInterval( xyl[0] , X , I , xi );
          if( xi >= 0 ){
            ii0 = xi;             ii1 = xi+1;
            X0  = X[ii0];         X1 = X[ii1];
          } else if( xi == - 2  ) {
            ii0 = I-1;            ii1 = 0;
            X0  = X[0] - ( X[1]-X[0] + X[I-1]-X[I-2] )/2;
            X1  = X[0];
          } else if( xi == -101 ) {
            ii0 = I-1;            ii1 = 0;
            X0  = X[I-1];
            X1  = X[I-1] + ( X[1]-X[0] + X[I-1]-X[I-2] )/2;
          }

          yj = GetInterval( xyl[1] , Y , J , yj );
          if( yj >= 0 ){
            jj0 = yj;             jj1 = yj+1;
            Y0  = Y[jj0];         Y1 = Y[jj1];
          } else if( yj == - 2  ) {
            jj0 = J-1;            jj1 = 0;
            Y0  = Y[0] - ( Y[1]-Y[0] + Y[J-1]-Y[J-2] )/2;
            Y1  = Y[0];
          } else if( yj == -101 ) {
            jj0 = J-1;            jj1 = 0;
            Y0  = Y[J-1];
            Y1  = Y[J-1] + ( Y[1]-Y[0] + Y[J-1]-Y[J-2] )/2;
          }

          break;

        case DECAY_TO_ZERO:
          if( XY[0] < OX || XY[0] > LX ){ isOutside= 1; break; } 
          if( XY[1] < OY || XY[1] > LY ){ isOutside= 1; break; } 
          isOutside = 0;
          
          xyl[0] = XY[0];
          xi = GetInterval( xyl[0] , X , I , xi );
          if( xi >= 0 ){
            ii0 = xi;             ii1 = xi+1;
            X0  = X[ii0];         X1 = X[ii1];
          } else if( xi == - 2  ) {
            ii0 = -1;             ii1 = 0;
            X0  = OX;             X1  = X[0];
          } else if( xi == -101 ) {
            ii0 = I-1;            ii1 = -1;
            X0  = X[I-1];         X1  = LX;
          }

          xyl[1] = XY[1];
          yj = GetInterval( xyl[1] , Y , J , yj );
          if( yj >= 0 ){
            jj0 = yj;             jj1 = yj+1;
            Y0  = Y[jj0];         Y1 = Y[jj1];
          } else if( yj == - 2  ) {
            jj0 = -1;             jj1 = 0;
            Y0  = OY;             Y1  = Y[0];
          } else if( yj == -101 ) {
            jj0 = J-1;            jj1 = -1;
            Y0  = Y[J-1];         Y1  = LY;
          }
          break;

      }
      if( isOutside ){ break; }
      
      D = dt/( ( X1 - X0 )*( Y1 - Y0 ) );
      
      VX00 = VX( ii0 , jj0 );      VX10 = VX( ii1 , jj0 );
      VX01 = VX( ii0 , jj1 );      VX11 = VX( ii1 , jj1 );

      VY00 = VY( ii0 , jj0 );      VY10 = VY( ii1 , jj0 );
      VY01 = VY( ii0 , jj1 );      VY11 = VY( ii1 , jj1 );

      while( xyl[0] >= X0 && xyl[0] <= X1 &&
             xyl[1] >= Y0 && xyl[1] <= Y1 &&
             m < N ) {

        U0 = X0 - xyl[0];          U1 = X1 - xyl[0];
        V0 = Y0 - xyl[1];          V1 = Y1 - xyl[1];

        v[0] = ( U1*V1*VX00 - U0*V1*VX10 - U1*V0*VX01 + U0*V0*VX11 )*D;
        v[1] = ( U1*V1*VY00 - U0*V1*VY10 - U1*V0*VY01 + U0*V0*VY11 )*D;

/* With this break is slower!!!
        if( v[0]==0 && v[1]==0 ){
          m = N;
          break;
        }
*/

        XY[0]  += v[0];
        XY[1]  += v[1];
        
        xyl[0] += v[0];
        xyl[1] += v[1];
        
        
        if( computeJAC ){
          Vx_x = ( ( VX01 - VX11 )*V0 + ( VX10 - VX00 )*V1 )*D;
          Vy_x = ( ( VY01 - VY11 )*V0 + ( VY10 - VY00 )*V1 )*D;
          Vx_y = ( ( VX10 - VX11 )*U0 + ( VX01 - VX00 )*U1 )*D;
          Vy_y = ( ( VY10 - VY11 )*U0 + ( VY01 - VY00 )*U1 )*D;
          
          deltaJx_x = Vx_x*Jx_x + Vx_y*Jy_x;
          deltaJy_x = Vy_x*Jx_x + Vy_y*Jy_x;
          deltaJx_y = Vx_x*Jx_y + Vx_y*Jy_y;
          deltaJy_y = Vy_x*Jx_y + Vy_y*Jy_y;
          
          Jx_x += deltaJx_x;
          Jx_y += deltaJx_y;
          Jy_x += deltaJy_x;
          Jy_y += deltaJy_y;
        }
        
        m++;
      }
    }
    Ox(p) = XY[0];     Oy(p) = XY[1];
    
    if( computeJAC ){
      OJx_x(p) = Jx_x;   OJx_y(p) = Jx_y;   
      OJy_x(p) = Jy_x;   OJy_y(p) = Jy_y;   
    }

  }

  if( nlhs > 2 ){
    for( p=0 ; p<nP ; p++ ){
      DetJ[p] = OJx_x(p)*OJy_y(p) - OJy_x(p)*OJx_y(p);
    }
  }

  EXIT: myFreeALLOCATES();
}
