/*
  phi = Euler4D( V , Gx , Gy , Gz , Gt , points , N , [t0 t1] ,
                       BOUNDARY_MODE , [ BOUNDARY_SIZE ]
                     );

      size(V) = [ numel(Gx) numel(Gy) numel(Gz) numel(Gt) 3 ];
 
 
  
    if N == []  set N = 500

      if [t0 t1] == [] return flow ( from t=Gt(1) to t=Gt(end) at all times Gt );

      only allowed  t1>t0 !!
 
 
      BOUNDARY_MODE:
           'value'                          stop the evolution outside the grid
           'closest'
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY+LZ)/3

 
 
    size( phi ) = size( points )                         if [t0 t1] is specified
    size( phi ) = [ numel(points)/3 , 3 , numel(Gt) ]    if [t0 t1]==[]  (flow case)


 
  
 
  [ phi , jac ] = Euler4D( V , Gx , Gy , Gz , Gt , points , N , [t0 t1] , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );
 
    size( jac ) = [ size( points )  , 3 ]                    if [t0 t1] is specified
    size( jac ) = [ numel(points)/3 , 3 , 3 , numel(Gt) ]    if [t0 t1]==[]  (flow case)
 


  [ phi , jac , DETjac ] = Euler4D( V , Gx , Gy , Gz , Gt , points , N , [t0 t1] , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );
 
    size( DETjac ) = [ numel(points)/3  , 1 ]            if [t0 t1] is specified
    size( DETjac ) = [ numel(points)/3  , numel(Gt) ]    if [t0 t1]==[]  (flow case)
 
 */

#include "myMEX.h"

#if !defined( real )
  #define     real       double
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

#define   N_def         500

#define VX(i,j,k,t)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VX[ (t)*IJK + (k)*IJ + (j)*I + (i) ] : 0 )
#define VY(i,j,k,t)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VY[ (t)*IJK + (k)*IJ + (j)*I + (i) ] : 0 )
#define VZ(i,j,k,t)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VZ[ (t)*IJK + (k)*IJ + (j)*I + (i) ] : 0 )

#define Px(p)          P[ (p)        ]
#define Py(p)          P[ (p) +   nP ]
#define Pz(p)          P[ (p) + 2*nP ]

#define Fx(p,t)        O[ (p)        + (t)*3*nP ]
#define Fy(p,t)        O[ (p) +   nP + (t)*3*nP ]
#define Fz(p,t)        O[ (p) + 2*nP + (t)*3*nP ]

#define FJx_x(p,t)    OJ[ (p)        + (t)*9*nP ]
#define FJy_x(p,t)    OJ[ (p) +   nP + (t)*9*nP ]
#define FJz_x(p,t)    OJ[ (p) + 2*nP + (t)*9*nP ]

#define FJx_y(p,t)    OJ[ (p) + 3*nP + (t)*9*nP ]
#define FJy_y(p,t)    OJ[ (p) + 4*nP + (t)*9*nP ]
#define FJz_y(p,t)    OJ[ (p) + 5*nP + (t)*9*nP ]

#define FJx_z(p,t)    OJ[ (p) + 6*nP + (t)*9*nP ]
#define FJy_z(p,t)    OJ[ (p) + 7*nP + (t)*9*nP ]
#define FJz_z(p,t)    OJ[ (p) + 8*nP + (t)*9*nP ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST };

  real    *X, *Y, *Z, OX, LX, OY, LY, OZ, LZ, T0, T1, *T;
  int     I,J,K,xi,yj,zk, ii0, ii1, jj0, jj1, kk0, kk1, tt0, tt1;
  int     IJ, IJK, S;

  real    *VX, *VY, *VZ, v[3], v0[3], v1[3];

  real    *P;
  int     nP, p;
  real    XYZ[3], xyzl[3];

  real    dt, tau, deltaT;
  int     N;
  real    D, U0, U1, V0, V1, W0, W1, X0, X1, Y0, Y1, Z0, Z1;
  real    VX0000, VX0010, VX0100, VX0110, VX1000, VX1010, VX1100, VX1110;
  real    VY0000, VY0010, VY0100, VY0110, VY1000, VY1010, VY1100, VY1110;
  real    VZ0000, VZ0010, VZ0100, VZ0110, VZ1000, VZ1010, VZ1100, VZ1110;

  real    VX0001, VX0011, VX0101, VX0111, VX1001, VX1011, VX1101, VX1111;
  real    VY0001, VY0011, VY0101, VY0111, VY1001, VY1011, VY1101, VY1111;
  real    VZ0001, VZ0011, VZ0101, VZ0111, VZ1001, VZ1011, VZ1101, VZ1111;
  
  real    Jx_x, Jx_y, Jx_z, Jy_x, Jy_y, Jy_z, Jz_x, Jz_y, Jz_z;
  real    deltaJx_x, deltaJx_y, deltaJx_z, deltaJy_x, deltaJy_y, deltaJy_z, deltaJz_x, deltaJz_y, deltaJz_z;
  real    Vx_x, Vx_y, Vx_z, Vy_x, Vy_y, Vy_z, Vz_x, Vz_y, Vz_z;
  
  real    boundary_size;
  enum    boundary_modes boundary_mode;
  
  real   *O, *OJ, *DetJ;
  char    STR[100];
  int     argN;
  int     isOutside, FLOW, computeJAC;
  int     Odims[3], ndims, d;
  
  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );
  K  = mySize( prhs[0] , 2 );
  S  = mySize( prhs[0] , 3 );
  IJ = I*J;
  IJK = IJ*K;

  
  if( S == 1 ){
    myErrMsgTxt("size(V,4) has to be at least 2. Else use Euler Trilinear.");
  }
  
  if( myNumel( prhs[0] ) != IJK*S*3 ){
    myErrMsgTxt("size(V) has to be [numel(Gx) numel(Gy) numel(Gz) numel(Gt) 3].");
  }
  
  VX = myGetPr( prhs[0] );
  VY = VX + IJK*S;
  VZ = VY + IJK*S;
  
  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(Gy) Coordinates do not coincide with size(V,2)."); }
  if( myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(Gz) Coordinates do not coincide with size(V,3)."); }
  if( myNumel( prhs[4] ) != S ){ myErrMsgTxt("numel(Gt) Coordinates do not coincide with size(V,4)."); }
  

  X  = myGetPr( prhs[1] );
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("Gx  Coordinates are not sorted."); }
  
  Y  = myGetPr( prhs[2] );
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("Gy  Coordinates are not sorted."); }
  
  Z  = myGetPr( prhs[3] );
  if( !checkIsSorted(Z,K) ){ myErrMsgTxt("Gz  Coordinates are not sorted."); }
  
  T  = myGetPr( prhs[4] );
  if( !checkIsSorted(T,S) ){ myErrMsgTxt("Gt  Coordinates are not sorted."); }



  nP = myNumel( prhs[5] );
  nP = nP/3;
  P  = myGetPr( prhs[5] );

  if( nrhs > 6 && myIsParameter( prhs[6] ) && !myIsEmpty( prhs[6] ) ){
    N = (int) myGetValue( prhs[6] ); 
  } else {
    N = N_def;
  }
  
  if( nrhs > 7 && !myIsEmpty( prhs[7] ) ){
    if( myNumel( prhs[7] ) != 2 ){ myErrMsgTxt("[t0 t1] bad specified."); }
    T0 = myGetValueIth( prhs[7] , 0 );
    T1 = myGetValueIth( prhs[7] , 1 );
    if( T0 <  T[ 0 ] ) { myErrMsgTxt("t0 cannot be smaller than Gt(1)."); }
    if( T1 >  T[S-1] ) { myErrMsgTxt("t1 cannot be greater than Gt(end)."); }
    FLOW = 0;
  }  else {
    FLOW = 1; T0 = T[0]; T1 = T[S-1];
  }
  if( T0 >= T1     ) { myErrMsgTxt("t1 has to be greater than t0."); }
  dt = (T1-T0)/N;
  
  /*Parsing arguments*/
  /*Defaults*/
  boundary_mode = VALUE;
  boundary_size = (X[I-1]-X[0]+Y[J-1]-Y[0]+Z[K-1]-Z[0])/3;
  
  argN = 8;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );

      if( ! myStrcmpi(STR,"circular")  || 
          ! myStrcmpi(STR,"periodic") 
           || ! myStrcmpi(STR,"circ") || ! myStrcmpi(STR,"per") ) { boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"value") || ! myStrcmpi(STR,"zero")   ) { boundary_mode = VALUE;    argN++; continue; }

      if( ! myStrcmpi( STR,"decay_to_zero") 
              || ! myStrcmpi( STR,"tozero") || ! myStrcmpi( STR,"decay") ){
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
  
  
  /*Creating output*/
  if( FLOW ){
    Odims[0]= nP;
    Odims[1]= 3;
    Odims[2]= S;
    plhs[0] = mxCreateNumericArray( 3 , Odims , mxREAL_CLASS , mxREAL );
  } else {
    plhs[0] = myDuplicateArrayWithClass( prhs[5] , mxREAL_CLASS , mxREAL );
  }
  O = (real *) mxGetData( plhs[0] );

  if( nlhs > 1 ){
    if( FLOW ){
      Odims[2] = 3;
      Odims[3] = S;
      plhs[1] = mxCreateNumericArray( 4 , Odims , mxREAL_CLASS , mxREAL );
    } else {
      ndims = myNDims( prhs[5] );
      for( d=0 ; d<ndims ; d++ ){
        Odims[d] = mySize( prhs[5] , d );
      }
      Odims[d++] = 3;
      plhs[1] = mxCreateNumericArray( d , Odims , mxREAL_CLASS , mxREAL );
    }
    OJ = (real *) mxGetData( plhs[1] );
    computeJAC = 1;
  } else {
    computeJAC = 0;
  }

  if( nlhs > 2 ){
    Odims[0] = nP;
    if( FLOW ){
      Odims[1] = S;
    } else {
      Odims[1] = 1;
    }
    plhs[2] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
    DetJ = (real *) mxGetData( plhs[2] );
  }
  /*END Creating output*/    
  
  switch( boundary_mode ){
    case VALUE:
      OX = X[0];  LX = X[I-1];
      OY = Y[0];  LY = Y[J-1];
      OZ = Z[0];  LZ = Z[K-1];
      break;
    case DECAY_TO_ZERO:
      OX = X[0]-boundary_size;  LX = X[I-1]+boundary_size;
      OY = Y[0]-boundary_size;  LY = Y[J-1]+boundary_size;
      OZ = Z[0]-boundary_size;  LZ = Z[K-1]+boundary_size;
      break;
    case CIRCULAR:
      OX = X[0] - (X[1]-X[0])/2;  LX = X[I-1] + ( X[I-1]-X[I-2] )/2 - OX;
      OY = Y[0] - (Y[1]-Y[0])/2;  LY = Y[J-1] + ( Y[J-1]-Y[J-2] )/2 - OY;
      OZ = Z[0] - (Z[1]-Z[0])/2;  LZ = Z[K-1] + ( Z[K-1]-Z[K-2] )/2 - OZ;
      break;
  }
  
  xi = -1;  yj = -1;   zk= -1;   tt0 = -1;
  isOutside = 0;
  for( p=0 ; p<nP ; p++ ){
    XYZ[0] = Px(p);    XYZ[1] = Py(p);    XYZ[2] = Pz(p);

    if( computeJAC ){
      Jx_x = 1;  Jx_y = 0;  Jx_z = 0;
      Jy_x = 0;  Jy_y = 1;  Jy_z = 0;
      Jz_x = 0;  Jz_y = 0;  Jz_z = 1;
    }
    
    tau = T0;
    if( FLOW ){
      Fx(p,0) = XYZ[0];   Fy(p,0) = XYZ[1];   Fz(p,0) = XYZ[2];
      if( computeJAC ){
        FJx_x(p,0) = Jx_x;    FJx_y(p,0) = Jx_y;    FJx_z(p,0) = Jx_z;
        FJy_x(p,0) = Jy_x;    FJy_y(p,0) = Jy_y;    FJy_z(p,0) = Jy_z;
        FJz_x(p,0) = Jz_x;    FJz_y(p,0) = Jz_y;    FJz_z(p,0) = Jz_z;
      }
    }

    while( tau < T1 ){
      tt0 = GetInterval( tau , T , S , tt0 );  tt1 = tt0+1;
      if( tt0 < 0 || tt1 < 0 || tt1 >= S ){
        mexPrintf("error, llamarme!!!  tau: %g  - tt0: %d  - tt1: %d\n",tau,tt0,tt1);
        myErrMsgTxt("tau out of range!!!");
      }

      switch( boundary_mode ){
        case VALUE:
          if( XYZ[0] < OX || XYZ[0] > LX ){ isOutside= 1; break; } 
          if( XYZ[1] < OY || XYZ[1] > LY ){ isOutside= 1; break; } 
          if( XYZ[2] < OZ || XYZ[2] > LZ ){ isOutside= 1; break; } 
          isOutside = 0;

          xyzl[0] = XYZ[0];
          xi = GetInterval( xyzl[0] , X , I , xi );
          ii0 = xi;             ii1 = xi+1;
          X0  = X[ xi ];         X1 = X[xi+1];
          
          xyzl[1] = XYZ[1];
          yj = GetInterval( xyzl[1] , Y , J , yj );
          jj0 = yj;             jj1 = yj+1;
          Y0  = Y[ yj ];         Y1 = Y[yj+1];
          
          xyzl[2] = XYZ[2];
          zk = GetInterval( xyzl[2] , Z , K , zk );
          kk0 = zk;             kk1 = zk+1;
          Z0  = Z[ zk ];         Z1 = Z[zk+1];
          break;

        case CLOSEST:
          xyzl[0] = XYZ[0];
          if( xyzl[0] < X[0] ){ xyzl[0] = X[0]; } else if( xyzl[0] > X[I-1] ){ xyzl[0] = X[I-1]; }

          xi = GetInterval( xyzl[0] , X , I , xi );
          ii0 = xi;             ii1 = xi+1;
          X0  = X[ xi ];         X1 = X[xi+1];
          
          xyzl[1] = XYZ[1];
          if( xyzl[1] < Y[0] ){ xyzl[1] = Y[0]; } else if( xyzl[1] > Y[J-1] ){ xyzl[1] = Y[J-1]; }
            
          yj = GetInterval( xyzl[1] , Y , J , yj );
          jj0 = yj;             jj1 = yj+1;
          Y0  = Y[ yj ];         Y1 = Y[yj+1];
          
          xyzl[2] = XYZ[2];
          if( xyzl[2] < Z[0] ){ xyzl[2] = Z[0]; } else if( xyzl[2] > Z[K-1] ){ xyzl[2] = Z[K-1]; }

          zk = GetInterval( xyzl[2] , Z , K , zk );
          kk0 = zk;             kk1 = zk+1;
          Z0  = Z[ zk ];         Z1 = Z[zk+1];
          break;

        case CIRCULAR:
//          if( XYZ[0] < OX ){             xyzl[0] = XYZ[0] - ( (int) ( (XYZ[0] - OX)/LX ) - 1 )*LX;
//          } else if( XYZ[0] > OX+LX ) {  xyzl[0] = XYZ[0] - ( (int) ( (XYZ[0] - OX)/LX )     )*LX;
//          } else {                       xyzl[0] = XYZ[0];
//          }
//  
//        if( yy < OY ){             yy = yy - ( (int) ( (yy - OY)/LY ) - 1 )*LY;
//        } else if( yy > OY+LY ) {  yy = yy - ( (int) ( (yy - OY)/LY )     )*LY;
//        }
//
//        if( zz < OZ ){             zz = zz - ( (int) ( (zz - OZ)/LZ ) - 1 )*LZ;
//        } else if( zz > OZ+LZ ) {  zz = zz - ( (int) ( (zz - OZ)/LZ )     )*LZ;
//        }
//        break;
          
          
          
          xyzl[0] = PutInside(XYZ[0],OX,LX);
          xyzl[1] = PutInside(XYZ[1],OY,LY);
          xyzl[2] = PutInside(XYZ[2],OZ,LZ);


          xi = GetInterval( xyzl[0] , X , I , xi );
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

          yj = GetInterval( xyzl[1] , Y , J , yj );
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

          zk = GetInterval( xyzl[2] , Z , K , zk );
          if( zk >= 0 ){
            kk0 = zk;             kk1 = zk+1;
            Z0  = Z[kk0];         Z1 = Z[kk1];
          } else if( zk == - 2  ) {
            kk0 = K-1;            kk1 = 0;
            Z0  = Z[0] - ( Z[1]-Z[0] + Z[K-1]-Z[K-2] )/2;
            Z1  = Z[0];
          } else if( zk == -101 ) {
            kk0 = K-1;            kk1 = 0;
            Z0  = Z[K-1];
            Z1  = Z[K-1] + ( Z[1]-Z[0] + Z[K-1]-Z[K-2] )/2;
          }

          break;

        case DECAY_TO_ZERO:
          if( XYZ[0] < OX || XYZ[0] > LX ){ isOutside= 1; break; } 
          if( XYZ[1] < OY || XYZ[1] > LY ){ isOutside= 1; break; } 
          if( XYZ[2] < OZ || XYZ[2] > LZ ){ isOutside= 1; break; } 
          isOutside = 0;
          
          xyzl[0] = XYZ[0];
          xi = GetInterval( xyzl[0] , X , I , xi );
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

          xyzl[1] = XYZ[1];
          yj = GetInterval( xyzl[1] , Y , J , yj );
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

          xyzl[2] = XYZ[2];
          zk = GetInterval( xyzl[2] , Z , K , zk );
          if( zk >= 0 ){
            kk0 = zk;             kk1 = zk+1;
            Z0  = Z[kk0];         Z1 = Z[kk1];
          } else if( zk == - 2  ) {
            kk0 = -1;      kk1 = 0;
            Z0  = OZ;      Z1  = Z[0];
          } else if( zk == -101 ) {
            kk0 = K-1;                   kk1 = -1;
            Z0  = Z[K-1];                Z1  = LZ;
          }

          break;
      }

      if( isOutside ){
        if( FLOW ){
          for( ; tt1 < S ; tt1++ ){
            Fx(p,tt1) = XYZ[0];   Fy(p,tt1) = XYZ[1];   Fz(p,tt1) = XYZ[2];
            if( computeJAC ){
              FJx_x(p,tt1) = Jx_x;    FJx_y(p,tt1) = Jx_y;    FJx_z(p,tt1) = Jx_z;
              FJy_x(p,tt1) = Jy_x;    FJy_y(p,tt1) = Jy_y;    FJy_z(p,tt1) = Jy_z;
              FJz_x(p,tt1) = Jz_x;    FJz_y(p,tt1) = Jz_y;    FJz_z(p,tt1) = Jz_z;
            }
          }
        }
        break;
      }

      D = 1/( ( X1 - X0 )*( Y1 - Y0 )*( Z1 - Z0 ) );
      
      VX0000 = VX( ii0 , jj0 , kk0 , tt0);
      VX0010 = VX( ii0 , jj0 , kk1 , tt0);
      VX0100 = VX( ii0 , jj1 , kk0 , tt0);
      VX0110 = VX( ii0 , jj1 , kk1 , tt0);
      VX1000 = VX( ii1 , jj0 , kk0 , tt0);
      VX1010 = VX( ii1 , jj0 , kk1 , tt0);
      VX1100 = VX( ii1 , jj1 , kk0 , tt0);
      VX1110 = VX( ii1 , jj1 , kk1 , tt0);

      VY0000 = VY( ii0 , jj0 , kk0 , tt0);
      VY0010 = VY( ii0 , jj0 , kk1 , tt0);
      VY0100 = VY( ii0 , jj1 , kk0 , tt0);
      VY0110 = VY( ii0 , jj1 , kk1 , tt0);
      VY1000 = VY( ii1 , jj0 , kk0 , tt0);
      VY1010 = VY( ii1 , jj0 , kk1 , tt0);
      VY1100 = VY( ii1 , jj1 , kk0 , tt0);
      VY1110 = VY( ii1 , jj1 , kk1 , tt0);

      VZ0000 = VZ( ii0 , jj0 , kk0 , tt0);
      VZ0010 = VZ( ii0 , jj0 , kk1 , tt0);
      VZ0100 = VZ( ii0 , jj1 , kk0 , tt0);
      VZ0110 = VZ( ii0 , jj1 , kk1 , tt0);
      VZ1000 = VZ( ii1 , jj0 , kk0 , tt0);
      VZ1010 = VZ( ii1 , jj0 , kk1 , tt0);
      VZ1100 = VZ( ii1 , jj1 , kk0 , tt0);
      VZ1110 = VZ( ii1 , jj1 , kk1 , tt0);
      

      VX0001 = VX( ii0 , jj0 , kk0 , tt1);
      VX0011 = VX( ii0 , jj0 , kk1 , tt1);
      VX0101 = VX( ii0 , jj1 , kk0 , tt1);
      VX0111 = VX( ii0 , jj1 , kk1 , tt1);
      VX1001 = VX( ii1 , jj0 , kk0 , tt1);
      VX1011 = VX( ii1 , jj0 , kk1 , tt1);
      VX1101 = VX( ii1 , jj1 , kk0 , tt1);
      VX1111 = VX( ii1 , jj1 , kk1 , tt1);

      VY0001 = VY( ii0 , jj0 , kk0 , tt1);
      VY0011 = VY( ii0 , jj0 , kk1 , tt1);
      VY0101 = VY( ii0 , jj1 , kk0 , tt1);
      VY0111 = VY( ii0 , jj1 , kk1 , tt1);
      VY1001 = VY( ii1 , jj0 , kk0 , tt1);
      VY1011 = VY( ii1 , jj0 , kk1 , tt1);
      VY1101 = VY( ii1 , jj1 , kk0 , tt1);
      VY1111 = VY( ii1 , jj1 , kk1 , tt1);

      VZ0001 = VZ( ii0 , jj0 , kk0 , tt1);
      VZ0011 = VZ( ii0 , jj0 , kk1 , tt1);
      VZ0101 = VZ( ii0 , jj1 , kk0 , tt1);
      VZ0111 = VZ( ii0 , jj1 , kk1 , tt1);
      VZ1001 = VZ( ii1 , jj0 , kk0 , tt1);
      VZ1011 = VZ( ii1 , jj0 , kk1 , tt1);
      VZ1101 = VZ( ii1 , jj1 , kk0 , tt1);
      VZ1111 = VZ( ii1 , jj1 , kk1 , tt1);
      
      while( xyzl[0] >= X0 && xyzl[0] <= X1 &&
             xyzl[1] >= Y0 && xyzl[1] <= Y1 &&
             xyzl[2] >= Z0 && xyzl[2] <= Z1 &&
             tau < T[tt1] && tau < T1 ) {

        U0 = X0 - xyzl[0];          U1 = X1 - xyzl[0];
        V0 = Y0 - xyzl[1];          V1 = Y1 - xyzl[1];
        W0 = Z0 - xyzl[2];          W1 = Z1 - xyzl[2];

        v0[0] = ( U1*V1*W1*VX0000 - U1*V1*W0*VX0010
                + U1*V0*W0*VX0110 - U1*V0*W1*VX0100
                + U0*V1*W0*VX1010 - U0*V1*W1*VX1000
                + U0*V0*W1*VX1100 - U0*V0*W0*VX1110 );

        v0[1] = ( U1*V1*W1*VY0000 - U1*V1*W0*VY0010
                + U1*V0*W0*VY0110 - U1*V0*W1*VY0100
                + U0*V1*W0*VY1010 - U0*V1*W1*VY1000
                + U0*V0*W1*VY1100 - U0*V0*W0*VY1110 );

        v0[2] = ( U1*V1*W1*VZ0000 - U1*V1*W0*VZ0010
                + U1*V0*W0*VZ0110 - U1*V0*W1*VZ0100
                + U0*V1*W0*VZ1010 - U0*V1*W1*VZ1000
                + U0*V0*W1*VZ1100 - U0*V0*W0*VZ1110 );

        
        v1[0] = ( U1*V1*W1*VX0001 - U1*V1*W0*VX0011
                + U1*V0*W0*VX0111 - U1*V0*W1*VX0101
                + U0*V1*W0*VX1011 - U0*V1*W1*VX1001
                + U0*V0*W1*VX1101 - U0*V0*W0*VX1111 );

        v1[1] = ( U1*V1*W1*VY0001 - U1*V1*W0*VY0011
                + U1*V0*W0*VY0111 - U1*V0*W1*VY0101
                + U0*V1*W0*VY1011 - U0*V1*W1*VY1001
                + U0*V0*W1*VY1101 - U0*V0*W0*VY1111 );

        v1[2] = ( U1*V1*W1*VZ0001 - U1*V1*W0*VZ0011
                + U1*V0*W0*VZ0111 - U1*V0*W1*VZ0101
                + U0*V1*W0*VZ1011 - U0*V1*W1*VZ1001
                + U0*V0*W1*VZ1101 - U0*V0*W0*VZ1111 );

        
        v[0] = ( ( T[tt1]-tau )*v0[0] + (tau-T[tt0])*v1[0] ) /( T[tt1]-T[tt0] )*D;
        v[1] = ( ( T[tt1]-tau )*v0[1] + (tau-T[tt0])*v1[1] ) /( T[tt1]-T[tt0] )*D;
        v[2] = ( ( T[tt1]-tau )*v0[2] + (tau-T[tt0])*v1[2] ) /( T[tt1]-T[tt0] )*D;
        
        
        if ( tau + dt > T[tt1] ){ 
          deltaT = T[tt1]-tau; 
        } else {
          deltaT = dt;
        }
        
        v[0] *= deltaT;
        v[1] *= deltaT;
        v[2] *= deltaT;
        
        XYZ[0] += v[0];
        XYZ[1] += v[1];
        XYZ[2] += v[2];

        xyzl[0] += v[0];
        xyzl[1] += v[1];
        xyzl[2] += v[2];

        if( computeJAC ){
          Vx_x = ( ( V1*( W0*(VX0010-VX1010) + W1*(VX1000-VX0000)) + V0*( W0*(VX1110-VX0110) + W1*(VX0100-VX1100)) )*( T[tt1] - tau ) +
                   ( V1*( W0*(VX0011-VX1011) + W1*(VX1001-VX0001)) + V0*( W0*(VX1111-VX0111) + W1*(VX0101-VX1101)) )*( tau - T[tt0] )
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vy_x = ( ( V1*( W0*(VY0010-VY1010) + W1*(VY1000-VY0000)) + V0*( W0*(VY1110-VY0110) + W1*(VY0100-VY1100)) )*( T[tt1] - tau ) +
                   ( V1*( W0*(VY0011-VY1011) + W1*(VY1001-VY0001)) + V0*( W0*(VY1111-VY0111) + W1*(VY0101-VY1101)) )*( tau - T[tt0] )
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vz_x = ( ( V1*( W0*(VZ0010-VZ1010) + W1*(VZ1000-VZ0000)) + V0*( W0*(VZ1110-VZ0110) + W1*(VZ0100-VZ1100)) )*( T[tt1] - tau ) +
                   ( V1*( W0*(VZ0011-VZ1011) + W1*(VZ1001-VZ0001)) + V0*( W0*(VZ1111-VZ0111) + W1*(VZ0101-VZ1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          

          Vx_y = ( ( U1*( W0*(VX0010-VX0110) + W1*(VX0100-VX0000)) + U0*( W0*(VX1110-VX1010) + W1*(VX1000-VX1100)) )*( T[tt1] - tau ) +
                   ( U1*( W0*(VX0011-VX0111) + W1*(VX0101-VX0001)) + U0*( W0*(VX1111-VX1011) + W1*(VX1001-VX1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vy_y = ( ( U1*( W0*(VY0010-VY0110) + W1*(VY0100-VY0000)) + U0*( W0*(VY1110-VY1010) + W1*(VY1000-VY1100)) )*( T[tt1] - tau ) +
                   ( U1*( W0*(VY0011-VY0111) + W1*(VY0101-VY0001)) + U0*( W0*(VY1111-VY1011) + W1*(VY1001-VY1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vz_y = ( ( U1*( W0*(VZ0010-VZ0110) + W1*(VZ0100-VZ0000)) + U0*( W0*(VZ1110-VZ1010) + W1*(VZ1000-VZ1100)) )*( T[tt1] - tau ) +
                   ( U1*( W0*(VZ0011-VZ0111) + W1*(VZ0101-VZ0001)) + U0*( W0*(VZ1111-VZ1011) + W1*(VZ1001-VZ1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;

          
          Vx_z = ( ( U1*( V0*(VX0100-VX0110) + V1*(VX0010-VX0000)) + U0*( V1*(VX1000-VX1010) + V0*(VX1110-VX1100)) )*( T[tt1] - tau ) +
                   ( U1*( V0*(VX0101-VX0111) + V1*(VX0011-VX0001)) + U0*( V1*(VX1001-VX1011) + V0*(VX1111-VX1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vy_z = ( ( U1*( V0*(VY0100-VY0110) + V1*(VY0010-VY0000)) + U0*( V1*(VY1000-VY1010) + V0*(VY1110-VY1100)) )*( T[tt1] - tau ) +
                   ( U1*( V0*(VY0101-VY0111) + V1*(VY0011-VY0001)) + U0*( V1*(VY1001-VY1011) + V0*(VY1111-VY1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;
          Vz_z = ( ( U1*( V0*(VZ0100-VZ0110) + V1*(VZ0010-VZ0000)) + U0*( V1*(VZ1000-VZ1010) + V0*(VZ1110-VZ1100)) )*( T[tt1] - tau ) +
                   ( U1*( V0*(VZ0101-VZ0111) + V1*(VZ0011-VZ0001)) + U0*( V1*(VZ1001-VZ1011) + V0*(VZ1111-VZ1101)) )*( tau - T[tt0] ) 
                 )/( T[tt1]-T[tt0] )*D*deltaT;

          deltaJx_x = Vx_x*Jx_x + Vx_y*Jy_x + Vx_z*Jz_x ;
          deltaJy_x = Vy_x*Jx_x + Vy_y*Jy_x + Vy_z*Jz_x ;
          deltaJz_x = Vz_x*Jx_x + Vz_y*Jy_x + Vz_z*Jz_x ;

          deltaJx_y = Vx_x*Jx_y + Vx_y*Jy_y + Vx_z*Jz_y ;
          deltaJy_y = Vy_x*Jx_y + Vy_y*Jy_y + Vy_z*Jz_y ;
          deltaJz_y = Vz_x*Jx_y + Vz_y*Jy_y + Vz_z*Jz_y ;

          deltaJx_z = Vx_x*Jx_z + Vx_y*Jy_z + Vx_z*Jz_z ;
          deltaJy_z = Vy_x*Jx_z + Vy_y*Jy_z + Vy_z*Jz_z ;
          deltaJz_z = Vz_x*Jx_z + Vz_y*Jy_z + Vz_z*Jz_z ;
          
          Jx_x += deltaJx_x;    Jx_y += deltaJx_y;   Jx_z += deltaJx_z;
          Jy_x += deltaJy_x;    Jy_y += deltaJy_y;   Jy_z += deltaJy_z;
          Jz_x += deltaJz_x;    Jz_y += deltaJz_y;   Jz_z += deltaJz_z;
        }        
        
        tau += deltaT;

        if( FLOW && tau >= T[tt1] ){
          Fx(p,tt1) = XYZ[0];   Fy(p,tt1) = XYZ[1];   Fz(p,tt1) = XYZ[2];
          if( computeJAC ){
            FJx_x(p,tt1) = Jx_x;    FJx_y(p,tt1) = Jx_y;    FJx_z(p,tt1) = Jx_z;
            FJy_x(p,tt1) = Jy_x;    FJy_y(p,tt1) = Jy_y;    FJy_z(p,tt1) = Jy_z;
            FJz_x(p,tt1) = Jz_x;    FJz_y(p,tt1) = Jz_y;    FJz_z(p,tt1) = Jz_z;
          }
        }
        
      }
    }
    
    if( !FLOW ){
      Fx(p,0) = XYZ[0];   Fy(p,0) = XYZ[1];   Fz(p,0) = XYZ[2];
      if( computeJAC ){
        FJx_x(p,0) = Jx_x;    FJx_y(p,0) = Jx_y;    FJx_z(p,0) = Jx_z;
        FJy_x(p,0) = Jy_x;    FJy_y(p,0) = Jy_y;    FJy_z(p,0) = Jy_z;
        FJz_x(p,0) = Jz_x;    FJz_y(p,0) = Jz_y;    FJz_z(p,0) = Jz_z;
      }
    }
  }
  
  if( ! FLOW ){ S = 1; }
  if( nlhs > 2 ){
    for( tt0 = 0 ; tt0 < S ; tt0++ ){
      for( p=0 ; p<nP ; p++ ){
        DetJ[ p + tt0*nP ] = FJx_x(p,tt0)*( FJy_y(p,tt0)*FJz_z(p,tt0) - FJy_z(p,tt0)*FJz_y(p,tt0) ) -
                             FJy_x(p,tt0)*( FJx_y(p,tt0)*FJz_z(p,tt0) - FJx_z(p,tt0)*FJz_y(p,tt0) ) +
                             FJz_x(p,tt0)*( FJx_y(p,tt0)*FJy_z(p,tt0) - FJx_z(p,tt0)*FJy_y(p,tt0) ) ;
      }
    }
  }

  EXIT: myFreeALLOCATES();
}
