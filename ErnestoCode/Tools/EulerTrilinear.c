/*
  phi = EulerTrilinear( V , Gx , Gy , Gz , points , N , 
                        BOUNDARY_MODE , [ BOUNDARY_SIZE ]
                      );

      if N == []  set N = 500
 
      BOUNDARY_MODE:
           'value'                          stop the evolution outside the grid
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY+LZ)/3


  [ phi , jac ] = EulerTrilinear( V , Gx , Gy , Gz , points , N , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );

    size(jac) = [ size( points ) , 3 ]


  [ phi , jac , DETjac ] = EulerTrilinear( V , Gx , Gy , Gz , points , N , BOUNDARY_MODE , [ BOUNDARY_SIZE ] );

    size(DETjac) = [ numel(points)/3  ,  1 ]

*/

#include "myMEX.h"

#if !defined( real )
  #define   real       double
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define   N_def         500

#define VX(i,j,k)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VX[ (k)*IJ + (j)*I + (i) ] : 0 )
#define VY(i,j,k)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VY[ (k)*IJ + (j)*I + (i) ] : 0 )
#define VZ(i,j,k)    ( ( (i)!=-1 && (j)!=-1 && (k)!=-1 ) ? VZ[ (k)*IJ + (j)*I + (i) ] : 0 )

#define Px(p)         P[ (p)        ]
#define Py(p)         P[ (p) +   nP ]
#define Pz(p)         P[ (p) + 2*nP ]

#define Ox(p)         O[ (p)        ]
#define Oy(p)         O[ (p) +   nP ]
#define Oz(p)         O[ (p) + 2*nP ]

#define OJx_x(p)      OJ[ (p)        ]
#define OJy_x(p)      OJ[ (p) +   nP ]
#define OJz_x(p)      OJ[ (p) + 2*nP ]

#define OJx_y(p)      OJ[ (p) + 3*nP ]
#define OJy_y(p)      OJ[ (p) + 4*nP ]
#define OJz_y(p)      OJ[ (p) + 5*nP ]

#define OJx_z(p)      OJ[ (p) + 6*nP ]
#define OJy_z(p)      OJ[ (p) + 7*nP ]
#define OJz_z(p)      OJ[ (p) + 8*nP ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST };

  real    *X, *Y, *Z, OX, LX, OY, LY, OZ, LZ;
  int     I,J,K,xi,yj,zk, ii0, ii1, jj0, jj1, kk0, kk1;
  int     IJ;

  real    *VX, *VY, *VZ, v[3];
  real    Jx_x, Jx_y, Jx_z, Jy_x, Jy_y, Jy_z, Jz_x, Jz_y, Jz_z;
  real    deltaJx_x, deltaJx_y, deltaJx_z, deltaJy_x, deltaJy_y, deltaJy_z, deltaJz_x, deltaJz_y, deltaJz_z;
  real    Vx_x, Vx_y, Vx_z, Vy_x, Vy_y, Vy_z, Vz_x, Vz_y, Vz_z;
  

  real    *P;
  int     nP, p;
  real    XYZ[3], xyzl[3];

  real    dt;
  int     m, N;
  real    D, U0, U1, V0, V1, W0, W1, X0, X1, Y0, Y1, Z0, Z1;
  real    VX000, VX001, VX010, VX011, VX100, VX101, VX110, VX111;
  real    VY000, VY001, VY010, VY011, VY100, VY101, VY110, VY111;
  real    VZ000, VZ001, VZ010, VZ011, VZ100, VZ101, VZ110, VZ111;

  real    boundary_size;
  enum    boundary_modes boundary_mode;
  
  real    *O, *OJ, *DetJ;
  char    STR[100];
  int     argN;
  int     isOutside;
  int     Odims[50], ndims, d, computeJAC;
  
  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );
  K  = mySize( prhs[0] , 2 );
  IJ = I*J;

  if( myNumel( prhs[0] ) != IJ*K*3 ){
    myErrMsgTxt("size(V) has to be [numel(Gx) numel(Gy) numel(Gz) 3].");
  }
  
  VX = myGetPr( prhs[0] );
  VY = VX + IJ*K;                    
  VZ = VY + IJ*K;                    
  
  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(Gy) Coordinates do not coincide with size(V,2)."); }
  if( myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(Gz) Coordinates do not coincide with size(V,3)."); }
  
  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("oGx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("oGy  Coordinates are not sorted."); }

  Z = myGetPr(prhs[3]);
  if( !checkIsSorted(Z,K) ){ myErrMsgTxt("oGz  Coordinates are not sorted."); }

  nP = myNumel( prhs[4] );
  nP = nP/3;
  P  = myGetPr( prhs[4] );

  if( nrhs > 5 && myIsParameter( prhs[5] ) && !myIsEmpty( prhs[5] ) ){ 
    N = (int) myGetValue( prhs[5] ); 
  } else {
    N = N_def;
  }
  dt = 1.0/N;
  
  
  /*Parsing arguments*/
  /*Defaults*/
  boundary_mode = VALUE;
  boundary_size = (X[I-1]-X[0]+Y[J-1]-Y[0]+Z[K-1]-Z[0])/3;
  
  argN = 6;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );

      if( ! myStrcmpi(STR,"circular") || ! myStrcmpi(STR,"periodic") || 
          ! myStrcmpi(STR,"circ")     || ! myStrcmpi(STR,"per")      ) { boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"zero")     ) { boundary_mode = VALUE;    argN++; continue; }
      if( ! myStrcmpi(STR,"closest")                                 ) { boundary_mode = CLOSEST;  argN++; continue; }

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
  plhs[0]  = myDuplicateArrayWithClass( prhs[4] , mxREAL_CLASS , mxREAL );
  O = (real *)mxGetData( plhs[0] );
  
  if( nlhs > 1 ){
    ndims = myNDims( prhs[4] );
    for( d=0 ; d<ndims ; d++ ){
      Odims[d] = mySize( prhs[4] , d );
    }
    Odims[d++] = 3;
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
  
  xi = -1;  yj = -1;   zk= -1;
  for( p=0 ; p<nP ; p++ ){
    XYZ[0] = Px(p);    XYZ[1] = Py(p);    XYZ[2] = Pz(p);

    if( computeJAC ){
      Jx_x = 1;  Jx_y = 0;  Jx_z = 0;
      Jy_x = 0;  Jy_y = 1;  Jy_z = 0;
      Jz_x = 0;  Jz_y = 0;  Jz_z = 1;
    }
    
    m = 0;
    while( m < N ){
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
          isOutside = 0;
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
      if( isOutside ){ break; }

      
      D = dt/( ( X1 - X0 )*( Y1 - Y0 )*( Z1 - Z0 ) );
      
      VX000 = VX( ii0 , jj0 , kk0 );
      VX001 = VX( ii0 , jj0 , kk1 );
      VX010 = VX( ii0 , jj1 , kk0 );
      VX011 = VX( ii0 , jj1 , kk1 );
      VX100 = VX( ii1 , jj0 , kk0 );
      VX101 = VX( ii1 , jj0 , kk1 );
      VX110 = VX( ii1 , jj1 , kk0 );
      VX111 = VX( ii1 , jj1 , kk1 );

      VY000 = VY( ii0 , jj0 , kk0 );
      VY001 = VY( ii0 , jj0 , kk1 );
      VY010 = VY( ii0 , jj1 , kk0 );
      VY011 = VY( ii0 , jj1 , kk1 );
      VY100 = VY( ii1 , jj0 , kk0 );
      VY101 = VY( ii1 , jj0 , kk1 );
      VY110 = VY( ii1 , jj1 , kk0 );
      VY111 = VY( ii1 , jj1 , kk1 );

      VZ000 = VZ( ii0 , jj0 , kk0 );
      VZ001 = VZ( ii0 , jj0 , kk1 );
      VZ010 = VZ( ii0 , jj1 , kk0 );
      VZ011 = VZ( ii0 , jj1 , kk1 );
      VZ100 = VZ( ii1 , jj0 , kk0 );
      VZ101 = VZ( ii1 , jj0 , kk1 );
      VZ110 = VZ( ii1 , jj1 , kk0 );
      VZ111 = VZ( ii1 , jj1 , kk1 );
      
      while( xyzl[0] >= X0 && xyzl[0] <= X1 &&
             xyzl[1] >= Y0 && xyzl[1] <= Y1 &&
             xyzl[2] >= Z0 && xyzl[2] <= Z1 &&
             m < N  ) {

        U0 = X0 - xyzl[0];          U1 = X1 - xyzl[0];
        V0 = Y0 - xyzl[1];          V1 = Y1 - xyzl[1];
        W0 = Z0 - xyzl[2];          W1 = Z1 - xyzl[2];

        v[0] = ( U1*V1*W1*VX000 - U1*V1*W0*VX001
               + U1*V0*W0*VX011 - U1*V0*W1*VX010
               + U0*V1*W0*VX101 - U0*V1*W1*VX100
               + U0*V0*W1*VX110 - U0*V0*W0*VX111 )*D;

        v[1] = ( U1*V1*W1*VY000 - U1*V1*W0*VY001
               + U1*V0*W0*VY011 - U1*V0*W1*VY010
               + U0*V1*W0*VY101 - U0*V1*W1*VY100
               + U0*V0*W1*VY110 - U0*V0*W0*VY111 )*D;

        v[2] = ( U1*V1*W1*VZ000 - U1*V1*W0*VZ001
               + U1*V0*W0*VZ011 - U1*V0*W1*VZ010
               + U0*V1*W0*VZ101 - U0*V1*W1*VZ100
               + U0*V0*W1*VZ110 - U0*V0*W0*VZ111 )*D;

/* With this break is slower!!!
        if( v[0]==0 && v[1]==0 && v[2]==0 ){
          m = N;
          break;
        }
*/

        XYZ[0] += v[0];
        XYZ[1] += v[1];
        XYZ[2] += v[2];

        xyzl[0] += v[0];
        xyzl[1] += v[1];
        xyzl[2] += v[2];
        
        
        if( computeJAC ){
          
          Vx_x =  ( V1*( W0*(VX001-VX101) + W1*(VX100-VX000)) + V0*( W0*(VX111-VX011) + W1*(VX010-VX110)) )*D;
          Vy_x =  ( V1*( W0*(VY001-VY101) + W1*(VY100-VY000)) + V0*( W0*(VY111-VY011) + W1*(VY010-VY110)) )*D;
          Vz_x =  ( V1*( W0*(VZ001-VZ101) + W1*(VZ100-VZ000)) + V0*( W0*(VZ111-VZ011) + W1*(VZ010-VZ110)) )*D;

          Vx_y =  ( U1*( W0*(VX001-VX011) + W1*(VX010-VX000)) + U0*( W0*(VX111-VX101) + W1*(VX100-VX110)) )*D;
          Vy_y =  ( U1*( W0*(VY001-VY011) + W1*(VY010-VY000)) + U0*( W0*(VY111-VY101) + W1*(VY100-VY110)) )*D;
          Vz_y =  ( U1*( W0*(VZ001-VZ011) + W1*(VZ010-VZ000)) + U0*( W0*(VZ111-VZ101) + W1*(VZ100-VZ110)) )*D;
          
          Vx_z =  ( U1*( V0*(VX010-VX011) + V1*(VX001-VX000)) + U0*( V1*(VX100-VX101) + V0*(VX111-VX110)) )*D;
          Vy_z =  ( U1*( V0*(VY010-VY011) + V1*(VY001-VY000)) + U0*( V1*(VY100-VY101) + V0*(VY111-VY110)) )*D;
          Vz_z =  ( U1*( V0*(VZ010-VZ011) + V1*(VZ001-VZ000)) + U0*( V1*(VZ100-VZ101) + V0*(VZ111-VZ110)) )*D;

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
        
        m++;
      }
    }
    Ox(p) = XYZ[0];    Oy(p) = XYZ[1];    Oz(p) = XYZ[2];

    if( computeJAC ){
      OJx_x(p) = Jx_x;   OJx_y(p) = Jx_y;   OJx_z(p) = Jx_z;
      OJy_x(p) = Jy_x;   OJy_y(p) = Jy_y;   OJy_z(p) = Jy_z;
      OJz_x(p) = Jz_x;   OJz_y(p) = Jz_y;   OJz_z(p) = Jz_z;
    }
  }
  
  if( nlhs > 2 ){
    for( p=0 ; p<nP ; p++ ){
      DetJ[p] = OJx_x(p)*( OJy_y(p)*OJz_z(p) - OJy_z(p)*OJz_y(p) ) -
                OJy_x(p)*( OJx_y(p)*OJz_z(p) - OJx_z(p)*OJz_y(p) ) +
                OJz_x(p)*( OJx_y(p)*OJy_z(p) - OJx_z(p)*OJy_y(p) ) ;
    }
  }

  EXIT: myFreeALLOCATES();
}

