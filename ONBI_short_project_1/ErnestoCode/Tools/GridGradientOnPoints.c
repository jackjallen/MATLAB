/*
  o = 3DGridGradientOnPoints( IM , oGx , oGy , oGz , POINTS , 
                ['om' 'omatrix'] , original_grid_matrix     ( DEF: eye(4)    )
                ['nm' 'nmatrix'] , points_matrix            ( DEF: eye(4)    )
                BOUNDARY_MODE , [ boundary_size ]           ( DEF: 'value'   )
 

          IM : a volume of data
          oGx, oGy, oGz : the grid where the data are defined (the original one)
          POINTS : the points at where interpolate the data
            size( IM ) = [ numel(oGx)  numel(oGy)  numel(oGz)    times   1 ]
            size(  o ) = [ numel(POINTS)/3  , 3 , times ]


          original_grid_matrix : a 4x4 Transformation Matrix

          points_matrix        : a 4x4 Transformation Matrix

  
          BOUNDARY_MODE : 'value' 'extrapolation_value'
                          'sym[metric]' 
                          'circ[ular]' , 'per[iodic]'
                          'closest'
                          'decay_to_zero' , 'zero' , 'tozero' , 'decay'
                          'decay_to_id'   , 'id'   , 'toid'     (not yet implemented)
                    the data do not decay along the singletons
              after 'decay' option, the boundary size could be specified
              by default... ( LX + LY + LZ )/3
 
*/

/*


I = I3D( zeros(10,20,5,4,2) );
 
d_vx_x = 1;
d_vx_y = 10;
d_vx_z = pi;
d_vy_x = 100;
d_vy_y = 200;
d_vy_z = -exp(1);
I.data(:,:,:,:,1) = repmat(  d_vx_x*I.XX + d_vx_y*I.YY + d_vx_z*I.ZZ , [1 1 1 4] );
I.data(:,:,:,:,2) = repmat(  d_vy_x*I.XX + d_vy_y*I.YY + d_vy_z*I.ZZ , [1 1 1 4] );

G = gradient( I );
squeeze( G(1,1,1,1,:,:) )
 
*/


#include "myMEX.h"

#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define IM(i,j,k,t)    IM[ (k)*IJ + (j)*I + (i) + (t)*IJK ]

// #define Ox(p,t)          O[ (p) + (t)*3*P       ]
// #define Oy(p,t)          O[ (p) + (t)*3*P +   P ]
// #define Oz(p,t)          O[ (p) + (t)*3*P + 2*P ]


#define Ox(p,t)          O[ (p) + (t)*P         ]
#define Oy(p,t)          O[ (p) + (t)*P +   P*T ]
#define Oz(p,t)          O[ (p) + (t)*P + 2*P*T ]


#define X(i)            X[ (i) ]
#define Y(j)            Y[ (j) ]
#define Z(k)            Z[ (k) ]

#define x(p)            xyz[ (p) ]
#define y(p)            xyz[ (p) +  P  ]
#define z(p)            xyz[ (p) + 2*P ]

#define nM(i,j)        nM[ 4*(j) + (i) ]
#define oM(i,j)        oM[ 4*(j) + (i) ]
#define ioM(i,j)      ioM[ 4*(j) + (i) ]
#define MAT(i,j)      MAT[ 4*(j) + (i) ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID , CLOSEST };
  enum    cases { ZERO , FORWARD , CENTER , TRILINEAR , DECAY0 , DECAY0_0 , DECAY0_1 , DECAY1 , DECAY1_0 , DECAY1_1 , BORDER };

  real  *IM, *O;

  real  *X, *Y, *Z, *xyz, *nM, *oM, ioM[16], MAT[16], DET;
  real  xx,yy,zz,xxnr,yynr,zznr;
  real  LX, LY, LZ, OX, OY, OZ;

  real  u,uu,v,vv,w,ww;
  real  x_1 , x0 , x1 , y_1 , y0 , y1 , z_1 , z0 , z1 ;
  
  int     ii0 , ii1 , ii_1 , jj0 , jj1 , jj_1 , kk0 , kk1 , kk_1 ;
  int     p , t ,ii , jj , kk;
  int     I ,J ,K ,IJ ,IJK ,
          P,
          T;
  int     Odims[3], ndims;

  enum    boundary_modes boundary_mode;
  enum    boundary_modes caseX, caseY, caseZ;

  int     isOut, flipX=0, flipY=0, flipZ=0;
  
  real    boundary_size;

  char    STR[100];
  int     argN;

  IM= myGetPr( prhs[0] );
  I = mySize( prhs[0] , 0 );
  J = mySize( prhs[0] , 1 );
  K = mySize( prhs[0] , 2 );
  IJ= I*J;
  IJK= IJ*K;

  ndims = myNDims( prhs[0] );
  T = 1;
  for( t=3 ; t<ndims ; t++ ){
    T *= *( mxGetDimensions(prhs[0]) + t );
  }  

  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(oGx) Coordinates do not coincide with size(IM,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(oGy) Coordinates do not coincide with size(IM,2)."); }
  if( myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(oGz) Coordinates do not coincide with size(IM,3)."); }

  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("oGx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("oGy  Coordinates are not sorted."); }

  Z = myGetPr(prhs[3]);
  if( !checkIsSorted(Z,K) ){ myErrMsgTxt("oGz  Coordinates are not sorted."); }

  if( I == 1 ){
    OX = X[0] - 0.5;                    LX = 1;
  } else {
    OX = X[0] - ( X[1]-X[0] )/2.0;        LX = X[I-1] + ( X[I-1]-X[I-2] )/2.0 - OX;
  }
  
  if( J == 1 ){
    OY = Y[0] - 0.5;                    LY = 1;
  } else {
    OY = Y[0] - ( Y[1]-Y[0] )/2.0;        LY = Y[J-1] + ( Y[J-1]-Y[J-2] )/2.0 - OY;
  }
  
  if( K == 1 ){
    OZ = Z[0] - 0.5;                    LZ = 1;
  } else {
    OZ = Z[0] - ( Z[1]-Z[0] )/2.0;        LZ = Z[K-1] + ( Z[K-1]-Z[K-2] )/2.0 - OZ;
  }

  xyz = myGetPr( prhs[4] );
  P   = myNumel( prhs[4] )/3;

  /*Parsing arguments*/
  /*Defaults*/
  boundary_mode = VALUE;
  boundary_size = (LX+LY+LZ)/3;
  oM = NULL;
  nM = NULL;
  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );
    if( ! myStrcmpi(STR,"closest")                             ){ boundary_mode = CLOSEST;   argN++; continue; }
    if( ! myStrcmpi(STR,"symmetric") || ! myStrcmpi(STR,"sym") ){ boundary_mode = SYMMETRIC; argN++; continue; }
    if( ! myStrcmpi(STR,"circular")  || 
        ! myStrcmpi(STR,"periodic") 
         || ! myStrcmpi(STR,"circ") || ! myStrcmpi(STR,"per") ) { boundary_mode = CIRCULAR; argN++; continue; }
    if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"extrapolation_value") )
                                                          { boundary_mode = VALUE;    argN++; continue; }

    if( ! myStrcmpi( STR,"decay_to_zero") || ! myStrcmpi( STR,"zero") 
            || ! myStrcmpi( STR,"tozero") || ! myStrcmpi( STR,"decay") ){
      boundary_mode = DECAY_TO_ZERO; argN++;
      if( nrhs > argN && myIsParameter( prhs[argN] ) ){
        if( !myIsEmpty( prhs[argN] ) ){ boundary_size = myGetValue( prhs[argN] ); }
        argN++; continue;
      }
      continue;
    }

    if( ! myStrcmpi( STR,"omatrix") || ! myStrcmpi( STR,"om") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) ){
        if( !myIsEmpty( prhs[argN] ) ){ oM = myGetPr( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word oMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }

    if( ! myStrcmpi( STR,"nmatrix") || ! myStrcmpi( STR,"nm") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) ){
        if( ~myIsEmpty( prhs[argN] ) ){ nM = myGetPr( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word nMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }
    
    argN++;
  }
  /*END Parsing arguments*/
  
  if( oM != NULL ){
    DET = oM(0,2)*oM(1,1)*oM(2,0) - oM(0,1)*oM(1,2)*oM(2,0) - oM(0,2)*oM(1,0)*oM(2,1) + oM(0,0)*oM(1,2)*oM(2,1) + oM(0,1)*oM(1,0)*oM(2,2) - oM(0,0)*oM(1,1)*oM(2,2);
    
    ioM(0,0)= (oM(1,2)*oM(2,1) - oM(1,1)*oM(2,2))/DET;
    ioM(0,1)= (oM(0,1)*oM(2,2) - oM(0,2)*oM(2,1))/DET;
    ioM(0,2)= (oM(0,2)*oM(1,1) - oM(0,1)*oM(1,2))/DET;
    ioM(0,3)= (oM(0,2)*oM(1,3)*oM(2,1) + oM(0,3)*oM(1,1)*oM(2,2) + oM(0,1)*oM(1,2)*oM(2,3) - oM(0,3)*oM(1,2)*oM(2,1) - oM(0,1)*oM(1,3)*oM(2,2) - oM(0,2)*oM(1,1)*oM(2,3))/DET;

    ioM(1,0)= (oM(1,0)*oM(2,2) - oM(1,2)*oM(2,0))/DET;
    ioM(1,1)= (oM(0,2)*oM(2,0) - oM(0,0)*oM(2,2))/DET;
    ioM(1,2)= (oM(0,0)*oM(1,2) - oM(0,2)*oM(1,0))/DET;
    ioM(1,3)= (oM(0,3)*oM(1,2)*oM(2,0) + oM(0,0)*oM(1,3)*oM(2,2) + oM(0,2)*oM(1,0)*oM(2,3) - oM(0,0)*oM(1,2)*oM(2,3) - oM(0,2)*oM(1,3)*oM(2,0) - oM(0,3)*oM(1,0)*oM(2,2))/DET;

    ioM(2,0)= (oM(1,1)*oM(2,0) - oM(1,0)*oM(2,1))/DET;
    ioM(2,1)= (oM(0,0)*oM(2,1) - oM(0,1)*oM(2,0))/DET;
    ioM(2,2)= (oM(0,1)*oM(1,0) - oM(0,0)*oM(1,1))/DET;
    ioM(2,3)= (oM(0,1)*oM(1,3)*oM(2,0) + oM(0,3)*oM(1,0)*oM(2,1) + oM(0,0)*oM(1,1)*oM(2,3) - oM(0,3)*oM(1,1)*oM(2,0) - oM(0,0)*oM(1,3)*oM(2,1) - oM(0,1)*oM(1,0)*oM(2,3))/DET;

    ioM(3,0)=0;    ioM(3,1)=0;    ioM(3,2)=0;    ioM(3,3)=1;
  }
  
  if( nM == NULL && oM == NULL ){
    MAT(3,3) = 0;
  } else if( nM == NULL && oM != NULL ){
    for(ii=0; ii<16; ii++){ MAT[ii] = ioM[ii]; }
  } else if( nM != NULL && oM == NULL ){
    for(ii=0; ii<16; ii++){ MAT[ii] = nM[ii]; }
  } else if( nM != NULL && oM != NULL ){
    for( ii = 0; ii<4 ; ii++ ){
      for( jj = 0; jj<4 ; jj++ ){
        MAT(ii,jj) = ioM(ii,0)*nM(0,jj) +
                     ioM(ii,1)*nM(1,jj) + 
                     ioM(ii,2)*nM(2,jj) + 
                     ioM(ii,3)*nM(3,jj);
      }
    }
  }
  
  
  
  
  /*Creating output*/
  Odims[0]= P;
  Odims[1]= 3;
  if( T > 1 ){
    Odims[2]= T;
    plhs[0] = mxCreateNumericArray( 3 , Odims , mxREAL_CLASS , mxREAL );
  } else {
    plhs[0] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
  }
  O = (real*)mxGetData( plhs[0] );
  /*END Creating output*/
  
  ii = -1; jj = -1; kk = -1;
  for( p=0 ; p<P ; p++ ){
    if( MAT(3,3) == 0 ) {
      xx= x(p);
      yy= y(p);
      zz= z(p);
    } else {
      xxnr= x(p);
      yynr= y(p);
      zznr= z(p);

      xx = xxnr*MAT(0,0) + yynr*MAT(0,1) + zznr*MAT(0,2) + MAT(0,3);
      yy = xxnr*MAT(1,0) + yynr*MAT(1,1) + zznr*MAT(1,2) + MAT(1,3);
      zz = xxnr*MAT(2,0) + yynr*MAT(2,1) + zznr*MAT(2,2) + MAT(2,3);
    }

    isOut = 0;
    switch( boundary_mode ){
      case VALUE:
        if( zz < Z(0) || zz > Z(K-1) ||
            xx < X(0) || xx > X(I-1) ||
            yy < Y(0) || yy > Y(J-1) ){
          isOut = 1; 
        }
        break;
        
      case CLOSEST:
        if( ( zz < Z(0) || zz > Z(K-1) ) &&
            ( xx < X(0) || xx > X(I-1) ) &&
            ( yy < Y(0) || yy > Y(J-1) ) ){
          isOut = 1; 
        }
        break;

      case SYMMETRIC:
        xx = PutInside(xx,OX,2.0*LX);
        if( xx > OX+LX ){   flipX = 1; xx = 2.0*(OX+LX) - xx;
        } else {            flipX = 0; 
        }

        yy = PutInside(yy,OY,2.0*LY);
        if( yy > OY+LY ){   flipY = 1; yy = 2.0*(OY+LY) - yy;
        } else {            flipY = 0; 
        }

        zz = PutInside(zz,OZ,2.0*LZ);
        if( zz > OZ+LZ ){   flipZ = 1; zz = 2.0*(OZ+LZ) - zz;
        } else {            flipZ = 0; 
        }

        break;

      case CIRCULAR:
        xx = PutInside(xx,OX,LX);
        yy = PutInside(yy,OY,LY);
        zz = PutInside(zz,OZ,LZ);

        break;

      case DECAY_TO_ZERO:
        if( K > 1 ){
          if( zz < Z(0)-boundary_size || zz > Z(K-1)+boundary_size ){ isOut = 1; break;  }
        } else {
          if( zz != Z(0)                                           ){ isOut = 1; break;  }
        }
          
        if( I > 1 ){
          if( xx < X(0)-boundary_size || xx > X(I-1)+boundary_size ){ isOut = 1; break;  }
        } else {
          if( xx != X(0)                                           ){ isOut = 1; break;  }
        }

        if( J > 1 ){
          if( yy < Y(0)-boundary_size || yy > Y(J-1)+boundary_size ){ isOut = 1; break;  }
        } else {
          if( yy != Y(0)                                           ){ isOut = 1; break;  }
        }
        break;
    }
    if( isOut ){
      for( t=0 ; t<T ; t++ ){
        Ox(p,t) = 0;
        Oy(p,t) = 0;
        Oz(p,t) = 0;
      } 
      continue;
    }

    ii = GetInterval( xx , X , I , ii );
    jj = GetInterval( yy , Y , J , jj );
    kk = GetInterval( zz , Z , K , kk );
    

//  ... -2   ) [.  0   ) [.  1   ) [.  2   ) ... [. I-3  ) [. I-2  ] (   -101  ...
//            *         *         *             *         *         *
//           G0        G1        G2            GI-3      GI-2      GI-1                 

    #define GRID          X
    #define LENGTH        I
    #define caseCOORD     caseX
    #define idx           ii
    #define idx0          ii0
    #define idx1          ii1
    #define idx_1         ii_1
    #define coord         xx
    #define coord0        x0
    #define coord1        x1
    #define coord_1       x_1
    #define local         u
    #define llocal        uu

    if( LENGTH == 1 ){

      idx0  = 0;   coord0 = GRID(0);
      idx1  = 0;   coord1 = GRID(0);
      local    = 1;   llocal = 0;
      caseCOORD = ZERO;

    } else if( boundary_mode == VALUE ){

      if( coord == GRID(0) ){
        idx0   = 0;    coord0  = GRID(0);
        idx1   = 1;    coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CIRCULAR ){

      if( coord < GRID(0) ){
        idx0   = LENGTH-1;  coord0  = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        idx1   = 0;    coord1  = GRID(0);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord > GRID(LENGTH-1) ){
        idx0   = LENGTH-1;  coord0  = GRID(LENGTH-1);
        idx1   = 0;    coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = LENGTH-1; coord_1 = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0;
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;  coord_1 = GRID(LENGTH-2);
        idx1   = 0;         coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == SYMMETRIC ){

      if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = 0;   coord_1 = GRID(0) - ( GRID(1)- GRID(0) );
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;   coord_1 = GRID(LENGTH-2);
        idx1   = LENGTH-1;   coord1  = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord < GRID(0) ){
        idx0  = 0;   coord0 = GRID(0);
        idx1  = 0;   coord1 = GRID(0) - ( GRID(1)- GRID(0) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord > GRID(LENGTH-1) ){
        idx0  = LENGTH-1;   coord0 = GRID(LENGTH-1);
        idx1  = LENGTH-1;   coord1 = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == DECAY_TO_ZERO ){

      if( coord == GRID(0) - boundary_size ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY0_0;
      } else if( coord == GRID(LENGTH-1) + boundary_size ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY1_1;
      } else if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = (coord-coord0)/boundary_size;
        llocal  = 0;
        caseCOORD = DECAY0;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 1 - (coord-coord0)/boundary_size;
        caseCOORD = DECAY1;
      } else if( coord == GRID(0) ){
        idx_1 = -1;     coord_1 = GRID(0) - boundary_size;
        idx0  = 0;      coord0  = GRID(0);
        idx1  = 1;      coord1  = GRID(1);
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY0_1;
      } else if( coord == GRID(LENGTH-1) ){
        idx_1 = LENGTH-2;    coord_1 = GRID(LENGTH-2);
        idx0  = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1  = 0;      coord1  = GRID(LENGTH-1) + boundary_size;
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY1_0;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CLOSEST ){

      if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 0;      coord0  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else  if( coord == GRID(0) ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 1;      coord1  = GRID(1);
        local   = 0;     llocal  = 1;
        caseCOORD = BORDER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = BORDER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    }
    idx = idx0;

        


    #define GRID          Y
    #define LENGTH        J
    #define caseCOORD     caseY
    #define idx           jj
    #define idx0          jj0
    #define idx1          jj1
    #define idx_1         jj_1
    #define coord         yy
    #define coord0        y0
    #define coord1        y1
    #define coord_1       y_1
    #define local         v
    #define llocal        vv

    if( LENGTH == 1 ){

      idx0  = 0;   coord0 = GRID(0);
      idx1  = 0;   coord1 = GRID(0);
      local    = 1;   llocal = 0;
      caseCOORD = ZERO;

    } else if( boundary_mode == VALUE ){

      if( coord == GRID(0) ){
        idx0   = 0;    coord0  = GRID(0);
        idx1   = 1;    coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CIRCULAR ){

      if( coord < GRID(0) ){
        idx0   = LENGTH-1;  coord0  = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        idx1   = 0;    coord1  = GRID(0);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord > GRID(LENGTH-1) ){
        idx0   = LENGTH-1;  coord0  = GRID(LENGTH-1);
        idx1   = 0;    coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = LENGTH-1; coord_1 = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0;
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;  coord_1 = GRID(LENGTH-2);
        idx1   = 0;         coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == SYMMETRIC ){

      if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = 0;   coord_1 = GRID(0) - ( GRID(1)- GRID(0) );
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;   coord_1 = GRID(LENGTH-2);
        idx1   = LENGTH-1;   coord1  = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord < GRID(0) ){
        idx0  = 0;   coord0 = GRID(0);
        idx1  = 0;   coord1 = GRID(0) - ( GRID(1)- GRID(0) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord > GRID(LENGTH-1) ){
        idx0  = LENGTH-1;   coord0 = GRID(LENGTH-1);
        idx1  = LENGTH-1;   coord1 = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == DECAY_TO_ZERO ){

      if( coord == GRID(0) - boundary_size ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY0_0;
      } else if( coord == GRID(LENGTH-1) + boundary_size ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY1_1;
      } else if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = (coord-coord0)/boundary_size;
        llocal  = 0;
        caseCOORD = DECAY0;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 1 - (coord-coord0)/boundary_size;
        caseCOORD = DECAY1;
      } else if( coord == GRID(0) ){
        idx_1 = -1;     coord_1 = GRID(0) - boundary_size;
        idx0  = 0;      coord0  = GRID(0);
        idx1  = 1;      coord1  = GRID(1);
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY0_1;
      } else if( coord == GRID(LENGTH-1) ){
        idx_1 = LENGTH-2;    coord_1 = GRID(LENGTH-2);
        idx0  = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1  = 0;      coord1  = GRID(LENGTH-1) + boundary_size;
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY1_0;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CLOSEST ){

      if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 0;      coord0  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else  if( coord == GRID(0) ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 1;      coord1  = GRID(1);
        local   = 0;     llocal  = 1;
        caseCOORD = BORDER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = BORDER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    }
    idx = idx0;
        


    #define GRID          Z
    #define LENGTH        K
    #define caseCOORD     caseZ
    #define idx           kk
    #define idx0          kk0
    #define idx1          kk1
    #define idx_1         kk_1
    #define coord         zz
    #define coord0        z0
    #define coord1        z1
    #define coord_1       z_1
    #define local         w
    #define llocal        ww

    if( LENGTH == 1 ){

      idx0  = 0;   coord0 = GRID(0);
      idx1  = 0;   coord1 = GRID(0);
      local    = 1;   llocal = 0;
      caseCOORD = ZERO;

    } else if( boundary_mode == VALUE ){

      if( coord == GRID(0) ){
        idx0   = 0;    coord0  = GRID(0);
        idx1   = 1;    coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CIRCULAR ){

      if( coord < GRID(0) ){
        idx0   = LENGTH-1;  coord0  = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        idx1   = 0;    coord1  = GRID(0);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord > GRID(LENGTH-1) ){
        idx0   = LENGTH-1;  coord0  = GRID(LENGTH-1);
        idx1   = 0;    coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      } else if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = LENGTH-1; coord_1 = GRID(0) - ( GRID(1)-GRID(0) )/2.0 - ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0;
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;  coord_1 = GRID(LENGTH-2);
        idx1   = 0;         coord1  = GRID(LENGTH-1) + ( GRID(1)-GRID(0) )/2.0 + ( GRID(LENGTH-1)-GRID(LENGTH-2) )/2.0 ;
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == SYMMETRIC ){

      if( coord == GRID(0) ){
        idx0   = 0;   coord0  = GRID(0);
        idx_1  = 0;   coord_1 = GRID(0) - ( GRID(1)- GRID(0) );
        idx1   = 1;   coord1  = GRID(1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-1;   coord0  = GRID(LENGTH-1);
        idx_1  = LENGTH-2;   coord_1 = GRID(LENGTH-2);
        idx1   = LENGTH-1;   coord1  = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else if( coord < GRID(0) ){
        idx0  = 0;   coord0 = GRID(0);
        idx1  = 0;   coord1 = GRID(0) - ( GRID(1)- GRID(0) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord > GRID(LENGTH-1) ){
        idx0  = LENGTH-1;   coord0 = GRID(LENGTH-1);
        idx1  = LENGTH-1;   coord1 = GRID(LENGTH-1) + ( GRID(LENGTH-1) - GRID(LENGTH-2) );
        local    = 1;   llocal = 0;
        caseCOORD = ZERO;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == DECAY_TO_ZERO ){

      if( coord == GRID(0) - boundary_size ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY0_0;
      } else if( coord == GRID(LENGTH-1) + boundary_size ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 0;
        caseCOORD = DECAY1_1;
      } else if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0) - boundary_size;
        idx1 = 0;      coord1  = GRID(0);
        local   = (coord-coord0)/boundary_size;
        llocal  = 0;
        caseCOORD = DECAY0;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1) + boundary_size;
        local   = 0;
        llocal  = 1 - (coord-coord0)/boundary_size;
        caseCOORD = DECAY1;
      } else if( coord == GRID(0) ){
        idx_1 = -1;     coord_1 = GRID(0) - boundary_size;
        idx0  = 0;      coord0  = GRID(0);
        idx1  = 1;      coord1  = GRID(1);
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY0_1;
      } else if( coord == GRID(LENGTH-1) ){
        idx_1 = LENGTH-2;    coord_1 = GRID(LENGTH-2);
        idx0  = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1  = 0;      coord1  = GRID(LENGTH-1) + boundary_size;
        local    = 0;
        llocal   = 1;
        caseCOORD = DECAY1_0;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    } else if( boundary_mode == CLOSEST ){

      if( idx == -2 ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 0;      coord0  = GRID(0);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else if( idx == -101 ){
        idx0 = LENGTH-1;    coord0  = GRID(LENGTH-1);
        idx1 = LENGTH-1;    coord1  = GRID(LENGTH-1);
        local   = 0;
        llocal  = 0;
        caseCOORD = ZERO;
      } else  if( coord == GRID(0) ){
        idx0 = 0;      coord0  = GRID(0);
        idx1 = 1;      coord1  = GRID(1);
        local   = 0;     llocal  = 1;
        caseCOORD = BORDER;
      } else if( coord == GRID(LENGTH-1) ){
        idx0   = LENGTH-2;  coord0  = GRID(LENGTH-2);
        idx1   = LENGTH-1;  coord1  = GRID(LENGTH-1);
        local     = 1;    llocal  = 0;
        caseCOORD = BORDER;
      } else if( coord == GRID(idx) ){
        idx0   = idx;   coord0  = GRID(idx0);
        idx_1  = idx-1; coord_1 = GRID(idx_1);
        idx1   = idx+1; coord1  = GRID(idx1);
        local     = 0;    llocal  = 1;
        caseCOORD = CENTER;
      } else {
        idx0   = idx;   coord0 = GRID(idx0);
        idx1   = idx+1; coord1 = GRID(idx1);
        local     = (coord-coord0)/(coord1-coord0); llocal = 1-local;
        caseCOORD = TRILINEAR;
      }

    }
    idx = idx0;
        


        



    
    for( t=0 ; t<T ; t++ ){          
      
      switch( caseX ){
        case ZERO:
          xx = 0;
          break;
        case CENTER:
          xx = (( IM( ii1  ,jj0,kk0,t) - IM( ii0  ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii1  ,jj1,kk0,t) - IM( ii0  ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii1  ,jj0,kk1,t) - IM( ii0  ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii1  ,jj1,kk1,t) - IM( ii0  ,jj1,kk1,t) ) *  v *  w  )/( x1  - x0  )/2.0 +
               (( IM( ii0  ,jj0,kk0,t) - IM( ii_1 ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii0  ,jj1,kk0,t) - IM( ii_1 ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii0  ,jj0,kk1,t) - IM( ii_1 ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii0  ,jj1,kk1,t) - IM( ii_1 ,jj1,kk1,t) ) *  v *  w  )/( x0  - x_1 )/2.0;
          break;
        case BORDER:                     
          xx = (( IM( ii1  ,jj0,kk0,t) - IM( ii0  ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii1  ,jj1,kk0,t) - IM( ii0  ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii1  ,jj0,kk1,t) - IM( ii0  ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii1  ,jj1,kk1,t) - IM( ii0  ,jj1,kk1,t) ) *  v *  w  )/( x1  - x0  )/2.0;
          break;
        case TRILINEAR:                     
          xx = (( IM( ii1  ,jj0,kk0,t) - IM( ii0  ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii1  ,jj1,kk0,t) - IM( ii0  ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii1  ,jj0,kk1,t) - IM( ii0  ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii1  ,jj1,kk1,t) - IM( ii0  ,jj1,kk1,t) ) *  v *  w  )/( x1  - x0  );
          break;
        case DECAY0:
          xx =  ( IM( ii1  ,jj0,kk0,t)  * vv * ww  +
                  IM( ii1  ,jj1,kk0,t)  *  v * ww  +
                  IM( ii1  ,jj0,kk1,t)  * vv *  w  +
                  IM( ii1  ,jj1,kk1,t)  *  v *  w  )/( x1  - x0  );
          break;
        case DECAY0_0:
          xx =  ( IM( ii1  ,jj0,kk0,t)  * vv * ww  +
                  IM( ii1  ,jj1,kk0,t)  *  v * ww  +
                  IM( ii1  ,jj0,kk1,t)  * vv *  w  +
                  IM( ii1  ,jj1,kk1,t)  *  v *  w  )/( x1  - x0  )/2.0;
          break;
        case DECAY0_1:
          xx = (( IM( ii1  ,jj0,kk0,t) - IM( ii0  ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii1  ,jj1,kk0,t) - IM( ii0  ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii1  ,jj0,kk1,t) - IM( ii0  ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii1  ,jj1,kk1,t) - IM( ii0  ,jj1,kk1,t) ) *  v *  w  )/( x1  - x0  )/2.0 +
                ( IM( ii0  ,jj0,kk0,t) * vv * ww  +
                  IM( ii0  ,jj1,kk0,t) *  v * ww  +
                  IM( ii0  ,jj0,kk1,t) * vv *  w  +
                  IM( ii0  ,jj1,kk1,t) *  v *  w  )/( x0  - x_1 )/2.0;
          break;
        case DECAY1:
          xx = -( IM( ii0  ,jj0,kk0,t)  * vv * ww  +
                  IM( ii0  ,jj1,kk0,t)  *  v * ww  +
                  IM( ii0  ,jj0,kk1,t)  * vv *  w  +
                  IM( ii0  ,jj1,kk1,t)  *  v *  w  )/( x1  - x0  );
          break;
        case DECAY1_1:
          xx = -( IM( ii0  ,jj0,kk0,t)  * vv * ww  +
                  IM( ii0  ,jj1,kk0,t)  *  v * ww  +
                  IM( ii0  ,jj0,kk1,t)  * vv *  w  +
                  IM( ii0  ,jj1,kk1,t)  *  v *  w  )/( x1  - x0  )/2.0;
          break;
        case DECAY1_0:
          xx = -( IM( ii0  ,jj0,kk0,t) * vv * ww  +
                  IM( ii0  ,jj1,kk0,t) *  v * ww  +
                  IM( ii0  ,jj0,kk1,t) * vv *  w  +
                  IM( ii0  ,jj1,kk1,t) *  v *  w  )/( x1  - x0  )/2.0 +
               (( IM( ii0  ,jj0,kk0,t) - IM( ii_1 ,jj0,kk0,t) ) * vv * ww  +
                ( IM( ii0  ,jj1,kk0,t) - IM( ii_1 ,jj1,kk0,t) ) *  v * ww  +
                ( IM( ii0  ,jj0,kk1,t) - IM( ii_1 ,jj0,kk1,t) ) * vv *  w  +
                ( IM( ii0  ,jj1,kk1,t) - IM( ii_1 ,jj1,kk1,t) ) *  v *  w   )/( x0  - x_1 )/2.0;
          break;
      }
      
      
      switch( caseY ){
        case ZERO:
          yy = 0;
          break;
        case CENTER:   
          yy = (( IM(ii0, jj1  ,kk0,t) - IM(ii0, jj0  ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj1  ,kk0,t) - IM(ii1, jj0  ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj1  ,kk1,t) - IM(ii0, jj0  ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj1  ,kk1,t) - IM(ii1, jj0  ,kk1,t) ) *  u *  w  )/( y1  - y0  )/2.0 +
               (( IM(ii0, jj0  ,kk0,t) - IM(ii0, jj_1 ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj0  ,kk0,t) - IM(ii1, jj_1 ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj0  ,kk1,t) - IM(ii0, jj_1 ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj0  ,kk1,t) - IM(ii1, jj_1 ,kk1,t) ) *  u *  w  )/( y0  - y_1 )/2.0;
          break;
        case BORDER:                     
          yy = (( IM(ii0, jj1  ,kk0,t) - IM(ii0, jj0  ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj1  ,kk0,t) - IM(ii1, jj0  ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj1  ,kk1,t) - IM(ii0, jj0  ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj1  ,kk1,t) - IM(ii1, jj0  ,kk1,t) ) *  u *  w  )/( y1  - y0  )/2.0;
          break;
        case TRILINEAR:
          yy = (( IM(ii0, jj1  ,kk0,t) - IM(ii0, jj0  ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj1  ,kk0,t) - IM(ii1, jj0  ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj1  ,kk1,t) - IM(ii0, jj0  ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj1  ,kk1,t) - IM(ii1, jj0  ,kk1,t) ) *  u *  w  )/( y1  - y0  );
          break;
        case DECAY0:
          yy =  ( IM(ii0, jj1  ,kk0,t)  * uu * ww  +
                  IM(ii1, jj1  ,kk0,t)  *  u * ww  +
                  IM(ii0, jj1  ,kk1,t)  * uu *  w  +
                  IM(ii1, jj1  ,kk1,t)  *  u *  w  )/( y1  - y0  );
          break;
        case DECAY0_0:
          yy =  ( IM(ii0, jj1  ,kk0,t)  * uu * ww  +
                  IM(ii1, jj1  ,kk0,t)  *  u * ww  +
                  IM(ii0, jj1  ,kk1,t)  * uu *  w  +
                  IM(ii1, jj1  ,kk1,t)  *  u *  w  )/( y1  - y0  )/2.0;
          break;
        case DECAY0_1:
          yy = (( IM(ii0, jj1  ,kk0,t) - IM(ii0, jj0  ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj1  ,kk0,t) - IM(ii1, jj0  ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj1  ,kk1,t) - IM(ii0, jj0  ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj1  ,kk1,t) - IM(ii1, jj0  ,kk1,t) ) *  u *  w  )/( y1  - y0  )/2.0 +
                ( IM(ii0, jj0  ,kk0,t) * uu * ww  +
                  IM(ii1, jj0  ,kk0,t) *  u * ww  +
                  IM(ii0, jj0  ,kk1,t) * uu *  w  +
                  IM(ii1, jj0  ,kk1,t) *  u *  w  )/( y0  - y_1 )/2.0;
          break;
        case DECAY1:
          yy = -( IM(ii0, jj0  ,kk0,t)  * uu * ww  +
                  IM(ii1, jj0  ,kk0,t)  *  u * ww  +
                  IM(ii0, jj0  ,kk1,t)  * uu *  w  +
                  IM(ii1, jj0  ,kk1,t)  *  u *  w  )/( y1  - y0  );
          break;
        case DECAY1_1:
          yy = -( IM(ii0, jj0  ,kk0,t)  * uu * ww  +
                  IM(ii1, jj0  ,kk0,t)  *  u * ww  +
                  IM(ii0, jj0  ,kk1,t)  * uu *  w  +
                  IM(ii1, jj0  ,kk1,t)  *  u *  w  )/( y1  - y0  )/2.0;
          break;
        case DECAY1_0:
          yy = -( IM(ii0, jj0  ,kk0,t) * uu * ww  +
                  IM(ii1, jj0  ,kk0,t) *  u * ww  +
                  IM(ii0, jj0  ,kk1,t) * uu *  w  +
                  IM(ii1, jj0  ,kk1,t) *  u *  w  )/( y1  - y0  )/2.0 +
               (( IM(ii0, jj0  ,kk0,t) - IM(ii0, jj_1 ,kk0,t) ) * uu * ww  +
                ( IM(ii1, jj0  ,kk0,t) - IM(ii1, jj_1 ,kk0,t) ) *  u * ww  +
                ( IM(ii0, jj0  ,kk1,t) - IM(ii0, jj_1 ,kk1,t) ) * uu *  w  +
                ( IM(ii1, jj0  ,kk1,t) - IM(ii1, jj_1 ,kk1,t) ) *  u *  w   )/( y0  - y_1 )/2.0;
          break;
      }

     
      switch( caseZ ){
        case ZERO:
          zz = 0;
          break;
        case CENTER:   
          zz = (( IM(ii0,jj0, kk1  ,t) - IM(ii0,jj0, kk0  ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk1  ,t) - IM(ii1,jj0, kk0  ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk1  ,t) - IM(ii0,jj1, kk0  ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk1  ,t) - IM(ii1,jj1, kk0  ,t) ) *  u *  v  )/( z1  - z0  )/2.0 +
               (( IM(ii0,jj0, kk0  ,t) - IM(ii0,jj0, kk_1 ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk0  ,t) - IM(ii1,jj0, kk_1 ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk0  ,t) - IM(ii0,jj1, kk_1 ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk0  ,t) - IM(ii1,jj1, kk_1 ,t) ) *  u *  v  )/( z0  - z_1 )/2.0;
          break;
        case BORDER:
          zz = (( IM(ii0,jj0, kk1  ,t) - IM(ii0,jj0, kk0  ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk1  ,t) - IM(ii1,jj0, kk0  ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk1  ,t) - IM(ii0,jj1, kk0  ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk1  ,t) - IM(ii1,jj1, kk0  ,t) ) *  u *  v  )/( z1  - z0  )/2.0;
          break;
        case TRILINEAR:
          zz = (( IM(ii0,jj0, kk1  ,t) - IM(ii0,jj0, kk0  ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk1  ,t) - IM(ii1,jj0, kk0  ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk1  ,t) - IM(ii0,jj1, kk0  ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk1  ,t) - IM(ii1,jj1, kk0  ,t) ) *  u *  v  )/( z1  - z0  );
          break;
        case DECAY0:
          zz =  ( IM(ii0,jj0, kk1  ,t)  * uu * vv  +
                  IM(ii1,jj0, kk1  ,t)  *  u * vv  +
                  IM(ii0,jj1, kk1  ,t)  * uu *  v  +
                  IM(ii1,jj1, kk1  ,t)  *  u *  v  )/( z1  - z0  );
          break;
        case DECAY0_0:
          zz =  ( IM(ii0,jj0, kk1  ,t)  * uu * vv  +
                  IM(ii1,jj0, kk1  ,t)  *  u * vv  +
                  IM(ii0,jj1, kk1  ,t)  * uu *  v  +
                  IM(ii1,jj1, kk1  ,t)  *  u *  v  )/( z1  - z0  )/2.0;
          break;
        case DECAY0_1:
          zz = (( IM(ii0,jj0, kk1  ,t) - IM(ii0,jj0, kk0  ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk1  ,t) - IM(ii1,jj0, kk0  ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk1  ,t) - IM(ii0,jj1, kk0  ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk1  ,t) - IM(ii1,jj1, kk0  ,t) ) *  u *  v  )/( z1  - z0  )/2.0 +
                ( IM(ii0,jj0, kk0  ,t) * uu * vv  +
                  IM(ii1,jj0, kk0  ,t) *  u * vv  +
                  IM(ii0,jj1, kk0  ,t) * uu *  v  +
                  IM(ii1,jj1, kk0  ,t) *  u *  v  )/( z0  - z_1 )/2.0;
          break;
        case DECAY1:
          zz = -( IM(ii0,jj0, kk0  ,t)  * uu * vv  +
                  IM(ii1,jj0, kk0  ,t)  *  u * vv  +
                  IM(ii0,jj1, kk0  ,t)  * uu *  v  +
                  IM(ii1,jj1, kk0  ,t)  *  u *  v  )/( z1  - z0  );
          break;
        case DECAY1_1:
          zz = -( IM(ii0,jj0, kk0  ,t)  * uu * vv  +
                  IM(ii1,jj0, kk0  ,t)  *  u * vv  +
                  IM(ii0,jj1, kk0  ,t)  * uu *  v  +
                  IM(ii1,jj1, kk0  ,t)  *  u *  v  )/( z1  - z0  )/2.0;
          break;
        case DECAY1_0:
          zz = -( IM(ii0,jj0, kk0  ,t) * uu * vv  +
                  IM(ii1,jj0, kk0  ,t) *  u * vv  +
                  IM(ii0,jj1, kk0  ,t) * uu *  v  +
                  IM(ii1,jj1, kk0  ,t) *  u *  v  )/( z1  - z0  )/2.0 +
               (( IM(ii0,jj0, kk0  ,t) - IM(ii0,jj0, kk_1 ,t) ) * uu * vv  +
                ( IM(ii1,jj0, kk0  ,t) - IM(ii1,jj0, kk_1 ,t) ) *  u * vv  +
                ( IM(ii0,jj1, kk0  ,t) - IM(ii0,jj1, kk_1 ,t) ) * uu *  v  +
                ( IM(ii1,jj1, kk0  ,t) - IM(ii1,jj1, kk_1 ,t) ) *  u *  v   )/( z0  - z_1 )/2.0;
          break;
      }

      if( flipX ){ xx *= -1; }
      if( flipY ){ yy *= -1; }
      if( flipZ ){ zz *= -1; }
      
      if( oM != NULL ){
        Ox(p,t) = xx*ioM(0,0) + yy*ioM(1,0) + zz*ioM(2,0);
        Oy(p,t) = xx*ioM(0,1) + yy*ioM(1,1) + zz*ioM(2,1);
        Oz(p,t) = xx*ioM(0,2) + yy*ioM(1,2) + zz*ioM(2,2);
      } else {
        Ox(p,t) = xx;
        Oy(p,t) = yy;
        Oz(p,t) = zz;
      }        
                  
    }
  }

  EXIT: myFreeALLOCATES();
}
