/*
  o = 3DGridGradientOnPoints( IM , oGx , oGy , oGz , POINTS , 
                ['om' 'omatrix'] , original_grid_matrix     ( DEF: eye(4)    )
                ['nm' 'nmatrix'] , points_matrix            ( DEF: eye(4)    )
                BOUNDARY_MODE , [ boundary_size ]        ( DEF: 'value'   )
 

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
                          'decay_to_zero' , 'zero' , 'tozero' , 'decay'
                          'decay_to_id'   , 'id'   , 'toid'     (not yet implemented)
                    the data do not decay along the singletons
              after 'decay' option, the boundary size could be specified
              by default... ( LX + LY + LZ )/3
 
*/

#include "myMEX.h"

#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define IM(i,j,k,t)    IM[ (k)*IJ + (j)*I + (i) + (t)*IJK ]

#define Ox(p,t)          O[ (p) + (t)*3*P       ]
#define Oy(p,t)          O[ (p) + (t)*3*P +   P ]
#define Oz(p,t)          O[ (p) + (t)*3*P + 2*P ]

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

real PutInsideSym( real , real );
real PutInside( real , real );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID };

  real  *IM, *O;

  real  *X, *Y, *Z, *xyz, *nM, *oM, ioM[16], MAT[16], DET;
  real  xx,yy,zz,dx,dy,dz,xxnr,yynr,zznr;
  real  LX, LY, LZ, OX, OY, OZ, iDX, iDY, iDZ;

  real  u,uu,v,vv,w,ww;
  
  int     ii0 , ii1 , jj0 , jj1 , kk0 , kk1;
  int     p , t ,ii , jj , kk;
  int     I ,J ,K ,IJ ,IJK ,
          P,
          T;
  int     Odims[3], ndims;

  enum    boundary_modes boundary_mode;

  int     isEqualSpaced;
  int     isOut, edgeX, edgeY, edgeZ;
  
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

  if( myNumel( prhs[1] ) != I ){ mexErrMsgTxt("numel(oGx) Coordinates do not coincide with size(IM,1)."); }
  if( myNumel( prhs[2] ) != J ){ mexErrMsgTxt("numel(oGy) Coordinates do not coincide with size(IM,2)."); }
  if( myNumel( prhs[3] ) != K ){ mexErrMsgTxt("numel(oGz) Coordinates do not coincide with size(IM,3)."); }
  
  X = myGetPr(prhs[1]); 
  if( I == 1 ){
    OX = X[0] - 0.5;                    LX = 1;
  } else {
    OX = X[0] - ( X[1]-X[0] )/2;        LX = X[I-1] + ( X[I-1]-X[I-2] )/2 - OX;
  }
  
  Y = myGetPr(prhs[2]);
  if( J == 1 ){
    OY = Y[0] - 0.5;                    LY = 1;
  } else {
    OY = Y[0] - ( Y[1]-Y[0] )/2;        LY = Y[J-1] + ( Y[J-1]-Y[J-2] )/2 - OY;
  }
  
  Z = myGetPr(prhs[3]);
  if( K == 1 ){
    OZ = Z[0] - 0.5;                    LZ = 1;
  } else {
    OZ = Z[0] - ( Z[1]-Z[0] )/2;        LZ = Z[K-1] + ( Z[K-1]-Z[K-2] )/2 - OZ;
  }

  if( I>1 ){ dx= X(1)-X(0); dx = 1.0/dx; } else { dx=1; }
  if( J>1 ){ dy= Y(1)-Y(0); dy = 1.0/dy; } else { dy=1; }
  if( K>1 ){ dz= Z(1)-Z(0); dz = 1.0/dz; } else { dz=1; }

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
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; mexErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );
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
      mexErrMsgTxt("After the word oMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }

    if( ! myStrcmpi( STR,"nmatrix") || ! myStrcmpi( STR,"nm") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) ){
        if( ~myIsEmpty( prhs[argN] ) ){ nM = myGetPr( prhs[argN] ); }
        argN++; continue;
      }
      mexErrMsgTxt("After the word nMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }
    
    argN++;
//     mexPrintf("%s - ",STR); mexErrMsgTxt("Invalid keyword");
  }
  
//   DISP( boundary_size );
  /*END Parsing arguments*/
  
  if( checkEqualSpaced( X , I ) && checkEqualSpaced( Y , J ) && checkEqualSpaced( Z , K ) ){
    isEqualSpaced = 1;
  } else {
    isEqualSpaced = 0;
  }

//   DISP( P );

  
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

//     DISP( xx );
//     DISP( yy );
//     DISP( zz );

    isOut = 0;
    switch( boundary_mode ){
      case VALUE:
        if( zz < OZ || zz > OZ+LZ ||
            xx < OX || xx > OX+LX ||
            yy < OY || yy > OY+LY ){
          isOut = 1; 
          
//           DISP( OX ); DISP( LX );
//           DISP( OY ); DISP( LY );
//           DISP( OZ ); DISP( LZ );
          
          
//           mexPrintf("aca\n");
        }
        break;
        
      case SYMMETRIC:
        xx = PutInsideSym( xx - OX , LX ) + OX;
        yy = PutInsideSym( yy - OY , LY ) + OY;
        zz = PutInsideSym( zz - OZ , LZ ) + OZ;
        break;
        
      case CIRCULAR:
        xx = PutInside( xx - OX , LX ) + OX;
        yy = PutInside( yy - OY , LY ) + OY;
        zz = PutInside( zz - OZ , LZ ) + OZ;
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

    }
    if( isOut ){
//       mexPrintf( "\n\nesta fuera\n\n" );
      for(t=0;t<T;t++){
        Ox(p,t) = 0;
        Oy(p,t) = 0;
        Oz(p,t) = 0;
      } 
      continue;
    }

//     mexPrintf( "\n\nesta adentro\n\n" );
  
    if( isEqualSpaced ){
      ii = (int) ( (xx-X(0))*dx ); if( ii < 0 || xx<X(0) ){ ii= -2; } else if( ii > I-2 ){ if( xx == X(I-1) ){ ii = I-2; } else { ii = -101; } }
      if( ii < I-2  &&  xx == X(ii+1) ){ ii++; }
      jj = (int) ( (yy-Y(0))*dy ); if( jj < 0 || yy<Y(0) ){ jj= -2; } else if( jj > J-2 ){ if( yy == Y(J-1) ){ jj = J-2; } else { jj = -101; } }
      if( jj < J-2  &&  yy == Y(jj+1) ){ jj++; }
      kk = (int) ( (zz-Z(0))*dz ); if( kk < 0 || zz<Z(0) ){ kk= -2; } else if( kk > K-2 ){ if( zz == Z(K-1) ){ kk = K-2; } else { kk = -101; } }
      if( kk < K-2  &&  zz == Z(kk+1) ){ kk++; }
    } else {
      ii = GetInterval( xx , X , I , ii );
      jj = GetInterval( yy , Y , J , jj );
      kk = GetInterval( zz , Z , K , kk );
    }
    
//     DISP(ii);
//     DISP(jj);
//     DISP(kk);
    
    
    edgeX = 0;
    edgeY = 0;
    edgeZ = 0;
    
    if( I == 1 ){
      u   = 1;  uu  = 0;  ii0 = 0;  ii1 = 0;
    } else if( ii >= 0 ){
      iDX = 1.0/( X(ii+1)-X(ii) );
      u   = ( xx-X(ii) ) * iDX;
      uu  = 1-u;  
      ii0 =  ii;  
      ii1 = ii+1;
      if( xx == X(ii) || xx == X(I-1) ){
        edgeX = 1;
      }
    } else if( ii == -2 ) {
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDX = 1.0/boundary_size;
          u   = ( xx - X(0) + boundary_size ) * iDX;
          uu  = 0; ii0 = 0;ii1 = 0;
          break;
        case SYMMETRIC: u = 0; uu = 1; ii0 = 0; ii1 = 0; break;
        case CIRCULAR:
          iDX = 2.0/( X(1)-X(0)+X(I-1)-X(I-2) );
          uu = ( X(0)-xx ) * iDX;
          u  = 1-uu;  ii0 = I-1;  ii1 = 0;
          break;
      }
    } else if( ii == -101 ){
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDX = 1.0/boundary_size;
          uu  = ( X(I-1) + boundary_size - xx ) * iDX;
          u   = 0; ii0 = I-1; ii1 = I-1;
          break;
        case SYMMETRIC: u = 0; uu = 1; ii0 = I-1; ii1 = I-1; break;
        case CIRCULAR:
          iDX = 2.0/( X(1)-X(0)+X(I-1)-X(I-2) );
          u = ( xx - X(I-1) )* iDX;
          uu = 1-u;  ii0 = I-1;  ii1 = 0;
          break;
      }
    }
    
    if( J == 1 ){
      v   = 1;  vv  = 0;  jj0 = 0;  jj1 = 0;
    } else if( jj >= 0 ){
      iDY = 1.0/( Y(jj+1)-Y(jj) );
      v   = ( yy-Y(jj) ) * iDY ;
      vv  = 1-v;  jj0 =  jj;  jj1 = jj+1;
      if( yy == Y(jj) || yy == Y(J-1) ){
        edgeY = 1;
      }
    } else if( jj == -2 ) {
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDY = 1.0/boundary_size;
          v   = ( yy - Y(0) + boundary_size ) * iDY;
          vv  = 0; jj0 = 0;jj1 = 0;
          break;
        case SYMMETRIC: v = 0; vv = 1; jj0 = 0; jj1 = 0; break;
        case CIRCULAR:
          iDY = 2.0/( Y(1)-Y(0)+Y(J-1)-Y(J-2) );
          vv = ( Y(0)-yy ) * iDY;
          v  = 1-vv;  jj0 = J-1;  jj1 = 0;
          break;
      }
    } else if( jj == -101 ){
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDY = 1.0/boundary_size;
          vv  = ( Y(J-1) + boundary_size - yy ) * iDY;
          v   = 0; jj0 = J-1; jj1 = J-1;
          break;
        case SYMMETRIC: v = 0; vv = 1; jj0 = J-1; jj1 = J-1; break;
        case CIRCULAR:
          iDY = 2.0/( Y(1)-Y(0)+Y(J-1)-Y(J-2) );
          v = ( yy - Y(J-1) ) * iDY;
          vv = 1-v;  jj0 = J-1;  jj1 = 0;
          break;
      }
    }
    
    if( K == 1 ){
      w   = 1;  ww  = 0;  kk0 = 0;  kk1 = 0;
    } else if( kk >= 0 ){
      iDZ = 1.0/( Z(kk+1)-Z(kk) );
      w   = ( zz-Z(kk) ) * iDZ ;
      ww  = 1-w;  kk0 =  kk;  kk1 = kk+1;
      if( zz == Z(kk) || zz == Z(K-1) ){
        edgeZ = 1;
      }
    } else if( kk == -2 ) {
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDZ = 1.0/boundary_size;
          w   = ( zz - Z(0) + boundary_size ) * iDZ;
          ww  = 0; kk0 = 0;kk1 = 0;
          break;
        case SYMMETRIC: w = 0; ww = 1; kk0 = 0; kk1 = 0; break;
        case CIRCULAR:
          iDZ = 2.0/( Z(1)-Z(0)+Z(K-1)-Z(K-2) );
          ww = ( Z(0)-zz ) * iDZ;
          w  = 1-ww;  kk0 = K-1;  kk1 = 0;
          break;
      }
    } else if( kk == -101 ){
      switch(boundary_mode){
        case VALUE: for(t=0;t<T;t++){Ox(p,t)=0;Oy(p,t)=0;Oz(p,t)=0;} continue;
        case DECAY_TO_ZERO:
          iDZ = 1.0/boundary_size;
          ww  = ( Z(K-1) + boundary_size - zz ) * iDZ;
          w   = 0; kk0 = K-1; kk1 = K-1;
          break;
        case SYMMETRIC: w = 0; ww = 1; kk0 = K-1; kk1 = K-1; break;
        case CIRCULAR:
          iDZ = 2.0/( Z(1)-Z(0)+Z(K-1)-Z(K-2) );
          w = ( zz - Z(K-1) ) * iDZ;
          ww = 1-w;  kk0 = K-1;  kk1 = 0;
          break;
      }
    }
    
//     DISP( u );
//     DISP( uu );
//     DISP( v );
//     DISP( vv );
//     DISP( w );
//     DISP( ww );
    
    for( t=0 ; t<T ; t++ ){
//       Ox(p,t)=  IM( ii0 , jj0 , kk0 ,t) * uu * vv * ww +
//                 IM( ii1 , jj0 , kk0 ,t) *  u * vv * ww +
//                 IM( ii0 , jj1 , kk0 ,t) * uu *  v * ww +
//                 IM( ii1 , jj1 , kk0 ,t) *  u *  v * ww +
//                 IM( ii0 , jj0 , kk1 ,t) * uu * vv *  w +
//                 IM( ii1 , jj0 , kk1 ,t) *  u * vv *  w +
//                 IM( ii0 , jj1 , kk1 ,t) * uu *  v *  w +
//                 IM( ii1 , jj1 , kk1 ,t) *  u *  v *  w ;
      
      if( edgeX ){
        if( ii == 0 ){
          xx = ( IM(  1   , jj0 , kk0 ,t) - IM(  0   , jj0 , kk0 ,t) )/( X( 1  ) - X( 0  ) );
        } else if( ii == I-2 ){
          xx = ( IM( I-1  , jj0 , kk0 ,t) - IM( I-2  , jj0 , kk0 ,t) )/( X(I-1 ) - X(I-2 ) );
        } else {
          xx = ( IM( ii+1 , jj0 , kk0 ,t) - IM( ii-1 , jj0 , kk0 ,t) )/( X(ii+1) - X(ii-1) );
        }
        
      } else {
        xx = ( -IM( ii0 , jj0 , kk0 ,t) * vv * ww * (uu==0?0:1) +
                IM( ii1 , jj0 , kk0 ,t) * vv * ww * ( u==0?0:1) +
               -IM( ii0 , jj1 , kk0 ,t) *  v * ww * (uu==0?0:1) +
                IM( ii1 , jj1 , kk0 ,t) *  v * ww * ( u==0?0:1) +
               -IM( ii0 , jj0 , kk1 ,t) * vv *  w * (uu==0?0:1) +
                IM( ii1 , jj0 , kk1 ,t) * vv *  w * ( u==0?0:1) +
               -IM( ii0 , jj1 , kk1 ,t) *  v *  w * (uu==0?0:1) +
                IM( ii1 , jj1 , kk1 ,t) *  v *  w * ( u==0?0:1) )*iDX;
      }


      if( edgeY ){
        if( jj == 0 ){
          yy = ( IM( ii0 ,  1   , kk0 ,t) - IM( ii0 ,  0   , kk0 ,t) )/( Y( 1  ) - Y( 0  ) );
        } else if( jj == J-2 ){
          yy = ( IM( ii0 , J-1  , kk0 ,t) - IM( ii0 , J-2  , kk0 ,t) )/( Y(J-1 ) - Y(J-2 ) );
        } else {
          yy = ( IM( ii0 , jj+1 , kk0 ,t) - IM( ii0 , jj-1 , kk0 ,t) )/( Y(jj+1) - Y(jj-1) );
        }
      } else {
        yy = ( -IM( ii0 , jj0 , kk0 ,t) * uu * ww * (vv==0?0:1) +
               -IM( ii1 , jj0 , kk0 ,t) *  u * ww * (vv==0?0:1) +
                IM( ii0 , jj1 , kk0 ,t) * uu * ww * ( v==0?0:1) +
                IM( ii1 , jj1 , kk0 ,t) *  u * ww * ( v==0?0:1) +
               -IM( ii0 , jj0 , kk1 ,t) * uu *  w * (vv==0?0:1) +
               -IM( ii1 , jj0 , kk1 ,t) *  u *  w * (vv==0?0:1) +
                IM( ii0 , jj1 , kk1 ,t) * uu *  w * ( v==0?0:1) +
                IM( ii1 , jj1 , kk1 ,t) *  u *  w * ( v==0?0:1) )*iDY;
      }
      
      if( edgeZ ){
        if( kk == 0 ){
          zz = ( IM( ii0 , jj0 ,  1   ,t) - IM( ii0 , jj0 ,  0   ,t) )/( Z( 1  ) - Z( 0  ) );
        } else if( kk == K-2 ){
          zz = ( IM( ii0 , jj0 , K-1  ,t) - IM( ii0 , jj0 , K-2  ,t) )/( Z(K-1 ) - Z(K-2 ) );
        } else {
          zz = ( IM( ii0 , jj0 , kk+1 ,t) - IM( ii0 , jj0 , kk-1 ,t) )/( Z(kk+1) - Z(kk-1) );
        }
      } else {
        zz = ( -IM( ii0 , jj0 , kk0 ,t) * uu * vv * (ww==0?0:1) +
               -IM( ii1 , jj0 , kk0 ,t) *  u * vv * (ww==0?0:1) +
               -IM( ii0 , jj1 , kk0 ,t) * uu *  v * (ww==0?0:1) +
               -IM( ii1 , jj1 , kk0 ,t) *  u *  v * (ww==0?0:1) +
                IM( ii0 , jj0 , kk1 ,t) * uu * vv * ( w==0?0:1) +
                IM( ii1 , jj0 , kk1 ,t) *  u * vv * ( w==0?0:1) +
                IM( ii0 , jj1 , kk1 ,t) * uu *  v * ( w==0?0:1) +
                IM( ii1 , jj1 , kk1 ,t) *  u *  v * ( w==0?0:1) )*iDZ;
      }
      
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
  myFreeALLOCATES();
}

real PutInsideSym( real x , real L ){
  int n = (int) (x/L);
  
  if( x > L ){
    if( ! (n%2) ){
      return( x - n*L );
    } else {
      return( -x + (n+1)*L );
    }
  }
  if( x < 0 ){
    if( ! (n%2) ){
      return( - x + n*L );
    } else {
      return( x - (n-1)*L );
    }
  }
  return(x);
}

real PutInside( real x , real L ){
  if( x > L ){
    return( x - ((int) (x/L))*L );
  }
  if( x < 0 ){
    return( x - ( ( (int) (x/L) ) - 1 )*L );
  }
  return(x);
}




