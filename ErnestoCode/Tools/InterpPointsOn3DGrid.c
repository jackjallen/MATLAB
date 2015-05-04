/*
  o = Interp3DGridOn3DGrid( IM , oGx , oGy , oGz , POINTS , 
                ['om' 'omatrix'] , original_grid_matrix     ( DEF: eye(4)    )
                ['nm' 'nmatrix'] , points_matrix            ( DEF: eye(4)    )
                INTERP_MODE                                 ( DEF: '?linear' )
                BOUNDARY_MODE , [ boundary_size ]           ( DEF: 'value'   )
                ['outside_value' 'outval'], outside_value   ( DEF: NaN       )
 

          IM : a volume of data
          oGx, oGy, oGz : the grid where the data are defined (the original one)
          POINTS : the points at where interpolate the data
            size( IM ) = [ numel(oGx)  numel(oGy)  numel(oGz)    times  c1  c2 ... ck  ]
            size(  o ) = [ numel(POINTS)/3  ,  (times*c1*c2* ... *ck)  ]

          original_grid_matrix : a 4x4 Transformation Matrix

          points_matrix        : a 4x4 Transformation Matrix

          INTERP_MODE :  'lin[ear]'  , '*lin[ear]'  , '?lin[ear]'
                         'nea[rest]' , '*nea[rest]' , '?nea[rest]'
                         'cub[ic]'
                         '*sinc'      (not yet implemented)
                 * -> supose grid equal spaced
                 ? -> control if the grid is equal spaced
  
          BOUNDARY_MODE : 'value' 'extrapolation_value'
                          'sym[metric]' 
                          'circ[ular]' , 'per[iodic]'
                          'closest'
                          'decay_to_zero' , 'zero' , 'tozero' , 'decay'
                          'decay_to_id'   , 'id'   , 'toid'     (not yet implemented)
                    the data do not decay along the singletons
              after 'decay' option, the boundary size could be specified
              by default... ( LX + LY + LZ )/3
 
          outside_value : (real) 
                    
*/

#include "myMEX.h"

#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define IM(i,j,k,t)    IM[ (k)*IJ + (j)*I + (i) + (t)*IJK ]
#define O(p,t)          O[ (p) + (t)*nP ]

#define X(i)            X[ (i) ]
#define Y(j)            Y[ (j) ]
#define Z(k)            Z[ (k) ]

#define x(p)            P[ (p) ]
#define y(p)            P[ (p) +  nP  ]
#define z(p)            P[ (p) + 2*nP ]

#define nM(i,j)        nM[ 4*(j) + (i) ]
#define oM(i,j)        oM[ 4*(j) + (i) ]
#define ioM(i,j)      ioM[ 4*(j) + (i) ]
#define MAT(i,j)      MAT[ 4*(j) + (i) ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)

real CINT( real , real , real , real , real , real , real , real , real );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    interp_modes { LINEAR , NEAREST , CUBIC , SINC };
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID , CLOSEST };

  real  *IM, *O;
  real  *X, *Y, *Z, *P, *nM, *oM, ioM[16], MAT[16], DET, *DX, *DY, *DZ;
  real  xx,yy,zz,dx,dy,dz,xxnr,yynr,zznr;
  real  LX, LY, LZ, OX, OY, OZ;

  real  u,uu,v,vv,w,ww;
  
  int     ii0 , ii1 , jj0 , jj1 , kk0 , kk1;
  int     p , t ,ii , jj , kk;
  int     I ,J ,K ,IJ ,IJK ,
          nP,
          T;
  int     Odims[2], ndims;

  enum    interp_modes interp_mode;
  enum    boundary_modes boundary_mode;

  int     isEqualSpaced;
  int     isOut;
  
  real  boundary_size;
  real  outside_value;

  char    STR[100];
  int     argN;

  IM= myGetPr( prhs[0] );
  I = mySize( prhs[0] , 0 );
  J = mySize( prhs[0] , 1 );
  K = mySize( prhs[0] , 2 );
  IJ  = I*J;
  IJK = IJ*K;

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

  DX = DualGrid(X,I); OX= DX[0]; LX= DX[I]-OX;
  DY = DualGrid(Y,J); OY= DY[0]; LY= DY[J]-OY;
  DZ = DualGrid(Z,K); OZ= DZ[0]; LZ= DZ[K]-OZ;

  if( I>1 ){ dx= X(1)-X(0); dx = 1.0/dx; } else { dx=1; }
  if( J>1 ){ dy= Y(1)-Y(0); dy = 1.0/dy; } else { dy=1; }
  if( K>1 ){ dz= Z(1)-Z(0); dz = 1.0/dz; } else { dz=1; }

  P = myGetPr( prhs[4] );
  nP   = myNumel( prhs[4] )/3;

  /*Parsing arguments*/
  /*Defaults*/
  interp_mode   = LINEAR;
  boundary_mode = VALUE;
  boundary_size = (LX+LY+LZ)/3;
  isEqualSpaced = -1;
  outside_value = 0.0; outside_value = outside_value/outside_value; /*NAN*/
  oM = NULL;
  nM = NULL;
  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );
    if( ! myStrcmpi(STR,"linear")   || ! myStrcmpi(STR,"lin")  ){ interp_mode = LINEAR;  isEqualSpaced =  0; argN++; continue; }
    if( ! myStrcmpi(STR,"*linear")  || ! myStrcmpi(STR,"*lin") ){ interp_mode = LINEAR;  isEqualSpaced =  1; argN++; continue; }
    if( ! myStrcmpi(STR,"?linear")  || ! myStrcmpi(STR,"?lin") ){ interp_mode = LINEAR;  isEqualSpaced = -1; argN++; continue; }
    if( ! myStrcmpi(STR,"nearest")  || ! myStrcmpi(STR,"nea")  ){ interp_mode = NEAREST; isEqualSpaced =  0; argN++; continue; }
    if( ! myStrcmpi(STR,"*nearest") || ! myStrcmpi(STR,"*nea") ){ interp_mode = NEAREST; isEqualSpaced =  1; argN++; continue; }
    if( ! myStrcmpi(STR,"?nearest") || ! myStrcmpi(STR,"?nea") ){ interp_mode = NEAREST; isEqualSpaced = -1; argN++; continue; }
    if( ! myStrcmpi(STR,"cubic")    || ! myStrcmpi(STR,"cub")  ){ interp_mode = CUBIC;   isEqualSpaced =  0; argN++; continue; }
    if( ! myStrcmpi(STR,"sinc")                              ){ interp_mode = SINC;    isEqualSpaced =  1; argN++; continue; }
    if( ! myStrcmpi(STR,"symmetric") || ! myStrcmpi(STR,"sym") ){ boundary_mode = SYMMETRIC; argN++; continue; }
    if( ! myStrcmpi(STR,"circular")  || 
        ! myStrcmpi(STR,"periodic") 
         || ! myStrcmpi(STR,"circ") || ! myStrcmpi(STR,"per") ) { boundary_mode = CIRCULAR; argN++; continue; }
    if( ! myStrcmpi(STR,"closest")                            ){ boundary_mode = CLOSEST; argN++; continue; }
    if( ! myStrcmpi(STR,"value") || ! myStrcmpi(STR,"extrapolation_value") )
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

    if( ! myStrcmpi( STR,"outside_value") || ! myStrcmpi( STR,"outval") || ! myStrcmpi( STR,"outvalue") ){
      argN++;
      if( nrhs > argN && myIsParameter( prhs[argN] ) ){
        if( !myIsEmpty( prhs[argN] ) ){ outside_value = myGetValue( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word OUTSIDE_VALUE a value has to be specified.\nIf empty value, default(NaN) is set.");
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

    mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
    
  }
  /*END Parsing arguments*/
  
  if( isEqualSpaced < 0 ){
    if( checkEqualSpaced( X , I ) && checkEqualSpaced( Y , J ) && checkEqualSpaced( Z , K ) ){
      isEqualSpaced = 1;
    } else {
      isEqualSpaced = 0;
    }
  }

  
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
  Odims[0]= nP;
  if( T > 1 ){
    Odims[1]= T;
    plhs[0] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
  } else {
    plhs[0] = mxCreateNumericArray( 1 , Odims , mxREAL_CLASS , mxREAL );
  }
  O = (real*)mxGetData( plhs[0] );
  /*END Creating output*/
  
  ii = -1; jj = -1; kk = -1;
  for( p=0 ; p<nP ; p++ ){
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
        if( zz < OZ || zz > OZ+LZ ||
            xx < OX || xx > OX+LX ||
            yy < OY || yy > OY+LY ){
          isOut = 1; 
        }
        break;
        
      case CLOSEST:
        if( xx < X[ 0 ] ){ xx = X[ 0 ]; } else  if( xx > X[I-1] ){ xx = X[I-1]; }
        if( yy < Y[ 0 ] ){ yy = Y[ 0 ]; } else  if( yy > Y[J-1] ){ yy = Y[J-1]; }
        if( zz < Z[ 0 ] ){ zz = Z[ 0 ]; } else  if( zz > Z[K-1] ){ zz = Z[K-1]; }
        break;

      case SYMMETRIC:
        xx = PutInside(xx,OX,2*LX);  if( xx > OX+LX ){   xx = 2*(OX+LX) - xx; }
        yy = PutInside(yy,OY,2*LY);  if( yy > OY+LY ){   yy = 2*(OY+LY) - yy; }
        zz = PutInside(zz,OZ,2*LZ);  if( zz > OZ+LZ ){   zz = 2*(OZ+LZ) - zz; }
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

/*
      case DECAY_TO_ID:
        if( zz < Z(0)-boundary_size || zz > Z(K-1)+boundary_size ||
            xx < X(0)-boundary_size || xx > X(I-1)+boundary_size ||
            yy < Y(0)-boundary_size || yy > Y(J-1)+boundary_size ){
          isOut = 1; 
        }
        break;
*/
        
    }
    if( isOut ){
      for(t=0;t<T;t++){
        O(p,t) = outside_value;
      } 
      continue;
    }
    
    
    switch( interp_mode ){
      case NEAREST:
        if( isEqualSpaced ){
          ii = (int) ( (xx-OX)*dx ); if( ii < 0 || xx<OX ){ ii= -2; } else if( ii > I-1 ){ if( xx == DX[I] ){ ii = I-1; } else { ii = -101; } }
          jj = (int) ( (yy-OY)*dy ); if( jj < 0 || yy<OY ){ jj= -2; } else if( jj > J-1 ){ if( yy == DY[J] ){ jj = J-1; } else { jj = -101; } }
          kk = (int) ( (zz-OZ)*dz ); if( kk < 0 || zz<OZ ){ kk= -2; } else if( kk > K-1 ){ if( zz == DZ[K] ){ kk = K-1; } else { kk = -101; } }
        } else {
          ii = GetInterval( xx , DX , I+1 , ii );
          jj = GetInterval( yy , DY , J+1 , jj );
          kk = GetInterval( zz , DZ , K+1 , kk );
        }
        
        if( boundary_mode == DECAY_TO_ZERO && ( ii < 0 || jj < 0 || kk < 0 ) ){
          for( t=0 ; t<T ; t++ ){
            O(p,t)       = 0;
          }
        } else {
          for( t=0 ; t<T ; t++ ){
            O(p,t)       = IM(ii,jj,kk,t);
          }
        }
        break;

      case LINEAR:
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

        if( I == 1 ){
          u   = 1;  uu  = 0;  ii0 = 0;  ii1 = 0;
        } else if( ii >= 0 ){
          u   = ( xx-X(ii) ) / ( X(ii+1)-X(ii) );
          uu  = 1-u;  ii0 =  ii;  ii1 = ii+1;  
        } else if( ii == -2 ) {
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              u   = ( xx - X(0) + boundary_size ) / boundary_size;
              uu  = 0; ii0 = 0;ii1 = 0;
              break;
            case SYMMETRIC: u = 0; uu = 1; ii0 = 0; ii1 = 0; break;
            case CIRCULAR:
              uu = ( X(0)-xx )/( X(1)-X(0)+X(I-1)-X(I-2) )*2;
              u  = 1-uu;  ii0 = I-1;  ii1 = 0;
              break;
          }
        } else if( ii == -101 ){
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              uu  = ( X(I-1) + boundary_size - xx ) / boundary_size;
              u   = 0; ii0 = I-1; ii1 = I-1;
              break;
            case SYMMETRIC: u = 0; uu = 1; ii0 = I-1; ii1 = I-1; break;
            case CIRCULAR:
              u = ( xx - X(I-1) )/( X(1)-X(0)+X(I-1)-X(I-2) )*2;
              uu = 1-u;  ii0 = I-1;  ii1 = 0;
              break;
          }
        }
          
        if( J == 1 ){
          v   = 1;  vv  = 0;  jj0 = 0;  jj1 = 0;
        } else if( jj >= 0 ){
          v   = ( yy-Y(jj) ) / ( Y(jj+1)-Y(jj) );
          vv  = 1-v;  jj0 =  jj;  jj1 = jj+1;  
        } else if( jj == -2 ) {
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              v   = ( yy - Y(0) + boundary_size ) / boundary_size;
              vv  = 0; jj0 = 0;jj1 = 0;
              break;
            case SYMMETRIC: v = 0; vv = 1; jj0 = 0; jj1 = 0; break;
            case CIRCULAR:
              vv = ( Y(0)-yy )/( Y(1)-Y(0)+Y(J-1)-Y(J-2) )*2;
              v  = 1-vv;  jj0 = J-1;  jj1 = 0;
              break;
          }
        } else if( jj == -101 ){
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              vv  = ( Y(J-1) + boundary_size - yy ) / boundary_size;
              v   = 0; jj0 = J-1; jj1 = J-1;
              break;
            case SYMMETRIC: v = 0; vv = 1; jj0 = J-1; jj1 = J-1; break;
            case CIRCULAR:
              v = ( yy - Y(J-1) )/( Y(1)-Y(0)+Y(J-1)-Y(J-2) )*2;
              vv = 1-v;  jj0 = J-1;  jj1 = 0;
              break;
          }
        }

        if( K == 1 ){
          w   = 1;  ww  = 0;  kk0 = 0;  kk1 = 0;
        } else if( kk >= 0 ){
          w   = ( zz-Z(kk) ) / ( Z(kk+1)-Z(kk) );
          ww  = 1-w;  kk0 =  kk;  kk1 = kk+1;  
        } else if( kk == -2 ) {
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              w   = ( zz - Z(0) + boundary_size ) / boundary_size;
              ww  = 0; kk0 = 0;kk1 = 0;
              break;
            case SYMMETRIC: w = 0; ww = 1; kk0 = 0; kk1 = 0; break;
            case CIRCULAR:
              ww = ( Z(0)-zz )/( Z(1)-Z(0)+Z(K-1)-Z(K-2) )*2;
              w  = 1-ww;  kk0 = K-1;  kk1 = 0;
              break;
          }
        } else if( kk == -101 ){
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              ww  = ( Z(K-1) + boundary_size - zz ) / boundary_size;
              w   = 0; kk0 = K-1; kk1 = K-1;
              break;
            case SYMMETRIC: w = 0; ww = 1; kk0 = K-1; kk1 = K-1; break;
            case CIRCULAR:
              w = ( zz - Z(K-1) )/( Z(1)-Z(0)+Z(K-1)-Z(K-2) )*2;
              ww = 1-w;  kk0 = K-1;  kk1 = 0;
              break;
          }
        }

        
        for( t=0 ; t<T ; t++ ){
          O(p,t)= 
              IM( ii0 , jj0 , kk0 ,t) * uu * vv * ww + 
              IM( ii1 , jj0 , kk0 ,t) *  u * vv * ww + 
              IM( ii0 , jj1 , kk0 ,t) * uu *  v * ww + 
              IM( ii1 , jj1 , kk0 ,t) *  u *  v * ww +
              IM( ii0 , jj0 , kk1 ,t) * uu * vv *  w + 
              IM( ii1 , jj0 , kk1 ,t) *  u * vv *  w +
              IM( ii0 , jj1 , kk1 ,t) * uu *  v *  w +
              IM( ii1 , jj1 , kk1 ,t) *  u *  v *  w ;
        }

        break;

      
      
      case CUBIC:
        ii = GetInterval( xx , X , I , ii );
        jj = GetInterval( yy , Y , J , jj );
        kk = GetInterval( zz , Z , K , kk );

        if( ii < 1 || ii > I-3 ){ continue; }
        if( jj < 1 || jj > J-3 ){ continue; }
        if( kk < 1 || kk > K-3 ){ continue; }
        
        for( t=0 ; t<T ; t++ ){
          O(p,t)= CINT( xx  
                       ,X(ii-1), CINT(yy,Y(jj-1),CINT(zz,Z(kk-1),IM(ii-1,jj-1,kk-1,t)
                                                        ,Z( kk ),IM(ii-1,jj-1, kk ,t)
                                                        ,Z(kk+1),IM(ii-1,jj-1,kk+1,t)
                                                        ,Z(kk+2),IM(ii-1,jj-1,kk+2,t)
                                                     )
                                        ,Y( jj ),CINT(zz,Z(kk-1),IM(ii-1, jj ,kk-1,t)
                                                        ,Z( kk ),IM(ii-1, jj , kk ,t)
                                                        ,Z(kk+1),IM(ii-1, jj ,kk+1,t)
                                                        ,Z(kk+2),IM(ii-1, jj ,kk+2,t)
                                                     )
                                        ,Y(jj+1),CINT(zz,Z(kk-1),IM(ii-1,jj+1,kk-1,t)
                                                        ,Z( kk ),IM(ii-1,jj+1, kk ,t)
                                                        ,Z(kk+1),IM(ii-1,jj+1,kk+1,t)
                                                        ,Z(kk+2),IM(ii-1,jj+1,kk+2,t)
                                                     )
                                        ,Y(jj+2),CINT(zz,Z(kk-1),IM(ii-1,jj+2,kk-1,t)
                                                        ,Z( kk ),IM(ii-1,jj+2, kk ,t)
                                                        ,Z(kk+1),IM(ii-1,jj+2,kk+1,t)
                                                        ,Z(kk+2),IM(ii-1,jj+2,kk+2,t)
                                                     ) 
                                     )
                       ,X( ii ), CINT(yy,Y(jj-1),CINT(zz,Z(kk-1),IM( ii ,jj-1,kk-1,t)
                                                        ,Z( kk ),IM( ii ,jj-1, kk ,t)
                                                        ,Z(kk+1),IM( ii ,jj-1,kk+1,t)
                                                        ,Z(kk+2),IM( ii ,jj-1,kk+2,t)
                                                     )
                                        ,Y( jj ),CINT(zz,Z(kk-1),IM( ii , jj ,kk-1,t)
                                                        ,Z( kk ),IM( ii , jj , kk ,t)
                                                        ,Z(kk+1),IM( ii , jj ,kk+1,t)
                                                        ,Z(kk+2),IM( ii , jj ,kk+2,t)
                                                     )
                                        ,Y(jj+1),CINT(zz,Z(kk-1),IM( ii ,jj+1,kk-1,t)
                                                        ,Z( kk ),IM( ii ,jj+1, kk ,t)
                                                        ,Z(kk+1),IM( ii ,jj+1,kk+1,t)
                                                        ,Z(kk+2),IM( ii ,jj+1,kk+2,t)
                                                     )
                                        ,Y(jj+2),CINT(zz,Z(kk-1),IM( ii ,jj+2,kk-1,t)
                                                        ,Z( kk ),IM( ii ,jj+2, kk ,t)
                                                        ,Z(kk+1),IM( ii ,jj+2,kk+1,t)
                                                        ,Z(kk+2),IM( ii ,jj+2,kk+2,t)
                                                     ) 
                                     )
                       ,X(ii+1), CINT(yy,Y(jj-1),CINT(zz,Z(kk-1),IM(ii+1,jj-1,kk-1,t)
                                                        ,Z( kk ),IM(ii+1,jj-1, kk ,t)
                                                        ,Z(kk+1),IM(ii+1,jj-1,kk+1,t)
                                                        ,Z(kk+2),IM(ii+1,jj-1,kk+2,t)
                                                     )
                                        ,Y( jj ),CINT(zz,Z(kk-1),IM(ii+1, jj ,kk-1,t)
                                                        ,Z( kk ),IM(ii+1, jj , kk ,t)
                                                        ,Z(kk+1),IM(ii+1, jj ,kk+1,t)
                                                        ,Z(kk+2),IM(ii+1, jj ,kk+2,t)
                                                     )
                                        ,Y(jj+1),CINT(zz,Z(kk-1),IM(ii+1,jj+1,kk-1,t)
                                                        ,Z( kk ),IM(ii+1,jj+1, kk ,t)
                                                        ,Z(kk+1),IM(ii+1,jj+1,kk+1,t)
                                                        ,Z(kk+2),IM(ii+1,jj+1,kk+2,t)
                                                     )
                                        ,Y(jj+2),CINT(zz,Z(kk-1),IM(ii+1,jj+2,kk-1,t)
                                                        ,Z( kk ),IM(ii+1,jj+2, kk ,t)
                                                        ,Z(kk+1),IM(ii+1,jj+2,kk+1,t)
                                                        ,Z(kk+2),IM(ii+1,jj+2,kk+2,t)
                                                     ) 
                                     )
                       ,X(ii+2), CINT(yy,Y(jj-1),CINT(zz,Z(kk-1),IM(ii+2,jj-1,kk-1,t)
                                                        ,Z( kk ),IM(ii+2,jj-1, kk ,t)
                                                        ,Z(kk+1),IM(ii+2,jj-1,kk+1,t)
                                                        ,Z(kk+2),IM(ii+2,jj-1,kk+2,t)
                                                     )
                                        ,Y( jj ),CINT(zz,Z(kk-1),IM(ii+2, jj ,kk-1,t)
                                                        ,Z( kk ),IM(ii+2, jj , kk ,t)
                                                        ,Z(kk+1),IM(ii+2, jj ,kk+1,t)
                                                        ,Z(kk+2),IM(ii+2, jj ,kk+2,t)
                                                     )
                                        ,Y(jj+1),CINT(zz,Z(kk-1),IM(ii+2,jj+1,kk-1,t)
                                                        ,Z( kk ),IM(ii+2,jj+1, kk ,t)
                                                        ,Z(kk+1),IM(ii+2,jj+1,kk+1,t)
                                                        ,Z(kk+2),IM(ii+2,jj+1,kk+2,t)
                                                     )
                                        ,Y(jj+2),CINT(zz,Z(kk-1),IM(ii+2,jj+2,kk-1,t)
                                                        ,Z( kk ),IM(ii+2,jj+2, kk ,t)
                                                        ,Z(kk+1),IM(ii+2,jj+2,kk+1,t)
                                                        ,Z(kk+2),IM(ii+2,jj+2,kk+2,t)
                                                     ) 
                                     ) );
        }
        
        
        
    }
  }
  
  EXIT: myFreeALLOCATES();
}


real CINT( real X , real x0 , real I0 , real x1 , real I1 , real x2 , real I2 , real x3 , real I3 ){
  real    D1, D2;
  real d12, dx2, dx1;
  
  d12 = x1-x2;
  dx2 =  X-x2;
  dx1 =  X-x1;
  
  D1 = (I0 - I1)/(x0 - x1) + (-I0 + I2)/(x0 - x2) + (I1 - I2)/(x1 - x2);
  D2 = (I1 - I2)/(x1 - x2) + (-I1 + I3)/(x1 - x3) + (I2 - I3)/(x2 - x3);
  
  return(
    (  
      dx1*I2*dx1*(2*X+x1-3*x2) + 
      dx1*( D2*dx1 + D1*dx2 )*dx2*d12 - 
      I1*dx2*dx2*(2*X-3*x1+x2)
    )
    /( d12*d12*d12 )
     
  );
  
  
}


/*


X = cumsum( 0.5 + 2*rand(10,1 ) ); X = X-X(1); X = X/X(end)*10;
DX = unique( [ X(1)-3  ;  X(1)-(X(2)-X(1))/2 ; (X(2:end)+X(1:end-1))/2 ; X(end)+(X(end)-X(end-1))/2  ; X(end)+3 ] );
I = cos( (X/10).^2*3*pi ); I(2) = 0.7;
I = repmat( I , [1 5 5] );

Xi = linspace( -6 , 16 , 5000 )'; Xi = [ Xi Xi*0 Xi*0 ];

subplot(621);
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Original');

subplot(623);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'nearest','value',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-Value');

subplot(624);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'linear','value',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Linear-Value');


subplot(625);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'nearest','decay',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-decay');

subplot(626);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'linear','decay',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Linear-decay');


subplot(627);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'nearest','circular',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-circular');

subplot(628);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'linear','circular',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Linear-circular');


subplot(629);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'nearest','symmetric',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-symmetric');

subplot(6,2,10);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'linear','symmetric',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Linear-symmetric');


 
 
 
subplot(6,2,11);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'nearest','closest',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-symmetric');

subplot(6,2,12);
plot( Xi(:,1) , ...
InterpPointsOn3DGrid(I,X,-2:2,-2:2,Xi,'outside_value',-0.2,'linear','closest',3) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Linear-symmetric');
 
 
*/
