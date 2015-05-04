/*
  o = Interp3DGridOn3DGrid( IM , X , Y , Z , Xn , Yn , Zn , T , mode , extrap_val )

*/

#include "mex.h"

#define IM(i,j,k,t)    IM[ (k)*IJ + (j)*I + (i) + (t)*IJK ]
#define O(i,j,k,t)      O[ (k)*IJn + (j)*In + (i) + (t)*IJKn ]
#define X(i)            X[ (i) ]
#define Y(j)            Y[ (j) ]
#define Z(k)            Z[ (k) ]
#define Xn(i)          Xn[ (i) ]
#define Yn(j)          Yn[ (j) ]
#define Zn(k)          Zn[ (k) ]
#define TM(i,j)        TM[ 4*(j) + (i) ]

int PrevIndex( double x , double *X , int I , int li ){
  int      i, ii, im, s;
  double   xm;
  
  s = 1;
  i = li;
  while( X(i) > x ){
    i = i - s;
    if( i < 0 ){
      i = 0;
      break;
    }
    s = s << 1;
  }
  
  s = 1;
  ii = i+1;
  while( X(ii) < x ){
    ii = ii + s;
    if( ii > I-1 ){
      ii = I-1;
      break;
    }
    s = s << 1;
  }

  while( ii-i > 1 ){
    im = (i+ii) >> 1;
    xm = X(im);
    
    if( x < xm ){
      ii = im;
    } else if( x > xm ){
      i = im;
    } else {
      return(im);
    }
  }
  return(i);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double  *IM, *O, *X, *Y, *Z, *Xn, *Yn, *Zn, *TM;
  double  xx,yy,zz,dx,dy,dz,u,v,w,xxnr,yynr,zznr;
  int     i,I,j,J,k,K,ii,jj,kk,uu,vv,ww,In,Jn,Kn,in,jn,kn,T,t,lii,ljj,lkk;
  long    IJ,IJn,IJK,IJKn;
  int     Odims[4], ndims;
  /*  1 :  nearest ,  2 :  linear */
  /* -1 : *nearest , -2 : *linear */
  int     mode;    
  double  extval;
  char    mode_str[50];

  IM= (double *)mxGetData(prhs[0]);
  I= (int) *( mxGetDimensions(prhs[0]) + 0 );
  J= (int) *( mxGetDimensions(prhs[0]) + 1 );
  IJ= I*J;
  ndims = mxGetNumberOfDimensions(prhs[0]);
  if( ndims > 2){
    K= (int) *( mxGetDimensions(prhs[0]) + 2 );
  } else {
    K=1;
  }
  IJK= IJ*K;
  T = 1;
  for( t=3 ; t<ndims ; t++ ){
    T = T * ( (int) *( mxGetDimensions(prhs[0]) + t ) );
  }
  
  X = (double *)mxGetData(prhs[1]);
  Y = (double *)mxGetData(prhs[2]);
  Z = (double *)mxGetData(prhs[3]);

  Xn = (double *)mxGetData(prhs[4]);
  Yn = (double *)mxGetData(prhs[5]);
  Zn = (double *)mxGetData(prhs[6]);
  In = mxGetNumberOfElements( prhs[4] );
  Jn = mxGetNumberOfElements( prhs[5] );
  Kn = mxGetNumberOfElements( prhs[6] );
  IJn  = In*Jn;
  IJKn = IJn*Kn;

  Odims[0]= In;
  Odims[1]= Jn;
  Odims[2]= Kn;
  Odims[3]= T;
  
  if( T > 1 ){
    plhs[0] = mxCreateNumericArray( 4 , Odims , mxDOUBLE_CLASS , mxREAL );  
  } else {
    plhs[0] = mxCreateNumericArray( 3 , Odims , mxDOUBLE_CLASS , mxREAL );  
  }
  O = (double *)mxGetData( plhs[0] );

  if( nrhs > 7 && mxGetNumberOfElements(prhs[7])>0 ) {
    TM = (double *)mxGetData( prhs[7] );
  } else {
    TM = NULL;
  }
  
  if( nrhs > 8 && mxGetNumberOfElements( prhs[8] )>0 ) {
    if( mxGetClassID(prhs[8]) != mxCHAR_CLASS ){
      mexErrMsgTxt("Incorrect mode type\n");
    }
    mxGetString( prhs[8], mode_str, 40 );  
    if( strcmp(mode_str,"linear") == 0 ){
      mode = 2;
    } else if( strcmp(mode_str,"*linear") == 0 ){
      mode = -2;
    } else if( strcmp(mode_str,"nearest") == 0 ){
      mode = 1;
    } else if( strcmp(mode_str,"*nearest") == 0 ){
      mode = -1;
    } else if( strcmp(mode_str,"") == 0 ){
      mode = 2;
    }  else { 
      mexErrMsgTxt("Incorrect mode\n");
    }  
  } else {
    mode = 2;
  }    

  if( nrhs > 9 ) {
    extval = *( mxGetPr(prhs[9]) );
  } else {
    extval = 0.0;  extval = 0.0/extval;
  }

  if( I>1 ){ dx= X(1)-X(0); } else { dx=1; }
  if( J>1 ){ dy= Y(1)-Y(0); } else { dy=1; }
  if( K>1 ){ dz= Z(1)-Z(0); } else { dz=1; }
  lii = 0; ljj = 0; lkk = 0;
  
  for( in=0 ; in<In ; in++ ){
  for( jn=0 ; jn<Jn ; jn++ ){
  for( kn=0 ; kn<Kn ; kn++ ){
    
    if( TM == NULL ) {
      zz= Zn(kn); if( zz < Z(0) || zz > Z(K-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }
      xx= Xn(in); if( xx < X(0) || xx > X(I-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }
      yy= Yn(jn); if( yy < Y(0) || yy > Y(J-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }
    } else {
      xxnr= Xn(in); 
      yynr= Yn(jn); 
      zznr= Zn(kn); 

      zz = xxnr*TM(2,0) + yynr*TM(2,1) + zznr*TM(2,2) + TM(2,3);
      if( zz < Z(0) || zz > Z(K-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }

      xx = xxnr*TM(0,0) + yynr*TM(0,1) + zznr*TM(0,2) + TM(0,3);
      if( xx < X(0) || xx > X(I-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }

      yy = xxnr*TM(1,0) + yynr*TM(1,1) + zznr*TM(1,2) + TM(1,3);
      if( yy < Y(0) || yy > Y(J-1) ){ for(t=0;t<T;t++){O(in,jn,kn,t)=extval;} continue; }
    }

    if( mode < 0 ){
      ii = (int) ( (xx-X(0))/dx ); if( ii == I-1 ){ ii= I-2; }
      jj = (int) ( (yy-Y(0))/dy ); if( jj == J-1 ){ jj= J-2; }
      kk = (int) ( (zz-Z(0))/dz ); if( kk == K-1 ){ kk= K-2; }
    } else {
      ii = PrevIndex( xx , X , I , lii ); lii= ii;
      jj = PrevIndex( yy , Y , J , ljj ); ljj= jj;
      kk = PrevIndex( zz , Z , K , lkk ); lkk= kk;
    }
    
    if( mode == 1 || mode == -1 ){

      if( I>1 && ii<I-1 && xx-X(ii) > X(ii+1)-xx ){ ii++; }
      if( J>1 && jj<J-1 && yy-Y(jj) > Y(jj+1)-yy ){ jj++; }
      if( K>1 && kk<K-1 && zz-Z(kk) > Z(kk+1)-zz ){ kk++; }
      for( t=0 ; t<T ; t++ ){
        O(in,jn,kn,t)= IM(ii,jj,kk,t);
      }

    } else if( mode == 2 || mode == -2 ){

      if( I>1 && J>1 && K>1 ){
        if(mode < 0 ){
          u = ( xx-X(ii) ) / dx ;
          v = ( yy-Y(jj) ) / dy ;
          w = ( zz-Z(kk) ) / dz ;
        } else {
          u = ( xx-X(ii) ) / ( X(ii+1)-X(ii) );
          v = ( yy-Y(jj) ) / ( Y(jj+1)-Y(jj) );
          w = ( zz-Z(kk) ) / ( Z(kk+1)-Z(kk) );
        }

        for( t=0 ; t<T ; t++ ){
          O(in,jn,kn,t)= 
                IM(ii  ,jj  ,kk  ,t) * (1-u)*(1-v)*(1-w) + 
                IM(ii+1,jj  ,kk  ,t) *   u  *(1-v)*(1-w) + 
                IM(ii  ,jj+1,kk  ,t) * (1-u)*  v  *(1-w) + 
                IM(ii+1,jj+1,kk  ,t) *   u  *  v  *(1-w) +
                IM(ii  ,jj  ,kk+1,t) * (1-u)*(1-v)*  w   + 
                IM(ii+1,jj  ,kk+1,t) *   u  *(1-v)*  w   +
                IM(ii  ,jj+1,kk+1,t) * (1-u)*  v  *  w   +
                IM(ii+1,jj+1,kk+1,t) *   u  *  v  *  w   ;
        }
        
      } else {
        if( I > 1 ){ 
          uu = 1; u = ( xx-X(ii) ) / ( X(ii+1)-X(ii) );  
        } else { 
          uu = 0; u = 0; 
        }
        if( J > 1 ){ 
          vv = 1; v = ( yy-Y(jj) ) / ( Y(jj+1)-Y(jj) ); 
        } else { 
          vv = 0; v = 0; 
        }
        if( K > 1 ){ 
          ww = 1; w = ( zz-Z(kk) ) / ( Z(kk+1)-Z(kk) ); 
        } else { 
          ww = 0; w = 0; 
        }

        for( t=0 ; t<T ; t++ ){
          O(in,jn,kn,t)= 
              IM( ii      , jj      ,kk      ,t) * (1-u)*(1-v)*(1-w) + 
              IM( ii+1*uu , jj      ,kk      ,t) *   u  *(1-v)*(1-w) + 
              IM( ii      , jj+1*vv ,kk      ,t) * (1-u)*  v  *(1-w) + 
              IM( ii+1*uu , jj+1*vv ,kk      ,t) *   u  *  v  *(1-w) +
              IM( ii      , jj      ,kk+1*ww ,t) * (1-u)*(1-v)*  w   + 
              IM( ii+1*uu , jj      ,kk+1*ww ,t) *   u  *(1-v)*  w   +
              IM( ii      , jj+1*vv ,kk+1*ww ,t) * (1-u)*  v  *  w   +
              IM( ii+1*uu , jj+1*vv ,kk+1*ww ,t) *   u  *  v  *  w   ;
        }
      }
    }
  }}}
}
