/*
I = randn( 100 , 200 ) < 0 ;

DX = (0:size(I,1)) + 1*( + 5 + rand(1,size(I,1)+1)*.5 );
DY = 0:size(I,2);


c = boundary( I , DX , DY,  - 0.2 );
figure
[DXX,DYY] = ndgrid(DX,DY);
surface( 'XData', DXX , 'YData', DYY , 'ZData', DXX*0 , 'CData', double(I) , 'facecolor', 'flat','edgecolor',[.6 .6 .6] )
line( 'XData', c(1,:) , 'YData', c(2,:) , 'linewidth' , 2 , 'Color',[0 0 1] , 'marker','.');
axis equal

% c = boundary( I , [] , [],  -0.2 );
% figure
% imagesc(I);
% line( 'XData', c(2,:) , 'YData', c(1,:) , 'linewidth' , 2 , 'Color',[0 0 1] , 'marker','.');
% axis equal
*/ 


#include "myMEX.h"

#define I(i,j)           I[ (j)*M + (i) ]
#define II(i,j)         II[ (j)*M + (i) ]

int           M,N,index,cierraagujero;
double        OFFSET = 0, *DX, *DY;
signed char   *I;
double        *xy;

void addxy( double x , double y){
  index++;
  if( ~myISNAN(x) ){
    if( index > 0){
      if( ~myISNAN(xy[2*index-2]) && ~myISNAN(xy[2*index]) ){
        if( xy[2*(index-2)+1] == xy[2*(index-1)+1] && xy[2*(index-1)+1] == x ){
          index--;
        } else if( xy[2*(index-2)] == xy[2*(index-1)] && xy[2*(index-1)] == y ){
          index--;
        }
      }
    }
  } else {
    if( cierraagujero > 0 ){
      xy[2*(index-1)+1]   = xy[2*cierraagujero+1];
      xy[2*cierraagujero] = xy[2*(index-1)];
      cierraagujero= -1;
    }
  }
  xy[2*index+1]=x; xy[2*index]=y;
}


int bits( signed char n , int b );
void uv( int i , int j , signed char n );
void tipo15(int i,int j);
void tipo15n(int i,int j);
void tipo01(int i,int j);
void tipo02(int i,int j);
void tipo04(int i,int j);
void tipo08(int i,int j);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  int           i,j,ii;
  signed char   *II;
  double        deltaX, deltaY;
  int           allocX, allocY, allocI = 0;
  int           MN;
  
  
  if( mxGetClassID(prhs[0]) != mxLOGICAL_CLASS ){
  	myErrMsgTxt("First argument have to be a logical array.");
  }
  
  M      = mySize( prhs[0] , 0 );
  N      = mySize( prhs[0] , 1 );
  MN     = M*N;
  
  if( nrhs > 1 ){
    if( myIsEmpty(prhs[1]) ){
      allocX = 1;
      deltaX = 1;
    } else if( myNumel(prhs[1]) == 1 ) {
      allocX = 1;
      deltaX = myGetValue( prhs[1] );
    } else {
      allocX = 0;
      if( myNumel(prhs[1]) != M+1 ){
        myErrMsgTxt("numel( DX ) have to be equal to size(I,1)+1\n");
      }
      DX = myGetPr(prhs[1]);
    }
  }
  if( allocX ){
//     mexPrintf("allocando X ( %d), %g\n",M,deltaX);
    DX = mxMalloc( (M+2)*sizeof( double ) );
    for( i = 0 ; i <= M ; i++ ){
      DX[i] = 0.5 + i*deltaX;
    }
//     mexPrintf("OK X\n"); myFlush();
  }

  
  if( nrhs > 2 ){
    if( myIsEmpty(prhs[2]) ){
      allocY = 1;
      deltaY = 1;
    } else if( myNumel(prhs[2]) == 1 ) {
      allocY = 1;
      deltaY = myGetValue( prhs[2] );
    } else {
      allocY = 0;
      if( myNumel(prhs[2]) != N+1 ){
        myErrMsgTxt("numel( DY ) have to be equal to size(I,2)+1\n");
      }
      DY = myGetPr(prhs[2]);
    }
  }
  if( allocY ){
//     mexPrintf("allocando Y ( %d), %g \n", N , deltaY);
    DY = mxMalloc( (N+2)*sizeof( double ) );
    for( i = 0 ; i <= N ; i++ ){
      DY[i] = 0.5 + i*deltaY;
    }
//     mexPrintf("OK Y\n"); myFlush();
  }
  

  if( nrhs > 3 && myNumel(prhs[3]) > 0 ){
    OFFSET = - myGetValue( prhs[3] );
  }
  
  II     = mxGetData( prhs[0] );

  for( j=0 ; j<N ; j++ ){
    for( i=0 ; i<M ; i++ ){
      if( !allocI && II(i,j) ){ 
        allocI = 1;
        I = mxMalloc( MN*sizeof( signed char ) );
        for( ii = 0 ; ii < MN ; ii++ ){ I[ii] = 0; }
      }
      if( !allocI ){ continue; }
      if( j > 0 ){
        if( II(i,j) != II(i,j-1) ){ I(i,j) += 1; }
      } else {       if( II(i,j) ){ I(i,j) += 1; } }
      if( j < N-1 ){
        if( II(i,j) != II(i,j+1) ){ I(i,j) += 4; }
      } else {       if( II(i,j) ){ I(i,j) += 4; } }
      if( i > 0 ){
        if( II(i,j) != II(i-1,j) ){ I(i,j) += 2; }
      } else {       if( II(i,j) ){ I(i,j) += 2; } }
      if( i < M-1 ){
        if( II(i,j) != II(i+1,j) ){ I(i,j) += 8; }
      } else {       if( II(i,j) ){ I(i,j) += 8; } }
      if( !II(i,j) ){ I(i,j) *= -1; }
    }
  }
  if( !allocI ){
    plhs[0] = mxCreateDoubleMatrix( 2 , 1 , mxREAL );
    xy      = mxGetPr( plhs[0] );
    xy[0]   = NAN;
    xy[1]   = NAN;
    goto EXIT;
  }

  plhs[0]= mxCreateDoubleMatrix( 2 , MN*10 , mxREAL );
  index  = -1;
  xy     = mxGetPr( plhs[0] );
  
  for( j=0 ; j<N ; j++ ){
    for( i=0 ; i<M ; i++ ){
      switch( I(i,j) ){
        case  15:   tipo15(i,j);   break;
        case   3:   tipo01(i,j);   break;
        case   7:   tipo01(i,j);   break;
        case  11:   tipo01(i,j);   break;
        case   4: 
          if( j<N-1 ){
            if(I(i,j+1) == 15){  tipo15(i,j+1); break; }
          }
          cierraagujero = index+1;
          tipo04(i,j);   
          break;
        case  -4:
          if( j<N-1 ){
            if(I(i,j+1) == 15){  tipo15(i,j+1); break; }
          } 
        case   8:
          if( i<M-1 ){
            if(I(i+1,j) == 15){  tipo15(i+1,j); break; }
          } 
        case  -8:
          if( i<M-1 ){
            if(I(i+1,j) == 15){  tipo15(i+1,j); break; }
          } 
      }
    }
  }
  if(index>=0){
    mxSetN( plhs[0] , index );
  } else {
    addxy( NAN , NAN );
    mxSetN( plhs[0] , 1 );
  }    

  EXIT:
    if( allocI ){ mxFree(I ); }
//     if( allocI ){ mexPrintf( "I free\n" ); }
    if( allocX ){ mxFree(DX); }
//     if( allocX ){ mexPrintf( "X free\n" ); }
    if( allocY ){ mxFree(DY); }
//     if( allocY ){ mexPrintf( "Y free\n" ); }
    myFreeALLOCATES();
}
  
void tipo01(i,j){
  if( I(i,j) < 0 ){ tipo04(i,j+1); return; }
  addxy( DY[  j ]+OFFSET , DX[  i ]+OFFSET );
  addxy( DY[  j ]+OFFSET , DX[ i+1]-OFFSET );
  uv( i , j   , 1 );
  uv( i , j-1 , 4 );

  if( bits( I(i,j) , 8 )) { 
    tipo08(i,j); 
  } else if(i<M-1 && bits( I(i+1, j ), 1 ) ){ 
    tipo01(i+1, j ); 
  } else if( j>0 && i<M-1 && bits(I(i+1,j-1),2)) {
    addxy( DY[  j ]+OFFSET , DX[ i+1]+OFFSET );
    tipo02(i+1,j-1); 
  } else {
    addxy( NAN , NAN );
  }
}

void tipo02(i,j){
  addxy( DY[ j+1]-OFFSET , DX[  i ]+OFFSET );
  addxy( DY[  j ]+OFFSET , DX[  i ]+OFFSET );
  uv( i ,j,2);
  uv(i-1,j,8);

  if( bits( I( i , j ), 1 )){ 
    tipo01( i , j );
  } else if( j>0 && bits( I( i ,j-1), 2 )){ 
    tipo02( i ,j-1); 
  } else if( i>0 && j>0 && bits(I(i-1,j-1),4)) {
    addxy( DY[  j ]-OFFSET , DX[  i ]+OFFSET );
    tipo04(i-1,j-1); 
  } else {
    addxy( NAN , NAN );
  }
}

void tipo04(i,j){
  addxy( DY[ j+1]-OFFSET , DX[ i+1]-OFFSET );
  addxy( DY[ j+1]-OFFSET , DX[  i ]+OFFSET );
  uv(i, j ,4);
  uv(i,j+1,1);

  if( bits( I( i , j ), 2 ) ){ 
    tipo02( i , j );
  } else if( i>0 && bits( I(i-1, j ), 4 ) ){ 
    tipo04(i-1, j );
  } else if( i>0 && j<N-1 && bits(I(i-1,j+1),8)) {
    addxy( DY[ j+1]-OFFSET , DX[  i ]-OFFSET );
    tipo08(i-1,j+1);
  } else {
    addxy( NAN , NAN );
  }
}

void tipo08(i,j){
  addxy( DY[  j ]+OFFSET , DX[ i+1]-OFFSET );
  addxy( DY[ j+1]-OFFSET , DX[ i+1]-OFFSET );
  uv( i ,j,8);
  uv(i+1,j,2);

  if( bits( I( i , j ), 4 ) ){ 
    tipo04( i , j ); 
  } else if( j<N-1 && bits( I( i ,j+1), 8 ) ){ 
    tipo08( i ,j+1);
  } else if( i<M-1 && j<N-1 && bits(I(i+1,j+1),1)) {
    addxy( DY[ j+1]+OFFSET , DX[ i+1]-OFFSET );
    tipo01(i+1,j+1);
  } else {
    addxy( NAN , NAN );
  }
}


void tipo15(int i,int j){
  addxy( DY[  j ]+OFFSET , DX[  i ]+OFFSET );
  addxy( DY[  j ]+OFFSET , DX[ i+1]-OFFSET );
  addxy( DY[ j+1]-OFFSET , DX[ i+1]-OFFSET );
  addxy( DY[ j+1]-OFFSET , DX[  i ]+OFFSET );
  addxy( DY[  j ]+OFFSET , DX[  i ]+OFFSET );

  I(i,j) = 0;
  uv( i-1 , j   , 8 );
  uv( i   , j-1 , 4 );
  uv( i+1 , j   , 2 );
  uv( i   , j+1 , 1 );
  addxy( NAN , NAN );
}

void tipo15n(int i,int j){
  addxy( j-1-OFFSET+1 , i+OFFSET  +1 );
  addxy( j-1-OFFSET+1 , i-1-OFFSET+1 );
  addxy( j+OFFSET  +1 , i-1-OFFSET+1 );
  addxy( j+OFFSET  +1 , i+OFFSET  +1 );
  addxy( j-1-OFFSET+1 , i+OFFSET  +1 );

  I(i,j) = 0;
  uv( i-1 , j   , 8 );
  uv( i   , j-1 , 4 );
  uv( i+1 , j   , 2 );
  uv( i   , j+1 , 1 );
  addxy( NAN , NAN );
}


int bits( signed char n , int b ){
  if( n<0 ){ n=-n; }
  switch( b ){
    case 1:
      if( n == 1 || n == 3 || n == 5 || n == 7 || n ==  9 || n == 11 || n == 13 || n==15 ){ return(1); }
      break;
    case 2:
      if( n == 2 || n == 3 || n == 6 || n == 7 || n == 10 || n == 11 || n == 14 || n==15 ){ return(1); }
      break;
    case 4:
      if( n == 4 || n == 5 || n == 6 || n == 7 || n == 12 || n == 13 || n == 14 || n==15 ){ return(1); }
      break;
    case 8:
      if( n>=8) { return(1); }
      break;
  }
  return(0);
}

void uv( int i , int j , signed char n ){
  if( i<M && j<N && i>=0 && j>= 0 ){
    if( I(i,j)<0 ){  I(i,j) += n;
    } else {         I(i,j) -= n;     }      
  }
  return;
}
