#include "myMEX.h"
#include "stdlib.h"

#define I(i,j)           I[ (j)*M + (i) ]
#define II(i,j)         II[ (j)*M + (i) ]
#define mxFlush()       mexEvalString("drawnow expose;");
#define X(j)            OX + (j)*DX
#define Y(i)            OY + (i)*DY
#define O               offset
#define isnan(x)        ~(x==x)

int           M,N,index,cierraagujero;
double        offset, OX, DX, OY, DY;
signed char   *I;
double        *xy;
double        _NAN_;


void addxy(double x,double y){
  index++;
  if( ~isnan(x) ){
    if( index > 0){
      if( ~isnan(xy[2*index-2]) && ~isnan(xy[2*index]) ){
        if( xy[2*(index-2)] == xy[2*(index-1)] && xy[2*(index-1)] == x ){
          index--;
        } else if( xy[2*(index-2)+1] == xy[2*(index-1)+1] && xy[2*(index-1)+1] == y ){
          index--;
        }
      }
    }
  } else {
    if( cierraagujero > 0 ){
      xy[2*(index-1)]       = xy[2*cierraagujero];
      xy[2*cierraagujero+1] = xy[2*(index-1)+1];
      cierraagujero= -1;
    }
  }
  xy[2*index]=x; xy[2*index+1]=y;
}



int bits( signed char n , int b );
void uv( int i , int j , signed char n );
void tipo15(int i,int j);
void tipo15n(int i,int j);
void tipo01(int i,int j);
void tipo02(int i,int j);
void tipo04(int i,int j);
void tipo08(int i,int j);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int           i,j;
  signed char   *II;
  
//  _NAN_= log(-1);         
//  _NAN_= 0xFFFFFFFF;
//   _NAN_= 0xFFC00000;
 _NAN_= 0.0; _NAN_ = _NAN_/_NAN_;
//   mexPrintf("NAN:  %X  -  %f\n" , 0xFFC00000, _NAN_ );
  
  if(nrhs>1){ offset = *(mxGetPr(prhs[1]));
  } else {    offset = 0.1;  }
  if(nrhs>2){ OX = *(mxGetPr(prhs[2]));
  } else    { OX = 0.5; }
  if(nrhs>3){ DX = *(mxGetPr(prhs[3]));
  } else    { DX = 1; }
  if(nrhs>4){ OY = *(mxGetPr(prhs[4]));
  } else    { OY = 0.5; }
  if(nrhs>5){ DY = *(mxGetPr(prhs[5]));
  } else    { DY = 1; }

  II     = mxGetData( prhs[0] );
  M      = (int) *( mxGetDimensions(prhs[0]) + 0 );
  N      = (int) *( mxGetDimensions(prhs[0]) + 1 );
  
  I = mxMalloc( M*N*sizeof( signed char ) );
  for( j=0 ; j<N ; j++ ){
    for( i=0 ; i<M ; i++ ){
      I(i,j) = 0;
    }
  }
  for( j=0 ; j<N ; j++ ){
    for( i=0 ; i<M ; i++ ){
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

//   for( i=0 ; i<M ; i++ ){
//     for( j=0 ; j<N ; j++ ){
//       mexPrintf("%2d ", I(i,j));
//     }
//     mexPrintf("\n");
//   }

  plhs[0]= mxCreateDoubleMatrix( 2 , M*N*2 , mxREAL );
  index  = -1;
  xy     = mxGetPr( plhs[0] );
  
  for( j=0 ; j<N ; j++ ){
    for( i=0 ; i<M ; i++ ){
      switch( I(i,j) ){
        case  15:   tipo15(i,j);   break;
//         case -15:  tipo15n(i,j);   break;
        case   3:   tipo01(i,j);   break;
        case   7:   tipo01(i,j);   break;
        case  11:   tipo01(i,j);   break;
        case   4: 
          if( j<N-1 ){
            if(I(i,j+1) == 15){  tipo15(i,j+1); break; }
//             if(I(i,j+1) ==-15){ tipo15n(i,j+1); break; }
          }
          cierraagujero = index+1;
          tipo04(i,j);   
//           cierraagujero = -1;
          break;
        case  -4:
          if( j<N-1 ){
            if(I(i,j+1) == 15){  tipo15(i,j+1); break; }
//             if(I(i,j+1) ==-15){ tipo15n(i,j+1); break; }
          } 
        case   8:
          if( i<M-1 ){
            if(I(i+1,j) == 15){  tipo15(i+1,j); break; }
//             if(I(i+1,j) ==-15){ tipo15n(i+1,j); break; }
          } 
        case  -8:
          if( i<M-1 ){
            if(I(i+1,j) == 15){  tipo15(i+1,j); break; }
//             if(I(i+1,j) ==-15){ tipo15n(i+1,j); break; }
          } 
      }
    }
  }
  mxFree(I);
  if(index>=0){
    mxSetN( plhs[0] , index );
  } else {
    mxSetN( plhs[0] , 0 );
    mxSetM( plhs[0] , 0 );
  }    
}
  
void tipo01(i,j){
  if( I(i,j)<0 ){ tipo04(i,j+1); return; }
  addxy( X( j )+O , Y( i )+O );
  addxy( X( j )+O , Y(i+1)-O );
  uv( i , j   , 1 );
  uv( i , j-1 , 4 );

  if( bits( I(i,j) , 8 )) { 
    tipo08(i,j); 
  } else if(i<M-1 && bits( I(i+1, j ), 1 ) ){ 
    tipo01(i+1, j ); 
  } else if( j>0 && i<M-1 && bits(I(i+1,j-1),2)) {
    addxy( X( j )+O , Y(i+1)+O );
    tipo02(i+1,j-1); 
  } else {
    addxy( _NAN_ , _NAN_ );
  }
}

void tipo02(i,j){
  addxy( X(j+1)-O , Y( i )+O );
  addxy( X( j )+O , Y( i )+O );
  uv( i ,j,2);
  uv(i-1,j,8);

  if( bits( I( i , j ), 1 )){ 
    tipo01( i , j );
  } else if( j>0 && bits( I( i ,j-1), 2 )){ 
    tipo02( i ,j-1); 
  } else if( i>0 && j>0 && bits(I(i-1,j-1),4)) {
    addxy( X( j )-O , Y( i )+O );
    tipo04(i-1,j-1); 
  } else {
    addxy( _NAN_ , _NAN_ );
  }
}

void tipo04(i,j){
  addxy( X(j+1)-O , Y(i+1)-O );
  addxy( X(j+1)-O , Y( i )+O );
  uv(i, j ,4);
  uv(i,j+1,1);

  if( bits( I( i , j ), 2 ) ){ 
    tipo02( i , j );
  } else if( i>0 && bits( I(i-1, j ), 4 ) ){ 
    tipo04(i-1, j );
  } else if( i>0 && j<N-1 && bits(I(i-1,j+1),8)) {
    addxy( X(j+1)-O , Y( i )-O );
    tipo08(i-1,j+1);
  } else {
    addxy( _NAN_ , _NAN_ );
  }
}

void tipo08(i,j){
  addxy( X( j )+O , Y(i+1)-O );
  addxy( X(j+1)-O , Y(i+1)-O );
  uv( i ,j,8);
  uv(i+1,j,2);

  if( bits( I( i , j ), 4 ) ){ 
    tipo04( i , j ); 
  } else if( j<N-1 && bits( I( i ,j+1), 8 ) ){ 
    tipo08( i ,j+1);
  } else if( i<M-1 && j<N-1 && bits(I(i+1,j+1),1)) {
    addxy( X(j+1)+O , Y(i+1)-O );
    tipo01(i+1,j+1);
  } else {
    addxy( _NAN_ , _NAN_ );
  }
}


void tipo15(int i,int j){
  addxy( X( j )+O , Y( i )+O );
  addxy( X( j )+O , Y(i+1)-O );
  addxy( X(j+1)-O , Y(i+1)-O );
  addxy( X(j+1)-O , Y( i )+O );
  addxy( X( j )+O , Y( i )+O );

  I(i,j) = 0;
  uv( i-1 , j   , 8 );
  uv( i   , j-1 , 4 );
  uv( i+1 , j   , 2 );
  uv( i   , j+1 , 1 );
  addxy( _NAN_ , _NAN_ );
}

void tipo15n(int i,int j){
  addxy( j-1-offset+1 , i+offset  +1 );
  addxy( j-1-offset+1 , i-1-offset+1 );
  addxy( j+offset  +1 , i-1-offset+1 );
  addxy( j+offset  +1 , i+offset  +1 );
  addxy( j-1-offset+1 , i+offset  +1 );

  I(i,j) = 0;
  uv( i-1 , j   , 8 );
  uv( i   , j-1 , 4 );
  uv( i+1 , j   , 2 );
  uv( i   , j+1 , 1 );
  addxy( _NAN_ , _NAN_ );
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

