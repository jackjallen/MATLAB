#include "mex.h"
#include "stdlib.h"
#include "math.h"

// #define RRAND(a,b)      ( rand()*1.0/RAND_MAX*((b)-(a)) + (a) )
#define mxFlush()           mexEvalString("drawnow expose;");

#define D(i,j)        D[ (j)*M + (i) ]
#define F(i,j)        F[ (j)*M + (i) ]
#define I(i,j)        I[ (j)*M + (i) ]
#define G(i,j)        G[ (j)*M + (i) ]
#define E(i,j)        E[ (j)*M + (i) ]

struct node {
    double  d;
    int     i;
    int     j;
    struct  node * next;
};
typedef struct node node;

void AddItem(node *Q, double d , int i, int j);
struct node getMin(node * Q );

double *I, *G, *E;
int    M,N;

double fun(int i, int j, int ni, int nj){
  double  diag;
  double  fD;

  fD = 0;
  if(ni==i || nj==j){ 
    diag= 1.4142135623730951;
  } else {
    diag= 1; 
  }
  
  return( E(ni,nj) + G(ni,nj)/diag + fD );
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int       i, j, ni, nj;
  double    *F, *D,  d, INF;
  node      Q,n;

  INF = -log(0);
  Q.d = INF; Q.next = NULL;
  
  M =  mxGetM( prhs[0] );
  N =  mxGetN( prhs[0] );
  I = mxGetPr( prhs[0] );
  E = mxGetPr( prhs[2] );
  G = mxGetPr( prhs[3] );
  
  plhs[0] = mxCreateDoubleMatrix( M , N , mxREAL );
  plhs[1] = mxCreateDoubleMatrix( M , N , mxREAL );
  
  F = mxGetPr( plhs[0] );
  D = mxGetPr( plhs[1] );
  
  for(j=0;j<N;j++){
    for(i=0;i<M;i++){
      D(i,j)= INF;
    }
  }
  
//   double    *p, *q;
//   mxArray   *func[3], *out[1];
// 
//   func[0] = prhs[2];
//   func[1] = mxCreateDoubleMatrix( 1 , 2 , mxREAL ); p = mxGetPr( func[1] );
//   func[2] = mxCreateDoubleMatrix( 1 , 2 , mxREAL ); q = mxGetPr( func[2] );
//   out[0] = mxCreateDoubleMatrix( 1 , 1 , mxREAL ); 
//   mexCallMATLAB( 1 , out , 3 , func , "feval" );
//   OUTPUT = mxGetPr( out[0] );
//         p[0] = i; p[1]= j;
//           q[0] = ni; q[1] = nj;
  
  i = (int) *( mxGetPr(prhs[1])+0 )-1;
  j = (int) *( mxGetPr(prhs[1])+1 )-1;
  
  D(i,j) = 0;
  F(i,j) = 0;
  AddItem( &Q , 0 , i , j );
  n = getMin( &Q );
  while( n.d < INF ){
    i  = n.i; j  = n.j;
    for( ni=i-1 ; ni<=i+1 ; ni++ ){
      for( nj=j-1 ; nj<=j+1 ; nj++ ){
        if( ni==i && nj==j ){ continue; }
        if( ni >= 0 && ni < M && nj >= 0 && nj < N ){
          d = D(i,j) + fun(i,j,ni,nj);
          if( d < D(ni,nj) ) { 
            F(ni,nj) = j*M + i + 1; 
            D(ni,nj) = d; 
            AddItem( &Q , d , ni , nj ); 
          }
        }
      }
    }
    n= getMin( &Q );
  }
}









void AddItem(node *Q, double d, int i , int j ) {
  node *n;
  while( Q->next != NULL ){
    n = Q->next;
    if( n->d  >= d ){   break;   }
    Q = n;
  }
  n = Q->next;
  Q->next = (struct node *) mxMalloc( sizeof(node) );
  Q = Q->next;
  Q->next = n;
  Q->d    = d;
  Q->i    = i;
  Q->j    = j;
}

struct node getMin(node * Q ){
  node *n;
  if( Q->next == NULL ){
    return(*Q);
  }
  n = Q->next;
  if( n != NULL ){
    Q->next = n->next;
  }
  return(*n);
}
