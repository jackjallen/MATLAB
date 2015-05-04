/*

 
 triplets = get_tripletes( ALL_edges , edges );

 
 */

#include "myMEX.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  int       nEDGES, f, nNODES;
  int       i,j,c1,c2,t,e, is_new;
  double    a,b,c;
  double    node1, node2, e1, e2, triplet;
  double    *TRIS, *EDGES, *E1, *E2, *NODES_1, *NODES_2;
  double    C1[1024], C2[1024];

  
  nEDGES = mxGetM( prhs[0] );
  E1     = mxGetPr( prhs[0] );
  
  E2 = E1 + nEDGES;
  
  
  plhs[0] = myCreateDoubleMatrix_E( 3 , 128 , mxREAL );
  TRIS = mxGetPr( plhs[0] );

  nNODES = mxGetM( prhs[1] );
  NODES_1 = mxGetPr( prhs[1] );
  NODES_2 = NODES_1 + nNODES;
  
  t = 0;
    
  for( f = 0 ; f < nNODES ; f++ ){

    node1 = NODES_1[f];
    node2 = NODES_2[f];
    
    if( node1 == node2 ){
      continue;
    }
    if( node1 > node2 ){
      a = node1;
      node1 = node2;
      node2 = a;
    }
    
    c1 = 0; c2 = 0;
    for( e = 0 ; e < nEDGES ; e++ ){

      e1 = E1[e];
      e2 = E2[e];

      if( e1 == node1 ){
        is_new = 1;
        for( i = 0 ; i < c1 ; i++ ){
          if( C1[i] == e2 ){ is_new = 0; break; }
        }
        if( is_new ){ C1[c1++] = e2; }
      } else if( e2 == node1 ){
        is_new = 1;
        for( i = 0 ; i < c1 ; i++ ){
          if( C1[i] == e1 ){ is_new = 0; break; }
        }
        if( is_new ){ C1[c1++] = e1; }
      }


      if( e2 == node2 ){
        C2[c2++] = e1;
      } else if( e1 == node2 ){
        C2[c2++] = e2;
      }

    }

    for( i = 0 ; i < c1 ; i++ ){
      triplet = 0;
      e1 = C1[i];

      for( j = 0 ; j < c2  ; j++ ){
        if( C2[j] == e1 ){
          triplet = e1;
          break;
        }
      }

      if( triplet ){
        if( triplet < node1 ){
          a = triplet;
          b = node1;
          c = node2;
        } else if( triplet < node2 ){
          a = node1;
          b = triplet;
          c = node2;
        } else {
          a = node1;
          b = node2;
          c = triplet;
        }
        
        
        if( myNumel( plhs[0] ) == t ){
          mxSetN( plhs[0] , mxGetN(plhs[0])*2 );
          mxSetPr( plhs[0] , mxRealloc( TRIS , myNumel( plhs[0] )*sizeof(double) ) );
          
/* /           memcpy( mxGetPr(plhs[0]) , TRIS , t*sizeof(double) ); */
          
          TRIS = mxGetPr( plhs[0] );
        }
        
        TRIS[ t     ] = a;
        TRIS[ t + 1 ] = b;
        TRIS[ t + 2 ] = c;
        t += 3;
      }

    }

  }
  
  mxSetN( plhs[0] , t/3 );
  
}

