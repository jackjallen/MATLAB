#include "mex.h"
#include "myMEX.h"

/*#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS*/
#undef real

#include "MESH2carvePolyhedron.h"
#include "carvePolyhedron2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  carve::poly::Polyhedron *MESH;

  MESH = MESH2carvePolyhedron( prhs[0] );
  


  plhs[0]= carvePolyhedron2MESH( MESH );

  EXIT:
    free(MESH);
    myFreeALLOCATES();
}

