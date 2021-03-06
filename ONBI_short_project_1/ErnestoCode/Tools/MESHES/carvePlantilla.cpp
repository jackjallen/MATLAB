#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define carveOBJ_TYPE      carve::__________

#include "_____________.h"
#include "MESH2carvePolyhedron.h"
#include "carvePolyhedron2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    /*Help Example	
    mexPrintf("vtkCleanPolyData( MESH ...\n");
    mexPrintf("\n");
    mexPrintf("SetTolerance             , real       ... Specify tolerance in terms of fraction of bounding box length\n");
    mexPrintf("ToleranceIsAbsoluteOn    , [] (*)     ...\n");
    mexPrintf("ToleranceIsAbsoluteOff   , []         ...\n");
    mexPrintf("\n");*/
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  carve::poly::Polyhedron *MESH;
  carveOBJ_TYPE	      *myVAR;

  MESH = MESH2carvePolyhedron( prhs[0] );
  
  /*Defaults*/
  myVAR->______________();
  myVAR->______________();
  myVAR->______________(___);
  myVAR->______________(___);
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
    /*   myErrMsgTxt( "No keywords." ); */
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( myVAR , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( myVAR , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( myVAR , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( myVAR , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
  

  plhs[0]= carvePolyhedron2MESH( myVAR->GetOutput() );

  EXIT:
    myVAR->Delete();
    MESH->Delete();
    myFreeALLOCATES();
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  /*[]-valued methods*/
  Call_0( ________________ );
  Call_0( ________________ );
  Call_0( ________________ );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  /*Real valued methods*/
  Call_1( ________________   , v );
  Call_1( ________________   , v );
  Call_1( ________________   , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
  /*String valued methods*/
  Call_1( ________________   , v );
  Call_1( ________________   , v );
  Call_1( ________________   , v );	
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
  /*Array valued methods*/
  Call_1( ________________   , v );
  Call_1( ________________   , v );
  Call_1( ________________   , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
