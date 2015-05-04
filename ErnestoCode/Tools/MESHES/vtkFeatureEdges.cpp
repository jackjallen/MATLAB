#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkFeatureEdges
#include "vtkFeatureEdges.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkFeatureEdges( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetFeatureAngle       , real (180) ... Specify the feature angle for extracting feature edges.\n");
    mexPrintf("\n");
    mexPrintf( "SetBoundaryEdges      , logical    ... Turn on/off the extraction of boundary edges.\n");
    mexPrintf( "BoundaryEdgesOn       , [] (*)     ...\n");
    mexPrintf( "BoundaryEdgesOff      , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetFeatureEdges       , logical    ... Turn on/off the extraction of feature edges.\n");
    mexPrintf( "FeatureEdgesOn        , []         ...\n");
    mexPrintf( "FeatureEdgesOff       , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetNonManifoldEdges   , logical    ... Specify  \n");
    mexPrintf( "NonManifoldEdgesOn    , []         ...\n");
    mexPrintf( "NonManifoldEdgesOff   , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetManifoldEdges      , logical    ... Turn on/off the extraction of non-manifold edges.\n");
    mexPrintf( "ManifoldEdgesOn       , []         ...\n");
    mexPrintf( "ManifoldEdgesOff      , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetColoring           , logical    ... Turn on/off the coloring of edges by type.\n");
    mexPrintf( "ColoringOn            , []         ...\n");
    mexPrintf( "ColoringOff           , [] (*)     ...\n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  vtkOBJ_TYPE         *FEED;

  MESH= MESH2vtkPolyData( prhs[0] );

  FEED= vtkOBJ_TYPE::New();
  FEED->SetInput( MESH );

  /*Defaults*/
  FEED->SetFeatureAngle( 180.0 );
  FEED->BoundaryEdgesOn();
  FEED->FeatureEdgesOff();
  FEED->NonManifoldEdgesOff();
  FEED->ManifoldEdgesOff();
  FEED->ColoringOff();
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
      mexErrMsgTxt( "No keywords." );
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( FEED , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( FEED , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( FEED , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( FEED , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  FEED->Update();
  plhs[0]= vtkPolyData2MESH( FEED->GetOutput() );

  EXIT:
    FEED->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( ManifoldEdgesOn       );
  Call_0( ManifoldEdgesOff      );
  Call_0( ColoringOn            );
  Call_0( ColoringOff           );
  Call_0( NonManifoldEdgesOn    );
  Call_0( NonManifoldEdgesOff   );
  Call_0( FeatureEdgesOn        );
  Call_0( FeatureEdgesOff       );
  Call_0( BoundaryEdgesOn       );
  Call_0( BoundaryEdgesOff      );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetBoundaryEdges      , v );
  Call_1( SetFeatureEdges       , v );
  Call_1( SetFeatureAngle       , v );
  Call_1( SetNonManifoldEdges   , v );
  Call_1( SetColoring           , v );
  Call_1( SetManifoldEdges      , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
