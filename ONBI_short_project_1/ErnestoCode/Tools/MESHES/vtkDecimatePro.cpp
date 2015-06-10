#include "mex.h"
#include "myMEX.h"

#define   real              double
#define   mxREAL_CLASS      mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkDecimatePro
#include "vtkDecimatePro.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkDecimate( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetPreserveTopology       , logical    ... Turn on/off whether to preserve the topology of the original mesh. If on, mesh splitting and hole elimination will not occur. This may limit the maximum reduction that may be achieved.\n");
    mexPrintf( "PreserveTopologyOn        , [] (*)     ...\n");
    mexPrintf( "PreserveTopologyOff       , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetSplitting              , logical    ... Turn on/off the splitting of the mesh at corners, along edges, at non-manifold points, or anywhere else a split is required. Turning splitting off will better preserve the original topology of the mesh, but you may not obtain the requested reduction.\n");
    mexPrintf( "SplittingOn               , [] (*)     ...\n");
    mexPrintf( "SplittingOff              , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetPreSplitMesh           , logical    ... In some cases you may wish to split the mesh prior to algorithm execution. This separates the mesh into semi-planar patches, which are disconnected from each other. This can give superior results in some cases. If the ivar PreSplitMesh ivar is enabled, the mesh is split with the specified SplitAngle. Otherwise mesh splitting is deferred as long as possible.\n");
    mexPrintf( "PreSplitMeshOn            , [] (*)     ...\n");
    mexPrintf( "PreSplitMeshOff           , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetAccumulateError        , logical    ... The computed error can either be computed directly from the mesh or the error may be accumulated as the mesh is modified. If the error is accumulated, then it represents a global error bounds, and the ivar MaximumError becomes a global bounds on mesh error. Accumulating the error requires extra memory proportional to the number of vertices in the mesh. If AccumulateError is off, then the error is not accumulated.\n");
    mexPrintf( "AccumulateErrorOn         , []         ...\n");
    mexPrintf( "AccumulateErrorOff        , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetBoundaryVertexDeletion , logical    ... Turn on/off the deletion of vertices on the boundary of a mesh. This may limit the maximum reduction that may be achieved.\n");
    mexPrintf( "BoundaryVertexDeletionOn  , []         ...\n");
    mexPrintf( "BoundaryVertexDeletionOff , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetErrorIsAbsolute        , logical    ... The MaximumError is normally defined as a fraction of the dataset bounding diagonal. By setting ErrorIsAbsolute to 1, the error is instead defined as that specified by AbsoluteError. By default ErrorIsAbsolute=0.\n");
    mexPrintf("\n");
    mexPrintf( "SetTargetReduction        , real (0.1) ... Specify the desired reduction in the total number of polygons (e.g., if TargetReduction is set to 0.9, this filter will try to reduce the data set to 10 per cent of its original size). Because of various constraints, this level of reduction may not be realized. If you want to guarantee a particular reduction, you must turn off PreserveTopology, turn on SplitEdges and BoundaryVertexDeletion, and set the MaximumError to VTK_DOUBLE_MAX (these ivars are initialized this way when the object is instantiated).\n");
    mexPrintf( "SetFeatureAngle           , real (180) ... Specify the mesh feature angle. This angle is used to define what an edge is (i.e., if the surface normal between two adjacent triangles is >= FeatureAngle, an edge exists).\n");
    mexPrintf( "SetSplitAngle             , real       ... Specify the mesh split angle. This angle is used to control the splitting of the mesh. A split line exists when the surface normals between two edge connected triangles are >= SplitAngle.\n");
    mexPrintf( "SetMaximumError           , real       ... Set the largest decimation error that is allowed during the decimation process. This may limit the maximum reduction that may be achieved. The maximum error is specified as a fraction of the maximum length of the input data bounding box.\n");
    mexPrintf( "SetAbsoluteError          , real       ... Same as MaximumError, but to be used when ErrorIsAbsolute is 1\n");
    mexPrintf( "SetDegree                 , real       ... If the number of triangles connected to a vertex exceeds Degree, then the vertex will be split. (NOTE: the complexity of the triangulation algorithm is proportional to Degree^2. Setting degree small can improve the performance of the algorithm.)\n");
    mexPrintf( "SetInflectionPointRatio   , real       ... Specify the inflection point ratio. An inflection point occurs when the ratio of reduction error between two iterations is greater than or equal to the InflectionPointRatio. \n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  vtkOBJ_TYPE         *DECI;

  MESH= MESH2vtkPolyData( prhs[0] );

  DECI= vtkOBJ_TYPE::New();
  DECI->SetInput( MESH );

  /*Defaults*/
  DECI->PreserveTopologyOn();
  DECI->SplittingOff();
  DECI->PreSplitMeshOn();
  DECI->SetTargetReduction( 0.1 );
  DECI->SetFeatureAngle( 180.0 );
  //     DECI->AccumulateErrorOn();
  //     DECI->BoundaryVertexDeletionOn();
  //     DECI->SetSplitAngle(double);
  //     DECI->SetMaximumError(double);
  //     DECI->SetAbsoluteError(double);
  //     DECI->SetInflectionPointRatio(double);
  //     DECI->SetDegree(int);
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
/*      myErrMsgTxt( "No keywords." ); */
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( DECI , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( DECI , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( DECI , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( DECI , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  DECI->Update();
  plhs[0]= vtkPolyData2MESH( DECI->GetOutput() );

  EXIT:
    DECI->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( PreserveTopologyOn        );
  Call_0( PreserveTopologyOff       );
  Call_0( SplittingOn               );
  Call_0( SplittingOff              );
  Call_0( PreSplitMeshOn            );
  Call_0( PreSplitMeshOff           );
  Call_0( AccumulateErrorOn         );
  Call_0( AccumulateErrorOff        );
  Call_0( BoundaryVertexDeletionOn  );
  Call_0( BoundaryVertexDeletionOff );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetTargetReduction        , v );
  Call_1( SetPreserveTopology       , v );
  Call_1( SetFeatureAngle           , v );
  Call_1( SetSplitting              , v );
  Call_1( SetSplitAngle             , v );
  Call_1( SetPreSplitMesh           , v );
  Call_1( SetMaximumError           , v );
  Call_1( SetAccumulateError        , v );
  Call_1( SetErrorIsAbsolute        , v );
  Call_1( SetAbsoluteError          , v );
  Call_1( SetBoundaryVertexDeletion , v );
  Call_1( SetDegree                 , v );
  Call_1( SetInflectionPointRatio   , v );
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
