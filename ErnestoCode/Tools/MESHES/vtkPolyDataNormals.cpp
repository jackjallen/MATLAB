#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkPolyDataNormals
#include "vtkPolyDataNormals.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkComputeNormals( MESH ...\n");
    mexPrintf("\n");
    mexPrintf("SetFeatureAngle          , real       ... Specify the angle that defines a sharp edge. If the difference in angle across neighboring polygons is greater than this value, the shared edge is considered 'sharp'.\n");
    mexPrintf("\n");
    mexPrintf("SetSplitting             , logical    ... Turn on/off the splitting of sharp edges.\n");
    mexPrintf("SplittingOn              , []         ...\n");
    mexPrintf("SplittingOff             , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf("SetConsistency           , logical    ... Turn on/off the enforcement of consistent polygon ordering.\n");
    mexPrintf("ConsistencyOn            , []         ...\n");
    mexPrintf("ConsistencyOff           , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf("SetAutoOrientNormals     , logical    ... Turn on/off the automatic determination of correct normal orientation. NOTE: This assumes a completely closed surface (i.e. no boundary edges) and no non-manifold edges. If these constraints do not hold, all bets are off. This option adds some computational complexity, and is useful if you don't want to have to inspect the rendered image to determine whether to turn on the FlipNormals flag. However, this flag can work with the FlipNormals flag, and if both are set, all the normals in the output will point 'inward'.\n");
    mexPrintf("AutoOrientNormalsOn      , [] (*)     ...\n"); 
    mexPrintf("AutoOrientNormalsOff     , []         ...\n"); 
    mexPrintf("\n");
    mexPrintf("SetComputePointNormals   , logical    ... Turn on/off the computation of point normals.\n");
    mexPrintf("ComputePointNormalsOn    , [] (*)     ...\n"); 
    mexPrintf("ComputePointNormalsOff   , []         ...\n"); 
    mexPrintf("\n");
    mexPrintf("SetComputeCellNormals    , logical    ... Turn on/off the computation of cell normals.\n");
    mexPrintf("ComputeCellNormalsOn     , []         ...\n"); 
    mexPrintf("ComputeCellNormalsOff    , [] (*)     ...\n"); 
    mexPrintf("\n");
    mexPrintf("SetFlipNormals           , logical    ... Turn on/off the global flipping of normal orientation. Flipping reverves the meaning of front and back for Frontface and Backface culling in vtkProperty. Flipping modifies both the normal direction and the order of a cell's points.\n");
    mexPrintf("FlipNormalsOn            , []         ...\n"); 
    mexPrintf("FlipNormalsOff           , [] (*)     ...\n"); 
    mexPrintf("\n");
    mexPrintf("SetNonManifoldTraversal  , logical    ... Turn on/off traversal across non-manifold edges. This will prevent problems where the consistency of polygonal ordering is corrupted due to topological loops.\n");
    mexPrintf("NonManifoldTraversalOn   , []         ...\n"); 
    mexPrintf("NonManifoldTraversalOff  , []         ...\n"); 
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();

  long                n_xyz, p;
  vtkPolyData         *MESH;
  vtkOBJ_TYPE         *N;

  MESH= MESH2vtkPolyData( prhs[0] );
  
  N= vtkOBJ_TYPE::New();
  N->SetInput( MESH );

  /*Defaults*/
  N->ComputeCellNormalsOff();
  N->ComputePointNormalsOn();
  N->ConsistencyOff();
  N->AutoOrientNormalsOn();
  N->SplittingOff();         
  N->FlipNormalsOff();
  /*END Defaults*/

  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
    /*  myErrMsgTxt( "No keywords." );  */
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( N , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( N , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( N , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( N , method , xyz );
    }
    argN++;
  }
  /*END Parsing arguments*/
  
  
  N->Update();  
  plhs[0]= vtkPolyData2MESH( N->GetOutput() );
  
  
  EXIT:
    N->Delete();
    MESH->Delete();
    myFreeALLOCATES();
}


void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( SplittingOn               );
  Call_0( SplittingOff              );
  Call_0( ConsistencyOn             );
  Call_0( ConsistencyOff            );
  Call_0( AutoOrientNormalsOn       );
  Call_0( AutoOrientNormalsOff      );
  Call_0( ComputePointNormalsOn     );
  Call_0( ComputePointNormalsOff    );
  Call_0( ComputeCellNormalsOn      );
  Call_0( ComputeCellNormalsOff     );
  Call_0( FlipNormalsOn             );
  Call_0( FlipNormalsOff            );
  Call_0( NonManifoldTraversalOn    );
  Call_0( NonManifoldTraversalOff   );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetFeatureAngle         , v );
  Call_1( SetSplitting            , v );
  Call_1( SetConsistency          , v );
  Call_1( SetAutoOrientNormals    , v );
  Call_1( SetComputePointNormals  , v );
  Call_1( SetComputeCellNormals   , v );
  Call_1( SetFlipNormals          , v );
  Call_1( SetNonManifoldTraversal , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
//   Call_1( SetFileName              , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName              , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
