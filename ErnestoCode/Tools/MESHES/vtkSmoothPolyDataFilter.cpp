#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkSmoothPolyDataFilter 
#include "vtkSmoothPolyDataFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

   if(!nrhs){
    mexPrintf("vtkSmoothPolyDataFilter( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetConvergence                    , real       ... Specify a convergence criterion for the iteration process. Smaller numbers result in more smoothing iterations.\n");  
    mexPrintf( "SetNumberOfIterations             , real (10)  ... Specify the number of iterations for Laplacian smoothing.\n");  
    mexPrintf( "SetRelaxationFactor               , real       ... Specify the relaxation factor for Laplacian smoothing. As in all iterative methods, the stability of the process is sensitive to this parameter. In general, small relaxation factors and large numbers of iterations are more stable than larger relaxation factors and smaller numbers of iterations.\n");  
    mexPrintf( "SetFeatureAngle                   , real       ... Specify the feature angle for sharp edge identification.\n");  
    mexPrintf( "SetEdgeAngle                      , real       ... Specify the edge angle to control smoothing along edges (either interior or boundary).\n");  
    mexPrintf("\n");
    mexPrintf( "SetFeatureEdgeSmoothing           , logical    ... Turn on/off smoothing along sharp interior edges.\n");  
    mexPrintf( "FeatureEdgeSmoothingOn            , [] (*)     ...\n");  
    mexPrintf( "FeatureEdgeSmoothingOff           , []         ...\n");  
    mexPrintf("\n");
    mexPrintf( "SetBoundarySmoothing              , logical    ... Turn on/off the smoothing of vertices on the boundary of the mesh.\n");  
    mexPrintf( "BoundarySmoothingOn               , [] (*)     ...\n");  
    mexPrintf( "BoundarySmoothingOff              , []         ...\n");  
    mexPrintf("\n");
    mexPrintf( "SetGenerateErrorScalars           , logical    ... Turn on/off the generation of scalar distance values.\n");  
    mexPrintf( "GenerateErrorScalarsOn            , []         ...\n");  
    mexPrintf( "GenerateErrorScalarsOff           , [] (*)     ...\n");  
    mexPrintf("\n");
    mexPrintf( "SetGenerateErrorVectors           , logical    ... Turn on/off the generation of error vectors.\n");  
    mexPrintf( "GenerateErrorVectorsOn            , []         ...\n");  
    mexPrintf( "GenerateErrorVectorsOff           , [] (*)     ...\n");  
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }




  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkSmoothPolyDataFilter      *SMOO;
  SMOO= vtkSmoothPolyDataFilter::New();
  SMOO->SetInput( MESH );

  /*Defaults*/             
  SMOO->GenerateErrorVectorsOff();     
  SMOO->GenerateErrorScalarsOff();   
  SMOO->BoundarySmoothingOn();     
  SMOO->FeatureEdgeSmoothingOn();   
  SMOO->SetNumberOfIterations(10);
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
      mexPrintf( "No keywords." );
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( SMOO , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( SMOO , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( SMOO , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( SMOO , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  SMOO->Update();
  plhs[0]= vtkPolyData2MESH( SMOO->GetOutput() );

  EXIT:
    SMOO->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( FeatureEdgeSmoothingOn       );
  Call_0( FeatureEdgeSmoothingOff      );
  Call_0( BoundarySmoothingOn          );
  Call_0( BoundarySmoothingOff         );
  Call_0( GenerateErrorScalarsOn       );
  Call_0( GenerateErrorScalarsOff      );
  Call_0( GenerateErrorVectorsOn       );
  Call_0( GenerateErrorVectorsOff      );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetGenerateErrorScalars   , v );
  Call_1( SetGenerateErrorVectors   , v );
  Call_1( SetFeatureAngle           , v );
  Call_1( SetEdgeAngle              , v );
  Call_1( SetBoundarySmoothing      , v );
  Call_1( SetFeatureEdgeSmoothing   , v );
  Call_1( SetConvergence            , v );
  Call_1( SetNumberOfIterations     , v );
  Call_1( SetRelaxationFactor       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
//   Call_1( SetClosestPoint                       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
