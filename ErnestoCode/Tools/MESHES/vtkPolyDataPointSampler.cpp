#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkPolyDataPointSampler 
#include "vtkPolyDataPointSampler.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

   if(!nrhs){
    mexPrintf("vtkPolyDataPointSampler( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetDistance                   , real( 0.01 )  ... Set/Get the approximate distance between points. This is an absolute distance measure.\n");
    mexPrintf("\n");
    mexPrintf( "SetGenerateVertexPoints       , int           ... Specify/retrieve a boolean flag indicating whether cell vertex points should be output.\n");
    mexPrintf( "GenerateVertexPointsOn        , []  (*)       ...\n");
    mexPrintf( "GenerateVertexPointsOff       , []            ...\n");
    mexPrintf("\n");
    mexPrintf( "SetGenerateEdgePoints         , int           ... Specify/retrieve a boolean flag indicating whether cell edges should be sampled to produce output points.\n");
    mexPrintf( "GenerateEdgePointsOn          , [] (*)        ...\n");
    mexPrintf( "GenerateEdgePointsOff         , []            ...\n");
    mexPrintf("\n"); 
    mexPrintf( "SetGenerateInteriorPoints     , int           ... Specify/retrieve a boolean flag indicating whether cell interiors should be sampled to produce output points.\n");
    mexPrintf( "GenerateInteriorPointsOn      , [] (*)        ...\n");
    mexPrintf( "GenerateInteriorPointsOff     , []            ...\n");
    mexPrintf("\n"); 
    mexPrintf( "SetGenerateVertices          , int            ... Specify/retrieve a boolean flag indicating whether cell vertices should be generated. Cell vertices are useful if you actually want to display the points (that is, for each point generated, a vertex is generated). Recall that VTK only renders vertices and not points.\n");
    mexPrintf( "GenerateVerticesOn           , []             ...\n");
    mexPrintf( "GenerateVerticesOff          , [] (*)         ...\n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }


  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkPolyDataPointSampler      *POIN;
  POIN= vtkPolyDataPointSampler::New();
  POIN->SetInput( MESH );

  /*Defaults*/             
  POIN->SetDistance(0.01);     
  POIN->GenerateVertexPointsOn();     
  POIN->GenerateEdgePointsOn();     
  POIN->GenerateInteriorPointsOn();     
  POIN->GenerateVerticesOff();     
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
      CallMethod( POIN , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( POIN , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( POIN , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( POIN , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  POIN->Update();
  plhs[0]= vtkPolyData2MESH( POIN->GetOutput() );

  EXIT:
    POIN->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( GenerateVertexPointsOn      );
  Call_0( GenerateVertexPointsOff     );
  Call_0( GenerateEdgePointsOn        );
  Call_0( GenerateEdgePointsOff       );
  Call_0( GenerateInteriorPointsOn    );
  Call_0( GenerateInteriorPointsOff   );
  Call_0( GenerateVerticesOn          );
  Call_0( GenerateVerticesOff         );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetDistance                 , v );
  Call_1( SetGenerateVertexPoints     , v );
  Call_1( SetGenerateEdgePoints       , v );
  Call_1( SetGenerateInteriorPoints   , v );
  Call_1( SetGenerateVertices         , v );
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
