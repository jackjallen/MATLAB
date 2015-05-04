#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkHull 
#include "vtkHull.h"
#include "vtkTriangleFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz, p[3], n[3];
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkHull( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "AddCubeVertexPlanes            , []         ... Add the 8 planes that represent the vertices of a cube - the combination of the three face planes connecting to a vertex - (1,1,1), (1,1,-1), (1,-1,1), (1,-1,1), (-1,1,1), (-1,1,-1), (-1,-1,1), (-1,-1-1).\n");  
    mexPrintf( "AddCubeEdgePlanes              , []         ... Add the 12 planes that represent the edges of a cube - halfway between the two connecting face planes - (1,1,0), (-1,-1,0), (-1,1,0), (1,-1,0), (0,1,1), (0,-1,-1), (0,1,-1), (0,-1,1), (1,0,1), (-1,0,-1), (1,0,-1), (-1,0,1).\n");  
    mexPrintf( "AddCubeFacePlanes              , [] (*)     ... Add the six planes that make up the faces of a cube - (1,0,0), (-1, 0, 0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1).\n");  
    mexPrintf( "AddRecursiveSpherePlanes       , int        ...\n");  
    mexPrintf("\n");
    mexPrintf( "AddPlane (double A, double B, double C) \n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkOBJ_TYPE      *HULL;
  HULL= vtkOBJ_TYPE::New();
  HULL->SetInput( MESH );

  /*Defaults*/
  HULL->AddCubeFacePlanes( );
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
      CallMethod( HULL , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( HULL , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( HULL , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( HULL , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  HULL->Update();
  
  vtkTriangleFilter   *TRI;
  TRI = vtkTriangleFilter::New();
  TRI->SetInput( HULL->GetOutput() );
  TRI->PassVertsOff();
  TRI->PassLinesOff();
  TRI->Update();
  
  plhs[0]= vtkPolyData2MESH( TRI->GetOutput() );

  EXIT:
    HULL->Delete();
    MESH->Delete();
    TRI->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( AddCubeVertexPlanes         );
  Call_0( AddCubeEdgePlanes           );
  Call_0( AddCubeFacePlanes           );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( AddRecursiveSpherePlanes       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
  Call_1( AddPlane                       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
