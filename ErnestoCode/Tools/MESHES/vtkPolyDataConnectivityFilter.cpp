#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkPolyDataConnectivityFilter
#include "vtkPolyDataConnectivityFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkPolyDataConnectivityFilter( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetExtractionMode                     , int        ... Control the extraction of connected surfaces.\n");
    mexPrintf("\n");
    mexPrintf( "SetExtractionModeToLargestRegion      , [] (*)     ...\n");
    mexPrintf( "SetExtractionModeToAllRegions         , []         ...\n");
    mexPrintf("\n");
    mexPrintf( "SetExtractionModeToPointSeededRegions , []         ...\n");
    mexPrintf( "SetExtractionModeToCellSeededRegions  , []         ...\n");
    mexPrintf( "InitializeSeedList                    , []         ... Initialize list of point ids/cell ids used to seed regions.\n");
    mexPrintf( "AddSeed                               , id         ... Add a seed id (point or cell id). Note: ids are 0-offset.\n");
    mexPrintf( "DeleteSeed                            , id         ... Delete a seed id (point or cell id). Note: ids are 0-offset.\n");
    mexPrintf("\n");
    mexPrintf( "SetExtractionModeToSpecifiedRegions   , []         ...\n");
    mexPrintf( "InitializeSpecifiedRegionList         , []         ... Initialize list of region ids to extract.\n");
    mexPrintf( "AddSpecifiedRegion                    , id         ... Add a region id to extract. Note: ids are 0-offset.\n");
    mexPrintf( "DeleteSpecifiedRegion                 , id         ... Delete a region id to extract. Note: ids are 0-offset.\n");
    mexPrintf("\n");
    mexPrintf( "SetExtractionModeToClosestPointRegion , []         ...\n");
    mexPrintf( "SetClosestPoint                       , xyz        ... Use to specify x-y-z point coordinates when extracting the region closest to a specified point.\n");
    mexPrintf("\n");
    mexPrintf( "SetColorRegions                       , logical    ... Turn on/off the coloring of connected regions.\n");
    mexPrintf( "ColorRegionsOn                        , []         ...\n");
    mexPrintf( "ColorRegionsOff                       , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetScalarConnectivity                 , int        ... Turn on/off connectivity based on scalar value. If on, cells are connected only if they share points AND one of the cells scalar values falls in the scalar range specified.\n");
    mexPrintf( "ScalarConnectivityOn                  , []         ... \n");
    mexPrintf( "ScalarConnectivityOff                 , [] (*)     ... \n");
    mexPrintf("\n");
    mexPrintf( "SetScalarRange                        , [range1,range2]    ... Set the scalar range to use to extract cells based on scalar connectivity.\n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkOBJ_TYPE         *CONN;
  CONN= vtkOBJ_TYPE::New();
  CONN->SetInput( MESH );

  /*Defaults*/
  CONN->SetExtractionModeToLargestRegion();
  CONN->ColorRegionsOff();
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
      CallMethod( CONN , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( CONN , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( CONN , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( CONN , method , xyz );
    }
    argN++;
  }
  /*END Parsing arguments*/
    
  CONN->Update();
  plhs[0]= vtkPolyData2MESH( CONN->GetOutput() );

  EXIT:
    CONN->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( SetExtractionModeToClosestPointRegion );
  Call_0( SetExtractionModeToLargestRegion      );
  Call_0( SetExtractionModeToAllRegions         );
  Call_0( SetExtractionModeToPointSeededRegions );
  Call_0( SetExtractionModeToCellSeededRegions  );
  Call_0( SetExtractionModeToSpecifiedRegions   );
  Call_0( InitializeSpecifiedRegionList         );
  Call_0( InitializeSeedList                    );
  Call_0( ColorRegionsOn                        );
  Call_0( ColorRegionsOff                       );
  Call_0( ScalarConnectivityOn                  );
  Call_0( ScalarConnectivityOff                 );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetExtractionMode                     , v );
  Call_1( AddSeed                               , v );
  Call_1( DeleteSeed                            , v );
  Call_1( AddSpecifiedRegion                    , v );
  Call_1( DeleteSpecifiedRegion                 , v );
  Call_1( SetColorRegions                       , v );
  Call_1( SetColorRegions                       , v );
  Call_1( SetScalarConnectivity                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
  Call_1( SetScalarRange                        , v );
  Call_1( SetClosestPoint                       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
