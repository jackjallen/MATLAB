#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkQuadricDecimation 
#include "vtkQuadricDecimation.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz, p[3], n[3];
  char                STR[2000], method[2000];

  if(!nrhs){
    mexPrintf("vtkQuadricDecimation( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetTargetReduction             , real (0.1)  ... Set/Get the desired reduction (expressed as a fraction of the original number of triangles). The actual reduction may be less depending on triangulation and topological constraints.\n");
    mexPrintf("\n");
    mexPrintf( "SetAttributeErrorMetric        , logical    ... Decide whether to include data attributes in the error metric. If off, then only geometric error is used to control the decimation. By default the attribute errors are off.\n");  
    mexPrintf( "AttributeErrorMetricOn         , []         ...\n");
    mexPrintf( "AttributeErrorMetricOff        , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetScalarsAttribute            , logical    ... If attribute errors are to be included in the metric (i.e., AttributeErrorMetric is on), then the following flags control which attributes are to be included in the error calculation. By default all of these are on.\n");  
    mexPrintf( "ScalarsAttributeOn             , []         ...\n");
    mexPrintf( "ScalarsAttributeOff            , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetVectorsAttribute            , logical    ... If attribute errors are to be included in the metric (i.e., AttributeErrorMetric is on), then the following flags control which attributes are to be included in the error calculation. By default all of these are on.\n");  
    mexPrintf( "VectorsAttributeOn             , []         ...\n");
    mexPrintf( "VectorsAttributeOff            , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetNormalsAttribute            , logical    ... If attribute errors are to be included in the metric (i.e., AttributeErrorMetric is on), then the following flags control which attributes are to be included in the error calculation. By default all of these are on.\n");  
    mexPrintf( "NormalsAttributeOn             , []         ...\n");
    mexPrintf( "NormalsAttributeOff            , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetTCoordsAttribute            , logical    ... If attribute errors are to be included in the metric (i.e., AttributeErrorMetric is on), then the following flags control which attributes are to be included in the error calculation. By default all of these are on.\n");  
    mexPrintf( "TCoordsAttributeOn             , []         ...\n");
    mexPrintf( "TCoordsAttributeOff            , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetTensorsAttribute            , logical    ... If attribute errors are to be included in the metric (i.e., AttributeErrorMetric is on), then the following flags control which attributes are to be included in the error calculation. By default all of these are on.\n");  
    mexPrintf( "TensorsAttributeOn             , []         ...\n");
    mexPrintf( "TensorsAttributeOff            , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf( "SetScalarsWeight               , real       ... Set/Get the scaling weight contribution of the attribute. These values are used to weight the contribution of the attributes towards the error metric.\n");
    mexPrintf( "SetVectorsWeight               , real       ... Set/Get the scaling weight contribution of the attribute. These values are used to weight the contribution of the attributes towards the error metric.\n");
    mexPrintf( "SetNormalsWeight               , real       ... Set/Get the scaling weight contribution of the attribute. These values are used to weight the contribution of the attributes towards the error metric.\n");
    mexPrintf( "SetTCoordsWeight               , real       ... Set/Get the scaling weight contribution of the attribute. These values are used to weight the contribution of the attributes towards the error metric.\n");
    mexPrintf( "SetTensorsWeight               , real       ... Set/Get the scaling weight contribution of the attribute. These values are used to weight the contribution of the attributes towards the error metric.\n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkOBJ_TYPE      *DECI;
  DECI= vtkOBJ_TYPE::New();
  DECI->SetInput( MESH );

  /*Defaults*/
  DECI->SetTargetReduction(0.1);
  DECI->AttributeErrorMetricOff();
  DECI->ScalarsAttributeOff();
  DECI->VectorsAttributeOff();
  DECI->NormalsAttributeOff();
  DECI->TCoordsAttributeOff();
  DECI->TensorsAttributeOff();
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
  Call_0( AttributeErrorMetricOn         );
  Call_0( AttributeErrorMetricOff        );
  Call_0( ScalarsAttributeOn             );
  Call_0( ScalarsAttributeOff            );
  Call_0( VectorsAttributeOn             );
  Call_0( VectorsAttributeOff            );
  Call_0( NormalsAttributeOn             );
  Call_0( NormalsAttributeOff            );
  Call_0( TCoordsAttributeOn             );
  Call_0( TCoordsAttributeOff            );
  Call_0( TensorsAttributeOn             );
  Call_0( TensorsAttributeOff            );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetTargetReduction             , v );
  Call_1( SetAttributeErrorMetric        , v );  
  Call_1( SetScalarsAttribute            , v );  
  Call_1( SetVectorsAttribute            , v );  
  Call_1( SetNormalsAttribute            , v );  
  Call_1( SetTCoordsAttribute            , v );  
  Call_1( SetTensorsAttribute            , v );
  Call_1( SetScalarsWeight               , v );
  Call_1( SetVectorsWeight               , v );
  Call_1( SetNormalsWeight               , v );
  Call_1( SetTCoordsWeight               , v );
  Call_1( SetTensorsWeight               , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
//   Call_1( SetFileName                 , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
//   Call_1( AddPlane                       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
