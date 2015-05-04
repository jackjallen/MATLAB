#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkShrinkPolyData 
#include "vtkShrinkPolyData.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

   if(!nrhs){
    mexPrintf("vtkShrinkPolyData( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetShrinkFactor            , real (0.8)    ... Set the fraction of shrink for each cell.\n");  
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }




  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkShrinkPolyData      *SHRI;
  SHRI= vtkShrinkPolyData::New();
  SHRI->SetInput( MESH );

  /*Defaults*/             
  SHRI->SetShrinkFactor( 0.8 );     
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
      CallMethod( SHRI , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( SHRI , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( SHRI , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( SHRI , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  SHRI->Update();
  plhs[0]= vtkPolyData2MESH( SHRI->GetOutput() );

  EXIT:
    SHRI->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetShrinkFactor     , v );
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
