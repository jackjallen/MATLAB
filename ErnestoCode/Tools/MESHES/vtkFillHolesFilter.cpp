#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkFillHolesFilter 
#include "vtkFillHolesFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

   if(!nrhs){
    mexPrintf("vtkFillHolesFilter( MESH ,\n");
    mexPrintf("\n");
    mexPrintf( "SetHoleSize      , real (1e-8)     ...Specify the maximum hole size to fill. This is represented as a radius to the bounding circumsphere containing the hole. Note that this is an approximate area; the actual area cannot be computed without first triangulating the hole.\n");
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }


  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkFillHolesFilter      *FILL;
  FILL= vtkFillHolesFilter::New();
  FILL->SetInput( MESH );

  /*Defaults*/             
  FILL->SetHoleSize( 1e8 );     
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
      CallMethod( FILL , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( FILL , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( FILL , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( FILL , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  FILL->Update();
  plhs[0]= vtkPolyData2MESH( FILL->GetOutput() );

  EXIT:
    FILL->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
//   Call_0( GenerateVerticesOff         );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetHoleSize          , v );
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
