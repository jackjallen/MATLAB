#include "mex.h"
#include "myMEX.h"

#define   real                double
#define   mxREAL_CLASS        mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkDensifyPolyData 
#include "vtkDensifyPolyData.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

   if(!nrhs){
    mexPrintf("vtkDensifyPolyData( MESH ...\n");
    mexPrintf("\n");
    mexPrintf( "SetNumberOfSubdivisions           , int (1)    ... Number of recursive subdivisions. Initial value is 1.\n");  
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }




  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkDensifyPolyData      *DENS;
  DENS= vtkDensifyPolyData::New();
  DENS->SetInput( MESH );

  /*Defaults*/             
  DENS->SetNumberOfSubdivisions( 1 );     
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
      CallMethod( DENS , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( DENS , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( DENS , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( DENS , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  DENS->Update();
  plhs[0]= vtkPolyData2MESH( DENS->GetOutput() );

  EXIT:
    DENS->Delete();
    MESH->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetNumberOfSubdivisions    , v );
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
