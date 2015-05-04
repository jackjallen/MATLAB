#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkPolyDataWriter
#include "vtkPolyDataWriter.h"
#include "MESH2vtkPolyData.h"

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];
  
  if(!nrhs){
    mexPrintf("vtkPolyDataWriter( MESH ...\n");
    mexPrintf("\n");
    mexPrintf("SetFileName            , str(kk.vtp)...\n");
    mexPrintf("\n");
    mexPrintf("SetFileType            , int        ...\n");
    mexPrintf("SetFileTypeToASCII     , []         ...\n");
    mexPrintf("SetFileTypeToBinary    , [] (*)     ...\n");
    mexPrintf("\n");
    mexPrintf("SetHeader              , str        ...\n");
    mexPrintf("SetScalarsName         , str        ...\n");
    mexPrintf("SetVectorsName         , str        ...\n");
    mexPrintf("SetTensorsName         , str        ...\n");
    mexPrintf("SetNormalsName         , str        ...\n");
    mexPrintf("SetTCoordsName         , str        ...\n");
    mexPrintf("SetLookupTableName     , str        ...\n");
    mexPrintf("SetFieldDataName       , str        ...\n");
    mexPrintf("\n");
    return;
  }
  
  ALLOCATES();
  vtkPolyData         *MESH;
  
  MESH= MESH2vtkPolyData( prhs[0] );
  
  vtkPolyDataWriter *W = vtkPolyDataWriter::New();
  W->SetInput( MESH );

  /*Defaults*/
  W->SetFileName("kk.vtp");     
  W->SetFileTypeToBinary();
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
      CallMethod( W , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( W , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( W , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( W , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/

  W->Write();
    
  EXIT:
    W->Delete();
    MESH->Delete();
    myFreeALLOCATES();

}



void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( SetFileTypeToASCII      );
  Call_0( SetFileTypeToBinary     );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real  v ){
  Call_1( SetFileType            , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
  Call_1( SetFileName            , v );
  Call_1( SetHeader              , v );
  Call_1( SetScalarsName         , v );
  Call_1( SetVectorsName         , v );
  Call_1( SetTensorsName         , v );
  Call_1( SetNormalsName         , v );
  Call_1( SetTCoordsName         , v );
  Call_1( SetLookupTableName     , v );
  Call_1( SetFieldDataName       , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
//  Call_1( SetFileType            , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
