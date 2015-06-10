#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkIntersectionPolyDataFilter

#include "vtkIntersectionPolyDataFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

#include "vtkSmartPointer.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
   	
    mexPrintf("vtkIntersectionPolyDataFilter( MESH1, MESH2 ...\n");
    mexPrintf("\n");
    mexPrintf("SetSplitFirstOutput      , logical    ... If on, the second output will be the first input mesh split by the intersection with the second input mesh. Defaults to on. \n");
    mexPrintf("SplitFirstOutputOn       , [] (*)     ...\n");
    mexPrintf("SplitFirstOutputOff      , []         ...\n");
    mexPrintf("\n");
    mexPrintf("SetSplitSecondOutput      , logical    ... If on, the third output will be the second input mesh split by the intersection with the first input mesh. Defaults to on.  \n");
    mexPrintf("SplitSecondOutputOn       , [] (*)     ...\n");
    mexPrintf("SplitSecondOutputOff      , []         ...\n");
    mexPrintf("\n");
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  vtkPolyData         *MESH1,*MESH2;
  vtkOBJ_TYPE	      *INTER;

  MESH1 = MESH2vtkPolyData( prhs[0] );
  MESH2 = MESH2vtkPolyData( prhs[1] );

  INTER = vtkOBJ_TYPE::New();
  INTER->SetInput( (int)0, MESH1 );
  INTER->SetInput( (int)1, MESH2 );
  
  /*Defaults*/
  INTER->SplitFirstOutputOn();
  INTER->SplitSecondOutputOn();
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 2;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
    /*   myErrMsgTxt( "No keywords." ); */
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( INTER , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( INTER , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( INTER , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( INTER , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
  

  INTER->Update();  
  plhs[0]= vtkPolyData2MESH( INTER->GetOutput(0) );
  if (nlhs>1){
  plhs[1]= vtkPolyData2MESH( INTER->GetOutput(1) );
  if (nlhs>2){
  plhs[2]= vtkPolyData2MESH( INTER->GetOutput(2) );
  }
  }

  EXIT:
    INTER->Delete();
    MESH1->Delete();
    MESH2->Delete();
    myFreeALLOCATES();
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  /*[]-valued methods*/
  Call_0( SplitFirstOutputOn   );
  Call_0( SplitFirstOutputOff  );
  Call_0( SplitSecondOutputOn  );
  Call_0( SplitSecondOutputOff );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  /*Real valued methods*/
  Call_1( SetSplitFirstOutput  , v );
  Call_1( SetSplitSecondOutput , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
  /*String valued methods*/
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
  /*Array valued methods*/
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
