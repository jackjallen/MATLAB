#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define  vtkOBJ_TYPE      vtkClipPolyData 
#include "vtkClipPolyData.h"
#include "vtkPlane.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int                 argN;
  double              v, *xyz, p[3], n[3];
  char                STR[2000], method[2000];

  if((!nrhs) || (nrhs < 2) ){
    mexPrintf("vtkClipPolyData( MESH , [ point ; normal ] ...\n");
    mexPrintf("\n");
    mexPrintf( "SetValue                       , real (0.0) ... Set the clipping value of the implicit function (if clipping with implicit function) or scalar value (if clipping with scalars). The default value is 0.0.\n");
    mexPrintf("\n");
    mexPrintf( "SetInsideOut                   , logical    ... Set/Get the InsideOut flag. When off, a vertex is considered inside the implicit function if its value is greater than the Value ivar. When InsideOutside is turned on, a vertex is considered inside the implicit function if its implicit function value is less than or equal to the Value ivar. InsideOut is off by default.\n");
    mexPrintf( "InsideOutOn                    , []         ...\n");  
    mexPrintf( "InsideOutOff                   , [] (*)     ...\n");  
    mexPrintf("\n");
    mexPrintf( "SetGenerateClipScalars         , logical    ... If this flag is enabled, then the output scalar values will be interpolated from the implicit function values, and not the input scalar data. If you enable this flag but do not provide an implicit function an error will be reported.\n");
    mexPrintf( "GenerateClipScalarsOn          , []         ...\n");  
    mexPrintf( "GenerateClipScalarsOff         , [] (*)     ...\n");  
    mexPrintf("\n");
    if( nlhs ){ plhs[0]= mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

  ALLOCATES();
  vtkPolyData         *MESH;
  MESH= MESH2vtkPolyData( prhs[0] );

  vtkOBJ_TYPE      *CLIP;
  CLIP= vtkOBJ_TYPE::New();
  CLIP->SetInput( MESH );

  
  vtkPlane *PLANE; PLANE = vtkPlane::New();
  if( mxIsNumeric( prhs[1] ) &&  mxGetNumberOfElements( prhs[1] ) == 6 && mxGetM( prhs[1] ) == 2 && mxGetN( prhs[1] ) == 3 ){
    p[0] = myGetValueIth( prhs[1] , 0 );
    p[1] = myGetValueIth( prhs[1] , 2 );
    p[2] = myGetValueIth( prhs[1] , 4 );
    n[0] = myGetValueIth( prhs[1] , 1 );
    n[1] = myGetValueIth( prhs[1] , 3 );
    n[2] = myGetValueIth( prhs[1] , 5 );
  } else if( mxIsNumeric( prhs[1] ) &&  mxGetNumberOfElements( prhs[1] ) == 9 && mxGetM( prhs[1] ) == 3 && mxGetN( prhs[1] ) == 3 ){
    xyz = myGetPr( prhs[1] );
    p[0] = xyz[0];
    p[1] = xyz[3];
    p[2] = xyz[6];
    
    n[0] = ( xyz[4] - xyz[3] )*( xyz[8] - xyz[6] ) - ( xyz[7] - xyz[6] )*( xyz[5] - xyz[3] );
    n[1] = ( xyz[7] - xyz[6] )*( xyz[2] - xyz[0] ) - ( xyz[8] - xyz[6] )*( xyz[1] - xyz[0] );
    n[2] = ( xyz[1] - xyz[0] )*( xyz[5] - xyz[3] ) - ( xyz[4] - xyz[3] )*( xyz[2] - xyz[0] );
  } else {
  /*   myErrMsgTxt( "Incorrect Second Argument. A matrix [ point;normal] or [point;point;point] is expected." ); */
  }
  
/* //   mexPrintf("SetOrigin  [ %.15g  %.15g  %.15g ] \n" , p[0],p[1],p[2] );
//   mexPrintf("SetNormal  [ %.15g  %.15g  %.15g ] \n" , n[0],n[1],n[2] ); */
  
  PLANE->SetOrigin(p); 
  PLANE->SetNormal(n);
  
  CLIP->SetClipFunction( PLANE );
  
  /*Defaults*/
  CLIP->SetValue( 0.0 );
  CLIP->InsideOutOff();
  CLIP->GenerateClipScalarsOff();
  /*END Defaults*/
  
  /*Parsing arguments*/
  argN = 2;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
     /*  myErrMsgTxt( "No keywords." ); */
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( CLIP , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( CLIP , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( CLIP , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( CLIP , method , xyz );
    }
    argN++;

  }
  /*END Parsing arguments*/
    
  CLIP->Update();
  plhs[0]= vtkPolyData2MESH( CLIP->GetOutput() );

  EXIT:
    CLIP->Delete();
    MESH->Delete();
    PLANE->Delete();
    myFreeALLOCATES();
  
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  Call_0( InsideOutOn               );
  Call_0( InsideOutOff              );
  Call_0( GenerateClipScalarsOn     );
  Call_0( GenerateClipScalarsOff    );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  Call_1( SetValue                , v );
  Call_1( SetInsideOut            , v );
  Call_1( SetGenerateClipScalars  , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
/* //   Call_1( SetFileName                 , v ); */
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
/* //   Call_1( SetClosestPoint                       , v ); */
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}
