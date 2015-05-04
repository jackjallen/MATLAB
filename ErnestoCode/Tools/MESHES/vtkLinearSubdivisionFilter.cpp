#include "mex.h"
#include "myMEX.h"

#define vtkOBJ_TYPE      vtkLinearSubdivisionFilter

#include "vtkLinearSubdivisionFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  if(!nrhs){
    mexPrintf("vtkLinearSubdivisionFilter( MESH )\n");
    mexPrintf("\n");
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  vtkPolyData         		*MESH;
  vtkOBJ_TYPE			*SUB;

  MESH = MESH2vtkPolyData( prhs[0] );
  
  SUB = vtkOBJ_TYPE::New();
  SUB->SetInput( MESH );
  
  /*Defaults*/
  /*END Defaults*/
  
  /*Parsing arguments*/
  /*END Parsing arguments*/
  
  SUB->Update();  
  plhs[0]= vtkPolyData2MESH( SUB->GetOutput() );

  EXIT:
    SUB->Delete();
    MESH->Delete();
    myFreeALLOCATES();
}

