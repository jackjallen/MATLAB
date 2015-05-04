#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "vtkPolyDataReader.h"
#include "vtkTriangleFilter.h"

#include "vtkPolyData2MESH.h"

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){
  char                STR[2000];
  
  if(!nrhs){
    mexPrintf("vtkPolyDataReader( FileName )\n");
    mexPrintf("\n");
    return;
  }
  
  ALLOCATES();

  mxGetString( prhs[0], STR, 1999 );
  
  
  vtkPolyDataReader *R = vtkPolyDataReader::New();
  R->SetFileName( STR );

  R->Update();
  

  vtkTriangleFilter   *T = vtkTriangleFilter::New();
  T->SetInput( R->GetOutput() );
  T->PassVertsOff();
  T->PassLinesOff();
  T->Update();  
  
  plhs[0] = vtkPolyData2MESH( T->GetOutput() );
  
  EXIT:
    R->Delete();
    myFreeALLOCATES();

}

