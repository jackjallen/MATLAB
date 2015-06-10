#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkBooleanOperationPolyDataFilter

#include "vtkBooleanOperationPolyDataFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"




void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
  int                 argN,N,NF,i,n;
  double              v, *xyz;
  char                STR[2000], method[2000];

  if(!nrhs){
    /*
    mexPrintf("vtkBooleanOperationPolyDataFilter( MESH ...\n");
    mexPrintf("\n");
    mexPrintf("SetOperation             	, int       	... Set the boolean operation to perform: Union (0), Intersection (1), Difference (2). Defaults to union. \n");
    mexPrintf("SetOperationToUnion 		, [] (*)    	...\n");
    mexPrintf("SetOperationToIntersection 	, []        	...\n");
    mexPrintf("SetOperationToDifference		, []        	...\n");
    mexPrintf("\n");
    mexPrintf("SetReorientDifferenceCells   	, int       	... Turn on/off cell reorientation of the intersection portion of the surface when the operation is set to DIFFERENCE. Defaults to on. \n");
    mexPrintf("ReorientDifferenceCellsOn 	, [] (*)    	...\n");
    mexPrintf("ReorientDifferenceCellsOff 	, []        	...\n");
    mexPrintf("\n");
    mexPrintf("SetTolerance             	, real (1e-6)   ... Set the tolerance used to determine when a point's absolute distance is considered to be zero. Defaults to 1e-6. \n");
    mexPrintf("\n");*/
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

 
  ALLOCATES();
  vtkPolyData         *MESH1,*MESH2,*S,*MESH_OUT,*MESH_IN;
  vtkOBJ_TYPE	      *BOOL;
  vtkDoubleArray      *ARRAY_ZERO,*ARRAY_ONE;
  vtkIdList           *IdPOINTS;
  
  
  IdPOINTS=vtkIdList::New();
  IdPOINTS->SetNumberOfIds(3);
  MESH1=MESH2vtkPolyData( prhs[0] );
  //if (MESH1->GetCellData()->GetArray("Normals",i)==NULL){mexErrMsgTxt("Error: triNormals field undefined in Mesh1");}
  N=MESH1->GetNumberOfCells();
  ARRAY_ZERO=vtkDoubleArray::New();
  ARRAY_ZERO->SetName( "labels" ); 
  ARRAY_ZERO->SetNumberOfComponents( 1 );
  ARRAY_ZERO->SetNumberOfTuples( N );
  for( n=0 ; n<N ; n++ ) {
      ARRAY_ZERO->SetComponent(n, 0, 0 );
  }
  MESH1->GetCellData()->AddArray( ARRAY_ZERO );
  
  
  MESH2=MESH2vtkPolyData( prhs[1] );
  //if (MESH2->GetCellData()->GetArray("Normals",i)==NULL){mexErrMsgTxt("Error: triNormals field undefined in Mesh2");}
  N= MESH2->GetNumberOfCells();
  ARRAY_ONE=vtkDoubleArray::New();
  ARRAY_ONE->SetName( "labels" ); 
  ARRAY_ONE->SetNumberOfComponents( 1 );
  ARRAY_ONE->SetNumberOfTuples( N );
  for( n=0 ; n<N ; n++ ) {
      ARRAY_ONE->SetComponent(n, 0, 1 );
  }
  MESH2->GetCellData()->AddArray( ARRAY_ONE );

  
  
  BOOL=vtkOBJ_TYPE::New();
  
  /*Defaults*/
  BOOL->ReorientDifferenceCellsOff();
  BOOL->SetTolerance(1e-6);
  BOOL->SetOperationToDifference();
  /*END Defaults*/

  // Direct Difference
  BOOL->SetInput( 0, MESH1 );
  BOOL->SetInput( 1, MESH2 );
  BOOL->Update();
  S=BOOL->GetOutput(0);	

   
  MESH_OUT=vtkPolyData::New();
  MESH_OUT->Allocate();
  MESH_IN=vtkPolyData::New();
  MESH_IN->Allocate();
  
  MESH_OUT->SetPoints(S->GetPoints());  
  MESH_IN->SetPoints(S->GetPoints());
    
  NF=S->GetNumberOfCells();

  for( n= 0; n<NF ; n++ ) 
      {
      if( S->GetCell(n)->GetNumberOfPoints() == 3 ) 
      {
            S->GetCellPoints(n,IdPOINTS);
            if (S->GetCellData()->GetArray("labels",i)->GetComponent(n,0)==0)
            {
                  MESH_OUT->InsertNextCell(5,IdPOINTS);
            }
            else
            {
                  MESH_IN->InsertNextCell(5,IdPOINTS);
            }
      }
      }
  
 
  MESH_OUT->Update();
  MESH_IN->Update();
  plhs[0]= vtkPolyData2MESH( MESH_OUT );
  plhs[3]= vtkPolyData2MESH( MESH_IN );


  MESH_OUT->Delete();
  MESH_IN->Delete();

  
 
//  Reverse Difference

  BOOL->SetInput( 0, MESH2 );
  BOOL->SetInput( 1, MESH1 );
  BOOL->Update();
  S=BOOL->GetOutput(0);	

   
  MESH_OUT=vtkPolyData::New();
  MESH_OUT->Allocate();
  MESH_IN=vtkPolyData::New();
  MESH_IN->Allocate();
  
  MESH_OUT->SetPoints(S->GetPoints());  
  MESH_IN->SetPoints(S->GetPoints());
    
  NF=S->GetNumberOfCells();

  for( n= 0; n<NF ; n++ ) 
      {
      if( S->GetCell(n)->GetNumberOfPoints() == 3 ) 
      {
            S->GetCellPoints(n,IdPOINTS);
            if (S->GetCellData()->GetArray("labels",i)->GetComponent(n,0)==1)
            {
                  MESH_OUT->InsertNextCell(5,IdPOINTS);
            }
            else
            {
                  MESH_IN->InsertNextCell(5,IdPOINTS);
            }
      }
      }
  
 
  MESH_OUT->Update();
  MESH_IN->Update();
  plhs[2]= vtkPolyData2MESH( MESH_OUT );
  plhs[1]= vtkPolyData2MESH( MESH_IN );


  MESH_OUT->Delete();
  MESH_IN->Delete();

  EXIT:
    IdPOINTS->Delete();  
    ARRAY_ZERO->Delete();
    ARRAY_ONE->Delete();
    MESH1->Delete();
    MESH2->Delete();
    BOOL->Delete();


    myFreeALLOCATES();
}



