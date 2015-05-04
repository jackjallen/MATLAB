#include "mex.h"

#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"
#include "vtkAreInside.cpp"

#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkPlane.h"
#include "vtkIdList.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  vtkPolyData       *M1;
  vtkPolyData       *M2;
  vtkPlane          *P1;
  vtkPlane          *P2;
  vtkClipPolyData   *C1;
  vtkClipPolyData   *C2;
  long              Noc, c;
  double            x0[3], x1[3], x2[3], N[3];
  vtkIdList         *NODES;
  vtkAppendPolyData *A;
  vtkCleanPolyData  *C;

  NODES = vtkIdList::New();

  M1= MESH2vtkPolyData( prhs[0] );
  M2= MESH2vtkPolyData( prhs[1] );
  
  P1= vtkPlane::New();
  C1= vtkClipPolyData::New();
  C1->SetInput( M1 );
  C1->SetClipFunction( P1 );
  C1->SetValue(0);

  P2= vtkPlane::New();
  C2= vtkClipPolyData::New();
  C2->SetInput( M1 );
  C2->SetClipFunction( P2 );
  C2->SetValue(0);

  A= vtkAppendPolyData::New();
  A->AddInput( C1->GetOutput() );
  A->AddInput( C2->GetOutput() );
  
  C= vtkCleanPolyData::New();
  C->SetInput( A->GetOutput() );
  
  Noc= M2->GetNumberOfCells();
  for(c=0; c<Noc; c++ ){
    M2->GetCellPoints( c , NODES );
    M2->GetPoint( NODES->GetId(0) , x0 );
    M2->GetPoint( NODES->GetId(1) , x1 );
    M2->GetPoint( NODES->GetId(2) , x2 );
    
    N[0]= ( x1[1] - x0[1] )*( x2[2] - x0[2] ) - ( x1[2] - x0[2] )*( x2[1] - x0[1] );
    N[1]= ( x1[2] - x0[2] )*( x2[0] - x0[0] ) - ( x1[0] - x0[0] )*( x2[2] - x0[2] );
    N[2]= ( x1[0] - x0[0] )*( x2[1] - x0[1] ) - ( x1[1] - x0[1] )*( x2[0] - x0[0] );

    P1->SetOrigin( x0[0], x0[1], x0[2] );
    P1->SetNormal(  N[0],  N[1],  N[2] );
    P2->SetOrigin( x0[0], x0[1], x0[2] );
    P2->SetNormal( -N[0], -N[1], -N[2] );
    C1->Update();
    C2->Update();

    A->Update();
    C->Update();
    
    M1->DeepCopy( C->GetOutput() );
  }
  
  NODES->Delete();
  P1->Delete();
  C1->Delete();
  P2->Delete();
  C2->Delete();
  A->Delete();
  C->Delete();
  
  M2->Delete();
 
  plhs[0]= vtkPolyData2MESH( M1 );
  
  M1->Delete();
}