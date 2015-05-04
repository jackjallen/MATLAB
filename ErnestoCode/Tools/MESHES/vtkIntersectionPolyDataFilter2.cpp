#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "MESH2vtkPolyData.h"

#include "vtkIntersectionPolyDataFilter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDelaunay2D.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkObjectFactory.h"
#include "vtkOBBTree.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkTransform.h"
#include "vtkTriangle.h"
#include "Impldef.h"

#include <map>
#include <queue>

// void CellsAreInside(vtkPolyData*, vtkPolyData*,  vtkIdList*, vtkIdList* );
// void CopyCells(vtkPolyData* , vtkPolyData* , vtkIdList*);
// void CopyCellsAndPoints(vtkPolyData*, vtkPolyData*, vtkIdList*);


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){


  if(!nrhs){
    /*
    mexPrintf("vtkBooleanOperationPolyDataFilter( MESH1, MESH2, [CleanOut]\n");
    mexPrintf("\n");
    mexPrintf("Cleanout             	, boolean       	... (*) If true the output meshes contain only cell points \n");
    mexPrintf("                                                     If false the output meshes contain all input points  \n");
    mexPrintf("\n");*/
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }
  
  ALLOCATES();
  // Set up new poly data for the inputs to build cells and links.
  
  vtkSmartPointer< vtkPolyData > mesh0 = vtkSmartPointer< vtkPolyData >::New();
  mesh0->DeepCopy(MESH2vtkPolyData( prhs[0] ));
  mesh0->SetSource(NULL);

  vtkSmartPointer< vtkPolyData > mesh1 = vtkSmartPointer< vtkPolyData >::New();
  mesh1->DeepCopy(MESH2vtkPolyData( prhs[1] ));
  mesh1->SetSource(NULL);

  vtkPolyData *outputIntersection = vtkPolyData::New();
  vtkSmartPointer< vtkPoints > outputIntersectionPoints =
  vtkSmartPointer< vtkPoints >::New();
  outputIntersection->SetPoints(outputIntersectionPoints);

  // Find the triangle-triangle intersections between mesh0 and mesh1
  vtkSmartPointer< vtkOBBTree > obbTree0 = vtkSmartPointer< vtkOBBTree >::New();
  obbTree0->SetDataSet(mesh0);
  obbTree0->SetNumberOfCellsPerNode(10);
  obbTree0->SetMaxLevel(1000000);
  obbTree0->SetTolerance(1e-6);
  obbTree0->AutomaticOn();
  obbTree0->BuildLocator();

  vtkSmartPointer< vtkOBBTree > obbTree1 = vtkSmartPointer< vtkOBBTree >::New();
  obbTree1->SetDataSet(mesh1);
  obbTree1->SetNumberOfCellsPerNode(10);
  obbTree1->SetMaxLevel(1000000);
  obbTree1->SetTolerance(1e-6);
  obbTree1->AutomaticOn();
  obbTree1->BuildLocator();

  // Set up the structure for determining exact triangle-triangle
  // intersections.
  Implem *impl = new Implem();
  impl->Mesh[0]  = mesh0;
  impl->Mesh[1]  = mesh1;
  impl->OBBTree1 = obbTree1;

  vtkSmartPointer< vtkCellArray > lines = vtkSmartPointer< vtkCellArray >::New();
  outputIntersection->SetLines(lines);
  impl->IntersectionLines = lines;

  // Add cell data arrays that map the intersection line to the cells
  // it splits.
  impl->CellIds[0] = vtkIdTypeArray::New();
  impl->CellIds[0]->SetName("Input0CellID");
  outputIntersection->GetCellData()->AddArray(impl->CellIds[0]);
  impl->CellIds[0]->Delete();
  impl->CellIds[1] = vtkIdTypeArray::New();
  impl->CellIds[1]->SetName("Input1CellID");
  outputIntersection->GetCellData()->AddArray(impl->CellIds[1]);
  impl->CellIds[1]->Delete();

  impl->PointCellIds[0] = vtkIdTypeArray::New();
  impl->PointCellIds[0]->SetName("PointCellsIDs");
  impl->PointCellIds[1] = vtkIdTypeArray::New();
  impl->PointCellIds[1]->SetName("PointCellsIDs");

  double bounds0[6], bounds1[6];
  mesh0->GetBounds(bounds0);
  mesh1->GetBounds(bounds1);
  for (int i = 0; i < 3; i++)
    {
    int minIdx = 2*i;
    int maxIdx = 2*i+1;
    if (bounds1[minIdx] < bounds0[minIdx])
      {
      bounds0[minIdx] = bounds1[minIdx];
      }
    if (bounds1[maxIdx] > bounds0[maxIdx])
      {
      bounds0[maxIdx] = bounds1[maxIdx];
      }
    }

  vtkSmartPointer< vtkPointLocator > pointMerger =
    vtkSmartPointer< vtkPointLocator >::New();
  pointMerger->SetTolerance(1e-6);
  pointMerger->InitPointInsertion(outputIntersection->GetPoints(), bounds0);
  impl->PointMerger = pointMerger;

  // This performs the triangle intersection search
  obbTree0->IntersectWithOBBTree
    (obbTree1, 0, Implem::FindTriangleIntersections,
     impl);

//   // Split the first output if so desired
//   if ( this->SplitFirstOutput )
//     {
//     mesh0->BuildLinks();
//     impl->SplitMesh(0, outputPolyData0, outputIntersection);
//     }
//   else
//     {
//     outputPolyData0->ShallowCopy( mesh0 );
//     }
// 
//   // Split the second output if desired
//   if ( this->SplitSecondOutput )
//     {
//     mesh1->BuildLinks();
//     impl->SplitMesh(1, outputPolyData1, outputIntersection);
//     }
//   else
//     {
//     outputPolyData1->ShallowCopy( mesh1 );
//     }
// 
//   impl->PointCellIds[0]->Delete();
//   impl->PointCellIds[1]->Delete();
//   delete impl;
// 
// 
// 
// 
  
    EXIT:

    myFreeALLOCATES();
}
