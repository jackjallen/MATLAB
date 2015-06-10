#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkIntersectionPolyDataFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"
#include "vtkDistancePolyDataFilter.h"

#include "vtkSmartPointer.h"

void SortPolyData(vtkPolyData* ,  vtkIdList* , vtkIdList* );
void CopyCells(vtkPolyData* , vtkPolyData* , vtkIdList*);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
  int                 argN,NF,i,n;
  double              v;


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
  
  vtkPolyData         *MESH1;
  MESH1=MESH2vtkPolyData( prhs[0] );

  vtkPolyData         *MESH2;
  MESH2=MESH2vtkPolyData( prhs[1] );
  
  vtkIntersectionPolyDataFilter *PolyDataIntersection = vtkIntersectionPolyDataFilter::New();
  PolyDataIntersection->SetInput(0, MESH1);
  PolyDataIntersection->SetInput(1, MESH2);
  PolyDataIntersection->SplitFirstOutputOn();
  PolyDataIntersection->SplitSecondOutputOn();
  PolyDataIntersection->Update();
  
  vtkDistancePolyDataFilter     *PolyDataDistance = vtkDistancePolyDataFilter::New();
  PolyDataDistance->ComputeSecondDistanceOn();
  PolyDataDistance->SetInputConnection(0, PolyDataIntersection->GetOutputPort( 1 ));
  PolyDataDistance->SetInputConnection(1, PolyDataIntersection->GetOutputPort( 2 )); 
  PolyDataDistance->Update();
  
  vtkPolyData* pd0 = PolyDataDistance->GetOutput();
  vtkPolyData* pd1 = PolyDataDistance->GetSecondDistanceOutput();
  
  
    
  vtkSmartPointer< vtkIdList > interList = vtkSmartPointer< vtkIdList >::New();
  vtkSmartPointer< vtkIdList > unionList = vtkSmartPointer< vtkIdList >::New();

  SortPolyData(pd0, interList, unionList);
  
  vtkPolyData         *outputSurface; 
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(pd0);
  outputSurface->SetPoints(pd0->GetPoints());
  CopyCells(pd0, outputSurface,  unionList);
  outputSurface->Squeeze();
  plhs[0]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(pd0);
  outputSurface->SetPoints(pd0->GetPoints());
  CopyCells(pd0, outputSurface,  interList);
  outputSurface->Squeeze();
  plhs[1]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
  interList->Reset();
  unionList->Reset();
  
  SortPolyData(pd1, interList, unionList);
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(pd1);
  outputSurface->SetPoints(pd1->GetPoints());
  CopyCells(pd1, outputSurface,  unionList);
  outputSurface->Squeeze();
  plhs[2]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(pd1);
  outputSurface->SetPoints(pd1->GetPoints());
  CopyCells(pd1, outputSurface,  interList);
  outputSurface->Squeeze();
  plhs[3]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();

  EXIT:

          PolyDataDistance->Delete();
          PolyDataIntersection->Delete();
          MESH2->Delete();
          MESH1->Delete();

    myFreeALLOCATES();
}


//-----------------------------------------------------------------------------
void SortPolyData(vtkPolyData* input,  vtkIdList* interList, vtkIdList* unionList){
  int numCells = input->GetNumberOfCells();
  
          
  double* dist = static_cast<double*>(input->GetPointData()->GetArray("Distance")->WriteVoidPointer(0, 0));
  vtkDoubleArray *distArray = vtkDoubleArray::SafeDownCast( input->GetCellData()->GetArray("Distance") );

  for (int cid = 0; cid < numCells; cid++)
    {

    if ( distArray->GetValue( cid ) > 0 )
      {
      unionList->InsertNextId( cid );
      }
    else
      {
      interList->InsertNextId( cid );
      }
    }
}

//-----------------------------------------------------------------------------
void CopyCells(vtkPolyData* in, vtkPolyData* out,
               vtkIdList* cellIds){
  
  vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell> ::New();
  vtkSmartPointer< vtkIdList > newCellPts = vtkSmartPointer< vtkIdList >::New();
  for ( vtkIdType cellId = 0; cellId < cellIds->GetNumberOfIds(); cellId++ )
    {
    in->GetCell( cellIds->GetId( cellId ), cell );
    vtkIdList *cellPts = cell->GetPointIds();
    vtkIdType numCellPts = cell->GetNumberOfPoints();

    for ( vtkIdType i = 0; i < numCellPts; i++ )
      {
      newCellPts->InsertId( i, cellPts->GetId( i ) );
      }

    
    vtkIdType newCellId = out->InsertNextCell( cell->GetCellType(), newCellPts );
    newCellPts->Reset();
    } // for all cells
  }

// //-----------------------------------------------------------------------------
// void CopyCells(vtkPolyData* in, vtkPolyData* out,
//                vtkIdList* cellIds){
//   
//    vtkIdType numPts = in->GetNumberOfPoints();
// 
//   if ( out->GetPoints() == NULL)
//     {
//     vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
//     out->SetPoints( points );
//     }
// 
//   vtkPoints *newPoints = out->GetPoints();
// 
//   vtkSmartPointer< vtkIdList > pointMap = vtkSmartPointer< vtkIdList >::New();
//   pointMap->SetNumberOfIds( numPts );
//   for ( vtkIdType i = 0; i < numPts; i++ )
//     {
//     pointMap->SetId(i, -1);
//     }
// 
//   // Filter the cells
//   vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell> ::New();
//   vtkSmartPointer< vtkIdList > newCellPts = vtkSmartPointer< vtkIdList >::New();
//   for ( vtkIdType cellId = 0; cellId < cellIds->GetNumberOfIds(); cellId++ )
//     {
//     in->GetCell( cellIds->GetId( cellId ), cell );
//     vtkIdList *cellPts = cell->GetPointIds();
//     vtkIdType numCellPts = cell->GetNumberOfPoints();
// 
//     for ( vtkIdType i = 0; i < numCellPts; i++ )
//       {
//       vtkIdType ptId = cellPts->GetId( i );
//       vtkIdType newId = pointMap->GetId( ptId );
//       if ( newId < 0 )
//         {
//         double x[3];
//         in->GetPoint( ptId, x );
//         newId = newPoints->InsertNextPoint( x );
//         pointMap->SetId( ptId, newId );
//         }
//       newCellPts->InsertId( i, newId );
//       }
// 
//     
//     vtkIdType newCellId = out->InsertNextCell( cell->GetCellType(), newCellPts );
//     newCellPts->Reset();
//     } // for all cells
// 
// }
