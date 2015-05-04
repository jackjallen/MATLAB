#include "mex.h"
#include "myMEX.h"
#include "Impldef.h"
#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "vtkOBBTree.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"
#include "vtkGenericCell.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkCellLocator.h"
#include <exception>
#include "vtkPolyDataWriter.h"


void CellsAreInside(vtkPolyData*, vtkPolyData*,  vtkIdList*, vtkIdList* );
void CopyCells(vtkPolyData* , vtkPolyData* , vtkIdList*);
void CopyCellsAndPoints(vtkPolyData*, vtkPolyData*, vtkIdList*);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
 
  bool                CO;

 
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
  
  if (nrhs<3) CO=true;
  else CO= myGetValue( prhs[2] ); 
  
  ALLOCATES();
  
  vtkPolyData         *MESH1;
  MESH1=MESH2vtkPolyData( prhs[0] );
  
  vtkPolyData         *MESH2;
  MESH2=MESH2vtkPolyData( prhs[1] );
 
//   vtkIntersectionPolyDataFilter *PolyDataIntersection = vtkIntersectionPolyDataFilter2::New();
//   PolyDataIntersection->SetInput(0, MESH1);
//   PolyDataIntersection->SetInput(1, MESH2);
//   PolyDataIntersection->SplitFirstOutputOff();
//   PolyDataIntersection->SplitSecondOutputOn();
//   
  
  
 // Set up new poly data for the inputs to build cells and links.
  
  vtkSmartPointer< vtkPolyData > mesh0 = vtkSmartPointer< vtkPolyData >::New();
  mesh0->DeepCopy(MESH1);
  mesh0->SetSource(NULL);

  vtkSmartPointer< vtkPolyData > mesh1 = vtkSmartPointer< vtkPolyData >::New();
  mesh1->DeepCopy(MESH2);
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
    (obbTree1, 0, Implem::FindTriangleIntersections,impl);
  
   vtkPolyDataWriter *W=vtkPolyDataWriter::New();
  W->SetInput(outputIntersection);
  W->SetFileName("InterLines.vtk");
  W->Write();
  W->Delete();
  
  
  // Split the first output 
  
  
  vtkPolyData         *MESH1_T;
  MESH1_T=vtkPolyData::New();
  mesh0->BuildLinks();
  
  impl->SplitMesh(0, MESH1_T, outputIntersection);
  
  
 return;     
  
  
  vtkSmartPointer< vtkIdList > interList = vtkSmartPointer< vtkIdList >::New();
  vtkSmartPointer< vtkIdList > unionList = vtkSmartPointer< vtkIdList >::New();
   
  vtkPolyData         *outputSurface; 

  CellsAreInside(MESH1_T, MESH2, interList,  unionList); 
         
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(MESH1_T);
  if (CO) CopyCellsAndPoints(MESH1_T, outputSurface,  unionList);
  else CopyCells(MESH1_T, outputSurface,  unionList);
  outputSurface->Squeeze();
  plhs[0]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
 
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(MESH1_T);
  if (CO) CopyCellsAndPoints(MESH1_T, outputSurface,  interList);
  else CopyCells(MESH1_T, outputSurface,  interList);
  outputSurface->Squeeze();
  plhs[1]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
  
  if (nlhs>2)
          
  {
  vtkPolyData         *MESH2_T;
  MESH2_T=vtkPolyData::New();
  mesh1->BuildLinks();
  impl->SplitMesh(1, MESH2_T, outputIntersection);
    
 
  interList->Reset();
  unionList->Reset();
  
  CellsAreInside(MESH2_T, MESH1, interList,  unionList); 
   
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(MESH2_T);
  if (CO) CopyCellsAndPoints(MESH2_T, outputSurface,  unionList);
  else CopyCells(MESH2_T, outputSurface,  unionList);
  outputSurface->Squeeze();
  plhs[2]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  
  outputSurface = vtkPolyData::New();
  outputSurface->Allocate(MESH2_T);
  if (CO) CopyCellsAndPoints(MESH2_T, outputSurface,  interList);
  else CopyCells(MESH2_T, outputSurface,  interList);
  outputSurface->Squeeze();
  plhs[3]= vtkPolyData2MESH( outputSurface );
  outputSurface->Delete();
  MESH2_T->Delete();
  }
  
  EXIT:
          impl->PointCellIds[0]->Delete();
          impl->PointCellIds[1]->Delete();
          delete impl;
          
          MESH1_T->Delete();
          MESH2->Delete();
          MESH1->Delete();

    myFreeALLOCATES();
}


// //-----------------------------------------------------------------------------
  void CellsAreInside(vtkPolyData* A, vtkPolyData* B,  vtkIdList* interList, vtkIdList* unionList)
  
  {
  
 
  vtkPolyDataNormals  *NORMALS;
  NORMALS= vtkPolyDataNormals::New();
  NORMALS->SetInput( B );
  NORMALS->ComputeCellNormalsOn();
  NORMALS->ComputePointNormalsOn();
  NORMALS->ConsistencyOff();
  NORMALS->AutoOrientNormalsOff();
  NORMALS->SplittingOff();
  NORMALS->Update();
  
  vtkPolyData* B_Normals;
  B_Normals=NORMALS->GetOutput();
  
  vtkCellLocator      *LOC;
  LOC= vtkCellLocator::New();
  LOC->SetDataSet( B_Normals );
  LOC->CacheCellBoundsOn();
  LOC->SetNumberOfCellsPerBucket( 2 );
  LOC->BuildLocator();
  LOC->Update();
  

  
  double dist, x[3],CellCent[3],closestPoint[3],N[3], *NN, pcoords[3], weigths[3], side,inumCellPts, bb[6];
  int j, sub;
  vtkIdList *cellPts, *CellPointIds;
  vtkIdType i, numCellPts, ptId, cellId, cell;
  vtkGenericCell  *Gcell;
  
  B_Normals->GetBounds(bb);
  Gcell = vtkGenericCell::New();
  CellPointIds = vtkIdList::New(); 
  cellPts=vtkIdList::New();


  
  for ( cellId = 0; cellId < A->GetNumberOfCells(); cellId++ )
  { 
    // Calculo del centro  
    A->GetCellPoints( cellId, cellPts );
    numCellPts=A->GetCell(cellId)->GetNumberOfPoints();
    inumCellPts  = 1/((double)numCellPts);
    CellCent[0]=0;CellCent[1]=0;CellCent[2]=0;
    for ( i = 0; i < numCellPts; i++ )
      {
        ptId = cellPts->GetId( i );
        A->GetPoint( ptId, x );
        for ( j = 0; j < 3; j++ )
            {
            CellCent[j]+=x[j];
            }
      }
      

    for (  j = 0; j < 3; j++ )
      {
        CellCent[j]=CellCent[j]*inumCellPts;
      }
    
    // Distancia del centro al otro MESH

    if( CellCent[0] < bb[0] || CellCent[0] > bb[1] ||
        CellCent[1] < bb[2] || CellCent[1] > bb[3] ||
        CellCent[2] < bb[4] || CellCent[2] > bb[5] ) 
    {
        unionList->InsertNextId( cellId );
        continue;
    }
    
     
    LOC->FindClosestPoint( CellCent , closestPoint , Gcell , cell, sub, dist);
            
    if( dist < 1e-10 ) 
    {
      // Sobre el poligono
      unionList->InsertNextId( cellId );
      interList->InsertNextId( cellId );
      continue;
    }
    
    B_Normals->GetCell(cell)->EvaluatePosition( closestPoint, NULL, sub, pcoords, dist, weigths );
    
    if( weigths[0]<1e-10 || weigths[1]<1e-10 || weigths[2]<1e-10 )
    {   
      // Sobre una arista
      B_Normals->GetCellPoints( cell , CellPointIds );
      NN = B_Normals->GetPointData()->GetNormals()->GetTuple3( CellPointIds->GetId(0) );
      for( j=0 ; j<3; j++ ){  N[j]  = NN[j]*weigths[0];  }
      NN = B_Normals->GetPointData()->GetNormals()->GetTuple3( CellPointIds->GetId(1) );
      for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[1];  }
      NN = B_Normals->GetPointData()->GetNormals()->GetTuple3( CellPointIds->GetId(2) );
      for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[2];  }
    } 
    else 
    {
      NN = B_Normals->GetCellData()->GetNormals()->GetTuple3(cell);
      for( j=0 ; j<3; j++ )
        {  
          N[j] = NN[j];  
        }
    }

    side= 0;
    for( j=0; j<3; j++ )
    {
      side += (closestPoint[j]-CellCent[j])*N[j];
    }
    
    if( side < 0 )
    {
      unionList->InsertNextId( cellId );
    } 
    else 
    {
      interList->InsertNextId( cellId );
    }
         
  }
  cellPts->Delete();
  CellPointIds->Delete();
  Gcell->Delete();
  LOC->Delete();
  NORMALS->Delete();
  
  
  }  

  
//-----------------------------------------------------------------------------
void CopyCells(vtkPolyData* in, vtkPolyData* out,
               vtkIdList* cellIds){
  
  vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell> ::New();
  vtkSmartPointer< vtkIdList > newCellPts = vtkSmartPointer< vtkIdList >::New();
  out->SetPoints(in->GetPoints());
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

//-----------------------------------------------------------------------------
void CopyCellsAndPoints(vtkPolyData* in, vtkPolyData* out,
               vtkIdList* cellIds)

{
  
   vtkIdType numPts = in->GetNumberOfPoints();

  if ( out->GetPoints() == NULL)
    {
    vtkSmartPointer< vtkPoints > points = vtkSmartPointer< vtkPoints >::New();
    out->SetPoints( points );
    }

  vtkPoints *newPoints = out->GetPoints();

  vtkSmartPointer< vtkIdList > pointMap = vtkSmartPointer< vtkIdList >::New();
  pointMap->SetNumberOfIds( numPts );
  for ( vtkIdType i = 0; i < numPts; i++ )
    {
    pointMap->SetId(i, -1);
    }

  // Filter the cells
  vtkSmartPointer< vtkGenericCell > cell = vtkSmartPointer< vtkGenericCell> ::New();
  vtkSmartPointer< vtkIdList > newCellPts = vtkSmartPointer< vtkIdList >::New();
  for ( vtkIdType cellId = 0; cellId < cellIds->GetNumberOfIds(); cellId++ )
    {
    in->GetCell( cellIds->GetId( cellId ), cell );
    vtkIdList *cellPts = cell->GetPointIds();
    vtkIdType numCellPts = cell->GetNumberOfPoints();

    for ( vtkIdType i = 0; i < numCellPts; i++ )
      {
      vtkIdType ptId = cellPts->GetId( i );
      vtkIdType newId = pointMap->GetId( ptId );
      if ( newId < 0 )
        {
        double x[3];
        in->GetPoint( ptId, x );
        newId = newPoints->InsertNextPoint( x );
        pointMap->SetId( ptId, newId );
        }
      newCellPts->InsertId( i, newId );
      }

    
    vtkIdType newCellId = out->InsertNextCell( cell->GetCellType(), newCellPts );
    newCellPts->Reset();
    } // for all cells

}
