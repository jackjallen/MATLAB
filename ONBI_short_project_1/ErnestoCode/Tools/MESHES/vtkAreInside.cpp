#include "vtkTriangleFilter.h"
#include "vtkPolyData.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkDataArray.h"
#include "vtkCellLocator.h"
#include "vtkCellData.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"

int IsInside( vtkPolyData *SURF, int cell, double *x, double *closestPoint){;
  double  N[3], *NN, pcoords[3], weigths[3], side, distance;
  int     sub, j;
  vtkIdList           *NODES;

  NODES = vtkIdList::New();

  side= 0;
  SURF->GetCell( cell )->EvaluatePosition( closestPoint, NULL, sub, pcoords, distance, weigths );

  if( weigths[0]<1e-12 || weigths[1]<1e-12 || weigths[2]<1e-12 ){
    SURF->GetCellPoints( cell , NODES );
    NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(0) );
    for( j=0 ; j<3; j++ ){  N[j]  = NN[j]*weigths[0];  }
    NN= SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(1) );
    for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[1];  }
    NN= SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(2) );
    for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[2];  }
    for( j=0; j<3; j++ ){
      side += (closestPoint[j]-x[j])*N[j];
    }
  } else {
    NN = SURF->GetCellData()->GetNormals()->GetTuple3(cell);
    for( j=0; j<3; j++ ){
      side += (closestPoint[j]-x[j])*NN[j];
    }
  }

  if( side < 0 ){
    return( -1 );
  } else {
    return( 1  );
  }
   NODES->Delete();
}

vtkPolyData *CleanAndNormals( vtkPolyData *MESH ){
  vtkPolyData         *SURF;
  vtkTriangleFilter   *TRIAN;
  vtkCleanPolyData    *CLEAN;
  vtkPolyDataNormals  *NORMALS;

  TRIAN= vtkTriangleFilter::New();
    TRIAN->SetInput( MESH );

  CLEAN = vtkCleanPolyData::New();
    CLEAN->SetInput( TRIAN->GetOutput() );
    CLEAN->PointMergingOn();

  NORMALS= vtkPolyDataNormals::New();
    NORMALS->SetInput( CLEAN->GetOutput() );
    NORMALS->ComputeCellNormalsOn();
    NORMALS->ComputePointNormalsOn();
    NORMALS->ConsistencyOn();
    NORMALS->AutoOrientNormalsOn();
    NORMALS->SplittingOff();
    //NORMALS->FlipNormalsOn();
    NORMALS->Update();

  SURF= vtkPolyData::New();
    SURF->DeepCopy(NORMALS->GetOutput());

  TRIAN->Delete();
  CLEAN->Delete();
  NORMALS->Delete();
    
  return( SURF );
}


// void vtkAreInside( vtkPolyData *MESH, vtkDataSet *POINTS , double narrow) {
//   long                Nop, i;
//   int                 cell, sub;
//   double              closestPoint[3], x[3], distance, bb[6];
// 
//   vtkDataArray        *SIDE;
//   vtkPolyData         *SURF;
//   vtkCellLocator      *LOC;
//   vtkGenericCell      *Gcell;
//   
//   SURF = CleanAndNormals( MESH );
//     SURF->GetBounds(bb);
//   
//   LOC= vtkCellLocator::New();
//     LOC->SetDataSet( SURF );
//  	  LOC->CacheCellBoundsOn();
//     LOC->SetNumberOfCellsPerBucket( 2 );
//     LOC->BuildLocator();
//   Gcell = vtkGenericCell::New();
// 
//   Nop = POINTS->GetNumberOfPoints();
//   SIDE = vtkIntArray::New();
//     SIDE->SetNumberOfTuples( Nop );
// 
//   for( i=0; i<Nop ; i++ ){
//     POINTS->GetPoint(i,x);
//     if( x[0] < bb[0] || x[0] > bb[1] ||
//         x[1] < bb[2] || x[1] > bb[3] ||
//         x[2] < bb[4] || x[2] > bb[5] ) {
//       SIDE->SetTuple1(i,-1);
//       continue;
//     }
//     //LOC->FindClosestPoint( x , closestPoint , Gcell ,cell, sub, distance);
//     if( ! LOC->FindClosestPointWithinRadius( x , narrow , closestPoint , Gcell , cell, sub, distance )	){
//       SIDE->SetTuple1(i,-2);
//       continue;
//     }
//     if( distance < 1e-10 ) {
//       SIDE->SetTuple1(i,0);
//       continue;
//     }
//     
//     SIDE->SetTuple1(i, IsInside( SURF , cell , x , closestPoint ) );
//   }
// 
//   for( i=0; i<Nop ; i++ ){
//     while( SIDE->GetTuple1(i)   !=  1 && i<Nop   ){
//       if( SIDE->GetTuple1(i) == -2 ){
//         SIDE->SetTuple1(i,-1);
//       }
//       i++;
//     }
//     while( SIDE->GetTuple1(i+1) != -1  && i<Nop-1 ){
//       if( SIDE->GetTuple1(i) == -2 ){
//         SIDE->SetTuple1(i,1);
//       }
//       i++;
//     }
//     SIDE->SetTuple1(i,-1);
//   }
// 
//   POINTS->GetPointData()->SetScalars(SIDE);
// 
//   SIDE->Delete();
//   LOC->Delete();
//   Gcell->Delete();
// }


void vtkAreInside( vtkPolyData *MESH, vtkDataSet *POINTS ) {
  long                Nop, i;
  int                 cell, sub;
  double              closestPoint[3], x[3], distance, bb[6];

  vtkDataArray        *SIDE;
  vtkPolyData         *SURF;
  vtkCellLocator      *LOC;
  vtkGenericCell      *Gcell;
  
  SURF = CleanAndNormals( MESH );
    SURF->GetBounds(bb);
  
  LOC= vtkCellLocator::New();
    LOC->SetDataSet( SURF );
 	  LOC->CacheCellBoundsOn();
    LOC->SetNumberOfCellsPerBucket( 2 );
    LOC->BuildLocator();
  Gcell = vtkGenericCell::New();

  Nop = POINTS->GetNumberOfPoints();
  SIDE = vtkIntArray::New();
    SIDE->SetNumberOfTuples( Nop );

  for( i=0; i<Nop ; i++ ){
    POINTS->GetPoint(i,x);
    if( x[0] < bb[0] || x[0] > bb[1] ||
        x[1] < bb[2] || x[1] > bb[3] ||
        x[2] < bb[4] || x[2] > bb[5] ) {
      SIDE->SetTuple1(i,-1);
      continue;
    }
    LOC->FindClosestPoint( x , closestPoint , Gcell ,cell, sub, distance);
    if( distance < 1e-10 ) {
      SIDE->SetTuple1(i,0);
      continue;
    }
    
    SIDE->SetTuple1(i, IsInside( SURF , cell , x , closestPoint ) );
  }
  POINTS->GetPointData()->SetScalars(SIDE);

  SIDE->Delete();
  LOC->Delete();
  Gcell->Delete();
}

// int vtkAreInside( vtkPolyData *MESH, double *x ) {
//   int                 cell, sub;
//   double              closestPoint[3], distance, bb[6];
// 
//   vtkPolyData         *SURF;
//   vtkCellLocator      *LOC;
//   
//   SURF = CleanAndNormals( MESH );
//     SURF->GetBounds(bb);
//   
//   LOC= vtkCellLocator::New();
//     LOC->SetDataSet( SURF );
//  	  LOC->CacheCellBoundsOn();
//     LOC->SetNumberOfCellsPerBucket( 2 );
//     LOC->BuildLocator();
// 
//   if( x[0] < bb[0] || x[0] > bb[1] ||
//       x[1] < bb[2] || x[1] > bb[3] ||
//       x[2] < bb[4] || x[2] > bb[5] ) {
//     LOC->Delete();
//     return(-1);
//   }
//   LOC->FindClosestPoint( x , closestPoint , cell, sub, distance);
//   LOC->Delete();
// 
//   if( distance < 1e-10 ) {
//     return(0);
//   }
// 
//   return( IsInside( SURF , cell , x , closestPoint ) );
// }




