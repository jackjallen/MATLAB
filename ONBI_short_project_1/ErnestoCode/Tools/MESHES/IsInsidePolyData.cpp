#include "mex.h"
#include "MESH2vtkPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkDataArray.h"
#include "vtkCellLocator.h"
#include "vtkCellData.h"
#include "vtkGenericCell.h"
#include "vtkIdList.h"



// #define mxFlush()           mexEvalString("drawnow expose;")

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  long                i, Nop;
  int                 cell, sub, j;
  double              closestPoint[3], x[3], distance, bb[6], *out, *points;
  double              N[3], *NN, pcoords[3], weigths[3], side;

  vtkPolyData         *SURF;
  SURF= MESH2vtkPolyData( prhs[0] );
//   mexPrintf("hay celdas SURF %d\n\n",SURF->GetNumberOfCells() );

  vtkCleanPolyData    *CLEAN;
  CLEAN = vtkCleanPolyData::New();
    CLEAN->SetInput( SURF );
    CLEAN->PointMergingOn();
    CLEAN->ConvertLinesToPointsOn();
    CLEAN->ConvertPolysToLinesOn();
    CLEAN->ConvertStripsToPolysOn();
//     CLEAN->Update();
//   mexPrintf("hay celdas CLEAN %d\n\n",CLEAN->GetOutput()->GetNumberOfCells() );

  vtkTriangleFilter   *TRIAN;
  TRIAN= vtkTriangleFilter::New();
    TRIAN->SetInput( CLEAN->GetOutput() );
    TRIAN->PassVertsOff();
    TRIAN->PassLinesOn();
//     TRIAN->Update();
//   mexPrintf("hay celdas TRIAN %d\n\n",TRIAN->GetOutput()->GetNumberOfCells() );

  vtkPolyDataNormals  *NORMALS;
  NORMALS= vtkPolyDataNormals::New();
    NORMALS->SetInput( TRIAN->GetOutput() );
    NORMALS->ComputeCellNormalsOn();
    NORMALS->ComputePointNormalsOn();
    NORMALS->ConsistencyOn();
    NORMALS->AutoOrientNormalsOn();
    NORMALS->SplittingOff();
//     NORMALS->NonManifoldTraversalOff();
    //NORMALS->FlipNormalsOn();
    NORMALS->Update();

//   mexPrintf("hay celdas en NormalsData %d\n\n",NORMALS->GetOutput()->GetCellData()->GetNormals()->GetNumberOfTuples() );
//   mexPrintf("hay celdas Normals %d\n\n",NORMALS->GetOutput()->GetNumberOfCells() );

  SURF->DeepCopy(NORMALS->GetOutput());
//   mexPrintf("hay celdas SURF %d\n\n",SURF->GetNumberOfCells() );
//   mexPrintf("hay celdas en SURFData %d\n\n",SURF->GetCellData()->GetNormals()->GetNumberOfTuples() );

  TRIAN->Delete();
  CLEAN->Delete();
  NORMALS->Delete();
  
  vtkCellLocator      *LOC;
  LOC= vtkCellLocator::New();
    LOC->SetDataSet( SURF );
 	  LOC->CacheCellBoundsOn();
    LOC->SetNumberOfCellsPerBucket( 2 );
    LOC->BuildLocator();

  vtkGenericCell      *Gcell;
  Gcell = vtkGenericCell::New();

  
  SURF->GetBounds(bb);
  
  Nop = mxGetM( prhs[1] );
//   Nop = 10;
//   mexPrintf("Nop= %d\n",Nop);

  plhs[0]= mxCreateDoubleMatrix( Nop , 1 , mxREAL );
  out = mxGetPr( plhs[0] );

  vtkIdList           *NODES;
  NODES = vtkIdList::New();
  
  points = mxGetPr( prhs[1] );
  for( i=0; i<Nop ; i++ ){
//     mexPrintf("i= %d",i);
    
    x[0]= *( points + i         );
    x[1]= *( points + i +   Nop );
    x[2]= *( points + i + 2*Nop );
    if( x[0] < bb[0] || x[0] > bb[1] ||
        x[1] < bb[2] || x[1] > bb[3] ||
        x[2] < bb[4] || x[2] > bb[5] ) {
      *(out+i)= -1;
      
//       mexPrintf("  fuera bb\n"); mxFlush();
      continue;
    }
    LOC->FindClosestPoint( x , closestPoint , Gcell ,((vtkIdType&)cell), sub, distance);
    if( distance < 1e-10 ) {
      *(out+i)= 0;
//       mexPrintf("  sobre superficie\n"); mxFlush();
      continue;
    }
    
    SURF->GetCell( cell )->EvaluatePosition( closestPoint, NULL, sub, pcoords, distance, weigths );
    if( weigths[0]<1e-10 || weigths[1]<1e-10 || weigths[2]<1e-10 ){
//       mexPrintf("  en borde "); mxFlush();

      SURF->GetCellPoints( cell , NODES );
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(0) );
      for( j=0 ; j<3; j++ ){  N[j]  = NN[j]*weigths[0];  }
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(1) );
      for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[1];  }
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(2) );
      for( j=0 ; j<3; j++ ){  N[j] += NN[j]*weigths[2];  }

//       mexPrintf("  done "); mxFlush();
    } else {
//       mexPrintf("  en celda %d: ", cell); mxFlush();
      NN = SURF->GetCellData()->GetNormals()->GetTuple3(cell);
      for( j=0 ; j<3; j++ ){  N[j] = NN[j];  }
//       mexPrintf("  done "); mxFlush();
    }

    side= 0;
    for( j=0; j<3; j++ ){
      side += (closestPoint[j]-x[j])*N[j];
    }
//     mexPrintf("  side: %f   ",side); mxFlush();
    
    if( side < 0 ){
      *(out+i) = -1;
    } else {
      *(out+i) = 1;
    }

//     mexPrintf(" OK \n"); mxFlush();
    
  }

  NODES->Delete();
  LOC->Delete();
  Gcell->Delete();
  SURF->Delete();
}

