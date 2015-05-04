#define NX  150
#define NY  150
#define NZ  150

#include "mex.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"
#include "vtkAreInside.cpp"

#include "vtkStructuredPoints.h"
#include "vtkContourFilter.h"

#define MAX(a,b) a>b?a:b
#define MIN(a,b) a<b?a:b

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  vtkPolyData     *M1;
  vtkPolyData     *M2;
  vtkStructuredPoints *I1;
  vtkStructuredPoints *I2;
  vtkContourFilter  *C;
  double          bb1[6], bb2[6], bb[6];
  double          dx, dy, dz;
  int             i1, i2;
  double          spacing[3], narrow;
  long            i, Nop;
  
  M1= MESH2vtkPolyData( prhs[0] );
  M2= MESH2vtkPolyData( prhs[1] );

  M1->GetBounds( bb1 );
  M2->GetBounds( bb2 );
  
  bb[0]= MIN( bb1[0],bb2[0] );
  bb[1]= MAX( bb1[1],bb2[1] );
  bb[2]= MIN( bb1[2],bb2[2] );
  bb[3]= MAX( bb1[3],bb2[3] );
  bb[4]= MIN( bb1[4],bb2[4] );
  bb[5]= MAX( bb1[5],bb2[5] );
  
 
  dx= bb[1]-bb[0];
  dy= bb[3]-bb[2];
  dz= bb[5]-bb[4];

  bb[0] -= 0.05*dx;
  bb[1] += 0.05*dx;
  bb[2] -= 0.05*dy;
  bb[3] += 0.05*dy;
  bb[4] -= 0.05*dz;
  bb[5] += 0.05*dz;

  I1 = vtkStructuredPoints::New();
    I1->SetOrigin(bb[0],bb[2],bb[4]);
    I1->SetSpacing( (bb[1]-bb[0])/(NX-1) ,
                    (bb[3]-bb[2])/(NY-1) ,
                    (bb[5]-bb[4])/(NZ-1) );
    I1->SetDimensions(NX,NY,NZ);

//     I1->SetSpacing( 1 , 1 , 1 );
//     I1->SetDimensions(NX,NY,NZ);
  
  I2 = vtkStructuredPoints::New();
    I2->DeepCopy(I1);

  I1->GetSpacing(spacing);
  if( spacing[0]>spacing[1] ){
    if(spacing[0]>spacing[2] ){
      narrow= spacing[0];
    } else {
      narrow= spacing[2];
    }
  } else {
    if(spacing[1]>spacing[2] ){
      narrow= spacing[1];
    } else {
      narrow= spacing[2];
    }
  }    
  narrow= narrow*3.0;
    
  vtkAreInside( M1 , I1 , narrow );  M1->Delete();
  vtkAreInside( M2 , I2 , narrow );  M2->Delete();

//   vtkAreInside( M1 , I1 );  M1->Delete();
//   vtkAreInside( M2 , I2 );  M2->Delete();
    
  Nop= I1->GetNumberOfPoints();
  
  for( i=0 ; i<Nop ; i++ ) {
    i1= I1->GetPointData()->GetScalars()->GetTuple1(i);
    if( i1 >= 0 ) {
      i2= I2->GetPointData()->GetScalars()->GetTuple1(i);
      if( i2 < 0) {
      	I1->GetPointData()->GetScalars()->SetTuple1(i,-1);
      }
    }
  }
  
  C= vtkContourFilter::New();
    C->SetInput(I1);
    C->SetValue(0,-0.9);
    C->ComputeNormalsOff(); 
    C->ComputeGradientsOff();
    C->ComputeScalarsOff();
    C->Update();

  plhs[0]= vtkPolyData2MESH( C->GetOutput() );

  C->Delete();
//  I1->Delete();
//  I2->Delete();
}



