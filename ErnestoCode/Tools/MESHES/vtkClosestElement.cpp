/*

[ element_id , xyz_closest_point , distance , barycentric_coordinates ] = vtkClosestElement( MESH , [ point1 ; point2 ; ... ] );


 
 vtkClosestElement( MESH )                                  crea el locator
 [outs] = vtkClosestElement( [ point1 ; point2 ; ...] )     sobre el ultimo locator creado, calcula los puntos
 vtkClosestElement( [] , [] )                               libera el locator
 
 
 */

#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "MESH2vtkPolyData.h"

#include "vtkPolyData.h"
#include "vtkCellLocator.h"
#include "vtkGenericCell.h"

struct LOCATOR {
    vtkCellLocator *LOC;
    vtkPolyData    *MESH;
    vtkGenericCell *CELL;
};
typedef struct LOCATOR LOCATOR;

static LOCATOR L;


void clean(){

  try{ 
    if( L.CELL != NULL ){
      L.CELL->Delete(); 
    }
    L.CELL = NULL;
/* //     mexPrintf("CELL  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING CELL\n");
  }

  try{ 
    if( L.LOC != NULL ){
      L.LOC->Delete(); 
    }
    L.LOC = NULL;
/* //     mexPrintf("LOC  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING LOC\n");
  }

  try{ 
    if( L.MESH != NULL ){
      L.MESH->Delete(); 
    }
    L.MESH = NULL;
/* //     mexPrintf("MESH  deleted\n"); */
  } catch( char * str ) {
    mexPrintf("error DELETING MESH\n");
  }

}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  double    *points;
  long      n_points, p;
  double    point[3];
  double    closest[3];
  int       subId, t;
  double    dist;
  double    *O_id, *O_x, *O_y, *O_z, *O_distance, *O_barycentric_coordinates;
  double    xyz[3], M[9], DET;
  
  
  if( nrhs == 1  &&  mxIsStruct( prhs[0] ) ){

   /*  //     vtkClosestElement( MESH ) */
    
    clean();

    L.LOC  = vtkCellLocator::New();  
    L.CELL = vtkGenericCell::New();  
    L.MESH = MESH2vtkPolyData( prhs[0] );

    L.LOC->SetDataSet( L.MESH );
    L.LOC->SetNumberOfCellsPerBucket( 5 );
    L.LOC->BuildLocator();
    
/* //     mexPrintf("\n\nCREADO\n\n"); */
    
  } else if( nrhs == 1  &&  mxIsNumeric(prhs[0]) ){
   /*  //     vtkClosestElement( [ point1 ; point2 ] ) */
    

    if( L.LOC == NULL  ||  L.CELL == NULL  || L.MESH == NULL ){
      mexPrintf("No hay creado un LOCATOR.");
    }

    point[0] = 0; point[1] = 0; point[2] = 0;
    try {
      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);
    } catch( char * str ) {
      mexPrintf("el LOCATOR no responde\n");
      clean();
      mexPrintf("el LOCATOR no responde.");
    }    
    
    n_points =  mxGetM( prhs[0] );
    points   = mxGetPr( prhs[0] );

   /*  //create the outputs */
    plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
    O_id   = mxGetPr( plhs[0] );
    if( nlhs > 1) {
      plhs[1] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_x   = mxGetPr( plhs[1] );
      O_y   = O_x + n_points;
      O_z   = O_y + n_points;
    }
    if( nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix( n_points,1,mxREAL );
      O_distance = mxGetPr( plhs[2] );
    }
    if( nlhs > 3) {
      plhs[3] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_barycentric_coordinates = mxGetPr( plhs[3] );
    }

    for( p=0 ; p<n_points ; p++ ) {
      point[0]= points[p];
      point[1]= points[p+n_points];
      point[2]= points[p+2*n_points];

      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);

      O_id[ p ]= t+1;         /*   //update the output id */

      if( nlhs > 1) {                   /*     //update the output point */
        O_x[ p ] = closest[0];
        O_y[ p ] = closest[1];
        O_z[ p ] = closest[2];
      }
      if( nlhs > 2) {                    /*    //update the output distance */
        O_distance[p] = sqrt( dist );
      } 
      if( nlhs > 3) {                     /*   //barycentric coordinates */
        #define M(i,j)  M[ i-1 + (j-1)*3 ]
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(0), xyz );
        M(1,1) = xyz[0];
        M(2,1) = xyz[1];
        M(3,1) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(1), xyz );
        M(1,2) = xyz[0];
        M(2,2) = xyz[1];
        M(3,2) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(2), xyz );
        M(1,3) = xyz[0];
        M(2,3) = xyz[1];
        M(3,3) = xyz[2];
        
        DET = 1.0/( M(3,1)*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + M(2,1)*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + M(1,1)*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) );
        
        O_barycentric_coordinates[ p              ] = ( closest[2]*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + closest[1]*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + closest[0]*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) )*DET;
        O_barycentric_coordinates[ p +   n_points ] = ( closest[2]*( M(1,1)*M(2,3) - M(1,3)*M(2,1) ) + closest[1]*( M(1,3)*M(3,1) - M(1,1)*M(3,3) ) + closest[0]*( M(2,1)*M(3,3) - M(2,3)*M(3,1) ) )*DET;
        O_barycentric_coordinates[ p + 2*n_points ] = ( closest[2]*( M(1,2)*M(2,1) - M(1,1)*M(2,2) ) + closest[1]*( M(1,1)*M(3,2) - M(1,2)*M(3,1) ) + closest[0]*( M(2,2)*M(3,1) - M(2,1)*M(3,2) ) )*DET;
      }

    }
    
  } else if( nrhs == 2  &&  mxIsStruct(prhs[0])  &&  mxIsNumeric(prhs[1]) ){
  /*   //     vtkClosestElement( MESH , [ point1 ; point2 ] ) */
    
    clean();    
    
    L.LOC  = vtkCellLocator::New();
    L.CELL = vtkGenericCell::New();
    L.MESH = MESH2vtkPolyData( prhs[0] );

    L.LOC->SetDataSet( L.MESH );
   /*  //     L.LOC->SetNumberOfCellsPerBucket( 5 ); */
    L.LOC->BuildLocator();

    
    n_points =  mxGetM( prhs[1] );
    points   = mxGetPr( prhs[1] );

   /*  //create the outputs */
    plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
    O_id   = mxGetPr( plhs[0] );
    if( nlhs > 1) {
      plhs[1] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_x   = mxGetPr( plhs[1] );
      O_y   = O_x + n_points;
      O_z   = O_y + n_points;
    }
    if( nlhs > 2) {
      plhs[2] = mxCreateDoubleMatrix( n_points,1,mxREAL );
      O_distance = mxGetPr( plhs[2] );
    }
    if( nlhs > 3) {
      plhs[3] = mxCreateDoubleMatrix( n_points,3,mxREAL );
      O_barycentric_coordinates = mxGetPr( plhs[3] );
    }

    
    for( p=0 ; p<n_points ; p++ ) {
      point[0]= points[p];
      point[1]= points[p+n_points];
      point[2]= points[p+2*n_points];

      L.LOC->FindClosestPoint( point , closest , L.CELL , ((vtkIdType&)t) , subId , dist);

      O_id[ p ]= t+1;         /*   //update the output id */

      if( nlhs > 1) {                /*       //update the output point */
        O_x[ p ] = closest[0];
        O_y[ p ] = closest[1];
        O_z[ p ] = closest[2];
      }
      if( nlhs > 2) {                      /*  //update the output distance */
        O_distance[p] = sqrt( dist );
      }
      
      if( nlhs > 3) {                     /*   //barycentric coordinates */
        #define M(i,j)  M[ i-1 + (j-1)*3 ]
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(0), xyz );
        M(1,1) = xyz[0];
        M(2,1) = xyz[1];
        M(3,1) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(1), xyz );
        M(1,2) = xyz[0];
        M(2,2) = xyz[1];
        M(3,2) = xyz[2];
        
        L.MESH->GetPoint( L.MESH->GetCell(t)->GetPointId(2), xyz );
        M(1,3) = xyz[0];
        M(2,3) = xyz[1];
        M(3,3) = xyz[2];
        
        DET = 1.0/( M(3,1)*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + M(2,1)*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + M(1,1)*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) );
        
        O_barycentric_coordinates[ p              ] = ( closest[2]*( M(1,3)*M(2,2) - M(1,2)*M(2,3) ) + closest[1]*( M(1,2)*M(3,3) - M(1,3)*M(3,2) ) + closest[0]*( M(2,3)*M(3,2) - M(2,2)*M(3,3) ) )*DET;
        O_barycentric_coordinates[ p +   n_points ] = ( closest[2]*( M(1,1)*M(2,3) - M(1,3)*M(2,1) ) + closest[1]*( M(1,3)*M(3,1) - M(1,1)*M(3,3) ) + closest[0]*( M(2,1)*M(3,3) - M(2,3)*M(3,1) ) )*DET;
        O_barycentric_coordinates[ p + 2*n_points ] = ( closest[2]*( M(1,2)*M(2,1) - M(1,1)*M(2,2) ) + closest[1]*( M(1,1)*M(3,2) - M(1,2)*M(3,1) ) + closest[0]*( M(2,2)*M(3,1) - M(2,1)*M(3,2) ) )*DET;
      }
      
      
    }

    clean();

  } else if(  nrhs == 2  &&  mxIsNumeric(prhs[0])  &&  mxIsNumeric(prhs[1])  ){
 /*   // vtkClosestElement( [] , [] ) */
    
    clean();
/* //     mexPrintf("LOCATOR borrado !!\n" ); */
    
  } else {
    mexPrintf("No entiendo la llamada.");
  }
  
}
