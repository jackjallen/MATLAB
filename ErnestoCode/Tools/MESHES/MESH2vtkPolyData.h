#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

vtkPolyData* MESH2vtkPolyData( const mxArray *m ){
  long        n_xyz, n_tri;
  double      *xyz, *tri;
  long        n , t;
  int         n_fields, f, n_cols, c, n_rows;
  double      *DATA;
  const char  *field_name;

  vtkPolyData     *POLY   = vtkPolyData::New();
  vtkPoints       *VERTS  = vtkPoints::New();
  vtkCellArray    *FACES  = vtkCellArray::New();
  
  //Queremos que todos los puntos esten en double
  VERTS->SetDataTypeToDouble ();


  //debug
  //mexPrintf("Capacity solo con un New()", VERTS->GetData().Capacity());

  //mexPrintf("Tamaño de memoria consumida solo con un New() %lf.\n", ((float)VERTS->GetActualMemorySize()));

  //VERTS->Allocate(500);
  //mexPrintf("Tamaño de memoria consumida despues de Allocate(500) %lf.\n", ((float)VERTS->GetActualMemorySize()));

  //filling the points
  if( mxGetField( m,0,"xyz" ) != NULL ){
    n_xyz    =  mxGetM( mxGetField( m, 0, "xyz") );
    xyz      = mxGetPr( mxGetField( m, 0, "xyz") );
  } else {
    n_xyz = 0;
  }

  //debug
  VERTS->SetNumberOfPoints(n_xyz);
  //mexPrintf("Capacity despues de un SetNumberOfPoints()", VERTS->GetData().Capacity());
  //mexPrintf("Tamaño de memoria consumida despues de un SetNumberOfPoints() %lf de %i datos.\n", ((float)VERTS->GetActualMemorySize()), n_xyz);


  for ( n=0 ; n<n_xyz ; n++ ) {
    VERTS->InsertPoint( n, xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz] );
    /*if( (n % 100 ) == 0)
    {
    	mexPrintf("Tamaño de memoria consumida despues de insertPoint() %lf de %i datos.\n", ((float)VERTS->GetActualMemorySize()), n);
        mexPrintf("Tamaño: %i datos.\n", n);
    }*/
  }

  //mexPrintf("Tamaño de memoria consumida despues de un SetNumberOfPoints() %lf de %i datos.\n", ((float)VERTS->GetActualMemorySize()), n_xyz);


  POLY->SetPoints( VERTS );

  //filling the cells
  if( mxGetField( m,0,"tri" ) != NULL ){
    n_tri    =  mxGetM( mxGetField( m, 0, "tri") );
    tri      = mxGetPr( mxGetField( m, 0, "tri") );
  } else {
    n_tri= 0;
  }

  for (t=0 ; t<n_tri; t++) {
    if( tri[t+2*n_tri] > 0 ) {
      FACES->InsertNextCell(3);
      FACES->InsertCellPoint( (int)tri[ t ] - 1 );
      FACES->InsertCellPoint( (int)tri[ t+n_tri ] - 1 );
      FACES->InsertCellPoint( (int)tri[ t+2*n_tri ] - 1 );
    } else {
      if( tri[t+n_tri] > 0 ) {
        FACES->InsertNextCell(2);
        FACES->InsertCellPoint( (int)tri[ t ] - 1 );
        FACES->InsertCellPoint( (int)tri[ t+n_tri ] - 1 );
      } else {
        if( tri[t+n_tri] > 0 ) {
          FACES->InsertNextCell(1);
          FACES->InsertCellPoint( (int)tri[ t ] - 1 );
        } 
      }
    }
  }

  POLY->SetPolys(  FACES );

  VERTS->Delete();
  FACES->Delete();


  //filling the attributes
  n_fields= mxGetNumberOfFields( m );
  for( f=0 ; f < n_fields ; f++ ){
    field_name= mxGetFieldNameByNumber( m , f );
    DATA  = mxGetPr( mxGetField( m, 0, field_name ) );
    n_rows=  mxGetM( mxGetField( m, 0, field_name ) );
    n_cols=  mxGetN( mxGetField( m, 0, field_name ) );

    if( !strncmp( field_name,"tri",3 ) && strcmp( field_name, "tri") ){
      vtkDoubleArray  *ARRAY  = vtkDoubleArray::New();
        ARRAY->SetName( field_name+3 ); 
        ARRAY->SetNumberOfComponents( n_cols );
        ARRAY->SetNumberOfTuples( n_rows );

        for( t=0 ; t<n_rows ; t++ ) {
          for( c=0 ; c<n_cols ; c++ ) {
            ARRAY->SetComponent(t,c, DATA[t+c*n_rows] );
          }
        }
        POLY->GetCellData()->AddArray( ARRAY );
      ARRAY->Delete();
      if( !strcmp(field_name,"triNORMALS")) {
        POLY->GetCellData()->SetActiveNormals("Normals");
      }
    } 
    if( !strncmp( field_name,"xyz",3 ) && strcmp( field_name, "xyz") ){
      vtkDoubleArray  *ARRAY  = vtkDoubleArray::New();
        ARRAY->SetName( field_name+3 ); 
        ARRAY->SetNumberOfComponents( n_cols );
        ARRAY->SetNumberOfTuples( n_rows );

        for( t=0 ; t<n_rows ; t++ ) {
          for( c=0 ; c<n_cols ; c++ ) {
            ARRAY->SetComponent(t,c, DATA[t+c*n_rows] );
          }
        }
        POLY->GetPointData()->AddArray( ARRAY );
      ARRAY->Delete();
      if( !strcmp(field_name,"xyzNORMALS")) {
        POLY->GetPointData()->SetActiveNormals("Normals");
      }
    }
    if( !strcmp( field_name,"uv" )  ){
      vtkDoubleArray  *ARRAY  = vtkDoubleArray::New();
        ARRAY->SetName( field_name );
        ARRAY->SetNumberOfComponents( n_cols );
        ARRAY->SetNumberOfTuples( n_rows );

        for( t=0 ; t<n_rows ; t++ ) {
          for( c=0 ; c<n_cols ; c++ ) {
            ARRAY->SetComponent(t,c, DATA[t+c*n_rows] );
          }
        }
        POLY->GetPointData()->SetTCoords( ARRAY );
      ARRAY->Delete();
      POLY->GetPointData()->SetActiveTCoords("UV");
    }
  }
  return POLY;
}


#if !defined( vtkOBJ_TYPE )
  #define   vtkOBJ_TYPE       vtkPolyData
#endif

#if !defined( real )
  #define   real       double
#endif


#define Call_0(m)    if( !strcmp(met,#m) ){ O->m( ); return; }
#define Call_1(m,v)  if( !strcmp(met,#m) ){ O->m(v); return; }
void CallMethod( vtkOBJ_TYPE * , char * );
void CallMethod( vtkOBJ_TYPE * , char * , char  * );
void CallMethod( vtkOBJ_TYPE * , char * , real    );
void CallMethod( vtkOBJ_TYPE * , char * , real  * );
