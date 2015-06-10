/*#ifndef real
#define real double
#endif


#define real_old real
#undef real*/

#ifdef _DEBUG_
#include "mmz.h"
#endif

#include <carve/input.hpp>
#include <carve/polyhedron_decl.hpp>
#include <carve/interpolator.hpp>
/*
#define real real_old
#undef real_old
*/

typedef carve::poly::Polyhedron poly_t;

carve::poly::Polyhedron *MESH2carvePolyhedron( const mxArray *m )
{
	//real mm;
	carve::input::PolyhedronData data;

	int        n_xyz, n_tri;
	int        n, t;
	double      *xyz, *tri;

	//filling the points

	if( mxGetField( m,0,"xyz" ) != NULL )
	{
		n_xyz    =  mxGetM( mxGetField( m, 0, "xyz") );
		xyz      =  mxGetPr( mxGetField( m, 0, "xyz") );
	} else
	{
		n_xyz = 0;
	}
	if((mxGetN( mxGetField( m, 0, "xyz")) != 3) && (mxGetN( mxGetField( m, 0, "xyz")) != 0))
	{
		mexErrMsgTxt("This point has not 3 dimensions and not empty.!! \n");
	}
	for ( n=0 ; n<n_xyz ; n++ )
	{
		data.addVertex(carve::geom::VECTOR(xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz]));
	}

	//filling the cells
	if( mxGetField( m,0,"tri" ) != NULL ){
		n_tri    =  mxGetM( mxGetField( m, 0, "tri") );
		tri      = mxGetPr( mxGetField( m, 0, "tri") );
	} else {
		n_tri= 0;
	}
	if( (mxGetN( mxGetField( m, 0, "tri")) != 3) && (mxGetN( mxGetField( m, 0, "tri")) != 0) )
	{
		mexErrMsgTxt("This face has not 3 vertex and not empty.!! \n");
	}
	for (t=0 ; t<n_tri; t++)
	{
		//carve recorre los indices de los puntos en sentido horario
		data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ), ( (int)tri[ t+2*n_tri ] - 1 ));
		//data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+2*n_tri ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));

		/*if( tri[t+2*n_tri] > 0 )
		{
			//carve recorre los indices de los puntos en sentido horario
			data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+2*n_tri ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));
		} else
		{
			if( tri[t+n_tri] > 0 )
			{
				//¿Esta bien?, es un vector pero no se si esta bien transformado a formato "carve"
				data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));
			} else
			{
				if( tri[t+n_tri] > 0 )
				{
					data.addFace(( (int)tri[ t ] - 1 ));
				}
			}
		}*/
	}
	//filling the attributes... ¿SE PUEDE?
	//return new carve::poly::Polyhedron(data.points, data.getFaceCount(), data.faceIndices);
#ifdef _DEBUG_
	mexPrintf("Numero de cara anyadidas = %i\n",data.getFaceCount());
	mexPrintf("Numero de n_xyz = %i\n",n_xyz);
	mexPrintf("Numero de n_xyz desde matlab = %i\n",mxGetM( mxGetField( m, 0, "xyz") ) );
	mexPrintf("Numero de n_tri = %i\n",n_tri);mmzFlush();
#endif
	return data.create();
}

carve::poly::Polyhedron *MESH2carvePolyhedron( const mxArray *m , int iIDFace)
{
	//real mm;
	carve::input::PolyhedronData data;

	long        n_xyz, n_tri;
	long        n, t;
	double      *xyz, *tri;

	//filling the points
	if( mxGetField( m,0,"xyz" ) != NULL )
	{
		n_xyz    =  mxGetM( mxGetField( m, 0, "xyz") );
		xyz      =  mxGetPr( mxGetField( m, 0, "xyz") );
	} else
	{
		n_xyz = 0;
	}
	if(mxGetN( mxGetField( m, 0, "xyz")) != 3)
	{
		mexErrMsgTxt("This point has not 3 dimensions.!! \n");
	}
	std::vector<poly_t::vertex_t> v;
	for ( n=0 ; n<n_xyz ; n++ )
	{
		v.push_back(
				poly_t::vertex_t(
						carve::geom::VECTOR(xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz])
								)
					);
		//data.addVertex(carve::geom::VECTOR(xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz]));
	}

	//filling the cells
	if( mxGetField( m,0,"tri" ) != NULL ){
		n_tri    =  mxGetM( mxGetField( m, 0, "tri") );
		tri      = mxGetPr( mxGetField( m, 0, "tri") );
	} else {
		n_tri= 0;
	}
	if(mxGetN( mxGetField( m, 0, "tri")) != 3)
	{
		mexErrMsgTxt("This face has not 3 vertex.!! \n");
	}


	 std::vector<poly_t::face_t> faces;
	 faces.reserve(n_tri);
	 carve::interpolate::FaceAttr<int> f_tex_num;

	for (t=0 ; t<n_tri; t++)
	{
		faces.push_back(poly_t::face_t( ( &v[tri[ t ] - 1] ), ( &v[tri[ t+n_tri ] - 1] ), ( &v[tri[ t+2*n_tri ] - 1] ) ) );

		f_tex_num.setAttribute(&faces[t], t);
		//data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+2*n_tri ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));

		/*if( tri[t+2*n_tri] > 0 )
		{
			//carve recorre los indices de los puntos en sentido horario
			data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+2*n_tri ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));
		} else
		{
			if( tri[t+n_tri] > 0 )
			{
				//¿Esta bien?, es un vector pero no se si esta bien transformado a formato "carve"
				data.addFace(( (int)tri[ t ] - 1 ), ( (int)tri[ t+n_tri ] - 1 ));
			} else
			{
				if( tri[t+n_tri] > 0 )
				{
					data.addFace(( (int)tri[ t ] - 1 ));
				}
			}
		}*/
	}
	//filling the attributes... ¿SE PUEDE?
	//return new carve::poly::Polyhedron(data.points, data.getFaceCount(), data.faceIndices);
	//return data.create();
	poly_t* pPoly = new poly_t(faces);


	poly_t::face_t &f = pPoly->faces[38];
	int tttt;
	tttt = f_tex_num.getAttribute(&f, 0);
	mexPrintf("Cara 38 %i\n", tttt);


	if(pPoly->faces[38].owner == pPoly)
		mexPrintf("Owner de la Cara 38 es pPoly\n");

	return pPoly;
}

#if !defined( carveOBJ_TYPE )
  #define   carveOBJ_TYPE       carve::csg::CSG::OP
#endif

#define Call_UNION(u1, u2, u3)    if( (!strcmp(met,u1)) || (!strcmp(met,u2)) || (!strcmp(met,u3)) ){ (*O) = carve::csg::CSG::UNION; return; }
#define Call_INTERSECTION(u1, u2, u3, u4)    if( (!strcmp(met,u1)) || (!strcmp(met,u2)) || (!strcmp(met,u3)) || (!strcmp(met,u4)) ){ (*O) = carve::csg::CSG::INTERSECTION;return; }
#define Call_SYMM_DIFF(u1, u2, u3)    if( (!strcmp(met,u1)) || (!strcmp(met,u2)) || (!strcmp(met,u3)) ){ (*O) = carve::csg::CSG::SYMMETRIC_DIFFERENCE; return; }
#define Call_A_MINUS_B(u1, u2, u3)    if( (!strcmp(met,u1)) || (!strcmp(met,u2)) || (!strcmp(met,u3)) ){ (*O) = carve::csg::CSG::A_MINUS_B; return; }
#define Call_B_MINUS_A(u1)    if( !strcmp(met,u1) ){ (*O) = carve::csg::CSG::B_MINUS_A; return; }
#define Call_ALL(u1)    if( !strcmp(met,u1) ){ (*O) = carve::csg::CSG::ALL; return; }
#define Call_SPLIT(u1)    if( !strcmp(met,u1) ){ (*O) = carve::csg::CSG::SPLIT; return; }
#define Call_SPLIT_FULL(u1)    if( !strcmp(met,u1) ){ (*O) = carve::csg::CSG::SPLIT_FULL; return; }

void CallMethod( carveOBJ_TYPE * , char * ) throw (int);


/*carve::poly::Polyhedron *MESH2carvePolyhedron( const mxArray *m )
{
	std::vector<carve::poly::Vertex<3>  > verts;
	std::vector<carve::poly::Face<3> *> faces;

	long        n_xyz, n_tri;
	long        n, t;
	double      *xyz, *tri;

	//filling the points
	if( mxGetField( m,0,"xyz" ) != NULL )
	{
		n_xyz    =  mxGetM( mxGetField( m, 0, "xyz") );
		xyz      =  mxGetPr( mxGetField( m, 0, "xyz") );
	} else
	{
		n_xyz = 0;
	}
	verts.reserve(n_xyz);

	for ( n=0 ; n<n_xyz ; n++ )
	{
		verts.push_back(carve::poly::Vertex<3>(carve::geom::VECTOR(xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz])));

	}

	//filling the cells
	if( mxGetField( m,0,"tri" ) != NULL ){
		n_tri    =  mxGetM( mxGetField( m, 0, "tri") );
		tri      = mxGetPr( mxGetField( m, 0, "tri") );
	} else {
		n_tri= 0;
	}
	faces.reserve(n_tri);
	for (t=0 ; t<n_tri; t++)
	{
		//carve recorre los indices de los puntos en sentido horario
		faces.push_back( new carve::poly::Face<3>(
				&verts[( (int)tri[ t ] - 1 )		],
				&verts[( (int)tri[ t+2*n_tri ] - 1 )],
				&verts[( (int)tri[ t+n_tri ] - 1 )	]
				) );
	}
	//filling the attributes... ¿SE PUEDE?
	return new carve::poly::Polyhedron(faces, verts);
}*/
/*carve::poly::Polyhedron *MESH2carvePolyhedron( const mxArray *m )
{
  long        n_xyz, n_tri;
  double      *xyz, *tri;
  long        n , t;
  int         n_fields, f, n_cols, c, n_rows;
  double      *DATA;
  const char  *field_name;

  vtkPolyData     *POLY   = vtkPolyData::New();
  vtkPoints       *VERTS  = vtkPoints::New();
  vtkCellArray    *FACES  = vtkCellArray::New();
  
  //filling the points
  if( mxGetField( m,0,"xyz" ) != NULL ){
    n_xyz    =  mxGetM( mxGetField( m, 0, "xyz") );
    xyz      = mxGetPr( mxGetField( m, 0, "xyz") );
  } else {
    n_xyz = 0;
  }
  for ( n=0 ; n<n_xyz ; n++ ) {
    VERTS->InsertPoint( n, xyz[n] , xyz[n+n_xyz] , xyz[n+2*n_xyz] );
  }
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
}*/

/*
#if !defined( real )
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
void CallMethod( vtkOBJ_TYPE * , char * , real  * );*/
