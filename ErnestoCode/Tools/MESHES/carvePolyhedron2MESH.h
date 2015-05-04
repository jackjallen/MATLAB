/*#ifndef real
#define real double
#endif

#define real_old real
#undef real*/


#include <stdlib.h>

#include <carve/input.hpp>
#include <carve/face_decl.hpp>
#include <carve/polyhedron_decl.hpp>

//#include "myMEX.h"
#include "mz.h"
#include "mmz.h"

/*
#define real real_old
#undef real_old*/


mxArray *carvePolyhedron2MESH( carve::poly::Polyhedron *POLY )
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *DATA;
	mxArray     *m;
	long        N, n, nVert;
	double      xyz[3];
	int         n_cols, c;
	int         n_fields, f;
	char        name[256];

	m   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	//NODES
	N= POLY->vertices.size();
	if( N> 0)
	{
		DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
		for( n= 0; n<N ; n++ ) {
			for( c=0 ; c<3 ; c++ )
			{
				xyz[c] = POLY->vertices[n].v.v[c];
				*( mxGetPr(DATA) + n + c*N )= xyz[c];
			}
		}
		mxAddField( m,    "xyz" );
		mxSetField( m, 0, "xyz", DATA );

		/*n_fields= POLY->GetPointData()->GetNumberOfArrays();
    for( f = 0 ; f < n_fields ; f++ ){
      sprintf( name , "xyz%s", POLY->GetPointData()->GetArray(f)->GetName() );
      if( !strcmp( name , "xyz(null)" ) ) {
        sprintf( name , "xyz%s", "NULL" );
      }
      n_cols= POLY->GetPointData()->GetArray(f)->GetNumberOfComponents();
      DATA = mxCreateDoubleMatrix( N , n_cols , mxREAL );
      for( n= 0; n<N ; n++ ) {
        for( c=0 ; c< n_cols ; c++ ){
		 *( mxGetPr(DATA) + n + c*N )= POLY->GetPointData()->GetArray(f)->GetComponent(n,c);
        }
      }
      mxAddField( m,    name );
      mxSetField( m, 0, name, DATA );
    }*/
	}

	//TRI
	mxAddField( m,    "tri" );
	N= POLY->faces.size();
	if( N > 0 ) {
		DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
		for( n= 0; n<N ; n++ )
		{
			nVert = POLY->faces[n].nVertices();
			if( nVert > 3)
			{
				mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
				mexErrMsgTxt("\n");
			}
			c = 0;
			for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
			{
				*( mxGetPr(DATA) + n + c*N ) = POLY->vertexToIndex((*j)) + 1;
				c++;
			}
			for( ; c<3 ; c++ ){
				*( mxGetPr(DATA) + n + c*N )= 0;
			}
		}
		mxSetField( m, 0, "tri", DATA );

		/*n_fields= POLY->GetCellData()->GetNumberOfArrays();
    for( f = 0 ; f < n_fields ; f++ ){
      sprintf( name , "tri%s", POLY->GetCellData()->GetArray(f)->GetName() );
      n_cols= POLY->GetCellData()->GetArray(f)->GetNumberOfComponents();
      DATA = mxCreateDoubleMatrix( N , n_cols , mxREAL );
      for( n= 0; n<N ; n++ ) {
        for( c=0 ; c< n_cols ; c++ ){
		 *( mxGetPr(DATA) + n + c*N )= POLY->GetCellData()->GetArray(f)->GetComponent(n,c);
        }
      }
      mxAddField( m,    name );
      mxSetField( m, 0, name, DATA );
    }*/
	}

	return(m);
}

mxArray *carvePolyhedron2MESH( carve::poly::Polyhedron *POLY , int iIDpoly)
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *DATA;
	mxArray     *m;
	long        N, n, nVert;
	double      xyz[3];
	int         n_cols, c;
	int         n_fields, f;
	char        name[256];
	int iNumFaces = 0, iContFaces = 0;

	m   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	//NODES
	N= POLY->vertices.size();
	if( N> 0)
	{
		DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
		for( n= 0; n<N ; n++ ) {
			for( c=0 ; c<3 ; c++ )
			{
				xyz[c] = POLY->vertices[n].v.v[c];
				*( mxGetPr(DATA) + n + c*N )= xyz[c];
			}
		}
		mxAddField( m,    "xyz" );
		mxSetField( m, 0, "xyz", DATA );
	}

	//TRI

	/*Cuenta cuantas caras hay del owner solicitado*/
	N= POLY->faces.size();

	for( n= 0; n<N ; n++ )
	{
		if(POLY->faces[n].iIDpoly == iIDpoly)
			iNumFaces++;
	}

	mxAddField( m,    "tri" );
	if( iNumFaces > 0 )
	{
		DATA = mxCreateDoubleMatrix( iNumFaces , 3 , mxREAL );
		iContFaces = 0;
		for( n= 0; n<N ; n++ )
		{
			if(POLY->faces[n].iIDpoly == iIDpoly)
			{
				nVert = POLY->faces[n].nVertices();
				if( nVert > 3)
				{
					mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
					mexErrMsgTxt("\n");
				}
				c = 0;
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA) + iContFaces + c*iNumFaces ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA) + iContFaces + c*iNumFaces )= 0;
				}
				iContFaces++;
			}
		}
		mxSetField( m, 0, "tri", DATA );
	}

	return(m);
}

mxArray *carvePolyhedron2MESH( carve::poly::Polyhedron *POLY , int iIDpoly, mxArray *copyFrom)
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *DATA;
	mxArray     *m;
	long        N, n, nVert;
	double      xyz[3];
	int         n_cols, c;
	int         n_fields, f;
	char        name[256];
	int iNumFaces = 0, iContFaces = 0;

#ifdef _DEBUG_MEX
	CreateTicToc_old();
	mzCreateTicToc(1);
	TIC_old();
	mzTIC(1);
#endif

	m   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	if( POLY == NULL )
	{
		mxAddField( m , "xyz" );
		mxAddField( m , "tri" );

		return( m );
	}

	//NODES
	if (copyFrom == NULL)
	{
		N= POLY->vertices.size();
		if( N> 0)
		{
			DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
			for( n= 0; n<N ; n++ ) {
				for( c=0 ; c<3 ; c++ )
				{
					xyz[c] = POLY->vertices[n].v.v[c];
					*( mxGetPr(DATA) + n + c*N )= xyz[c];
				}
			}
			mxAddField( m,    "xyz" );
			mxSetField( m, 0, "xyz", DATA );
		}
	}
	else
	{
		funcCreate pFuncCreate;
		pFuncCreate = GetCreatePointer();

		mxAddField( m,    "xyz" );
		DATA = mxGetField( copyFrom, 0, "xyz");
		mxSetField( m , 0 , "xyz" , (  pFuncCreate((void *)DATA) ) );
		    //make the SharedCopy
	}

	//TRI

	/*Cuenta cuantas caras hay del owner solicitado*/
	N= POLY->faces.size();

	for( n= 0; n<N ; n++ )
	{
		if(POLY->faces[n].iIDpoly == iIDpoly)
			iNumFaces++;
	}

	mxAddField( m,    "tri" );
	if( iNumFaces > 0 )
	{
		DATA = mxCreateDoubleMatrix( iNumFaces , 3 , mxREAL );
		iContFaces = 0;
		for( n= 0; n<N ; n++ )
		{
			if(POLY->faces[n].iIDpoly == iIDpoly)
			{
				nVert = POLY->faces[n].nVertices();
				if( nVert > 3)
				{
					mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
					mexErrMsgTxt("\n");
				}
				c = 0;
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA) + iContFaces + c*iNumFaces ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA) + iContFaces + c*iNumFaces )= 0;
				}
				iContFaces++;
			}
		}
		mxSetField( m, 0, "tri", DATA );
	}
#ifdef _DEBUG_MEX
	if(copyFrom != NULL)
	{
		mexPrintf("xyz Original %p, tri Original %p\n", mxGetPr(mxGetField( copyFrom, 0, "xyz")), mxGetPr(mxGetField( copyFrom, 0, "tri")));
		mexPrintf("xyz Copia %p, tri Copia %p\n", mxGetPr(mxGetField( m, 0, "xyz")), mxGetPr(mxGetField( m, 0, "tri")));
	}

	TOC_old("TIC TOC, OK: ");
	mzTOC(1, "carvePolyhedron2MESH. mi TIC TOC.");
#endif

	return(m);
}

int carvePolyhedron2MESH( carve::poly::Polyhedron *POLY , mxArray * &pOUT1, mxArray * &pOUT2)
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *DATA1, *DATA2;
	mxArray     *m1, *m2;
	long        N, n, nVert;
	double      xyz[3];
	int         n_cols, c;
	int         n_fields, f;
	char        name[256];
	int iNumFaces1 = 0, iContFaces1 = 0, iNumFaces2 = 0, iContFaces2 = 0;
	int iret = -1;

	m1   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);
	m2   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	//NODES
	N= POLY->vertices.size();
	if( N> 0)
	{
		DATA1 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		DATA2 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		for( n= 0; n<N ; n++ ) {
			for( c=0 ; c<3 ; c++ )
			{
				xyz[c] = POLY->vertices[n].v.v[c];
				*( mxGetPr(DATA1) + n + c*N )= xyz[c];
				*( mxGetPr(DATA2) + n + c*N )= xyz[c];
			}
		}
		mxAddField( m1,    "xyz" );
		mxSetField( m1, 0, "xyz", DATA1 );
		mxAddField( m2,    "xyz" );
		mxSetField( m2, 0, "xyz", DATA2 );
	}

	//TRI

	/*Cuenta cuantas caras hay del owner solicitado*/
	N= POLY->faces.size();

	for( n= 0; n<N ; n++ )
	{
		if(POLY->faces[n].iIDpoly == 0)
			iNumFaces1++;
		if(POLY->faces[n].iIDpoly == 1)
			iNumFaces2++;
	}

	if( ( iNumFaces1 > 0 ) || ( iNumFaces2 > 0) )
	{
		if( iNumFaces1 > 0 )
		{
			DATA1 = mxCreateDoubleMatrix( iNumFaces1 , 3 , mxREAL );
			iContFaces1 = 0;
		}
		if( iNumFaces2 > 0 )
		{
			DATA2 = mxCreateDoubleMatrix( iNumFaces2 , 3 , mxREAL );
			iContFaces2 = 0;
		}
		for( n= 0; n<N ; n++ )
		{
			nVert = POLY->faces[n].nVertices();
			if( nVert > 3)
			{
				mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
				mexErrMsgTxt("\n");
			}
			c = 0;
			if(POLY->faces[n].iIDpoly == 0)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA1) + iContFaces1 + c*iNumFaces1 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA1) + iContFaces1 + c*iNumFaces1 )= 0;
				}
				iContFaces1++;
			}
			if(POLY->faces[n].iIDpoly == 1)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA2) + iContFaces2 + c*iNumFaces2 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA2) + iContFaces2 + c*iNumFaces2 )= 0;
				}
				iContFaces2++;
			}
		}
		mxAddField( m1,    "tri" );
		mxSetField( m1, 0, "tri", DATA1 );
		mxAddField( m2,    "tri" );
		mxSetField( m2, 0, "tri", DATA2 );

		pOUT1 = m1;
		pOUT2 = m2;
	}
	iret = 0;
	return(iret);
}
int carvePolyhedron2MESH( carve::poly::Polyhedron *POLY , mxArray * &pOUT1, mxArray * &pOUT2, mxArray * &pOUT3, mxArray * &pOUT4)
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *DATA1, *DATA2, *DATA3, *DATA4;
	mxArray     *m1, *m2, *m3, *m4;
	long        N, n, nVert;
	double      xyz[3];
	int         n_cols, c;
	int         n_fields, f;
	char        name[256];
	int iNumFaces1 = 0, iContFaces1 = 0, iNumFaces2 = 0, iContFaces2 = 0, iNumFaces3 = 0, iContFaces3 = 0, iNumFaces4 = 0, iContFaces4 = 0;
	int iret = -1;

	m1   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);
	m2   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);
	m3   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);
	m4   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	//NODES
	N= POLY->vertices.size();
	if( N> 0)
	{
		DATA1 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		DATA2 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		DATA3 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		DATA4 = mxCreateDoubleMatrix( N , 3 , mxREAL );
		for( n= 0; n<N ; n++ ) {
			for( c=0 ; c<3 ; c++ )
			{
				xyz[c] = POLY->vertices[n].v.v[c];
				*( mxGetPr(DATA1) + n + c*N )= xyz[c];
				*( mxGetPr(DATA2) + n + c*N )= xyz[c];
				*( mxGetPr(DATA3) + n + c*N )= xyz[c];
				*( mxGetPr(DATA4) + n + c*N )= xyz[c];
			}
		}
		mxAddField( m1,    "xyz" );
		mxSetField( m1, 0, "xyz", DATA1 );
		mxAddField( m2,    "xyz" );
		mxSetField( m2, 0, "xyz", DATA2 );
		mxAddField( m3,    "xyz" );
		mxSetField( m3, 0, "xyz", DATA3 );
		mxAddField( m4,    "xyz" );
		mxSetField( m4, 0, "xyz", DATA4 );
	}

	//TRI

	/*Cuenta cuantas caras hay del owner solicitado*/
	N= POLY->faces.size();

	for( n= 0; n<N ; n++ )
	{
		if(POLY->faces[n].iIDpoly == 0)
			iNumFaces1++;
		if(POLY->faces[n].iIDpoly == 1)
			iNumFaces2++;
		if(POLY->faces[n].iIDpoly == 2)
			iNumFaces3++;
		if(POLY->faces[n].iIDpoly == 3)
			iNumFaces4++;
	}

	if( ( iNumFaces1 > 0 ) || ( iNumFaces2 > 0) || ( iNumFaces3 > 0)|| ( iNumFaces4 > 0) )
	{
		if( iNumFaces1 > 0 )
		{
			DATA1 = mxCreateDoubleMatrix( iNumFaces1 , 3 , mxREAL );
			iContFaces1 = 0;
		}
		if( iNumFaces2 > 0 )
		{
			DATA2 = mxCreateDoubleMatrix( iNumFaces2 , 3 , mxREAL );
			iContFaces2 = 0;
		}
		if( iNumFaces3 > 0 )
		{
			DATA3 = mxCreateDoubleMatrix( iNumFaces3 , 3 , mxREAL );
			iContFaces3 = 0;
		}
		if( iNumFaces4 > 0 )
		{
			DATA4 = mxCreateDoubleMatrix( iNumFaces4 , 3 , mxREAL );
			iContFaces4 = 0;
		}
		for( n= 0; n<N ; n++ )
		{
			nVert = POLY->faces[n].nVertices();
			if( nVert > 3)
			{
				mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
				mexErrMsgTxt("\n");
			}
			c = 0;
			if(POLY->faces[n].iIDpoly == 0)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA1) + iContFaces1 + c*iNumFaces1 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA1) + iContFaces1 + c*iNumFaces1 )= 0;
				}
				iContFaces1++;
			}
			if(POLY->faces[n].iIDpoly == 1)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA2) + iContFaces2 + c*iNumFaces2 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA2) + iContFaces2 + c*iNumFaces2 )= 0;
				}
				iContFaces2++;
			}
			if(POLY->faces[n].iIDpoly == 2)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA3) + iContFaces3 + c*iNumFaces3 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA3) + iContFaces3 + c*iNumFaces3 )= 0;
				}
				iContFaces3++;
			}
			if(POLY->faces[n].iIDpoly == 3)
			{
				for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
				{
					*( mxGetPr(DATA4) + iContFaces4 + c*iNumFaces4 ) = POLY->vertexToIndex((*j)) + 1;
					c++;
				}
				for( ; c<3 ; c++ ){
					*( mxGetPr(DATA4) + iContFaces4 + c*iNumFaces4 )= 0;
				}
				iContFaces4++;
			}
		}
		mxAddField( m1,    "tri" );
		mxSetField( m1, 0, "tri", DATA1 );
		mxAddField( m2,    "tri" );
		mxSetField( m2, 0, "tri", DATA2 );
		mxAddField( m3,    "tri" );
		mxSetField( m3, 0, "tri", DATA3 );
		mxAddField( m4,    "tri" );
		mxSetField( m4, 0, "tri", DATA4 );

		pOUT1 = m1;
		pOUT2 = m2;
		pOUT3 = m3;
		pOUT4 = m4;
	}
	iret = 0;
	return(iret);
}

//mxArray *carvePolyhedron2MESH( carve::poly::Polyhedron *POLY , int iIDFace){
//  const char  *names[] = {""};
//  const int   dims[1] = {1};
//  mxArray     *DATA;
//  mxArray     *IDFace;
//  mxArray     *m;
//  long        N, n, nVert;
//  double      xyz[3];
//  int         n_cols, c;
//  int         n_fields, f;
//  char        name[256];
//
//  m   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);
//
////NODES
//  N= POLY->vertices.size();
//  if( N> 0)
//  {
//    DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
//    for( n= 0; n<N ; n++ ) {
//      for( c=0 ; c<3 ; c++ )
//      {
//    	  xyz[c] = POLY->vertices[n].v.v[c];
//    	  *( mxGetPr(DATA) + n + c*N )= xyz[c];
//      }
//    }
//    mxAddField( m,    "xyz" );
//    mxSetField( m, 0, "xyz", DATA );
//
//  }
//
////TRI
//  carve::interpolate::FaceAttr<int> f_tex_num;
//
//  N= POLY->faces.size();
//  if( N > 0 ) {
//    DATA = mxCreateDoubleMatrix( N , 3 , mxREAL );
//    IDFace = mxCreateDoubleMatrix( N , 1 , mxREAL );
//    for( n= 0; n<N ; n++ )
//    {
//    	poly_t::face_t &f = POLY->faces[n];
//    	int t;
//    	/*if (f_tex_num.hasAttribute(&f))
//    	{
//    		t = f_tex_num.getAttribute(&f, 0);
//    	}*/
//    	t = f_tex_num.getAttribute(&f, 0);
//    	mexPrintf("Face % i ; De la t: %i\n", n, t);
//
//    	nVert = POLY->faces[n].nVertices();
//    	if( nVert > 3)
//    	{
//    		mexPrintf("Warning!!!: This face has %d vertex.\n", nVert);
//    		mexErrMsgTxt("\n");
//    	}
//    	c = 0;
//    	*(mxGetPr(IDFace) + n) = t;
//    	for (carve::poly::Face<3>::vertex_iter_t j = POLY->faces[n].vbegin(); j != POLY->faces[n].vend(); ++j)
//    	{
//    		*( mxGetPr(DATA) + n + c*N ) = POLY->vertexToIndex((*j)) + 1;
//    		c++;
//    	}
//    	for( ; c<3 ; c++ ){
//    		*( mxGetPr(DATA) + n + c*N )= 0;
//    	}
//    }
//    mxAddField( m,    "tri" );
//    mxSetField( m, 0, "tri", DATA );
//
//    mxAddField( m,    "triLabel" );
//    mxSetField( m, 0, "triLabel", IDFace );
//
//    /*n_fields= POLY->GetCellData()->GetNumberOfArrays();
//    for( f = 0 ; f < n_fields ; f++ ){
//      sprintf( name , "tri%s", POLY->GetCellData()->GetArray(f)->GetName() );
//      n_cols= POLY->GetCellData()->GetArray(f)->GetNumberOfComponents();
//      DATA = mxCreateDoubleMatrix( N , n_cols , mxREAL );
//      for( n= 0; n<N ; n++ ) {
//        for( c=0 ; c< n_cols ; c++ ){
//     *( mxGetPr(DATA) + n + c*N )= POLY->GetCellData()->GetArray(f)->GetComponent(n,c);
//        }
//      }
//      mxAddField( m,    name );
//      mxSetField( m, 0, name, DATA );
//    }*/
//  }
//
//  return(m);
//}
