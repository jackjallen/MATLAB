#include "libParseargs.h"

//#include "mex.h"
//#include "myMEX.h"
#include "mmz.h"
#include "mz.h"

/*#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS*/
#undef real
#undef min
#undef max

#define carveOBJ_TYPE      carve::csg::CSG::OP

#include <carve/csg.hpp>
#include <carve/tree.hpp>
#include <carve/csg_triangulator.hpp>
#include <carve/csg_collector.hpp>

#include "MESH2carvePolyhedron.h"
#include "carvePolyhedron2MESH.h"

struct carve_options
{
	double epsilon;
	bool canonicalize;
	bool rescale;
	bool improve;
	carve::csg::CSG::CLASSIFY_TYPE classifier;

	carve_options(): epsilon(carve::EPSILON), canonicalize(false), rescale(false), improve(false), classifier(carve::csg::CSG::CLASSIFY_NORMAL){}
};


#include <stdlib.h>
void exit (int& status){throw;};

static void clearMEX(void)
{
	/*if(MESH1 != NULL);{ free(MESH1); MESH1 = NULL;}
	if(MESH2 != NULL);{ free(MESH2); MESH2 = NULL;}
	if(MESH_OUT != NULL);{ free(MESH_OUT); MESH_OUT = NULL;}*/
	//myFreeALLOCATES();
	mexPrintf("Borrado de MESH1, MESH2 y MESH_OUT en el mexAtExit.\n");
	mexPrintf("\n\nClearing mex...\n\n");
	return;
}

carve::poly::Polyhedron *compute_direct( carve::poly::Polyhedron  *a, carve::poly::Polyhedron  *b, carve::csg::CSG::OP op);
carve::poly::Polyhedron *compute_by_pipeline( carve::poly::Polyhedron  *a, carve::poly::Polyhedron  *b,
		carve::csg::CSG::OP op, carve::csg::CSG_TreeNode **TREE, carve_options options);

mxArray* createOutputEmpty();

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	char                STR[2048];

	//if( (nrhs != 3) || ( !mxIsStruct(prhs[0]) ) || ( !mxIsStruct(prhs[1]) ) || ( !mxIsChar(prhs[2]) ))
	if( (nrhs == 0) || (nrhs < 3) )
	{
		/*Help Example	*/
		mexPrintf("carveIntersect( IN_MESH1, IN_MESH2, OPERATION)\n");
		mexPrintf("\n");
		mexPrintf("OPERATION\t, string       ... Specify boolean operation ==\n");
		mexPrintf("\t\t{'u','union','+'}, \n"
				"\t\t{'i','int','intersect','intersection'},  \n"
				"\t\t{'sd','symdiff','symmetricdifference'},  \n"
				"\t\t{'minus','-','a_minus_b'},  \n"
				"\t\t{'b_minus_a'}\n\n"
                "\t\t{'all'}, \n"
				"\t\t{'split'}, \n"
				"\t\t{'splitFull'}\n\n");
		mexPrintf("\t\tAdditional options: \n"
				"\t\t{'Canonicalize' (false), 'Rescale' (false), 'Edge' (NORMAL), IMprove (false), 'EPsilon' (sqrt(EPS(1)))},  \n");
		mexPrintf("\n");
		return;
	}
	if( ( !mxIsStruct(prhs[0]) ) || ( !mxIsStruct(prhs[1]) ) || ( !mxIsChar(prhs[2]) ))
	{
		/*Help Example	*/
		mexPrintf("carveIntersect( IN_MESH1, IN_MESH2, OPERATION)\n");
		mexPrintf("\n");
		mexPrintf("OPERATION\t, string       ... Specify boolean operation ==\n");
		mexPrintf("\t\t{'u','union','+'}, \n"
				"\t\t{'i','int','intersect','intersection'},  \n"
				"\t\t{'sd','symdiff','symmetricdifference'},  \n"
				"\t\t{'minus','-','a_minus_b'},  \n"
				"\t\t{'b_minus_a'}, \n"
				"\t\t{'all'}, \n"
				"\t\t{'split'}, \n"
				"\t\t{'splitFull'}\n\n");

		mexPrintf("\t\tAdditional options: \n"
				"\t\t{'Canonicalize' (false), 'Rescale' (false), 'Edge' (NORMAL), IMprove (false), 'EPsilon' (sqrt(EPS(1)))},  \n");
		mexPrintf("\n");

		if(( !mxIsStruct(prhs[0]) ) || ( !mxIsStruct(prhs[1]) ))
		{
			mexPrintf("MESH Data have to be both Matlab structure types.\n");
		}
		if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
		mexErrMsgTxt("\n");
		//return;
	}

	//Operacion
	carve::csg::CSG::OP Op;
	mxGetString( prhs[2], STR, 2047 );
	try{
		CallMethod( &Op , STR );
	}catch(int i)
	{
		plhs[0] = createOutputEmpty();
		mexErrMsgTxt("Carve select boolean operation fail.\n");
	}
	switch(Op)
	{
		case carve::csg::CSG::UNION:
		case carve::csg::CSG::INTERSECTION:
		case carve::csg::CSG::A_MINUS_B:
		case carve::csg::CSG::B_MINUS_A:
		case carve::csg::CSG::SYMMETRIC_DIFFERENCE:
		case carve::csg::CSG::ALL:
		{
			if(nlhs > 1)
			{
				mexPrintf("Operation %s has only one output!!.\n", STR);
				mexErrMsgTxt("\n");
			}
			break;
		}
		case carve::csg::CSG::SPLIT:
		{
			if(nlhs != 3)
			{
				mexPrintf("Operation %s has three outputs!!.\n", STR);
				mexErrMsgTxt("\n");
			}
			break;
		}
		case carve::csg::CSG::SPLIT_FULL:
		{
			if( nlhs != 6 )
			{
				mexPrintf("Operation %s has six outputs!!.\n", STR);
				mexErrMsgTxt("\n");
			}
			break;
		}
		default:
		{
			mexPrintf("Operation %s no valid!!.\n", STR);
			mexErrMsgTxt("\n");
		}
	}
	mexAtExit(clearMEX);

	ALLOCATES();
	carve::poly::Polyhedron *MESH1 = NULL;
	carve::poly::Polyhedron *MESH2 = NULL;
	carve::poly::Polyhedron *MESH_OUT = NULL;
	carve::csg::CSG_TreeNode *TREE = NULL;

	//variables de configuracion del proceso
	carve_options options;

	int iIdxFirstArg = 3;
	if(nrhs > iIdxFirstArg)
	{
		/*Parseado de argumentos*/
		parseargs myParseargs;
		int iNumArgsOut, iRet, iNumKeys = 2;
		mxArray **ppArgsOut = NULL;
		mxArray **ppArgsKeysOut;

		//Debemos inicializar el parser para que todas las lineas a continuacion sean iguales.
		//Como solo nos llegaran las claves y no el CELL con el varargin, nos lo tenemos que fabricar:
		//varargin = {'Canonicalize, Rescale, Edge, EPsilon, carve::EPSILON'}

		/*Creacion de las keys a comprobar*/
		mxArray** varargin;
		varargin = new mxArray*[nrhs - iIdxFirstArg + 1];	//numero de parametros de entrada menos los 3 de los MESH1,MESH2 y operacion y mas el CELL a construir

		//los mxArray de entrada que son las claves a buscar, + el CELL con todos los posibles en donde los buscaremos:
		//CELL con todas las claves
#define iNA 5
		char str[iNA][13] = {"Canonicalize", "Rescale", "Edge", "EPsilon", "IMprove"};
		double dEpsilon = (carve::EPSILON);

		varargin[0] = myParseargs.createCELLArgs((iNA+1), NULL);
		for(int j = 0; j < iNA; j++)
		{
			varargin[0] = myParseargs.InsertArgInCELL(j, str[j], mxCHAR_CLASS);
		}
		varargin[0] = myParseargs.InsertArgInCELL(iNA, &dEpsilon, mxDOUBLE_CLASS);

#ifdef _DEBUG_
		myParseargs.printCELLArgs();
#endif


		/*
	  mwSize ndim = 2;
	  mwSize dims[] = {1,5};
	  varargin[0] = mxCreateCellArray(ndim,dims);

	  for(int j = 0; j < 4; j++)
	  {
		  mxSetCell(((mxArray*)varargin[0]), j, mxCreateString(str[j]) );
	  }
	  mxSetCell(((mxArray*)varargin[0]), 4, mxCreateDoubleScalar(carve::EPSILON));
		 */

		//Claves a buscar
		for(int i = 0; i < (nrhs - iIdxFirstArg); i++)
		{
			varargin[i+1] = ((mxArray*)prhs[iIdxFirstArg + i]);
		}
		//mexPrintf("Clave %i = %g.\n", 4 , (*((double*)mxGetData(varargin[5]))) );

		/*Fin Creacion de las keys a comprobar*/

		/*Parser*/
		//De esta forma, se le pasa solo los argumentos que son keys
		iRet = myParseargs.InitParser(&iNumArgsOut, ppArgsOut, (nrhs - iIdxFirstArg + 1), varargin);

#ifdef _DEBUG_
		mexPrintf("return init parser %i, iNumArgsOut=%i.\n", iRet,iNumArgsOut );
#endif

		//buscar en varargin ( "c", "canonicalize" )
		if(myParseargs.parse(0, &iNumArgsOut, ppArgsOut, iNumArgsOut, (mxArray**)ppArgsOut, ppArgsKeysOut, false, iNumKeys, "canonicalize", "c") > 0)
		{
			options.canonicalize = true;
#ifdef _DEBUG_
			mexPrintf("Encontrado canonicalize().\n");
#endif

		}

		//buscar en varargin ( "rescale", "r" )
		if(myParseargs.parse(0, &iNumArgsOut, ppArgsOut, iNumArgsOut, (mxArray**)ppArgsOut, ppArgsKeysOut, false, iNumKeys, "rescale", "r") > 0)
		{
			options.rescale = true;
#ifdef _DEBUG_
			mexPrintf("Encontrado rescale().\n");
#endif
		}

		//buscar en varargin ( "edge", "e" )
		if(myParseargs.parse(0, &iNumArgsOut, ppArgsOut, iNumArgsOut, (mxArray**)ppArgsOut, ppArgsKeysOut, false, iNumKeys, "edge", "e") > 0)
		{
			options.classifier = carve::csg::CSG::CLASSIFY_EDGE;
#ifdef _DEBUG_
			mexPrintf("Encontrado edge.\n");
#endif
		}

		//buscar en varargin ( "IMprove", "IM" )
		if(myParseargs.parse(0, &iNumArgsOut, ppArgsOut, iNumArgsOut, (mxArray**)ppArgsOut, ppArgsKeysOut, false, iNumKeys, "improve", "im") > 0)
		{
			options.improve = true;
#ifdef _DEBUG_
			mexPrintf("Encontrado improve.\n");
#endif
		}

		//buscar en varargin ( "epsilon", "p" )
		if(myParseargs.parse(1, &iNumArgsOut, ppArgsOut, iNumArgsOut, (mxArray**)ppArgsOut, ppArgsKeysOut, false, iNumKeys, "epsilon", "ep") > 0)
		{
			options.epsilon = *((double*)mxGetData(ppArgsKeysOut[0]));
#ifdef _DEBUG_
			mexPrintf("Encontrado epsion = %g.\n", options.epsilon);
#endif
		}

		/*Fin Parser*/

		/*destroy argumentos artificiales*/
		mxDestroyArray(varargin[0]);
		delete [] varargin;

	}
	/*fin de parseado de argumentos*/

	//Procesamiento
	MESH1 = MESH2carvePolyhedron( prhs[0] );
	MESH2 = MESH2carvePolyhedron( prhs[1] );

	if(MESH1 == NULL)
		mexErrMsgTxt("Error in MESH1 building...\n");
	if(MESH2 == NULL)
		mexErrMsgTxt("Error in MESH2 building...\n");
#ifdef _DEBUG_MEX

	//Prueba de busqueda de propietarios
	int iNumFOut, iNumFA = 0, iNumFB = 0;
	iNumFOut = MESH1->faces.size();
	for(int iNF = 0; iNF < iNumFOut; iNF++)
	{
		if(MESH1->faces[iNF].owner == MESH1)
		{
			iNumFA++;
		}
		if(MESH1->faces[iNF].owner == MESH2)
		{
			iNumFB++;
		}
	}
	mexPrintf("Fin del MESH2carvePolyhedron. Numero de caras totales de MESH_1 = %d, Num Caras que pretenecen a MESH1 = %i, Num Caras que pretenecen a MESH2 = %i.\n", iNumFOut, iNumFA, iNumFB);
	iNumFA = 0, iNumFB = 0;
	iNumFOut = MESH2->faces.size();
	for(int iNF = 0; iNF < iNumFOut; iNF++)
	{
		if(MESH2->faces[iNF].owner == MESH1)
		{
			iNumFA++;
		}
		if(MESH2->faces[iNF].owner == MESH2)
		{
			iNumFB++;
		}
	}
	mexPrintf("Fin del MESH2carvePolyhedron. Numero de caras totales de MESH_2 = %d, Num Caras que pretenecen a MESH1 = %i, Num Caras que pretenecen a MESH2 = %i.\n", iNumFOut, iNumFA, iNumFB);
#endif

	/*//Operacion
  mxGetString( prhs[2], STR, 2047 );
  try{
	  CallMethod( &Op , STR );
  }catch(int i)
  {
	  plhs[0] = createOutputEmpty();
	  delete MESH1;
	  delete MESH2;
	  myFreeALLOCATES();
	  mexErrMsgTxt("Carve select boolean operation fail.\n");
  }*/
	//Intersection
	try
	{
		//MESH_OUT = compute_direct(MESH1,MESH2,carve::csg::CSG::A_MINUS_B);
#ifdef _DEBUG_MEX
		mexPrintf("Iniciando Interseccion.\n");
#endif

		MESH_OUT = compute_by_pipeline(MESH1,MESH2, Op, &TREE, options);

#ifdef _DEBUG_MEX
		mexPrintf("Fin de Interseccion.\n");
#endif
	}catch(int i)
	{
		plhs[0] = createOutputEmpty();
		//delete MESH1;
		//delete MESH2;
		if(TREE != NULL) delete TREE;
		if(MESH_OUT != NULL) delete MESH_OUT;
		mmzFreeALLOCATES();
		mexErrMsgTxt("Carve exec boolean operation fail.\n");
	}
	catch(...)
	{
		plhs[0] = createOutputEmpty();
		if(TREE != NULL) delete TREE;
		if(MESH_OUT != NULL) delete MESH_OUT;
		mmzFreeALLOCATES();
		mexErrMsgTxt("Carve exec boolean operation results in unknown catastrofic fail.\n");
	}

	if(MESH_OUT != NULL)
	{
#ifdef _DEBUG_MEX

		//Prueba de busqueda de propietarios
		iNumFOut = MESH1->faces.size();
		iNumFA = 0, iNumFB = 0;
		for(int iNF = 0; iNF < iNumFOut; iNF++)
		{
			if(MESH1->faces[iNF].owner == MESH1)
			{
				iNumFA++;
			}
			if(MESH1->faces[iNF].owner == MESH2)
			{
				iNumFB++;
			}
		}
		mexPrintf("Fin del MESH2carvePolyhedron. Numero de caras totales de MESH_1 = %d, Num Caras que pretenecen a MESH1 = %i, Num Caras que pretenecen a MESH2 = %i.\n", iNumFOut, iNumFA, iNumFB);


		iNumFA = 0, iNumFB = 0;
		iNumFOut = MESH_OUT->faces.size();
		for(int iNF = 0; iNF < iNumFOut; iNF++)
		{
			if(MESH_OUT->faces[iNF].owner_orig == MESH1)
			{
				iNumFA++;
			}
			if(MESH_OUT->faces[iNF].owner_orig == MESH2)
			{
				iNumFB++;
			}
		}
		mexPrintf("Fin del procesamiento. Numero de caras totales de MESH_OUT = %d, Num Caras que pretenecen a MESH1 = %i, Num Caras que pretenecen a MESH2 = %i.\n", iNumFOut, iNumFA, iNumFB);


		mexPrintf("Owners: MESH1 = %p, MESH1_orig = %p, MESH2 = %p, MESH2_orig = %p, MESH_OUT = %p, MESH_OUT_orig = %p\n",
				MESH1->faces[0].owner, MESH1->faces[0].owner_orig, MESH2->faces[0].owner, MESH2->faces[0].owner_orig, MESH_OUT->faces[0].owner, MESH_OUT->faces[0].owner_orig);
#endif

		switch(Op)
		{
		case carve::csg::CSG::UNION:
		case carve::csg::CSG::INTERSECTION:
		case carve::csg::CSG::A_MINUS_B:
		case carve::csg::CSG::B_MINUS_A:
		case carve::csg::CSG::SYMMETRIC_DIFFERENCE:
		case carve::csg::CSG::ALL:
		{
			plhs[0] = carvePolyhedron2MESH( MESH_OUT);
			break;
		}
		case carve::csg::CSG::SPLIT:
		{
			/*plhs[0] = carvePolyhedron2MESH( MESH_OUT, 0);
			plhs[1] = carvePolyhedron2MESH( MESH_OUT, 1);*/
			/*carvePolyhedron2MESH( MESH_OUT, plhs[0], plhs[1]);*/
			plhs[0] = carvePolyhedron2MESH( MESH_OUT, 0, NULL);
			plhs[1] = carvePolyhedron2MESH( MESH_OUT, 1, plhs[0]);
			plhs[2] = carvePolyhedron2MESH( MESH_OUT, 2, plhs[0]);


			/*funcCreate pFuncCreate;
			pFuncCreate = GetCreatePointer();

			plhs[1] = pFuncCreate((void *)prhs[0]);*/
			break;
		}
		case carve::csg::CSG::SPLIT_FULL:
		{
			/*plhs[0] = carvePolyhedron2MESH( MESH_OUT, 0);
			plhs[1] = carvePolyhedron2MESH( MESH_OUT, 1);
			plhs[2] = carvePolyhedron2MESH( MESH_OUT, 2);
			plhs[3] = carvePolyhedron2MESH( MESH_OUT, 3);*/
			/*carvePolyhedron2MESH( MESH_OUT, plhs[0], plhs[1], plhs[2], plhs[3]);*/

			plhs[0] = carvePolyhedron2MESH( MESH_OUT, 0, NULL);
			plhs[1] = carvePolyhedron2MESH( MESH_OUT, 1, plhs[0]);
			plhs[2] = carvePolyhedron2MESH( MESH_OUT, 2, plhs[0]);
			plhs[3] = carvePolyhedron2MESH( MESH_OUT, 3, NULL);
			plhs[4] = carvePolyhedron2MESH( MESH_OUT, 4, plhs[3]);
			plhs[5] = carvePolyhedron2MESH( MESH_OUT, 5, plhs[3]);
			break;
		}
		}
	}
	else
	{
		for(int i = 0; i < nlhs; i++)
			plhs[i] = createOutputEmpty();
	}

	EXIT:
	for(int i = 0; i < nlhs; i++)
		if(plhs[i] == NULL){plhs[i] = createOutputEmpty(); mexPrintf("Building output MESh fail.\n");}
	//delete MESH1;
	//delete MESH2;
	if(TREE != NULL) delete TREE;
	if(MESH_OUT != NULL) delete MESH_OUT;
	mmzFreeALLOCATES();
}

carve::poly::Polyhedron *compute_direct( carve::poly::Polyhedron  *a, carve::poly::Polyhedron  *b, carve::csg::CSG::OP op)
{
	carve::csg::CSG myCSG;

	myCSG.hooks.registerHook(new carve::csg::CarveTriangulator, carve::csg::CSG::Hooks::PROCESS_OUTPUT_FACE_BIT);

	//Interseccion
	return myCSG.compute(a, b, carve::csg::CSG::A_MINUS_B);
}
carve::poly::Polyhedron *compute_by_pipeline(
		carve::poly::Polyhedron *a,
		carve::poly::Polyhedron *b,
		carve::csg::CSG::OP op,
		carve::csg::CSG_TreeNode **TREE, carve_options options
)
{
	try
	{
		carve::poly::Polyhedron  *result  = NULL;
		carve::csg::CSG myCSG;
		carve::csg::CSG_TreeNode *lhs  = NULL;
		carve::csg::CSG_TreeNode *rhs  = NULL;
		//carve::csg::CSG_TreeNode *TREE  = NULL;

		lhs = new carve::csg::CSG_PolyNode(a, true);
		if(lhs == NULL){ mexPrintf(" Warning in create treeNodes\n" ); return NULL; }

		rhs = new carve::csg::CSG_PolyNode(b, true);
		if(rhs == NULL) { mexPrintf(" Warning in create treeNodes\n" ); delete lhs; return NULL;}

		/*TREE = new carve::csg::CSG_OPNode(
		 * lhs,
		 * rhs,
		 * operation,
		 * options.rescale,
		 * options.classifier);*/
		(*TREE) = new carve::csg::CSG_OPNode(lhs, rhs,
				op,
				options.rescale,
				options.classifier);

		if( ((*TREE) == NULL) )
		{
			mexPrintf(" Warning in create TREE\n" );
			if(lhs != NULL) delete lhs;
			if(rhs != NULL) delete rhs;
			return NULL;
		}
		/*else
			mexPrintf("TREE created\n" );*/

		if (options.improve)
		{
			myCSG.hooks.registerHook(new carve::csg::CarveTriangulatorWithImprovement, carve::csg::CSG::Hooks::PROCESS_OUTPUT_FACE_BIT);
#ifdef _DEBUG_MEX
			mexPrintf("Processing With Improve...\n");
#endif
		}
		else
			myCSG.hooks.registerHook(new carve::csg::CarveTriangulator, carve::csg::CSG::Hooks::PROCESS_OUTPUT_FACE_BIT);

		//mexPrintf("Hook registered\n" );

		/*carve::setEpsilon(1.7);
		carve::setEpsilon(1e-8);*/

		//Interseccion
		carve::setEpsilon(options.epsilon);
		result = (*TREE)->eval(myCSG);

		if (options.canonicalize) result->canonicalize();

		if( (result == NULL) )
		{
			mexPrintf(" Warning in create results\n" );
			if((*TREE) != NULL) delete (*TREE);
			return NULL;
		}
		//delete TREE;

		//throw ("pp");
		return result;
	}
	catch(int i)
	{
		mexPrintf("Catastrofic Fail!! \n");
		throw (-1);
	}
	catch(...)
	{
		mexPrintf("Unknown Catastrofic Fail!! \n");
		throw;
	}
}

void CallMethod( carveOBJ_TYPE *O , char *met ) throw (int)
{
	try{
		Call_UNION( "u", "union", "+" );
		Call_INTERSECTION( "i", "int", "intersect", "intersection" );
		Call_SYMM_DIFF( "sd", "symdiff", "symmetricdifference" );
		Call_A_MINUS_B( "minus" , "-", "a_minus_b" );
		Call_B_MINUS_A( "b_minus_a" );
		Call_ALL( "all" );
		Call_SPLIT( "split" );
		Call_SPLIT_FULL( "splitFull" );

		mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); mmzFlush();
		//mexErrMsgTxt("\n");
		throw(-1);
	}
	catch(int i)
	{
		throw (-1);
	}
}

mxArray* createOutputEmpty()
{
	const char  *names[] = {""};
	const int   dims[1] = {1};
	mxArray     *m;

	m   = mxCreateStructArray(1, (const mwSize*)(dims), 0, names);

	mxAddField( m,    "xyz" );
	mxSetField( m, 0, "xyz", NULL );

	mxAddField( m,    "tri" );
	mxSetField( m, 0, "tri", NULL );

	return m;
}
