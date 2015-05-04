
#include "newmyMEX.h"

#include "inv4x4.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ALLOCATES();

	real *M, *O;
	int NX, NY;
	int Odims[2];

	if(nrhs < 1)
		mexErrMsgTxt("no input matrix!!.");
	if(nrhs > 1)
		mexPrintf("only first argument wil be inverted!!.");

	if(mxIsComplex(prhs[0]))
		mexErrMsgTxt("input matrix is a complex data!!.");
	if(mxIsEmpty(prhs[0]))
		mexErrMsgTxt("input matrix is empty!!.");
	if(!mxIsNumeric(prhs[0]))
		mexErrMsgTxt("input matrix is not numeric!!.");

	NX = mySize( prhs[0] , 0 );
	NY = mySize( prhs[0] , 1 );

	if( NX != NY ){ myErrMsgTxt("no square matrix!!."); }
	
	M = (real*)myGetPr( prhs[0] );

	/*Creating output*/
	Odims[0] = NX;
	Odims[1] = NY;
	plhs[0] = mxCreateNumericArray( 2 , Odims , mxDOUBLE_CLASS , mxREAL );
	O = (real*)myGetPr( plhs[0] );
	/*END Creating output*/
	
	switch (NX)
	{
	case 1:
		O[0] = 1/M[0];
		break;
	case 2:
		invt2x2m(O, M);
		break;
	case 3:
		invt3x3m(O, M);
		break;
	case 4:
		if( (M[3] == 0) && (M[7] == 0) && (M[11] == 0) && (M[15] == 1) )
		{
			invt4x4mh(O, M);
		}
		else
		{
			invt4x4m(O, M);
			//myErrMsgTxt("inverting matrix 4 x 4 not homogeneous is not implemented yet.");
		}
			break;
	default:
		myErrMsgTxt("max size matrix for inverting: 4 x 4.");
		break;
	}

	  EXIT: myFreeALLOCATES();
}

