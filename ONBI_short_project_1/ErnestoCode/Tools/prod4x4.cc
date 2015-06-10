
#include "newmyMEX.h"

#include "inv4x4.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ALLOCATES();

	real *Ndata, *Mdata, *O;
	int NX, NY, MX, MY;
	int Odims[2];

	if(nrhs < 2)
		mexErrMsgTxt("no enough input matrix!!.");
	if(nrhs > 2)
		mexPrintf("only two first matrix wil be multiplied!!.");

	if( (mxIsComplex(prhs[0])) || (mxIsComplex(prhs[1])) )
		mexErrMsgTxt("any input matrix is a complex data!!.");
	if( (mxIsEmpty(prhs[0])) || (mxIsEmpty(prhs[1])) )
		mexErrMsgTxt("any input matrix is empty!!.");
	if( (!mxIsNumeric(prhs[0])) || (!mxIsNumeric(prhs[1])) )
		mexErrMsgTxt("any input matrix is not numeric!!.");

	NX = mySize( prhs[0] , 0 );
	NY = mySize( prhs[0] , 1 );
	Ndata = (real*)myGetPr( prhs[0] );

	MX = mySize( prhs[1] , 0 );
	MY = mySize( prhs[1] , 1 );

	if( NY != MX ){ myErrMsgTxt("it is not able to multiply this matrix, error in dimensions!!."); }

	Mdata = (real*)myGetPr( prhs[1] );

	/*Creating output*/
	Odims[0] = NX;
	Odims[1] = MY;
	plhs[0] = mxCreateNumericArray( 2 , Odims , mxDOUBLE_CLASS , mxREAL );
	O = (real*)myGetPr( plhs[0] );
	/*END Creating output*/
	

	/* c=a*b (a=mxn; b=nxp; c=mxp) */
	prodm_matlab(O, Ndata, Mdata, NX, NY, MY);


	EXIT: myFreeALLOCATES();
}

