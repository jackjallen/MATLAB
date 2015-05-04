#include "myMEX.h"

#if !defined( real )
  #define   real       double
#endif
#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif

#define InterpMACRO( x , y )              \
  if( x <= ITX[0] ){                      \
    y = ITY[0];                           \
  } else if( x >= ITX[N-1] ){             \
    y = ITY[N-1];                         \
  } else {                                \
    ii = 1;                               \
    while( x > ITX[ii] ){ii++;} ii--;     \
    y = ITY[ii] + ( ITY[ii+1]-ITY[ii] )/( ITX[ii+1]-ITX[ii] )*( x-ITX[ii] );  \
  }


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  double            *O;
  long              id=0, IJK;
  unsigned char     *dataUC;
  char              *dataC;
  unsigned short    *dataUS;
  short             *dataS;
  unsigned int      *dataUI;
  int               *dataI;
  float             *dataF;
  double            *dataD;
  mxLogical         *dataL;
  double             x;
  double            *ITX, *ITY;
  int               N, ii;
  
  ITX = myGetPr( prhs[1] );
  N   = myNumel( prhs[1] )/2;
  ITY = ITX + N;
  
  IJK = myNumel( prhs[0] );
  
  plhs[0] = myDuplicateArrayWithClass( prhs[0] , mxREAL_CLASS , mxREAL );
  O = (double *) mxGetData( plhs[0] );
  
  switch( mxGetClassID( prhs[0] ) ){
    case mxLOGICAL_CLASS:
      dataL = (mxLogical *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){
        if( dataL[id] ){
          O[id] = 1;
        } else {
          O[id] = 0;
        }
      }
      break;

    case mxUINT8_CLASS:
      dataUC = (unsigned char *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataUC[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxINT8_CLASS:
      dataC = (char *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataC[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxUINT16_CLASS:
      dataUS = (unsigned short *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataUS[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxINT16_CLASS:
      dataS = (short *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataS[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxUINT32_CLASS:
      dataUI = (unsigned int *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataUI[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxINT32_CLASS:
      dataI = (int *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataI[id];
        InterpMACRO( x , O[id] );
      }
      break;

    case mxSINGLE_CLASS:
      dataF = (float *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = (double) dataF[id];
        InterpMACRO( x , O[id] );
      }
      break;
    
    case mxDOUBLE_CLASS:
      dataD = (double *) mxGetData( prhs[0] );
      for( ; id < IJK ; id++ ){ 
        x = dataD[id];
        InterpMACRO( x , O[id] );
      }
      break;
  }

  EXIT: myFreeALLOCATES();
}
