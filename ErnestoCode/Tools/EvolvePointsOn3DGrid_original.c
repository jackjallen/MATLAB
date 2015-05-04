/*
   phi = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , points , 
                           BOUNDARY_MODE , [ BOUNDARY_SIZE ]
                          'maxSTEPS'      , 50                    dt_min = 1/maxSTEPS
                          'minStepsPerVoxel'      , 1             dt_max = S/V/minStepsPerVoxel   (S, min voxel size ) (V, max V(x) )
                          'relTOL'        , 1/100                 at each step control the error
                          'flow'          , [ 0  1 ]              a vector of Ts where evaluate the output FLOW (exept at the first t)
                          { 'rk23'  'euler' 's_euler' 'rk45'  }
                          { 'jac_plus' 'jac_times' }
                          'mask' , logical( numel(points) )
                          'initJACs'  ,  numeric( [npoints x 3 x 3 ] )
                          'matrix' , eye(4)
                          'verbose'                               print percentage every 10%
                       );

      BOUNDARY_MODE:
           'value'   (default)              stop the evolution outside the grid
           'closest'
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY)/2

 
  phi = EvolvePointsOn3DGrid( V , Gx , Gy , Gz   , { grillaX , grillaY , grillaZ } , 'matrix' , [4x4] )
  phi = EvolvePointsOn3DGrid( V , Gx , Gy , '2d' , { grillaX , grillaY } , 'matrix' , [4x4] )

  
  if FLOW is a scalar, integrate between t0 = 0   and  t_end = FLOW
  if FLOW is a vector, integrate with t0 = F(1)  and  t_end = F(end)
      and return, coordinates, jacobians and determinants at times F(2:end)
  if you also want the initial coords, jacs and dets use for example
      F = [ 0 , 0 , 0.1 , 0.2 , 0.3 , 0.5 , 1.0 , 2.0 ]  (replicate t0!!)
 


   phi = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , points , 'initJACs' , init_jacs )

   [phi,jac] = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , points , 'initJACs' , init_jacs )

   [phi,jac,det] = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , points , 'initJACs' , init_jacs )

 
    numel( points )  have to be multiple of NSD=3 or of NSD=2
 
    size( initJACS ) have to be equal to [ size(points)  NSD ]
 
 
 
    size( phi )  =  [ size( points , 1:end-1 )  ,  numel(FLOW)-1  , NSD ]
    size( jac )  =  [ size( points , 1:end-1 )  ,  numel(FLOW)-1  , NSD , NSD ]
    size( det )  =  [ size( points , 1:end-1 )  ,  numel(FLOW)-1  ]
 
*/

#include "myMEX.h"
#if !defined( memcpy )
  #include "string.h"
#endif


#define    integration_method_def   RK23
#define    jacobian_method_def      PLUS
#define    boundary_mode_def        VALUE
#define    maxSTEPS_def             50
#define    minStepsPerVoxel_def     1
#define    relTOL_def               1.0/100


#define times( o , v , s )        o[0] = v[0]*s; o[1] = v[1]*s; o[2] = v[2]*s;

#define Mtimes3x3( o , a , b)     o[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];  \
                                  o[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];  \
                                  o[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];  \
                                  o[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];  \
                                  o[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];  \
                                  o[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];  \
                                  o[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];  \
                                  o[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];  \
                                  o[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8]


#define EVOLVE_JAC          if( nlhs > 1 ){                                                \
                              getDV( XYZT , Dv );                                            \
                              switch( jacobian_method ){                                     \
                                case PLUS:                                                   \
                                  Mtimes3x3( dJAC , Dv , JAC );                              \
                                  for( d = 0 ; d<9 ; d++ ){  JAC[d] += h*dJAC[d]; }          \
                                  break;                                                     \
                                case TIMES:                                                  \
                                  for( d = 0 ; d<9 ; d++ ){ Dv[d] *= h; }                    \
                                  MatrixExp3x3( Dv );                                        \
                                  Mtimes3x3( dJAC , Dv , JAC );                              \
                                  memcpy( JAC , dJAC , 9*sizeof( real ) );                   \
                                  break;                                                     \
                              }                                                              \
                            }

#define MAT(i,j)      MAT[ 4*(j) + (i) ]


#define     readBOUNDS                                                            \
      if( nrhs > argN && ! mxIsChar(prhs[argN]) ){                                \
        if( !myIsEmpty( prhs[argN] ) ){                                           \
          if( myNumel( prhs[argN] ) == 1 ){                                       \
            boundary_size[0] = myGetValue(prhs[argN]);                            \
            boundary_size[1] = boundary_size[0];                                  \
            boundary_size[2] = boundary_size[0];                                  \
            boundary_size[3] = boundary_size[0];                                  \
            boundary_size[4] = boundary_size[0];                                  \
            boundary_size[5] = boundary_size[0];                                  \
          } else if( myNumel( prhs[argN] ) == 6 ){                                \
            boundary_size[0] = myGetValueIth( prhs[argN] , 0 );                   \
            boundary_size[1] = myGetValueIth( prhs[argN] , 1 );                   \
            boundary_size[2] = myGetValueIth( prhs[argN] , 2 );                   \
            boundary_size[3] = myGetValueIth( prhs[argN] , 3 );                   \
            boundary_size[4] = myGetValueIth( prhs[argN] , 4 );                   \
            boundary_size[5] = myGetValueIth( prhs[argN] , 5 );                   \
          } else {                                                                \
            myErrMsgTxt("1 or 6 values expected as BoundarySize");                \
          }                                                                       \
        }                                                                         \
        argN++; continue;                                                         \
      }



/*
for( d = 0 ; d<9 ; d++ ){ Dv[d] *= h; }                                                                  \
mexPrintf("\n maxnorm( expm( [ %0.15e  %0.15e  %0.15e ; %0.15e  %0.15e  %0.15e ; %0.15e  %0.15e %0.15e ] ) - ",Dv[0],Dv[3],Dv[6],Dv[1],Dv[4],Dv[7],Dv[2],Dv[5],Dv[8] );
MatrixExp3( Dv );                                                                  \
mexPrintf("[ %0.15e  %0.15e  %0.15e ; %0.15e  %0.15e  %0.15e ; %0.15e  %0.15e %0.15e ] )\n\n", Dv[0],Dv[3],Dv[6],Dv[1],Dv[4],Dv[7],Dv[2],Dv[5],Dv[8] );
*/

/*
// real  Ct45[8] = { 0 , 0 , 0.20000000000000001 , 0.29999999999999999 , 0.80000000000000004 , 0.88888888888888884 , 1 , 1 };
*/
#define  Ct45_1         0.0   
#define  Ct45_2         0.20000000000000001
#define  Ct45_3         0.29999999999999999
#define  Ct45_4         0.80000000000000004
#define  Ct45_5         0.88888888888888884
#define  Ct45_6         1.0
#define  Ct45_7         1.0
/*
// real  A45[8][8] = { { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
//                   { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
//                   { 0 , 0.20000000000000001 , 0 , 0 , 0 , 0 , 0 , 0} ,
//                   { 0 , 0.074999999999999997 , 0.22500000000000001 , 0 , 0 , 0 , 0 , 0} ,
//                   { 0 , 0.97777777777777775 , -3.7333333333333334 , 3.5555555555555554 , 0 , 0 , 0 , 0 } ,
//                   { 0 , 2.9525986892242035 , -11.595793324188385 , 9.8228928516994358 , -0.29080932784636487 , 0 , 0 , 0 } ,
//                   { 0 , 2.8462752525252526 , -10.757575757575758 , 8.9064227177434727 , 0.27840909090909088 , -0.2735313036020583 , 0 , 0 } ,
//                   { 0 , 0.091145833333333329 , 0 , 0.44923629829290207 , 0.65104166666666663 , -0.322376179245283 , 0.13095238095238096 , 0 } 
//                 };
*/
#define  A45_21   0.20000000000000001
#define  A45_31   0.074999999999999997
#define  A45_32   0.22500000000000001
#define  A45_41   0.97777777777777775
#define  A45_42   -3.7333333333333334
#define  A45_43   3.5555555555555554
#define  A45_51   2.9525986892242035
#define  A45_52   -11.595793324188385
#define  A45_53   9.8228928516994358
#define  A45_54   -0.29080932784636487
#define  A45_61   2.8462752525252526
#define  A45_62   -10.757575757575758
#define  A45_63   8.9064227177434727
#define  A45_64   0.27840909090909088
#define  A45_65   -0.2735313036020583
#define  A45_71   0.091145833333333329
#define  A45_73   0.44923629829290207
#define  A45_74   0.65104166666666663
#define  A45_75   -0.322376179245283
#define  A45_76   0.13095238095238096                
/*
// real  B45[8] = { 0 , 0.091145833333333329  , 0 , 0.44923629829290207 , 0.65104166666666663 , -0.322376179245283 , 0.13095238095238096 , 0 };
*/
#define  B45_1         0.091145833333333329   
#define  B45_2         0.0
#define  B45_3         0.44923629829290207
#define  B45_4         0.65104166666666663
#define  B45_5        -0.322376179245283
#define  B45_6         0.13095238095238096
#define  B45_7         0.0
/*
// real  D45[8] = { 0 , 0.0012326388888888888 , 0 , -0.0042527702905061394 , 0.036979166666666667 , -0.05086379716981132 , 0.041904761904761903 , -0.025000000000000001 };
*/
#define  D45_1         0.0012326388888888888   
#define  D45_2         0.0
#define  D45_3        -0.0042527702905061394
#define  D45_4         0.036979166666666667
#define  D45_5        -0.05086379716981132
#define  D45_6         0.041904761904761903
#define  D45_7        -0.025000000000000001






/*
// real  Ct23[5] = { 0 , 0 , 0.5 , 0.75 , 1 };
*/
#define  Ct23_1         0.0
#define  Ct23_2         0.5
#define  Ct23_3         0.75
#define  Ct23_4         1.0
/*
// real  A23[5][5] = { {  0 , 0 , 0 , 0 , 0 } , 
//                     {  0 , 0 , 0 , 0 , 0 } , 
//                     {  0 , 0.5 , 0 , 0 , 0 } , 
//                     {  0 , 0 , 0.75 , 0 , 0} , 
//                     {  0 , 0.22222222222222221 , 0.33333333333333331 , 0.44444444444444442 , 0}
//                   };
*/
#define  A23_21   0.5
#define  A23_32   0.75
#define  A23_41   0.22222222222222221
#define  A23_42   0.33333333333333331
#define  A23_43   0.44444444444444442
/*
// real  B23[5] = {  0 , 0.22222222222222221 , 0.33333333333333331 , 0.44444444444444442 , 0};
*/
#define  B23_1   0.22222222222222221
#define  B23_2   0.33333333333333331
#define  B23_3   0.44444444444444442
/*
// real  D23[5] = { 0 , -0.069444444444444448 , 0.083333333333333329 , 0.1111111111111111 , -0.125 };
*/
#define  D23_1   -0.069444444444444448
#define  D23_2   0.083333333333333329
#define  D23_3   0.1111111111111111
#define  D23_4   -0.125



enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST } boundary_mode;
real    *VX, *VY, *VZ, *X, *Y, *Z, OX, OY, OZ, LX, LY, LZ;
int     I , J , K , IJ, IJK;
real    boundary_size[6];
int     NSD;

typedef struct CELL {
  int  ii, jj, kk;
  real X0, X1, Y0, Y1, Z0, Z1;
  real VX0000, VX1000, VX0100, VX1100, VX0010, VX1010, VX0110, VX1110;
  real VY0000, VY1000, VY0100, VY1100, VY0010, VY1010, VY0110, VY1110;
  real VZ0000, VZ1000, VZ0100, VZ1100, VZ0010, VZ1010, VZ0110, VZ1110;
  real DX, DY, DZ, D, iD;
  int  isOutside;
} CELL; 
CELL C = { -1 , -1 , -1 , 0 , -1 , 0 , -1 , 0 , -1 };


_inline void setCELL( real x , real y , real z );
_inline void getV( real *xyzt , real *v );
_inline void getDV( real *xyzt , real *Dv );
        void MatrixExp3x3( real * );
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    integration_methods  { EULER , S_EULER , RK45 , RK23 }    integration_method;
  enum    jacobian_methods     { PLUS   , TIMES }                   jacobian_method;

  real    XYZT[4], xk[4], v[3], vp[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], err[3], Dv[9];
  real    xxnr, yynr, zznr;
  real    maxV, minD, hmax, hmin, tol, h;
  real    minStepsPerVoxel, relTOL;

  real    prct , last_prct=0;
  unsigned char         *MASK, VERBOSE;
  real    JAC[9], dJAC[9];
  
  real    *XYZ, *F, *Xn, *Yn, *Zn;
  int     Odims[50], ndims, d, nP, p, nF, ff, in, jn, kn, In, Jn, Kn;

  int     maxSTEPS, Nadvice = 0;
  real    MAT[16];


  real    *O;
  char    STR[100];
  int     argN;
  
  real    *initJAC, *OJAC, *DET;
  
  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );
  
  if( mxIsChar( prhs[3] ) ){
    mxGetString( prhs[3], STR, 100 );
    if( myStrcmpi(STR,"2d") ) { myErrMsgTxt("The only valid keyword is 2D."); }
    K = 1;
    NSD = 2;
  } else if( myIsEmpty( prhs[3] ) ){
    K = 1;
    NSD = 2;
  } else {
    K  = mySize( prhs[0] , 2 );
    NSD = 3;
  }
  IJ = I*J;
  IJK = IJ*K;

  if( myNumel( prhs[0] ) != IJK*NSD ){
    myErrMsgTxt("size(V) has to be [numel(Gx) numel(Gy) numel(Gz) 3]   or  [numel(Gx) numel(Gy) 2] in 2D case.");
  }
  
  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(Gy) Coordinates do not coincide with size(V,2)."); }
  if( NSD == 3 && myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(Gz) Coordinates do not coincide with size(V,3)."); }

  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("Gx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("Gy  Coordinates are not sorted."); }

  if( NSD == 3 ){
    Z = myGetPr(prhs[3]);
    if( !checkIsSorted(Z,K) ){ myErrMsgTxt("Gz  Coordinates are not sorted."); }
  } else {
    Z = mxMalloc( 1*sizeof( real ) );
    Z[0] = 0;
  }

  
  if( mxIsCell( prhs[4] ) ){
    
    XYZ = NULL;
    
    if( myNumel( prhs[4] ) != NSD ){
      myErrMsgTxt("In cell provided case, %d  vectors are expected.", NSD);
    }
    
    Xn = myGetPr( mxGetCell( prhs[4] , 0 ) );
    In = myNumel( mxGetCell( prhs[4] , 0 ) );

    Yn = myGetPr( mxGetCell( prhs[4] , 1 ) );
    Jn = myNumel( mxGetCell( prhs[4] , 1 ) );

    if( NSD == 3 ){
      Zn = myGetPr( mxGetCell( prhs[4] , 2 ) );
      Kn = myNumel( mxGetCell( prhs[4] , 2 ) );
    } else {
      Zn  = mxMalloc( 1*sizeof( real ) );
      Zn[0] = 0;
      Kn = 1;
    }

    nP = In * Jn * Kn;
    
  } else {
    
    Xn = NULL;
    
    nP  = myNumel( prhs[4] );

    if( ( nP%NSD ) ){ myErrMsgTxt("Number of XYZ have to be multiple of %d.", NSD); }
    nP = nP/NSD;

    In = nP;
    Jn = 1;
    Kn = 1;
    
    XYZ  = myGetPr( prhs[4] );

  }
  

  VX = myGetPr( prhs[0] );
  VY = VX + IJK;
  if( NSD == 3 ){
    VZ = VY + IJK;
  } else {
    VZ = VY;      /* para evitar problemas!!!! */
  }

  
  /*Parsing arguments*/
  /*Defaults*/
  boundary_size[0]    = ( I > 1 ) ? ( X[ 1 ] - X[ 0 ])/2 : 0.5;
  boundary_size[1]    = ( I > 1 ) ? ( X[I-1] - X[I-2])/2 : 0.5;
  boundary_size[2]    = ( J > 1 ) ? ( Y[ 1 ] - Y[ 0 ])/2 : 0.5;
  boundary_size[3]    = ( J > 1 ) ? ( Y[J-1] - Y[J-2])/2 : 0.5;
  boundary_size[4]    = ( K > 1 ) ? ( Z[ 1 ] - Z[ 0 ])/2 : 0.5;
  boundary_size[5]    = ( K > 1 ) ? ( Z[K-1] - Z[K-2])/2 : 0.5;
  
  integration_method  = integration_method_def;
  jacobian_method     = jacobian_method_def;
  boundary_mode       = boundary_mode_def;
  relTOL              = relTOL_def;
  minStepsPerVoxel    = minStepsPerVoxel_def;
  maxSTEPS            = maxSTEPS_def;
  F                   = NULL;
  initJAC             = NULL;
  MASK                = NULL;
  MAT(3,3)            = 0;
  VERBOSE             = 0;

  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );

    if( ! myStrcmpi(STR,"verbose")                                 ) { VERBOSE = 1;                   argN++; continue; }


    if( ! myStrcmpi(STR,"s_euler")                                 ) { integration_method = S_EULER;  argN++; continue; }
    if( ! myStrcmpi(STR,"euler")                                   ) { integration_method = EULER;    argN++; continue; }
    if( ! myStrcmpi(STR,"rk23")                                    ) { integration_method = RK23;     argN++; continue; }
    if( ! myStrcmpi(STR,"rk45")                                    ) { integration_method = RK45;     argN++; continue; }


    if( !myStrcmpi(STR,"jac_times") || !myStrcmpi(STR,"jactimes")  ) { jacobian_method = TIMES;       argN++; continue; }
    if( !myStrcmpi(STR,"jac_plus")  || !myStrcmpi(STR,"jacplus")   ) { jacobian_method = PLUS;        argN++; continue; }


    if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"zero")     ) { boundary_mode = VALUE;    argN++; continue; }
    if( ! myStrcmpi(STR,"closest")                                 ) { boundary_mode = CLOSEST;  argN++; continue; }

    if( ! myStrcmpi(STR,"circular") || ! myStrcmpi(STR,"periodic") || 
        ! myStrcmpi(STR,"circ")     || ! myStrcmpi(STR,"per")      ) {
      boundary_mode = CIRCULAR; argN++;  readBOUNDS;  continue;
    }

    if( ! myStrcmpi( STR,"decay_to_zero") || 
        ! myStrcmpi( STR,"tozero")  || ! myStrcmpi( STR,"decay")   ) {
          boundary_mode = DECAY_TO_ZERO; argN++;  readBOUNDS;  continue;
    }

    if( ! myStrcmpi( STR,"max") || ! myStrcmpi( STR,"maxsteps") ){
      argN++;
      if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
        if( !myIsEmpty( prhs[argN] ) ){ maxSTEPS = myGetValue(prhs[argN]); }
        argN++; continue;
      } else {
        myErrMsgTxt("After maxSTEPS a number is expected.");
      }
      continue;
    }

    if( ! myStrcmpi( STR,"min") || ! myStrcmpi( STR,"minStepsPerVoxel") ){
      argN++;
      if( nrhs > argN  && !mxIsChar(prhs[argN]) ){
        if( !myIsEmpty( prhs[argN] ) ){ minStepsPerVoxel = myGetValue(prhs[argN]); }
        argN++; continue;
      } else {
        myErrMsgTxt("After minStepsPerVoxel a number is expected.");
      }
      continue;
    }


    if( ! myStrcmpi( STR,"reltol") ){
      argN++;
      if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
        if( !myIsEmpty( prhs[argN] ) ){ relTOL = myGetValue(prhs[argN]); }
        argN++; continue;
      } else {
        myErrMsgTxt("After relTOL a number is expected.");
      }
      continue;
    }


    if( !myStrcmpi( STR,"f")  ||  !myStrcmpi( STR,"flow") ){
      argN++;
      if( nrhs > argN  &&  !mxIsChar(prhs[argN])  &&  !myIsEmpty( prhs[argN] ) ){
        if( myNumel( prhs[argN] ) == 1 ){
          nF = 2;
          F = mxMalloc( 2*sizeof( real ) );
          F[0] = 0;
          F[1] = myGetValue( prhs[argN] );
        } else {
          nF = myNumel( prhs[argN] );
          F  = myGetPr( prhs[argN] );
        }
        argN++;
      } else {
        myErrMsgTxt("After FLOW a number or a vector is expected.");
      }
      continue;
    }


    if( ! myStrcmpi( STR,"matrix") || ! myStrcmpi( STR,"mat") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) &&  mySize( prhs[argN] , 0 ) == 4  &&  mySize( prhs[argN] , 1 ) == 4  && myNDims( prhs[argN] ) == 2 ){
        if( !myIsEmpty( prhs[argN] ) ){
          MAT(0,0) = myGetValueIth( prhs[argN] ,      0 );
          MAT(1,0) = myGetValueIth( prhs[argN] ,      1 );
          MAT(2,0) = myGetValueIth( prhs[argN] ,      2 );
          MAT(3,0) = myGetValueIth( prhs[argN] ,      3 );
          MAT(0,1) = myGetValueIth( prhs[argN] ,  4 + 0 );
          MAT(1,1) = myGetValueIth( prhs[argN] ,  4 + 1 );
          MAT(2,1) = myGetValueIth( prhs[argN] ,  4 + 2 );
          MAT(3,1) = myGetValueIth( prhs[argN] ,  4 + 3 );
          MAT(0,2) = myGetValueIth( prhs[argN] ,  8 + 0 );
          MAT(1,2) = myGetValueIth( prhs[argN] ,  8 + 1 );
          MAT(2,2) = myGetValueIth( prhs[argN] ,  8 + 2 );
          MAT(3,2) = myGetValueIth( prhs[argN] ,  8 + 3 );
          MAT(0,3) = myGetValueIth( prhs[argN] , 12 + 0 );
          MAT(1,3) = myGetValueIth( prhs[argN] , 12 + 1 );
          MAT(2,3) = myGetValueIth( prhs[argN] , 12 + 2 );
          MAT(3,3) = myGetValueIth( prhs[argN] , 12 + 3 );

          if( MAT(3,0) != 0 || MAT(3,1) != 0 || MAT(3,2) != 0 || MAT(3,3) != 1 ){
            myErrMsgTxt("The matrix must be an homogeneous 4x4 matrix.");
          }

        }
        argN++; continue;
      }
      myErrMsgTxt("After the word MATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }


    if( ! myStrcmpi( STR,"mask") ){
      argN++;
      if( nrhs > argN && mxIsLogical( prhs[argN] ) &&  ( myNumel( prhs[argN] ) == nP  ||  myIsEmpty( prhs[argN] ) )  ){
        if( !myIsEmpty( prhs[argN] ) ){
          MASK = (unsigned char *) mxGetData( prhs[argN] );
        }
        argN++; continue;
      }
      myErrMsgTxt("After the word MASK a nP logical is expected.");
    }


    if( ! myStrcmpi( STR,"initJACs") ){
      argN++;
      if( nrhs > argN && ( myNumel( prhs[argN] ) == nP*NSD*NSD  ||  myIsEmpty( prhs[argN] ) )  ){
        if( !myIsEmpty( prhs[argN] ) ){

          ndims = myNDims( prhs[argN] );

          if( mySize( prhs[argN] , ndims-2 ) != NSD || mySize( prhs[argN] , ndims-1 ) != NSD ){
            myErrMsgTxt("After the word initJACs  nPxNSDxNSD array is expected.");
          }

          initJAC = myGetPr( prhs[argN] );

        }
        argN++; continue;
      }
      myErrMsgTxt("After the word initJACs  nPxNSDxNSD array is expected.");
    }

    mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
  }
  /*END Parsing arguments*/
  
  if( F == NULL ){
    nF = 2;
    F = mxMalloc( 2*sizeof( real ) );
    F[0] = 0;
    F[1] = 1;
  }
  if( !checkIsSorted(F,nF) ){ myErrMsgTxt("FLOW times are not sorted."); }

  
  if( MAT(3,3) != 0 ){
    if( NSD == 3 ){
      if( MAT(3,0) != 0 || MAT(3,1) != 0 || MAT(3,2) != 0 || MAT(3,3) != 1 ){
        myErrMsgTxt("Invalid matrices");
      }
    } else if( NSD == 2 ){
      if( MAT(3,0) != 0 || MAT(3,1) != 0 || MAT(3,2) != 0 || MAT(3,3) != 1 || 
          MAT(2,0) != 0 || MAT(2,1) != 0 || MAT(2,2) != 1 || MAT(2,3) != 0 ||
          MAT(0,2) != 0 || MAT(1,2) != 0 ){
        myErrMsgTxt("Invalid 2D matrices.");
      }
    }
  }
  
  
  
  /*Checking sizes*/
  if( Xn == NULL ){

    ndims = myNDims( prhs[4] );
    for( d=0 ; d<ndims ; d++ ){
      Odims[d] = mySize( prhs[4] , d );
    }
    if( Odims[ndims-1] != NSD ){
      myErrMsgTxt("algun error... la ultima dimension no tiene size NSD .");
    }
    Odims[ndims-1] = nF-1;
    Odims[ndims  ] = NSD;
    Odims[ndims+1] = NSD;
    
  } else {
    
    ndims = 4;
    Odims[0] = In;
    Odims[1] = Jn;
    Odims[2] = Kn;
    Odims[3] = nF-1;
    Odims[4] = NSD;
    Odims[5] = NSD;
    
  }
  /*END Checking sizes*/
  
  /*Creating output*/
  plhs[0] = mxCreateNumericArray( ndims+1 , Odims , mxREAL_CLASS , mxREAL );
  O = (real *) mxGetData( plhs[0] );

  if( nlhs > 1 ){
    plhs[1] = mxCreateNumericArray( ndims+2 , Odims , mxREAL_CLASS , mxREAL );
    OJAC = (real *) mxGetData( plhs[1] );
  }

  if( nlhs > 2 ){
    plhs[2] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
    DET = (real *) mxGetData( plhs[2] );
  }
  /*END Creating output*/


  switch( boundary_mode ){
    case VALUE:
      OX = ( I>1 ) ? X[ 0 ] : X[ 0 ] - 0.5;   LX = ( I>1 ) ? X[I-1] : X[I-1] + 0.5;
      OY = ( J>1 ) ? Y[ 0 ] : Y[ 0 ] - 0.5;   LY = ( J>1 ) ? Y[J-1] : Y[J-1] + 0.5;
      OZ = ( K>1 ) ? Z[ 0 ] : Z[ 0 ] - 0.5;   LZ = ( K>1 ) ? Z[K-1] : Z[K-1] + 0.5;
      break;
    case CLOSEST:
      OX = X[ 0 ];  LX = X[I-1];
      OY = Y[ 0 ];  LY = Y[J-1];
      OZ = Z[ 0 ];  LZ = Z[K-1];
      break;
    case DECAY_TO_ZERO:
      OX = X[ 0 ];  LX = X[I-1];
      OY = Y[ 0 ];  LY = Y[J-1];
      OZ = Z[ 0 ];  LZ = Z[K-1];
      break;
    case CIRCULAR:
      OX = X[ 0 ] - boundary_size[0];     LX = X[I-1] + boundary_size[1];
      OY = Y[ 0 ] - boundary_size[2];     LY = Y[J-1] + boundary_size[3];
      OZ = Z[ 0 ] - boundary_size[4];     LZ = Z[K-1] + boundary_size[5];
      
      break;
  }

  
//   if( integration_method == S_EULER  ||  integration_method == RK23  || integration_method == RK45  ){
    minD = 1.0;
    if( I > 1 ){ minD = MIN( minD , (LX-OX)/I ); }
    if( J > 1 ){ minD = MIN( minD , (LY-OY)/J ); }
    if( K > 1 ){ minD = MIN( minD , (LZ-OZ)/K ); }
    tol  = minD*relTOL;
//   }
  
  hmin     = ( F[nF-1] - F[0] )/maxSTEPS;
  minStepsPerVoxel = 1.0/minStepsPerVoxel;

  // reset the CELL
  C.ii = C.jj =  C.kk = -1;
  C.X0 = C.Y0 = C.Z0 =  0;
  C.X1 = C.Y1 = C.Z1 = -1;
  C.VX0000 = C.VX1000 = C.VX0100 = C.VX1100 = C.VX0010 = C.VX1010 = C.VX0110 = C.VX1110 = 0;
  C.VY0000 = C.VY1000 = C.VY0100 = C.VY1100 = C.VY0010 = C.VY1010 = C.VY0110 = C.VY1110 = 0;
  C.VZ0000 = C.VZ1000 = C.VZ0100 = C.VZ1100 = C.VZ0010 = C.VZ1010 = C.VZ0110 = C.VZ1110 = 0;
  C.DX = C.DY = C.DZ = 1; 
  C.D = C.iD = 1;
  C.isOutside = 1;
  
  
  
  
  XYZT[0] = XYZT[1] = XYZT[2] = XYZT[3] = 0;
  v[0] = v[1] = v[2] = 0;

  p = -1;
  for( kn=0 ; kn < Kn ; kn++ ){
  for( jn=0 ; jn < Jn ; jn++ ){
  for( in=0 ; in < In ; in++ ){
    if( utIsInterruptPending() ){
      mexPrintf("USER INTERRUP in  EvolvePointsOn3DGrid  !!!\n");
      myErrMsgTxt("USER INTERRUP in  EvolvePointsOn3DGrid  !!!");
    }
    p++;
    if( VERBOSE ){
      prct = ( (double)(p+1) )/((double)nP)*100;
      if( prct >= last_prct ){
        if( last_prct != 0 ){
          mexPrintf("%3d%% done  ( %10d   of   %10d  points )\n" , (int) prct , p+1, nP ); myFlush();
        }
        last_prct += 10;
      }
    }

    
    if( MASK != NULL  &&  MASK[p] != 1 ){
      XYZT[0] = NAN;
    } else if( Xn == NULL ){
                      XYZT[0] = XYZ[ p        ];
                      XYZT[1] = XYZ[ p +   nP ];
      if( NSD == 3 ){ XYZT[2] = XYZ[ p + 2*nP ]; }
    } else {
                      XYZT[0] = Xn[ in ];
                      XYZT[1] = Yn[ jn ];
      if( NSD == 3 ){ XYZT[2] = Zn[ kn ]; }
    }
    
//     mexPrintf("XYZ:  %g,%g,%g\n" , XYZT[0],XYZT[1],XYZT[2]);
    
    if( myISNAN( XYZT[0] ) || myISNAN( XYZT[1] ) || myISNAN( XYZT[2] ) ){
      for( ff = 1 ; ff < nF ; ff++ ){

                        O[ p + ( (ff-1)            )*nP ] = NAN;
                        O[ p + ( (ff-1) +   (nF-1) )*nP ] = NAN;
        if( NSD == 3 ){ O[ p + ( (ff-1) + 2*(nF-1) )*nP ] = NAN; }
        
        if( nlhs > 1 ){
            OJAC[ p + ( (ff-1)            )*nP ] = NAN;
            OJAC[ p + ( (ff-1) +   (nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 2*(nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 3*(nF-1) )*nP ] = NAN;

          if( NSD == 3 ){
            OJAC[ p + ( (ff-1) + 4*(nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 5*(nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 6*(nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 7*(nF-1) )*nP ] = NAN;
            OJAC[ p + ( (ff-1) + 8*(nF-1) )*nP ] = NAN;
          }
        }
                        
        if( nlhs > 2 ){
          DET[ p + (ff-1)*nP ] = NAN;
        }
        
      }
//       mexPrintf("continue 1\n");
      continue;
    }
    

    
    if( MAT(3,3) == 1 ) {
      xxnr= XYZT[0];
      yynr= XYZT[1];
      zznr= XYZT[2];

      XYZT[0] = xxnr*MAT(0,0) + yynr*MAT(0,1) + zznr*MAT(0,2) + MAT(0,3);
      XYZT[1] = xxnr*MAT(1,0) + yynr*MAT(1,1) + zznr*MAT(1,2) + MAT(1,3);
      XYZT[2] = xxnr*MAT(2,0) + yynr*MAT(2,1) + zznr*MAT(2,2) + MAT(2,3);
    }
    
//     mexPrintf("XYZ after rotation:  %g,%g,%g\n" , XYZT[0],XYZT[1],XYZT[2]);
    
    
    XYZT[3] = F[0];
    
    if( nlhs > 1 ){
      if( initJAC != NULL ){
        if( NSD == 3 ){
          JAC[0] = initJAC[ p        ];
          JAC[1] = initJAC[ p +   nP ];
          JAC[2] = initJAC[ p + 2*nP ];
          JAC[3] = initJAC[ p + 3*nP ];
          JAC[4] = initJAC[ p + 4*nP ];
          JAC[5] = initJAC[ p + 5*nP ];
          JAC[6] = initJAC[ p + 6*nP ];
          JAC[7] = initJAC[ p + 7*nP ];
          JAC[8] = initJAC[ p + 8*nP ];
        } else {
          JAC[0] = initJAC[ p        ];
          JAC[1] = initJAC[ p +   nP ];
          JAC[3] = initJAC[ p + 2*nP ];
          JAC[4] = initJAC[ p + 3*nP ];
          JAC[2] = JAC[5] = JAC[6] = JAC[7] = 0.0;
          JAC[8] = 1.0;
        }
      } else {
        JAC[0] = JAC[4] = JAC[8] = 1.0;
        JAC[1] = JAC[2] = JAC[3] = JAC[5] = JAC[6] = JAC[7] = 0.0;
      }
    }
    
    for( ff = 1 ; ff < nF ; ff++ ){

      switch( integration_method ){
        case EULER:

          h    = hmin;
//     mexPrintf("h:  %g\n" , h );
          while( XYZT[3] < F[ff] ){
            getV( XYZT , v );
//     mexPrintf("v:  %g,%g,%g\n" , v[0],v[1],v[2]);

            if( !v[0] && !v[1] && !v[2] ){ break; }
            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }

            EVOLVE_JAC

            XYZT[0] += h*v[0];
            XYZT[1] += h*v[1];
            XYZT[2] += h*v[2];
            XYZT[3] += h;
            
//     mexPrintf("XYZT:  %g,%g,%g,%g\n" , XYZT[0],XYZT[1],XYZT[2],XYZT[3]);
            
          }
          break;

        case S_EULER:

          h = 1;
          while( XYZT[3] < F[ff] ){
            getV( XYZT , v );

            if( !v[0] && !v[1] && !v[2] ){ break; }
            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }

            maxV = MAX3( ABS( v[0] ) , ABS( v[1] ) , ABS( v[2] ) );
            hmax = MIN( MAX( minD/maxV*minStepsPerVoxel , hmin ) , minStepsPerVoxel );
            h    = MIN( h , hmax );

            EVOLVE_JAC

            XYZT[0] += h*v[0];
            XYZT[1] += h*v[1];
            XYZT[2] += h*v[2];
            XYZT[3] += h;
          }
          break;

        case RK45:

          getV( XYZT , v );
          h = 1;
          while( XYZT[3] < F[ff] ){
            if( !v[0] && !v[1] && !v[2] ){ break; }

            maxV = MAX3( ABS( v[0] ) , ABS( v[1] ) , ABS( v[2] ) );
            hmax = MIN( MAX( minD/maxV*minStepsPerVoxel , hmin ) , minStepsPerVoxel );
            h    = MIN( h , hmax );

            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }

            times( k1 , v , h );

            xk[0] = XYZT[0] + A45_21*k1[0];
            xk[1] = XYZT[1] + A45_21*k1[1];
            xk[2] = XYZT[2] + A45_21*k1[2];
            getV( xk , k2 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );


            xk[0] = XYZT[0] + A45_31*k1[0] + A45_32*k2[0];
            xk[1] = XYZT[1] + A45_31*k1[1] + A45_32*k2[1];
            xk[2] = XYZT[2] + A45_31*k1[2] + A45_32*k2[2];
            getV( xk , k3 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k3 , k3 , h );

            xk[0] = XYZT[0] + A45_41*k1[0] + A45_42*k2[0] + A45_43*k3[0];
            xk[1] = XYZT[1] + A45_41*k1[1] + A45_42*k2[1] + A45_43*k3[1];
            xk[2] = XYZT[2] + A45_41*k1[2] + A45_42*k2[2] + A45_43*k3[2];
            getV( xk , k4 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k4 , k4 , h );

            xk[0] = XYZT[0] + A45_51*k1[0] + A45_52*k2[0] + A45_53*k3[0] + A45_54*k4[0];
            xk[1] = XYZT[1] + A45_51*k1[1] + A45_52*k2[1] + A45_53*k3[1] + A45_54*k4[1];
            xk[2] = XYZT[2] + A45_51*k1[2] + A45_52*k2[2] + A45_53*k3[2] + A45_54*k4[2];
            getV( xk , k5 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k5 , k5 , h );

            xk[0] = XYZT[0] + A45_61*k1[0] + A45_62*k2[0] + A45_63*k3[0] + A45_64*k4[0] + A45_65*k5[0];
            xk[1] = XYZT[1] + A45_61*k1[1] + A45_62*k2[1] + A45_63*k3[1] + A45_64*k4[1] + A45_65*k5[1];
            xk[2] = XYZT[2] + A45_61*k1[2] + A45_62*k2[2] + A45_63*k3[2] + A45_64*k4[2] + A45_65*k5[2];
            getV( xk , k6 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k6 , k6 , h );

            memcpy( vp , v , 3*sizeof( real ) );

            xk[0] = XYZT[0] + B45_1*k1[0] + B45_3*k3[0] + B45_4*k4[0] + B45_5*k5[0] + B45_6*k6[0];
            xk[1] = XYZT[1] + B45_1*k1[1] + B45_3*k3[1] + B45_4*k4[1] + B45_5*k5[1] + B45_6*k6[1];
            xk[2] = XYZT[2] + B45_1*k1[2] + B45_3*k3[2] + B45_4*k4[2] + B45_5*k5[2] + B45_6*k6[2];
            xk[3] = XYZT[3] + h;
            getV( xk , v );
            if( C.isOutside  &&  h > hmin ){ 
              memcpy( v , vp , 3*sizeof( real ) );
              h = hmin; 
              continue; 
            }
            times( k7 , v , h );

            err[0] = D45_1*k1[0] + D45_3*k3[0] + D45_4*k4[0] + D45_5*k5[0] + D45_6*k6[0] + D45_7*k7[0];
            err[1] = D45_1*k1[1] + D45_3*k3[1] + D45_4*k4[1] + D45_5*k5[1] + D45_6*k6[1] + D45_7*k7[1];
            err[2] = D45_1*k1[2] + D45_3*k3[2] + D45_4*k4[2] + D45_5*k5[2] + D45_6*k6[2] + D45_7*k7[2];
            err[0] = MAX3( ABS( err[0] ) , ABS( err[1] ) , ABS( err[2] ) ) ;

            if( err[0] > tol ){
              if( h > hmin ){
                h = h/( 1.1*sqrt( sqrt( err[0]/tol ) ) );
                h = MAX( MIN( h , hmax ) , hmin );
                memcpy( v , vp , 3*sizeof( real ) );
                continue;
              } else {
                Nadvice = MAX( Nadvice , ceil( pow( err[0]/tol , 1.0/5.0 ) / h ) );
              }
            }


            EVOLVE_JAC

            memcpy( XYZT , xk , 4*sizeof( real ) );

            h = h/MAX( 0.1 , 1.1*sqrt( sqrt( err[0]/tol ) ) );
            h = MAX( MIN( h , hmax ) , hmin );

          }
          break;


        case RK23:

          getV( XYZT , v );
          h = 1;
          while( XYZT[3] < F[ff] ){
            if( !v[0] && !v[1] && !v[2] ){ break; }

            maxV = MAX3( ABS( v[0] ) , ABS( v[1] ) , ABS( v[2] ) );
            hmax = MIN( MAX( minD/maxV*minStepsPerVoxel , hmin ) , minStepsPerVoxel );
            h    = MIN( h , hmax );

            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }

            times( k1 , v , h );

            xk[0] = XYZT[0] + A23_21*k1[0];
            xk[1] = XYZT[1] + A23_21*k1[1];
            xk[2] = XYZT[2] + A23_21*k1[2];
            getV( xk , k2 );  if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );

            xk[0] = XYZT[0] + A23_32*k2[0];
            xk[1] = XYZT[1] + A23_32*k2[1];
            xk[2] = XYZT[2] + A23_32*k2[2];
            getV( xk , k3 );  if( C.isOutside  &&  h > hmin ){ h = hmin; continue; } 
            times( k3 , k3 , h );

            memcpy( vp , v , 3*sizeof( real ) );

            xk[0] = XYZT[0] + B23_1*k1[0] + B23_2*k2[0] + B23_3*k3[0];
            xk[1] = XYZT[1] + B23_1*k1[1] + B23_2*k2[1] + B23_3*k3[1];
            xk[2] = XYZT[2] + B23_1*k1[2] + B23_2*k2[2] + B23_3*k3[2];
            xk[3] = XYZT[3] + h;
            getV( xk , v );
            if( C.isOutside  &&  h > hmin ){ 
              memcpy( v , vp , 3*sizeof( real ) );
              h = hmin; 
              continue; 
            }
            times( k4 , v , h );

            err[0] = D23_1*k1[0] + D23_2*k1[0] + D23_3*k3[0] + D23_4*k4[0];
            err[1] = D23_1*k1[1] + D23_2*k1[1] + D23_3*k3[1] + D23_4*k4[1];
            err[2] = D23_1*k1[2] + D23_2*k1[2] + D23_3*k3[2] + D23_4*k4[2];
            err[0] = MAX3( ABS( err[0] ) , ABS( err[1] ) , ABS( err[2] ) ) ;

            if( err[0] > tol ){
              if( h > hmin ){
                h = h/( 1.1*sqrt( err[0]/tol ) );
                h = MAX( MIN( h , hmax ) , hmin );
                memcpy( v , vp , 3*sizeof( real ) );
                continue;
              } else {
                Nadvice = MAX( Nadvice , ceil( pow( err[0]/tol , 1.0/3.0 ) / h ) );
              }
            }

            EVOLVE_JAC

            memcpy( XYZT , xk , 4*sizeof( real ) );

            h = h/MAX( 0.1 , 1.1*sqrt( err[0]/tol ) );
            h = MAX( MIN( h , hmax ) , hmin );

          }
          break;

      }

      if( NSD == 3 ){

        O[ p + ( (ff-1)            )*nP ] = XYZT[0];
        O[ p + ( (ff-1) +   (nF-1) )*nP ] = XYZT[1];
        O[ p + ( (ff-1) + 2*(nF-1) )*nP ] = XYZT[2];

        if( nlhs > 1 ){
          OJAC[ p + ( (ff-1)            )*nP ] = JAC[0];
          OJAC[ p + ( (ff-1) +   (nF-1) )*nP ] = JAC[1];
          OJAC[ p + ( (ff-1) + 2*(nF-1) )*nP ] = JAC[2];
          OJAC[ p + ( (ff-1) + 3*(nF-1) )*nP ] = JAC[3];
          OJAC[ p + ( (ff-1) + 4*(nF-1) )*nP ] = JAC[4];
          OJAC[ p + ( (ff-1) + 5*(nF-1) )*nP ] = JAC[5];
          OJAC[ p + ( (ff-1) + 6*(nF-1) )*nP ] = JAC[6];
          OJAC[ p + ( (ff-1) + 7*(nF-1) )*nP ] = JAC[7];
          OJAC[ p + ( (ff-1) + 8*(nF-1) )*nP ] = JAC[8];
        }

        if( nlhs > 2 ){
          DET[ p + (ff-1)*nP ] =   JAC[0]*( JAC[4]*JAC[8] - JAC[5]*JAC[7] )
                                 - JAC[1]*( JAC[3]*JAC[8] - JAC[5]*JAC[6] )
                                 + JAC[2]*( JAC[3]*JAC[7] - JAC[4]*JAC[6] );
        }


      } else {

        O[ p + ( (ff-1)            )*nP ] = XYZT[0];
        O[ p + ( (ff-1) +   (nF-1) )*nP ] = XYZT[1];

        if( nlhs > 1 ){
          OJAC[ p + ( (ff-1)            )*nP ] = JAC[0];
          OJAC[ p + ( (ff-1) +   (nF-1) )*nP ] = JAC[1];
          OJAC[ p + ( (ff-1) + 2*(nF-1) )*nP ] = JAC[3];
          OJAC[ p + ( (ff-1) + 3*(nF-1) )*nP ] = JAC[4];
        }

        if( nlhs > 2 ){
          DET[ p + (ff-1)*nP ] = JAC[0]*JAC[4] - JAC[3]*JAC[1];
        }
        
      }

    }
  }}}

  if( Nadvice > 0 ){
    sprintf( STR , "Adviced maxSTEPS: %d" , Nadvice );
    mexWarnMsgIdAndTxt( "EXP:V:IncreaseMaxSteps" , STR );
  }


  EXIT: myFreeALLOCATES();
}


_inline void setCELL( real x , real y , real z ){
  int   ii0 , ii1, jj0 , jj1, kk0, kk1;
  int   isSingularCase;
  int   n;
  
  C.isOutside = 0;
  isSingularCase = 0;
  switch( boundary_mode ){
    case VALUE:
      if( x < OX || x > LX  ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( y < OY || y > LY  ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( z < OZ || z > LZ  ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      
      if( I > 1 ){
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0];
        ii1  = ii0 + 1;     C.X1 = X[ii1];
      } else {
        ii0 = 0; C.X0 = OX;
        ii1 = 0; C.X1 = LX;
      }

      if( J > 1 ){
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0];
        jj1  = jj0 + 1;     C.Y1 = Y[jj1];
      } else {
        jj0 = 0; C.Y0 = OY;
        jj1 = 0; C.Y1 = LY;
      }

      if( K > 1 ){
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0];
        kk1  = kk0 + 1;     C.Z1 = Z[kk1];
      } else {
        kk0 = 0; C.Z0 = OZ;
        kk1 = 0; C.Z1 = LZ;
      }
      
      break;

    case CLOSEST:
      if( I == 1 ){
        ii0 = 0;   C.X0 = x - 0.5;
        ii1 = 0;   C.X1 = x + 0.5;
      } else if( x < OX ){
        C.ii = 0;
        ii0  = 0;  C.X0 = x - 0.5;
        ii1  = 0;  C.X1 = x + 0.5;
      } else if( x > LX ){
        C.ii = MAX( 0 , I-2 );
        ii0  = I-1;  C.X0 = x - 0.5;;
        ii1  = I-1;  C.X1 = x + 0.5;;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0];
        ii1  = ii0 + 1;     C.X1 = X[ii1];
      }
      
      if( J == 1 ){
        jj0 = 0;   C.Y0 = y - 0.5;
        jj1 = 0;   C.Y1 = y + 0.5;
      } else if( y < OY ){
        C.jj = 0;
        jj0  = 0;  C.Y0 = y - 0.5;
        jj1  = 0;  C.Y1 = y + 0.5;
      } else if( y > LY ){
        C.jj = MAX( 0 , J-2 );
        jj0  = J-1;  C.Y0 = y - 0.5;
        jj1  = J-1;  C.Y1 = y + 0.5;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0];
        jj1  = jj0 + 1;     C.Y1 = Y[jj1];
      }

      if( K == 1 ){
        kk0 = 0;   C.Z0 = z - 0.5;
        kk1 = 0;   C.Z1 = z + 0.5;
      } else if( z < OZ ){
        C.kk = 0;
        kk0  = 0;  C.Z0 = z - 0.5;
        kk1  = 0;  C.Z1 = z + 0.5;
      } else if( z > LZ ){
        C.kk = MAX( 0 , K-2 );
        kk0  = K-1;  C.Z0 = z - 0.5;
        kk1  = K-1;  C.Z1 = z + 0.5;
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0];
        kk1  = kk0 + 1;     C.Z1 = Z[kk1];
      }      

      break;


    case DECAY_TO_ZERO:
      
      if( x < OX-boundary_size[0] || x > LX+boundary_size[1] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( y < OY-boundary_size[2] || y > LY+boundary_size[3] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( z < OZ-boundary_size[3] || z > LZ+boundary_size[5] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      
      if( x < OX ){
        isSingularCase = 1;
        C.ii = 0;
        ii0  = -1;        C.X0 = OX-boundary_size[0];
        ii1  = 0;         C.X1 = OX;
      } else if( x > LX ) {
        isSingularCase = 1;
        C.ii = MAX( 0 , I-2 );
        ii0  = I-1;       C.X0 = LX;
        ii1  = -1;        C.X1 = LX+boundary_size[1];
      } else if( I == 1 ){
        C.ii = 0;
        ii0  = 0;         C.X0 = -0.5;
        ii1  = 0;         C.X1 =  0.5;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;      C.X0 = X[ii0];
        ii1  = ii0 + 1;   C.X1 = X[ii1];
      }
      
      if( y < OY ){
        isSingularCase = 1;
        C.jj = 0;
        jj0  = -1;        C.Y0 = OY-boundary_size[2];
        jj1  = 0;         C.Y1 = OY;
      } else if( y > LY ) {
        isSingularCase = 1;
        C.jj = MAX( 0 , J-2 );
        jj0  = J-1;       C.Y0 = LY;
        jj1  = -1;        C.Y1 = LY+boundary_size[3];
      } else if( J == 1 ){
        C.jj = 0;
        jj0  = 0;         C.Y0 = -0.5;
        jj1  = 0;         C.Y1 =  0.5;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;      C.Y0 = Y[jj0];
        jj1  = jj0 + 1;   C.Y1 = Y[jj1];
      }

      if( z < OZ ){
        isSingularCase = 1;
        C.kk = 0;
        kk0  = -1;        C.Z0 = OZ-boundary_size[4];
        kk1  = 0;         C.Z1 = OZ;
      } else if( z > LZ ) {
        isSingularCase = 1;
        C.kk = MAX( 0 , K-2 );
        kk0  = K-1;       C.Z0 = LZ;
        kk1  = -1;        C.Z1 = LZ+boundary_size[5];
      } else if( K == 1 ){
        C.kk = 0;
        kk0  = 0;         C.Z0 = -0.5;
        kk1  = 0;         C.Z1 =  0.5;
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;      C.Z0 = Z[kk0];
        kk1  = kk0 + 1;   C.Z1 = Z[kk1];
      }

      break;

    case CIRCULAR:

      if( x < OX || x > LX ){
        n = floor( ( x - OX )/( LX-OX ) );
        x = x - n*( LX-OX );
      } else { n = 0; }
      if( I == 1 ) {
        ii0  = 0;           C.X0 = x - 0.5 + n*( LX-OX );
        ii1  = 0;           C.X1 = x + 0.5 + n*( LX-OX );
      } else if( x < X[0] ){
        ii0  = I-1;         C.X0 = X[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LX-OX );
        ii1  = 0;           C.X1 = X[ 0 ] + n*( LX-OX );
      } else if( x > X[I-1] ){
        ii0  = I-1;         C.X0 = X[ii0] + n*( LX-OX );
        ii1  = 0;           C.X1 = X[I-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LX-OX );
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0] + n*( LX-OX );
        ii1  = ii0 + 1;     C.X1 = X[ii1] + n*( LX-OX );
      }
      

      if( y < OY || y > LY ){
        n = floor( ( y - OY )/( LY-OY ) );
        y = y - n*( LY-OY );
      } else { n = 0; }
      if( J == 1 ) {
        jj0  = 0;           C.Y0 = y - 0.5 + n*( LY-OY );
        jj1  = 0;           C.Y1 = y + 0.5 + n*( LY-OY );
      } else if( y < Y[0] ){
        jj0  = J-1;         C.Y0 = Y[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LY-OY );
        jj1  = 0;           C.Y1 = Y[ 0 ] + n*( LY-OY );
      } else if( y > Y[J-1] ){
        jj0  = J-1;         C.Y0 = Y[jj0] + n*( LY-OY );
        jj1  = 0;           C.Y1 = Y[J-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LY-OY );
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0] + n*( LY-OY );
        jj1  = jj0 + 1;     C.Y1 = Y[jj1] + n*( LY-OY );
      }
      

      if( z < OZ || z > LZ ){
        n = floor( ( z - OZ )/( LZ-OZ ) );
        z = z - n*( LZ-OZ );
      } else { n = 0; }
      if( K == 1 ) {
        kk0  = 0;           C.Z0 = z - 0.5 + n*( LZ-OZ );
        kk1  = 0;           C.Z1 = z + 0.5 + n*( LZ-OZ );
      } else if( z < Z[0] ){
        kk0  = K-1;         C.Z0 = Z[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LZ-OZ );
        kk1  = 0;           C.Z1 = Z[ 0 ] + n*( LZ-OZ );
      } else if( z > Z[K-1] ){
        kk0  = K-1;         C.Z0 = Z[kk0] + n*( LZ-OZ );
        kk1  = 0;           C.Z1 = Z[K-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LZ-OZ );
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0] + n*( LZ-OZ );
        kk1  = kk0 + 1;     C.Z1 = Z[kk1] + n*( LZ-OZ );
      }
      
      break;
  }

  
  C.DX = C.X1 - C.X0;
  C.DY = C.Y1 - C.Y0;
  C.DZ = C.Z1 - C.Z0;
  
  C.D  = C.DX * C.DY * C.DZ;
  C.iD = 1.0/C.D;
  
  
  if( !isSingularCase ){
    #define  VV(i,j,k)   VX[ (i)  +  (j)*I  +  (k)*IJ  ]
    
//     mexPrintf("ii0 , jj0 , kk0:  %d,%d,%d\n" , ii0 , jj0 , kk0 );
    
    C.VX0000 = VV(  ii0  ,  jj0  ,  kk0  );
//     mexPrintf("set VX0000:  %g\n" , C.VX0000 );
    
    C.VX1000 = VV(  ii1  ,  jj0  ,  kk0  );
    C.VX0100 = VV(  ii0  ,  jj1  ,  kk0  );
    C.VX1100 = VV(  ii1  ,  jj1  ,  kk0  );
    C.VX0010 = VV(  ii0  ,  jj0  ,  kk1  );
    C.VX1010 = VV(  ii1  ,  jj0  ,  kk1  );
    C.VX0110 = VV(  ii0  ,  jj1  ,  kk1  );
    C.VX1110 = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    #define  VV(i,j,k)   VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VY0000 = VV(  ii0  ,  jj0  ,  kk0  );
    C.VY1000 = VV(  ii1  ,  jj0  ,  kk0  );
    C.VY0100 = VV(  ii0  ,  jj1  ,  kk0  );
    C.VY1100 = VV(  ii1  ,  jj1  ,  kk0  );
    C.VY0010 = VV(  ii0  ,  jj0  ,  kk1  );
    C.VY1010 = VV(  ii1  ,  jj0  ,  kk1  );
    C.VY0110 = VV(  ii0  ,  jj1  ,  kk1  );
    C.VY1110 = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VZ0000 = VV(  ii0  ,  jj0  ,  kk0  );
      C.VZ1000 = VV(  ii1  ,  jj0  ,  kk0  );
      C.VZ0100 = VV(  ii0  ,  jj1  ,  kk0  );
      C.VZ1100 = VV(  ii1  ,  jj1  ,  kk0  );
      C.VZ0010 = VV(  ii0  ,  jj0  ,  kk1  );
      C.VZ1010 = VV(  ii1  ,  jj0  ,  kk1  );
      C.VZ0110 = VV(  ii0  ,  jj1  ,  kk1  );
      C.VZ1110 = VV(  ii1  ,  jj1  ,  kk1  );
      #undef VV
    }

  } else {
    
    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VX[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VX0000 = VV(  ii0  ,  jj0  ,  kk0  );
    C.VX1000 = VV(  ii1  ,  jj0  ,  kk0  );
    C.VX0100 = VV(  ii0  ,  jj1  ,  kk0  );
    C.VX1100 = VV(  ii1  ,  jj1  ,  kk0  );
    C.VX0010 = VV(  ii0  ,  jj0  ,  kk1  );
    C.VX1010 = VV(  ii1  ,  jj0  ,  kk1  );
    C.VX0110 = VV(  ii0  ,  jj1  ,  kk1  );
    C.VX1110 = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VY0000 = VV(  ii0  ,  jj0  ,  kk0  );
    C.VY1000 = VV(  ii1  ,  jj0  ,  kk0  );
    C.VY0100 = VV(  ii0  ,  jj1  ,  kk0  );
    C.VY1100 = VV(  ii1  ,  jj1  ,  kk0  );
    C.VY0010 = VV(  ii0  ,  jj0  ,  kk1  );
    C.VY1010 = VV(  ii1  ,  jj0  ,  kk1  );
    C.VY0110 = VV(  ii0  ,  jj1  ,  kk1  );
    C.VY1110 = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VZ0000 = VV(  ii0  ,  jj0  ,  kk0  );
      C.VZ1000 = VV(  ii1  ,  jj0  ,  kk0  );
      C.VZ0100 = VV(  ii0  ,  jj1  ,  kk0  );
      C.VZ1100 = VV(  ii1  ,  jj1  ,  kk0  );
      C.VZ0010 = VV(  ii0  ,  jj0  ,  kk1  );
      C.VZ1010 = VV(  ii1  ,  jj0  ,  kk1  );
      C.VZ0110 = VV(  ii0  ,  jj1  ,  kk1  );
      C.VZ1110 = VV(  ii1  ,  jj1  ,  kk1  );
      #undef VV
    }

  }

  return;
}





_inline void getV( real *xyzt , real *v ){
  real  u0, u1, v0, v1, w0, w1;
  real  u0v0w0 ,
        u1v0w0 ,
        u0v1w0 ,
        u1v1w0 ,
        u0v0w1 ,
        u1v0w1 ,
        u0v1w1 ,
        u1v1w1 ;
  
  #define x xyzt[0]
  #define y xyzt[1]
  #define z xyzt[2]
  #define t xyzt[3]
  if( 
      ( x < C.X0 )  ||  ( x > C.X1 )   ||
      ( y < C.Y0 )  ||  ( y > C.Y1 )   ||
      ( z < C.Z0 )  ||  ( z > C.Z1 )   
  ){
//     mexPrintf("setCELL\n");
    setCELL( x , y , z );
//     mexPrintf("[ii,C.X0,x,C.X1]:%d,%g,%g,%g  ", C.ii , C.X0 , x , C.X1 );
//     mexPrintf("[jj,C.Y0,y,C.Y1]:%d,%g,%g,%g  ", C.jj , C.Y0 , y , C.Y1 );
//     mexPrintf("[kk,C.Z0,z,C.Z1]:%d,%g,%g,%g  ", C.kk , C.Z0 , z , C.Z1 );
//     mexPrintf("[C.iD,C.D]:%g,%g\n",C.iD,C.D);
  }
  
  if( C.isOutside ){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return;
  }
  
  
  if( NSD == 3 ){
    
    u0 = x - C.X0;   u1 = C.X1 - x;
    v0 = y - C.Y0;   v1 = C.Y1 - y;
    w0 = z - C.Z0;   w1 = C.Z1 - z;

    u0v0w0 = u0*v0*w0;
    u1v0w0 = u1*v0*w0;
    u0v1w0 = u0*v1*w0;
    u1v1w0 = u1*v1*w0;
    u0v0w1 = u0*v0*w1;
    u1v0w1 = u1*v0*w1;
    u0v1w1 = u0*v1*w1;
    u1v1w1 = u1*v1*w1;

    v[0] = (  C.VX0000*u1v1w1 +
              C.VX1000*u0v1w1 +
              C.VX0100*u1v0w1 +
              C.VX1100*u0v0w1 +
              C.VX0010*u1v1w0 +
              C.VX1010*u0v1w0 +
              C.VX0110*u1v0w0 +
              C.VX1110*u0v0w0
           )* C.iD;

    v[1] = (  C.VY0000*u1v1w1 +
              C.VY1000*u0v1w1 +
              C.VY0100*u1v0w1 +
              C.VY1100*u0v0w1 +
              C.VY0010*u1v1w0 +
              C.VY1010*u0v1w0 +
              C.VY0110*u1v0w0 +
              C.VY1110*u0v0w0
           )* C.iD;
  
    v[2] = (  C.VZ0000*u1v1w1 +
              C.VZ1000*u0v1w1 +
              C.VZ0100*u1v0w1 +
              C.VZ1100*u0v0w1 +
              C.VZ0010*u1v1w0 +
              C.VZ1010*u0v1w0 +
              C.VZ0110*u1v0w0 +
              C.VZ1110*u0v0w0
           )* C.iD;

  } else { 
    
    u0 = x - C.X0;   u1 = C.X1 - x;
    v0 = y - C.Y0;   v1 = C.Y1 - y;

//     mexPrintf("u0,u1:  %g,%g\n",u0,u1);
//     mexPrintf("v0,v1:  %g,%g\n",v0,v1);
    
    u0v0w1 = u0*v0;
    u1v0w1 = u1*v0;
    u0v1w1 = u0*v1;
    u1v1w1 = u1*v1;

//     mexPrintf("C.VX0000: %g\n",C.VX0000);
//     mexPrintf("C.VX1000: %g\n",C.VX1000);
//     mexPrintf("C.VX0100: %g\n",C.VX0100);
//     mexPrintf("C.VX1100: %g\n",C.VX1100);

//     mexPrintf("C.VY0000: %g\n",C.VY0000);
//     mexPrintf("C.VY1000: %g\n",C.VY1000);
//     mexPrintf("C.VY0100: %g\n",C.VY0100);
//     mexPrintf("C.VY1100: %g\n",C.VY1100);
    
    
    
    v[0] = (  C.VX0000*u1v1w1 +
              C.VX1000*u0v1w1 +
              C.VX0100*u1v0w1 +
              C.VX1100*u0v0w1 
           )* C.iD;

    v[1] = (  C.VY0000*u1v1w1 +
              C.VY1000*u0v1w1 +
              C.VY0100*u1v0w1 +
              C.VY1100*u0v0w1 
           )* C.iD;
  
    v[2] = 0; 

//     mexPrintf("v[0]:  %g\n",v[0]);
//     mexPrintf("v[1]:  %g\n",v[1]);
//     mexPrintf("v[2]:  %g\n",v[2]);
  }

  return;
  #undef  x
  #undef  y
  #undef  z
  #undef  t
}





_inline void getDV( real *xyzt , real *Dv ){
  real  u0, u1, v0, v1, w0, w1;
  
  #define x xyzt[0]
  #define y xyzt[1]
  #define z xyzt[2]
  #define t xyzt[3]
  if( 
      ( x < C.X0 )  ||  ( x > C.X1 )   ||
      ( y < C.Y0 )  ||  ( y > C.Y1 )   ||
      ( z < C.Z0 )  ||  ( z > C.Z1 )   
  ){
    setCELL( x , y , z );
  }
  
  if( C.isOutside ){
    Dv[0] = Dv[1] = Dv[2] = Dv[3] = Dv[4] = Dv[5] = Dv[6] = Dv[7] = Dv[8] = 0.0;
    return;
  }
  
  u0 = x - C.X0;   u1 = C.X1 - x;
  v0 = y - C.Y0;   v1 = C.Y1 - y;
  w0 = z - C.Z0;   w1 = C.Z1 - z;
  
  Dv[0] = ( - C.VX0000*v1*w1
            + C.VX1000*v1*w1
            - C.VX0100*v0*w1
            + C.VX1100*v0*w1
            - C.VX0010*v1*w0
            + C.VX1010*v1*w0
            - C.VX0110*v0*w0
            + C.VX1110*v0*w0
          )* C.iD;

  Dv[1] = ( - C.VY0000*v1*w1
            + C.VY1000*v1*w1
            - C.VY0100*v0*w1
            + C.VY1100*v0*w1
            - C.VY0010*v1*w0
            + C.VY1010*v1*w0
            - C.VY0110*v0*w0
            + C.VY1110*v0*w0
          )* C.iD;

  Dv[2] = ( - C.VZ0000*v1*w1
            + C.VZ1000*v1*w1
            - C.VZ0100*v0*w1
            + C.VZ1100*v0*w1
            - C.VZ0010*v1*w0
            + C.VZ1010*v1*w0
            - C.VZ0110*v0*w0
            + C.VZ1110*v0*w0
          )* C.iD;

  
  Dv[3] = ( - C.VX0000*u1*w1
            - C.VX1000*u0*w1
            + C.VX0100*u1*w1
            + C.VX1100*u0*w1
            - C.VX0010*u1*w0
            - C.VX1010*u0*w0
            + C.VX0110*u1*w0
            + C.VX1110*u0*w0
          )* C.iD;

  Dv[4] = ( - C.VY0000*u1*w1
            - C.VY1000*u0*w1
            + C.VY0100*u1*w1
            + C.VY1100*u0*w1
            - C.VY0010*u1*w0
            - C.VY1010*u0*w0
            + C.VY0110*u1*w0
            + C.VY1110*u0*w0
          )* C.iD;

  Dv[5] = ( - C.VZ0000*u1*w1
            - C.VZ1000*u0*w1
            + C.VZ0100*u1*w1
            + C.VZ1100*u0*w1
            - C.VZ0010*u1*w0
            - C.VZ1010*u0*w0
            + C.VZ0110*u1*w0
            + C.VZ1110*u0*w0
          )* C.iD;
  
  Dv[6] = ( - C.VX0000*u1*v1
            - C.VX1000*u0*v1
            - C.VX0100*u1*v0
            - C.VX1100*u0*v0
            + C.VX0010*u1*v1
            + C.VX1010*u0*v1
            + C.VX0110*u1*v0
            + C.VX1110*u0*v0
           )* C.iD;

  Dv[7] = ( - C.VY0000*u1*v1
            - C.VY1000*u0*v1
            - C.VY0100*u1*v0
            - C.VY1100*u0*v0
            + C.VY0010*u1*v1
            + C.VY1010*u0*v1
            + C.VY0110*u1*v0
            + C.VY1110*u0*v0
           )* C.iD;
  
  Dv[8] = ( - C.VZ0000*u1*v1
            - C.VZ1000*u0*v1
            - C.VZ0100*u1*v0
            - C.VZ1100*u0*v0
            + C.VZ0010*u1*v1
            + C.VZ1010*u0*v1
            + C.VZ0110*u1*v0
            + C.VZ1110*u0*v0
          )* C.iD;
  
  return;
  #undef  x
  #undef  y
  #undef  z
  #undef  t
}


void MatrixExp3x3( real *A ){
  real    normA;
  real    A1[9] = { 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 };
  real    A2[9];
  real    A4[9];
  real    A6[9];
  real    A8[9];
  real    U[9] ;
  real    V[9] ;
  real    aux[9];
  int     i, s;
  real    scale = 2.0;
  real    det;
  
  normA = MAX3(  ABS( A[0] ) + ABS( A[1] ) + ABS( A[2] ) ,
                 ABS( A[3] ) + ABS( A[4] ) + ABS( A[5] ) ,
                 ABS( A[6] ) + ABS( A[7] ) + ABS( A[8] ) );

  
  if( normA <= 1.495585217958292e-2 ){

    Mtimes3x3( A2 , A , A );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A2[i] +  60.0*A1[i]; 
      V[i] = 12.0*A2[i] + 120.0*A1[i];
    }

  } else if( normA <= 0.2539398330063230 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A4[i] +  420.0*A2[i] + 15120.0*A1[i]; 
      V[i] = 30.0*A4[i] + 3360.0*A2[i] + 30240.0*A1[i];
    }

  } else if( normA <= 0.9504178996162932 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A6[i] +  1512.0*A4[i] +  277200.0*A2[i] +  8648640.0*A1[i]; 
      V[i] = 56.0*A6[i] + 25200.0*A4[i] + 1995840.0*A2[i] + 17297280.0*A1[i];
    }

  } else if( normA <= 2.097847961257068 ){
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    Mtimes3x3( A8 , A6 , A2 );
    for(i=0;i<9;i++){ 
      U[i] =  1.0*A8[i] +   3960.0*A6[i] +  2162160.0*A4[i] +  302702400.0*A2[i] +  8821612800.0*A1[i]; 
      V[i] = 90.0*A8[i] + 110880.0*A6[i] + 30270240.0*A4[i] + 2075673600.0*A2[i] + 17643225600.0*A1[i];
    }
    
  } else if( normA <= 5.371920351148152 ){
    
    
    Mtimes3x3( A2 , A  , A  );
    Mtimes3x3( A4 , A2 , A2 );
    Mtimes3x3( A6 , A4 , A2 );
    
    for(i=0;i<9;i++){ 
      U[i] =   1.0*A6[i] +  16380.0*A4[i] +   40840800.0*A2[i];
      V[i] = 182.0*A6[i] + 960960.0*A4[i] + 1323241920.0*A2[i];
    }
    
    Mtimes3x3( aux , A6 , U ); memcpy( U , aux , 9*sizeof( real ) );
    Mtimes3x3( aux , A6 , V ); memcpy( V , aux , 9*sizeof( real ) );
        
    for(i=0;i<9;i++){ 
      U[i] = U[i] +  33522128640.0*A6[i] +  10559470521600.0*A4[i] + 1187353796428800.0*A2[i] + 32382376266240000.0*A1[i];
      V[i] = V[i] + 670442572800.0*A6[i] + 129060195264000.0*A4[i] + 7771770303897600.0*A2[i] + 64764752532480000.0*A1[i];
    }
    
  } else {

    #define log2   0.69314718055994529
    s = (int) ( log( normA/5.371920351148152 )/log2  +  1 );
    for( i = 1; i<s ; i++ ){
      scale = scale*2.0;
    }
    scale = 1.0/scale;
    
    aux[0] = A[0]*scale;
    aux[1] = A[1]*scale;
    aux[2] = A[2]*scale;
    aux[3] = A[3]*scale;
    aux[4] = A[4]*scale;
    aux[5] = A[5]*scale;
    aux[6] = A[6]*scale;
    aux[7] = A[7]*scale;
    aux[8] = A[8]*scale;
    
    MatrixExp3x3( aux );
    memcpy( A , aux , 9*sizeof( real ) );
    
    for( i=0 ; i<s ; i++ ){
      Mtimes3x3( aux , A , A );
      memcpy( A , aux , 9*sizeof( real ) );
    }
    
    return;

  }
  
  Mtimes3x3( aux , A , U ); memcpy( U , aux , 9*sizeof( real ) );
  for(i=0;i<9;i++){ 
    A2[i] = V[i]-U[i];
    A4[i] = V[i]+U[i];
  }

  det =   A2[0]*( A2[4]*A2[8] - A2[5]*A2[7] )
        - A2[1]*( A2[3]*A2[8] - A2[5]*A2[6] )
        + A2[2]*( A2[3]*A2[7] - A2[4]*A2[6] );
  det = 1.0/det;
  
  aux[0] = ( A2[4]*A2[8] - A2[5]*A2[7] )*det;
  aux[1] = ( A2[2]*A2[7] - A2[1]*A2[8] )*det;
  aux[2] = ( A2[1]*A2[5] - A2[2]*A2[4] )*det;
  aux[3] = ( A2[5]*A2[6] - A2[3]*A2[8] )*det;
  aux[4] = ( A2[0]*A2[8] - A2[2]*A2[6] )*det;
  aux[5] = ( A2[2]*A2[3] - A2[0]*A2[5] )*det;
  aux[6] = ( A2[3]*A2[7] - A2[4]*A2[6] )*det;
  aux[7] = ( A2[1]*A2[6] - A2[0]*A2[7] )*det;
  aux[8] = ( A2[0]*A2[4] - A2[1]*A2[3] )*det;
  
  Mtimes3x3( A , aux , A4 );
  
  return;
}
