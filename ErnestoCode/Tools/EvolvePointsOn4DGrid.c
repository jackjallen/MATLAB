/*
   phi = EvolvePointsOn4DGrid( V , Gx , Gy , Gz , Gt , points ,
                           TIME_INTERP_MODE 
                           BOUNDARY_MODE  , [ BOUNDARY_SIZE ]
                          'maxSTEPS'      , 50                 dt_min = 1/maxSTEPS
                          'minStepsPerVoxel'      , 1                  dt_max = S/V/minStepsPerVoxel   (S, min voxel size ) (V, max V(x) )
                          'relTOL'        , 1/100              at each step control the error
                          'flow'          , [ 0   1 ]          a vector of Ts where evaluate the output FLOW (exept at the first t)
                          { 'rk23'  'euler' 's_euler' 'rk45'  }
                          { 'jac_plus' 'jac_times' }
                       );

      TIME_INTERP_MODE:
           'lineal'                         by default
           'constant'

 
      BOUNDARY_MODE:
           'value'   (default)              stop the evolution outside the grid
           'closest'
           'per[iodic]','circ[ular]'        periodic boundary conditions
           'decay[_to_zero]'                decay to zero at distance BOUNDARY_SIZE
                              after 'decay' specify BOUNDARY_SIZE 
                                        by default BOUNDARY_SIZE = (LX+LY+LZ)/2

  phi = EvolvePointsOn3DGrid( V , Gx , Gy , '2d' , Gt , { points , initJACS } )

  

 
  if FLOW is a scalar, integrate between t0 = 0   and  t_end = FLOW
  if FLOW is a vector, integrate with t0 = F(1)  and  t_end = F(end)
      and return, coordinates, jacobians and determinants at times F(2:end)
  if you also want the initial coords, jacs and dets use for example
      F = [ 0 , 0 , 0.1 , 0.2 , 0.3 , 0.5 , 1.0 , 2.0 ]  (replicate t0!!)
 

   phi = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , Gt , { points , initJACS } )

   [phi,jac] = EvolvePointsOn3DGrid( V , Gx , Gy , Gz , Gt , { points , initJACS } )

   [phi,jac,det_jac] = EvolvePointsOn3DGrid( V , Gx , Gy , Gz ,  Gt , { points , initJACS } )

 
    numel( points )  have to be multiple of 3
 
    size( initJACS ) have to be equal to [ size(points)  3 ]
 
 
 
    size( phi )  =  [ size( points )  numel(FLOW)-1  ]
    size( jac )  =  [ size( points )  3   numel(FLOW)-1  ]
    size( det_jac )  =  [ numel( points )/3   numel(FLOW)-1  ]
 
*/

#include "myMEX.h"
#if !defined( memcpy )
  #include "string.h"
#endif

#if !defined( real )
  #define   real       double
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define    integration_method_def   RK23
#define    jacobian_method_def      PLUS
#define    boundary_mode_def        VALUE
#define    tinterp_mode_def         LINEAL
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


#define EVOLVE_JAC          if( computeJAC ){                                                \
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
enum    tinterp_modes  { LINEAL , CONSTANT } tinterp_mode;
real    *VX, *VY, *VZ, *X, *Y, *Z, *T, OX, OY, OZ, LX, LY, LZ;
int     I , J , K , S , nT , IJ , IJK , IJKS, dim3D;
real    boundary_size;

typedef struct CELL {
  int  ii, jj, kk, tt;
  real X0, X1, Y0, Y1, Z0, Z1, T0, T1;
  real VX0000, VX1000, VX0100, VX1100, VX0010, VX1010, VX0110, VX1110, VX0001, VX1001, VX0101, VX1101, VX0011, VX1011, VX0111, VX1111;
  real VY0000, VY1000, VY0100, VY1100, VY0010, VY1010, VY0110, VY1110, VY0001, VY1001, VY0101, VY1101, VY0011, VY1011, VY0111, VY1111;
  real VZ0000, VZ1000, VZ0100, VZ1100, VZ0010, VZ1010, VZ0110, VZ1110, VZ0001, VZ1001, VZ0101, VZ1101, VZ0011, VZ1011, VZ0111, VZ1111;
  real DX, DY, DZ, DT, D, iD;
  int  isOutside;
} CELL; 
CELL C =  { -1 , -1 , -1 , -1, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };


_inline void setCELL( real x , real y , real z , real t );
_inline void getV( real *xyzt , real *v );
_inline void getDV( real *xyzt , real *Dv );
        void MatrixExp3x3( real * );
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    integration_methods  { EULER , S_EULER , RK45 , RK23 }    integration_method;
  enum    jacobian_methods     { PLUS   , TIMES }                   jacobian_method;

  real    XYZT[4], xk[4], v[3], vp[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], err[3], Dv[9];
  real    maxV, minD, hmax, hmin, tol, h;
  real    minStepsPerVoxel, relTOL;
  
  real    JAC[9], dJAC[9];
  
  real    *P, *F;
  int     Odims[50], nD, d, nP, p, nF, ff;

  int     maxSTEPS, Nadvice = 0;


  real    *O;
  char    STR[100];
  int     argN;
  
  int     NSD  , provideINIT_JAC = 0, computeJAC = 0, computeDET = 0;
  real    *initJAC, *OJAC, *DET;
  

  I  = mySize( prhs[0] , 0 );
  J  = mySize( prhs[0] , 1 );

  if( mxIsChar( prhs[3] ) ){
    mxGetString( prhs[3], STR, 100 );
    if( myStrcmpi(STR,"2d") ) { myErrMsgTxt("The only valid keyword is 2D."); }
    K = 1;
    dim3D = 0;
    NSD = 2;
  } else {
    K  = mySize( prhs[0] , 2 );
    dim3D = 1;
    NSD = 3;
  }
  S  = mySize( prhs[0] , 3 );
  IJ = I*J;
  IJK = IJ*K;
  IJKS = IJK*S;

  if( myNumel( prhs[0] ) != IJKS*NSD ){
    myErrMsgTxt("size(V) has to be [numel(Gx) numel(Gy) numel(Gz)  nTimes  3]   or  [numel(Gx) numel(Gy) nTimes 2] in 2D case.");
  }
  
  VX = myGetPr( prhs[0] );
  VY = VX + IJKS;
  VZ = VY + IJKS;

  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(Gy) Coordinates do not coincide with size(V,2)."); }
  if( dim3D && myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(Gz) Coordinates do not coincide with size(V,3)."); }

  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("Gx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("Gy  Coordinates are not sorted."); }
  
  if( dim3D ){
    Z = myGetPr(prhs[3]);
    if( !checkIsSorted(Z,K) ){ myErrMsgTxt("Gz  Coordinates are not sorted."); }
  } else {
    Z = mxMalloc( 1*sizeof( real ) );
    Z[0] = 0;
  }
   T = myGetPr(prhs[4]);
  nT = myNumel(prhs[4]);
  if( nT < 2 ){
    myErrMsgTxt("Gt  must define a non-empty interval!!.");
  }
  
  if( !checkIsSorted(T,nT) ){ myErrMsgTxt("Gt  Coordinates are not sorted."); }

  
  if( mxIsCell( prhs[5] ) ){
    provideINIT_JAC = 1;
    
    /* checking the sizes of the initial Points and Jacobians */
    nD = myNDims( mxGetCell( prhs[5] , 0 ));

    if( ( nD + 1 ) != myNDims( mxGetCell( prhs[5] , 1 ) ) ){
      myErrMsgTxt("init_JAC must have size equal to [ size(P) %d ].", NSD);
    }
    
    for( d = 0 ; d < nD ; d++ ){
      if( mySize( mxGetCell( prhs[5] , 0 ) , d ) != mySize( mxGetCell( prhs[5] , 1 ) , d ) ) {
        myErrMsgTxt("init_JAC must have size equal to [ size(P) %d ].", NSD);
      }
    }

    if( mySize( mxGetCell( prhs[5] , 1 ) , d+1 ) ) {
      myErrMsgTxt("init_JAC must have size equal to [ size(P) %d ].", NSD);
    }
    
    nP = myNumel( mxGetCell( prhs[5] , 0 ) );
    if( ( nP%NSD ) ){ myErrMsgTxt("Number of P have to be multiple of %d.", NSD); }
    nP = nP/NSD;
    /*END checking the sizes of the initial Points and Jacobians*/

    P  = myGetPr( mxGetCell( prhs[5] , 0 ) );
    initJAC = myGetPr( mxGetCell( prhs[5] , 1 ) );

  } else {
    
    nP = myNumel( prhs[5] );
    if( ( nP%NSD ) ){ myErrMsgTxt("Number of P have to be multiple of %d.", NSD); }
    nP = nP/NSD;
    P  = myGetPr( prhs[5] );

  }
  
  /*Parsing arguments*/
  /*Defaults*/
  boundary_size = (X[I-1]-X[0]+Y[J-1]-Y[0]+Z[K-1]-Z[0])/NSD;
  integration_method  = integration_method_def;
  jacobian_method     = jacobian_method_def;
  boundary_mode       = boundary_mode_def;
  tinterp_mode        = tinterp_mode_def;
  relTOL              = relTOL_def;
  minStepsPerVoxel    = minStepsPerVoxel_def;
  maxSTEPS            = maxSTEPS_def;
  F                   = NULL;
  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );

      if( ! myStrcmpi(STR,"lineal")                                  ) { tinterp_mode = LINEAL;         argN++; continue; }
      if( ! myStrcmpi(STR,"constant")                                ) { tinterp_mode = CONSTANT;       argN++; continue; }

      if( ! myStrcmpi(STR,"s_euler")                                 ) { integration_method = S_EULER;  argN++; continue; }
      if( ! myStrcmpi(STR,"euler")                                   ) { integration_method = EULER;    argN++; continue; }
      if( ! myStrcmpi(STR,"rk23")                                    ) { integration_method = RK23;     argN++; continue; }
      if( ! myStrcmpi(STR,"rk45")                                    ) { integration_method = RK45;     argN++; continue; }

      
      if( !myStrcmpi(STR,"jac_times") || !myStrcmpi(STR,"jactimes")  ) { jacobian_method = TIMES;       argN++; continue; }
      if( !myStrcmpi(STR,"jac_plus")  || !myStrcmpi(STR,"jacplus")   ) { jacobian_method = PLUS;        argN++; continue; }
      
      
      if( ! myStrcmpi(STR,"circular") || ! myStrcmpi(STR,"periodic") || 
          ! myStrcmpi(STR,"circ")     || ! myStrcmpi(STR,"per")      ) { boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"value")    || ! myStrcmpi(STR,"zero")     ) { boundary_mode = VALUE;    argN++; continue; }
      if( ! myStrcmpi(STR,"closest")                                 ) { boundary_mode = CLOSEST;  argN++; continue; }

      if( ! myStrcmpi( STR,"decay_to_zero") || 
          ! myStrcmpi( STR,"tozero")  || ! myStrcmpi( STR,"decay")   ) { boundary_mode = DECAY_TO_ZERO; argN++;
        if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
          if( !myIsEmpty( prhs[argN] ) ){ boundary_size = myGetValue(prhs[argN]); }
          argN++; continue;
        }
        continue;
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
      
      
      mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
    }
  }
  /*END Parsing arguments*/
  
  if( F == NULL ){
    nF = 2;
    F = mxMalloc( 2*sizeof( real ) );
    F[0] = 0;
    F[1] = 1;
  }
  if( !checkIsSorted(F,nF) ){ myErrMsgTxt("FLOW times are not sorted."); }


  if( F[0] < T[0]  ||  F[nF-1] > T[nT-1] ){
    myErrMsgTxt("FLOW interval must be inside T interval.");
  }
  
  if( tinterp_mode == LINEAL   &&   S != nT ){
    myErrMsgTxt("Using  LINEAL  time_interpolation_mode, size(V,4) have to be equal to numel(Gt) .");
  }
  if( tinterp_mode == CONSTANT   &&   S != nT-1 ){
    myErrMsgTxt("Using  LINEAL  time_interpolation_mode, size(V,4) have to be equal to numel(Gt)-1 .");
  }
  
  
  /*Creating output*/
  nD = myNDims( prhs[5] );
  for( d=0 ; d<nD ; d++ ){
    Odims[d] = mySize( prhs[5] , d );
  }
  Odims[d++] = nF-1;
  plhs[0] = mxCreateNumericArray( d , Odims , mxREAL_CLASS , mxREAL );
  O = (real *) mxGetData( plhs[0] );


  if( nlhs > 1 ){
    nD = myNDims( prhs[5] );
    for( d=0 ; d<nD ; d++ ){
      Odims[d] = mySize( prhs[5] , d );
    }
    Odims[d++] = NSD;
    Odims[d++] = nF-1;
    plhs[1] = mxCreateNumericArray( d , Odims , mxREAL_CLASS , mxREAL );
    OJAC = (real *) mxGetData( plhs[1] );
    computeJAC = 1;
  }
  

  if( nlhs > 2 ){
    Odims[0] = nP;
    Odims[1] = nF-1;
    plhs[2] = mxCreateNumericArray( 2 , Odims , mxREAL_CLASS , mxREAL );
    DET = (real *) mxGetData( plhs[2] );
    computeDET = 1;
  }
  /*END Creating output*/


  switch( boundary_mode ){
    case VALUE:
      OX = X[0];  LX = X[I-1];
      OY = Y[0];  LY = Y[J-1];
      OZ = Z[0];  LZ = Z[K-1];
      break;
    case CLOSEST:
      OX = X[0];  LX = X[I-1];
      OY = Y[0];  LY = Y[J-1];
      OZ = Z[0];  LZ = Z[K-1];
      break;
    case DECAY_TO_ZERO:
      OX = X[0];  LX = X[I-1];
      OY = Y[0];  LY = Y[J-1];
      OZ = Z[0];  LZ = Z[K-1];
      break;
    case CIRCULAR:
      if( I == 1 ){
        OX = X[0] - 0.5;  LX = X[0] + 0.5;
      } else {
        OX = X[0] - ( X[1]-X[0] )/2.0;  LX = X[I-1] + ( X[I-1]-X[I-2] )/2.0;
      }
      
      if( J == 1 ){
        OY = Y[0] - 0.5;  LY = Y[0] + 0.5;
      } else {
        OY = Y[0] - ( Y[1]-Y[0] )/2.0;  LY = Y[J-1] + ( Y[J-1]-Y[J-2] )/2.0;
      }

      if( K == 1 ){
        OZ = Z[0] - 0.5;  LZ = Z[0] + 0.5;
      } else {
        OZ = Z[0] - ( Z[1]-Z[0] )/2.0;  LZ = Z[K-1] + ( Z[K-1]-Z[K-2] )/2.0;
      }
      
      break;
  }

  
  if( integration_method == S_EULER  ||  integration_method == RK23  || integration_method == RK45  ){
    minD = 1.0;
    if( I > 1 ){ minD = MIN( minD , (LX-OX)/I ); }
    if( J > 1 ){ minD = MIN( minD , (LY-OY)/J ); }
    if( K > 1 ){ minD = MIN( minD , (LZ-OZ)/K ); }
    tol  = minD*relTOL;
  }
  
  hmin     = ( F[nF-1] - F[0] )/maxSTEPS;
  minStepsPerVoxel = 1.0/minStepsPerVoxel;

  XYZT[2] = 0;
  for( p=0 ; p < nP ; p++ ){
    XYZT[0] = P[ p        ];
    XYZT[1] = P[ p +   nP ];
    if( dim3D ){ XYZT[2] = P[ p + 2*nP ]; }
    XYZT[3] = F[0];
    
    if( computeJAC ){
      if( provideINIT_JAC ){
        if( dim3D ){
          for( d = 0 ; d < 9 ; d++ ){ JAC[d] = initJAC[ p + d*nP ]; }
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
          while( XYZT[3] < F[ff] ){
            getV( XYZT , v );

            if( !v[0] && !v[1] && !v[2] ){ break; }
            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }

            EVOLVE_JAC

            XYZT[0] += h*v[0];
            XYZT[1] += h*v[1];
            XYZT[2] += h*v[2];
            XYZT[3] += h;
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
            xk[3] = XYZT[3] + h*Ct45_2;
            getV( xk , k2 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );


            xk[0] = XYZT[0] + A45_31*k1[0] + A45_32*k2[0];
            xk[1] = XYZT[1] + A45_31*k1[1] + A45_32*k2[1];
            xk[2] = XYZT[2] + A45_31*k1[2] + A45_32*k2[2];
            xk[3] = XYZT[3] + h*Ct45_3;
            getV( xk , k3 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k3 , k3 , h );

            xk[0] = XYZT[0] + A45_41*k1[0] + A45_42*k2[0] + A45_43*k3[0];
            xk[1] = XYZT[1] + A45_41*k1[1] + A45_42*k2[1] + A45_43*k3[1];
            xk[2] = XYZT[2] + A45_41*k1[2] + A45_42*k2[2] + A45_43*k3[2];
            xk[3] = XYZT[3] + h*Ct45_4;
            getV( xk , k4 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k4 , k4 , h );

            xk[0] = XYZT[0] + A45_51*k1[0] + A45_52*k2[0] + A45_53*k3[0] + A45_54*k4[0];
            xk[1] = XYZT[1] + A45_51*k1[1] + A45_52*k2[1] + A45_53*k3[1] + A45_54*k4[1];
            xk[2] = XYZT[2] + A45_51*k1[2] + A45_52*k2[2] + A45_53*k3[2] + A45_54*k4[2];
            xk[3] = XYZT[3] + h*Ct45_5;
            getV( xk , k5 );  
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k5 , k5 , h );

            xk[0] = XYZT[0] + A45_61*k1[0] + A45_62*k2[0] + A45_63*k3[0] + A45_64*k4[0] + A45_65*k5[0];
            xk[1] = XYZT[1] + A45_61*k1[1] + A45_62*k2[1] + A45_63*k3[1] + A45_64*k4[1] + A45_65*k5[1];
            xk[2] = XYZT[2] + A45_61*k1[2] + A45_62*k2[2] + A45_63*k3[2] + A45_64*k4[2] + A45_65*k5[2];
            xk[3] = XYZT[3] + h*Ct45_6;
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
            xk[3] = XYZT[3] + h*Ct23_2;
            getV( xk , k2 );  if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );

            xk[0] = XYZT[0] + A23_32*k2[0];
            xk[1] = XYZT[1] + A23_32*k2[1];
            xk[2] = XYZT[2] + A23_32*k2[2];
            xk[3] = XYZT[3] + h*Ct23_3;
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
      if( dim3D ){

        O[ p        + (ff-1)*3*nP ] = XYZT[0];     
        O[ p +   nP + (ff-1)*3*nP ] = XYZT[1];
        O[ p + 2*nP + (ff-1)*3*nP ] = XYZT[2];

        if( computeJAC ){
          for( d = 0 ; d < 9 ; d++ ){
            OJAC[ p + d*nP + (ff-1)*9*nP ] = JAC[d];
          }

          
          if( computeDET ){
            DET[ p + (ff-1)*nP ] =   JAC[0]*( JAC[4]*JAC[8] - JAC[5]*JAC[7] )
                                   - JAC[1]*( JAC[3]*JAC[8] - JAC[5]*JAC[6] )
                                   + JAC[2]*( JAC[3]*JAC[7] - JAC[4]*JAC[6] );
          }
        }


      } else {

        O[ p        + (ff-1)*2*nP ] = XYZT[0];
        O[ p +   nP + (ff-1)*2*nP ] = XYZT[1];

        if( computeJAC ){
          OJAC[ p +        (ff-1)*4*nP ] = JAC[0];
          OJAC[ p +   nP + (ff-1)*4*nP ] = JAC[1];
          OJAC[ p + 2*nP + (ff-1)*4*nP ] = JAC[3];
          OJAC[ p + 3*nP + (ff-1)*4*nP ] = JAC[4];
          
          if( computeDET ){
            DET[ p + (ff-1)*nP ] = JAC[0]*JAC[4] - JAC[3]*JAC[1];
          }
        }

      }

    }
  }

  if( Nadvice > 0 ){
    sprintf( STR , "Adviced maxSTEPS: %d" , Nadvice );
    mexWarnMsgIdAndTxt( "EXP:V:IncreaseMaxSteps" , STR );
  }


  EXIT: myFreeALLOCATES();
}





_inline void setCELL( real x , real y , real z , real t ){
  int   ii0 , ii1, jj0 , jj1, kk0, kk1, tt0 , tt1;
  int   isSingularCase;
  int   n;
  
  C.isOutside = 0;
  isSingularCase = 0;
  switch( boundary_mode ){
    case VALUE:
      if( ( I>1 && ( x < OX || x > LX ) ) || ( I==1 && ( x < OX-0.5 || x > LX+0.5 ) ) ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( ( J>1 && ( y < OY || y > LY ) ) || ( J==1 && ( y < OY-0.5 || y > LY+0.5 ) ) ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( ( K>1 && ( z < OZ || z > LZ ) ) || ( K==1 && ( z < OZ-0.5 || z > LZ+0.5 ) ) ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      
      if( I > 1 ){
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0];
        ii1  = ii0 + 1;     C.X1 = X[ii1];
      } else {
        ii0 = 0; C.X0 = OX - 0.5;
        ii1 = 0; C.X1 = LX + 0.5;
      }

      if( J > 1 ){
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0];
        jj1  = jj0 + 1;     C.Y1 = Y[jj1];
      } else {
        jj0 = 0; C.Y0 = OY - 0.5;
        jj1 = 0; C.Y1 = LY + 0.5;
      }

      if( K > 1 ){
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0];
        kk1  = kk0 + 1;     C.Z1 = Z[kk1];
      } else {
        kk0 = 0; C.Z0 = OZ - 0.5;
        kk1 = 0; C.Z1 = LZ + 0.5;
      }

      break;

    case CLOSEST:

      if( x < OX ){
        C.ii = 0;
        ii0  = 0;  C.X0 = x - 100;
        ii1  = 0;  C.X1 = OX;
      } else if( x > LX ){
        C.ii = MAX( 0 , I-2 );
        ii0  = I-1;  C.X0 = LX;
        ii1  = I-1;  C.X1 = x + 100;
      } else if( I == 1 ){
        ii0 = 0;   C.X0 = OX - 1;
        ii1 = 0;   C.X1 = LX + 1;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0];
        ii1  = ii0 + 1;     C.X1 = X[ii1];
      }
      
      if( y < OY ){
        C.jj = 0;
        jj0  = 0;  C.Y0 = y - 100;
        jj1  = 0;  C.Y1 = OY;
      } else if( y > LY ){
        C.jj = MAX( 0 , J-2 );
        jj0  = J-1;  C.Y0 = LY;
        jj1  = J-1;  C.Y1 = y + 100;
      } else if( J == 1 ){
        jj0 = 0;   C.Y0 = OY - 1;
        jj1 = 0;   C.Y1 = LY + 1;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0];
        jj1  = jj0 + 1;     C.Y1 = Y[jj1];
      }

      if( z < OZ ){
        C.kk = 0;
        kk0  = 0;  C.Z0 = z - 100;
        kk1  = 0;  C.Z1 = OZ;
      } else if( z > LZ ){
        C.kk = MAX( 0 , K-2 );
        kk0  = K-1;  C.Z0 = LZ;
        kk1  = K-1;  C.Z1 = z + 100;
      } else if( K == 1 ){
        kk0 = 0;   C.Z0 = OZ - 1;  
        kk1 = 0;   C.Z1 = LZ + 1;  
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0];
        kk1  = kk0 + 1;     C.Z1 = Z[kk1];
      }      

      break;


    case DECAY_TO_ZERO:
      
      if( x < OX-boundary_size || x > LX+boundary_size ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( y < OY-boundary_size || y > LY+boundary_size ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( z < OZ-boundary_size || z > LZ+boundary_size ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      
      if( x < OX ){
        isSingularCase = 1;
        C.ii = 0;
        ii0  = -1;        C.X0 = OX-boundary_size;
        ii1  = 0;         C.X1 = OX;
      } else if( x > LX ) {
        isSingularCase = 1;
        C.ii = MAX( 0 , I-2 );
        ii0  = I-1;       C.X0 = LX;
        ii1  = -1;        C.X1 = LX+boundary_size;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;      C.X0 = X[ii0];
        ii1  = ii0 + 1;   C.X1 = X[ii1];
      }
      
      if( y < OY ){
        isSingularCase = 1;
        C.jj = 0;
        jj0  = -1;        C.Y0 = OY-boundary_size;
        jj1  = 0;         C.Y1 = OY;
      } else if( y > LY ) {
        isSingularCase = 1;
        C.jj = MAX( 0 , J-2 );
        jj0  = J-1;       C.Y0 = LX;
        jj1  = -1;        C.Y1 = LX+boundary_size;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;      C.Y0 = Y[jj0];
        jj1  = jj0 + 1;   C.Y1 = Y[jj1];
      }

      if( z < OZ ){
        isSingularCase = 1;
        C.kk = 0;
        kk0  = -1;        C.Z0 = OZ-boundary_size;
        kk1  = 0;         C.Z1 = OZ;
      } else if( z > LZ ) {
        isSingularCase = 1;
        C.kk = MAX( 0 , K-2 );
        kk0  = K-1;       C.Z0 = LZ;
        kk1  = -1;        C.Z1 = LZ+boundary_size;
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

      if( x < X[0] ){
        ii0  = I-1;         C.X0 = OX - ( X[I-1]-X[I-2] )/2.0 + n*( LX-OX );
        ii1  = 0;           C.X1 = X[0] + n*( LX-OX );
      } else if( x > X[I-1] ){
        ii0  = I-1;         C.X0 = X[ii0] + n*( LX-OX );
        ii1  = 0;           C.X1 = LX + ( X[1]-X[0] )/2.0 + n*( LX-OX );
      } else if( I == 1) {
        ii0  = 0;           C.X0 = OX + n*(LX-OX);
        ii1  = 0;           C.X1 = LX + n*(LX-OX);
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        ii0  = C.ii;        C.X0 = X[ii0] + n*( LX-OX );
        ii1  = ii0 + 1;     C.X1 = X[ii1] + n*( LX-OX );
      }
      
      
      if( y < OY || y > LY ){
        n = floor( ( y - OY )/( LY-OY ) );
        y = y - n*( LY-OY );
      } else { n = 0; }

      if( y < Y[0] ){
        jj0  = J-1;         C.Y0 = OY - ( Y[J-1]-Y[J-2] )/2.0 + n*( LY-OY );
        jj1  = 0;           C.Y1 = Y[0] + n*( LY-OY );
      } else if( y > Y[J-1] ){
        jj0  = J-1;         C.Y0 = Y[jj0] + n*( LY-OY );
        jj1  = 0;           C.Y1 = LY + ( Y[1]-Y[0] )/2.0 + n*( LY-OY );
      } else if( J == 1) {
        jj0  = 0;           C.Y0 = OY + n*(LY-OY);
        jj1  = 0;           C.Y1 = LY + n*(LY-OY);
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        jj0  = C.jj;        C.Y0 = Y[jj0] + n*( LY-OY );
        jj1  = jj0 + 1;     C.Y1 = Y[jj1] + n*( LY-OY );
      }


      if( z < OZ || z > LZ ){
        n = floor( ( z - OZ )/( LZ-OZ ) );
        z = z - n*( LZ-OZ );
      } else { n = 0; }

      if( z < Z[0] ){
        kk0  = K-1;         C.Z0 = OZ - ( Z[K-1]-Z[K-2] )/2.0 + n*( LZ-OZ );
        kk1  = 0;           C.Z1 = Z[0] + n*( LZ-OZ );
      } else if( z > Z[K-1] ){
        kk0  = K-1;         C.Z0 = Z[kk0] + n*( LZ-OZ );
        kk1  = 0;           C.Z1 = LZ + ( Z[1]-Z[0] )/2.0 + n*( LZ-OZ );
      } else if( K == 1) {
        kk0  = 0;           C.Z0 = OZ + n*(LZ-OZ);
        kk1  = 0;           C.Z1 = LZ + n*(LZ-OZ);
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        kk0  = C.kk;        C.Z0 = Z[kk0] + n*( LZ-OZ );
        kk1  = kk0 + 1;     C.Z1 = Z[kk1] + n*( LZ-OZ );
      }

      break;
  }
  
  switch( tinterp_mode ){
    case LINEAL:
      C.tt  = GetInterval( t , T , nT , C.tt );
      tt0   = C.tt;         C.T0 = tt0;
      tt1   = tt0+1;        C.T1 = tt1;
      break;
    
    case CONSTANT: 
      C.tt  = GetInterval( t , T , nT , C.tt );
      tt0   = C.tt;         C.T0 = tt0;
      tt1   = tt0;          C.T1 = tt0+1;
      break;
  }
  
  C.DX = C.X1 - C.X0;
  C.DY = C.Y1 - C.Y0;
  C.DZ = C.Z1 - C.Z0;
  C.DT = C.T1 - C.T0;
  
  C.D  = C.DX * C.DY * C.DZ * C.DT;
  C.iD = 1.0/C.D;
  
  
  if( !isSingularCase ){
    #define  VV(i,j,k,t)   VX[ (i)  +  (j)*I  +  (k)*IJ  + (t)*IJK ]
    C.VX0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
    C.VX1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
    C.VX0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
    C.VX1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
    C.VX0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
    C.VX1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
    C.VX0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
    C.VX1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

    if( tinterp_mode == LINEAL ){
      C.VX0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
      C.VX1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
      C.VX0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
      C.VX1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
      C.VX0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
      C.VX1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
      C.VX0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
      C.VX1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
    }
    #undef VV

    #define  VV(i,j,k,t)   VY[ (i)  +  (j)*I  +  (k)*IJ  + (t)*IJK ]
    C.VY0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
    C.VY1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
    C.VY0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
    C.VY1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
    C.VY0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
    C.VY1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
    C.VY0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
    C.VY1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

    if( tinterp_mode == LINEAL ){
      C.VY0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
      C.VY1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
      C.VY0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
      C.VY1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
      C.VY0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
      C.VY1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
      C.VY0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
      C.VY1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
    }
    #undef VV

    if( dim3D ){
      #define  VV(i,j,k,t)   VZ[ (i)  +  (j)*I  +  (k)*IJ  + (t)*IJK ]
      C.VZ0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
      C.VZ1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
      C.VZ0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
      C.VZ1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
      C.VZ0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
      C.VZ1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
      C.VZ0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
      C.VZ1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

      if( tinterp_mode == LINEAL ){
        C.VZ0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
        C.VZ1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
        C.VZ0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
        C.VZ1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
        C.VZ0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
        C.VZ1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
        C.VZ0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
        C.VZ1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
      }
      #undef VV
    }


  } else {
    
    #define  VV(i,j,k,t)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VX[ (i)  +  (j)*I  +  (k)*IJ + (t)*IJK ]
    C.VX0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
    C.VX1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
    C.VX0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
    C.VX1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
    C.VX0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
    C.VX1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
    C.VX0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
    C.VX1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

    if( tinterp_mode == LINEAL ){
      C.VX0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
      C.VX1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
      C.VX0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
      C.VX1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
      C.VX0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
      C.VX1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
      C.VX0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
      C.VX1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
    }
    #undef VV

    #define  VV(i,j,k,t)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VY[ (i)  +  (j)*I  +  (k)*IJ  + (t)*IJK ]
    C.VY0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
    C.VY1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
    C.VY0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
    C.VY1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
    C.VY0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
    C.VY1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
    C.VY0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
    C.VY1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

    if( tinterp_mode == LINEAL ){
      C.VY0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
      C.VY1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
      C.VY0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
      C.VY1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
      C.VY0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
      C.VY1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
      C.VY0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
      C.VY1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
    }
    #undef VV

    if( dim3D ){
      #define  VV(i,j,k,t)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VZ[ (i)  +  (j)*I  +  (k)*IJ  + (t)*IJK ]
      C.VZ0000 = VV(  ii0  ,  jj0  ,  kk0  , tt0 );
      C.VZ1000 = VV(  ii1  ,  jj0  ,  kk0  , tt0 );
      C.VZ0100 = VV(  ii0  ,  jj1  ,  kk0  , tt0 );
      C.VZ1100 = VV(  ii1  ,  jj1  ,  kk0  , tt0 );
      C.VZ0010 = VV(  ii0  ,  jj0  ,  kk1  , tt0 );
      C.VZ1010 = VV(  ii1  ,  jj0  ,  kk1  , tt0 );
      C.VZ0110 = VV(  ii0  ,  jj1  ,  kk1  , tt0 );
      C.VZ1110 = VV(  ii1  ,  jj1  ,  kk1  , tt0 );

      if( tinterp_mode == LINEAL ){
        C.VZ0001 = VV(  ii0  ,  jj0  ,  kk0  , tt1 );
        C.VZ1001 = VV(  ii1  ,  jj0  ,  kk0  , tt1 );
        C.VZ0101 = VV(  ii0  ,  jj1  ,  kk0  , tt1 );
        C.VZ1101 = VV(  ii1  ,  jj1  ,  kk0  , tt1 );
        C.VZ0011 = VV(  ii0  ,  jj0  ,  kk1  , tt1 );
        C.VZ1011 = VV(  ii1  ,  jj0  ,  kk1  , tt1 );
        C.VZ0111 = VV(  ii0  ,  jj1  ,  kk1  , tt1 );
        C.VZ1111 = VV(  ii1  ,  jj1  ,  kk1  , tt1 );
      }
      #undef VV
    }

  }

  return;
}



_inline void getV( real *xyzt , real *v ){
  real  u0, u1, v0, v1, w0, w1, s0, s1;
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
      ( z < C.Z0 )  ||  ( z > C.Z1 )   || 
      ( t < C.T0 )  ||  ( t > C.T1 )
  ){
    setCELL( x , y , z , t );
  }
  
  if( C.isOutside ){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return;
  }
  
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
  
  if( tinterp_mode == LINEAL ){
    s0 = t - C.T0;   s1 = C.T1 - t;
    
    v[0] = (  C.VX0000*u1v1w1*s1 +    C.VX0001*u1v1w1*s0 +
              C.VX1000*u0v1w1*s1 +    C.VX1001*u0v1w1*s0 +
              C.VX0100*u1v0w1*s1 +    C.VX0101*u1v0w1*s0 +
              C.VX1100*u0v0w1*s1 +    C.VX1101*u0v0w1*s0 +
              C.VX0010*u1v1w0*s1 +    C.VX0011*u1v1w0*s0 +
              C.VX1010*u0v1w0*s1 +    C.VX1011*u0v1w0*s0 +
              C.VX0110*u1v0w0*s1 +    C.VX0111*u1v0w0*s0 +
              C.VX1110*u0v0w0*s1 +    C.VX1111*u0v0w0*s0  
           )* C.iD;
  
    v[1] = (  C.VY0000*u1v1w1*s1 +    C.VY0001*u1v1w1*s0 +
              C.VY1000*u0v1w1*s1 +    C.VY1001*u0v1w1*s0 +
              C.VY0100*u1v0w1*s1 +    C.VY0101*u1v0w1*s0 +
              C.VY1100*u0v0w1*s1 +    C.VY1101*u0v0w1*s0 +
              C.VY0010*u1v1w0*s1 +    C.VY0011*u1v1w0*s0 +
              C.VY1010*u0v1w0*s1 +    C.VY1011*u0v1w0*s0 +
              C.VY0110*u1v0w0*s1 +    C.VY0111*u1v0w0*s0 +
              C.VY1110*u0v0w0*s1 +    C.VY1111*u0v0w0*s0  
           )* C.iD;
    
    v[2] = (  C.VZ0000*u1v1w1*s1 +    C.VZ0001*u1v1w1*s0 +
              C.VZ1000*u0v1w1*s1 +    C.VZ1001*u0v1w1*s0 +
              C.VZ0100*u1v0w1*s1 +    C.VZ0101*u1v0w1*s0 +
              C.VZ1100*u0v0w1*s1 +    C.VZ1101*u0v0w1*s0 +
              C.VZ0010*u1v1w0*s1 +    C.VZ0011*u1v1w0*s0 +
              C.VZ1010*u0v1w0*s1 +    C.VZ1011*u0v1w0*s0 +
              C.VZ0110*u1v0w0*s1 +    C.VZ0111*u1v0w0*s0 +
              C.VZ1110*u0v0w0*s1 +    C.VZ1111*u0v0w0*s0  
           )* C.iD;

  } else {
    
    v[0] = (  C.VX0000*u1v1w1 +
              C.VX1000*u0v1w1 +
              C.VX0100*u1v0w1 +
              C.VX1100*u0v0w1 +
              C.VX0010*u1v1w0 +
              C.VX1010*u0v1w0 +
              C.VX0110*u1v0w0 +
              C.VX1110*u0v0w0 
           )* C.iD * C.DT;
  
    v[1] = (  C.VY0000*u1v1w1 +
              C.VY1000*u0v1w1 +
              C.VY0100*u1v0w1 +
              C.VY1100*u0v0w1 +
              C.VY0010*u1v1w0 +
              C.VY1010*u0v1w0 +
              C.VY0110*u1v0w0 +
              C.VY1110*u0v0w0 
           )* C.iD * C.DT;
    
    v[2] = (  C.VZ0000*u1v1w1 +
              C.VZ1000*u0v1w1 +
              C.VZ0100*u1v0w1 +
              C.VZ1100*u0v0w1 +
              C.VZ0010*u1v1w0 +
              C.VZ1010*u0v1w0 +
              C.VZ0110*u1v0w0 +
              C.VZ1110*u0v0w0 
           )* C.iD * C.DT;

  }

  return;
  #undef  x
  #undef  y
  #undef  z
  #undef  t
}





_inline void getDV( real *xyzt , real *Dv ){
  real  u0, u1, v0, v1, w0, w1, s0, s1;
  
  #define x xyzt[0]
  #define y xyzt[1]
  #define z xyzt[2]
  #define t xyzt[3]
  if( 
      ( x < C.X0 )  ||  ( x > C.X1 )   ||
      ( y < C.Y0 )  ||  ( y > C.Y1 )   ||
      ( z < C.Z0 )  ||  ( z > C.Z1 )   ||
      ( t < C.T0 )  ||  ( t > C.T1 )
  ){
    setCELL( x , y , z , t );
  }
  
  if( C.isOutside ){
    Dv[0] = Dv[1] = Dv[2] = Dv[3] = Dv[4] = Dv[5] = Dv[6] = Dv[7] = Dv[8] = 0.0;
    return;
  }
  
  u0 = x - C.X0;   u1 = C.X1 - x;
  v0 = y - C.Y0;   v1 = C.Y1 - y;
  w0 = z - C.Z0;   w1 = C.Z1 - z;
  
  if( tinterp_mode == LINEAL ){
    s0 = t - C.T0;   s1 = C.T1 - t;

    Dv[0] = ( - C.VX0000*v1*w1*s1      - C.VX0001*v1*w1*s0
              + C.VX1000*v1*w1*s1      + C.VX1001*v1*w1*s0
              - C.VX0100*v0*w1*s1      - C.VX0101*v0*w1*s0
              + C.VX1100*v0*w1*s1      + C.VX1101*v0*w1*s0
              - C.VX0010*v1*w0*s1      - C.VX0011*v1*w0*s0
              + C.VX1010*v1*w0*s1      + C.VX1011*v1*w0*s0
              - C.VX0110*v0*w0*s1      - C.VX0111*v0*w0*s0
              + C.VX1110*v0*w0*s1      + C.VX1111*v0*w0*s0
            )* C.iD;

    Dv[1] = ( - C.VY0000*v1*w1*s1      - C.VY0001*v1*w1*s0
              + C.VY1000*v1*w1*s1      + C.VY1001*v1*w1*s0
              - C.VY0100*v0*w1*s1      - C.VY0101*v0*w1*s0
              + C.VY1100*v0*w1*s1      + C.VY1101*v0*w1*s0
              - C.VY0010*v1*w0*s1      - C.VY0011*v1*w0*s0
              + C.VY1010*v1*w0*s1      + C.VY1011*v1*w0*s0
              - C.VY0110*v0*w0*s1      - C.VY0111*v0*w0*s0
              + C.VY1110*v0*w0*s1      + C.VY1111*v0*w0*s0
            )* C.iD;

    Dv[2] = ( - C.VZ0000*v1*w1*s1      - C.VZ0001*v1*w1*s0
              + C.VZ1000*v1*w1*s1      + C.VZ1001*v1*w1*s0
              - C.VZ0100*v0*w1*s1      - C.VZ0101*v0*w1*s0
              + C.VZ1100*v0*w1*s1      + C.VZ1101*v0*w1*s0
              - C.VZ0010*v1*w0*s1      - C.VZ0011*v1*w0*s0
              + C.VZ1010*v1*w0*s1      + C.VZ1011*v1*w0*s0
              - C.VZ0110*v0*w0*s1      - C.VZ0111*v0*w0*s0
              + C.VZ1110*v0*w0*s1      + C.VZ1111*v0*w0*s0
            )* C.iD;


    Dv[3] = ( - C.VX0000*u1*w1*s1      - C.VX0001*u1*w1*s0
              - C.VX1000*u0*w1*s1      - C.VX1001*u0*w1*s0
              + C.VX0100*u1*w1*s1      + C.VX0101*u1*w1*s0
              + C.VX1100*u0*w1*s1      + C.VX1101*u0*w1*s0
              - C.VX0010*u1*w0*s1      - C.VX0011*u1*w0*s0
              - C.VX1010*u0*w0*s1      - C.VX1011*u0*w0*s0
              + C.VX0110*u1*w0*s1      + C.VX0111*u1*w0*s0
              + C.VX1110*u0*w0*s1      + C.VX1111*u0*w0*s0
            )* C.iD;

    Dv[4] = ( - C.VY0000*u1*w1*s1      - C.VY0001*u1*w1*s0
              - C.VY1000*u0*w1*s1      - C.VY1001*u0*w1*s0
              + C.VY0100*u1*w1*s1      + C.VY0101*u1*w1*s0
              + C.VY1100*u0*w1*s1      + C.VY1101*u0*w1*s0
              - C.VY0010*u1*w0*s1      - C.VY0011*u1*w0*s0
              - C.VY1010*u0*w0*s1      - C.VY1011*u0*w0*s0
              + C.VY0110*u1*w0*s1      + C.VY0111*u1*w0*s0
              + C.VY1110*u0*w0*s1      + C.VY1111*u0*w0*s0
            )* C.iD;

    Dv[5] = ( - C.VZ0000*u1*w1*s1      - C.VZ0001*u1*w1*s0
              - C.VZ1000*u0*w1*s1      - C.VZ1001*u0*w1*s0
              + C.VZ0100*u1*w1*s1      + C.VZ0101*u1*w1*s0
              + C.VZ1100*u0*w1*s1      + C.VZ1101*u0*w1*s0
              - C.VZ0010*u1*w0*s1      - C.VZ0011*u1*w0*s0
              - C.VZ1010*u0*w0*s1      - C.VZ1011*u0*w0*s0
              + C.VZ0110*u1*w0*s1      + C.VZ0111*u1*w0*s0
              + C.VZ1110*u0*w0*s1      + C.VZ1111*u0*w0*s0
            )* C.iD;

    Dv[6] = ( - C.VX0000*u1*v1*s1      - C.VX0001*u1*v1*s0
              - C.VX1000*u0*v1*s1      - C.VX1001*u0*v1*s0
              - C.VX0100*u1*v0*s1      - C.VX0101*u1*v0*s0
              - C.VX1100*u0*v0*s1      - C.VX1101*u0*v0*s0
              + C.VX0010*u1*v1*s1      + C.VX0011*u1*v1*s0
              + C.VX1010*u0*v1*s1      + C.VX1011*u0*v1*s0
              + C.VX0110*u1*v0*s1      + C.VX0111*u1*v0*s0
              + C.VX1110*u0*v0*s1      + C.VX1111*u0*v0*s0
             )* C.iD;

    Dv[7] = ( - C.VY0000*u1*v1*s1      - C.VY0001*u1*v1*s0
              - C.VY1000*u0*v1*s1      - C.VY1001*u0*v1*s0
              - C.VY0100*u1*v0*s1      - C.VY0101*u1*v0*s0
              - C.VY1100*u0*v0*s1      - C.VY1101*u0*v0*s0
              + C.VY0010*u1*v1*s1      + C.VY0011*u1*v1*s0
              + C.VY1010*u0*v1*s1      + C.VY1011*u0*v1*s0
              + C.VY0110*u1*v0*s1      + C.VY0111*u1*v0*s0
              + C.VY1110*u0*v0*s1      + C.VY1111*u0*v0*s0
             )* C.iD;

    Dv[8] = ( - C.VZ0000*u1*v1*s1      - C.VZ0001*u1*v1*s0
              - C.VZ1000*u0*v1*s1      - C.VZ1001*u0*v1*s0
              - C.VZ0100*u1*v0*s1      - C.VZ0101*u1*v0*s0
              - C.VZ1100*u0*v0*s1      - C.VZ1101*u0*v0*s0
              + C.VZ0010*u1*v1*s1      + C.VZ0011*u1*v1*s0
              + C.VZ1010*u0*v1*s1      + C.VZ1011*u0*v1*s0
              + C.VZ0110*u1*v0*s1      + C.VZ0111*u1*v0*s0
              + C.VZ1110*u0*v0*s1      + C.VZ1111*u0*v0*s0
             )* C.iD;    

  } else {
    
    Dv[0] = ( - C.VX0000*v1*w1
              + C.VX1000*v1*w1
              - C.VX0100*v0*w1
              + C.VX1100*v0*w1
              - C.VX0010*v1*w0
              + C.VX1010*v1*w0
              - C.VX0110*v0*w0
              + C.VX1110*v0*w0
            )* C.iD * C.DT;

    Dv[1] = ( - C.VY0000*v1*w1
              + C.VY1000*v1*w1
              - C.VY0100*v0*w1
              + C.VY1100*v0*w1
              - C.VY0010*v1*w0
              + C.VY1010*v1*w0
              - C.VY0110*v0*w0
              + C.VY1110*v0*w0
            )* C.iD * C.DT;

    Dv[2] = ( - C.VZ0000*v1*w1
              + C.VZ1000*v1*w1
              - C.VZ0100*v0*w1
              + C.VZ1100*v0*w1
              - C.VZ0010*v1*w0
              + C.VZ1010*v1*w0
              - C.VZ0110*v0*w0
              + C.VZ1110*v0*w0
            )* C.iD * C.DT;


    Dv[3] = ( - C.VX0000*u1*w1
              - C.VX1000*u0*w1
              + C.VX0100*u1*w1
              + C.VX1100*u0*w1
              - C.VX0010*u1*w0
              - C.VX1010*u0*w0
              + C.VX0110*u1*w0
              + C.VX1110*u0*w0
            )* C.iD * C.DT;

    Dv[4] = ( - C.VY0000*u1*w1
              - C.VY1000*u0*w1
              + C.VY0100*u1*w1
              + C.VY1100*u0*w1
              - C.VY0010*u1*w0
              - C.VY1010*u0*w0
              + C.VY0110*u1*w0
              + C.VY1110*u0*w0
            )* C.iD * C.DT;

    Dv[5] = ( - C.VZ0000*u1*w1
              - C.VZ1000*u0*w1
              + C.VZ0100*u1*w1
              + C.VZ1100*u0*w1
              - C.VZ0010*u1*w0
              - C.VZ1010*u0*w0
              + C.VZ0110*u1*w0
              + C.VZ1110*u0*w0
            )* C.iD * C.DT;

    Dv[6] = ( - C.VX0000*u1*v1
              - C.VX1000*u0*v1
              - C.VX0100*u1*v0
              - C.VX1100*u0*v0
              + C.VX0010*u1*v1
              + C.VX1010*u0*v1
              + C.VX0110*u1*v0
              + C.VX1110*u0*v0
             )* C.iD * C.DT;

    Dv[7] = ( - C.VY0000*u1*v1
              - C.VY1000*u0*v1
              - C.VY0100*u1*v0
              - C.VY1100*u0*v0
              + C.VY0010*u1*v1
              + C.VY1010*u0*v1
              + C.VY0110*u1*v0
              + C.VY1110*u0*v0
             )* C.iD * C.DT;

    Dv[8] = ( - C.VZ0000*u1*v1
              - C.VZ1000*u0*v1
              - C.VZ0100*u1*v0
              - C.VZ1100*u0*v0
              + C.VZ0010*u1*v1
              + C.VZ1010*u0*v1
              + C.VZ0110*u1*v0
              + C.VZ1110*u0*v0
             )* C.iD * C.DT;
    
  }
    
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
