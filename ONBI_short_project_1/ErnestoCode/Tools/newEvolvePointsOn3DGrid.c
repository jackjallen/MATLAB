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

#include "newmyMEX.h" /* gav */
#include "newlib.h"

#if !defined( memcpy )
  #include "string.h"
#endif


#define    integration_method_def   EULER
/* #define    integration_method_def   RK23 */
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


#define EVOLVE_JAC          if( nlhs > 1 ){                                                  \
                              newgetDV( XYZT , Dv );                                         \
                              switch( jacobian_method ){                                     \
                                case PLUS:                                                   \
                                  Mtimes3x3( dJAC , Dv , JAC );                              \
                                  for( d = 0 ; d<9 ; d++ ){  JAC[d] += h*dJAC[d]; }          \
                                  /* print_matrixt("JAC in EVOLVE_JAC", JAC, 3, 3); */       \
                                  break;                                                     \
                                case TIMES:                                                  \
                                  for( d = 0 ; d<9 ; d++ ){ Dv[d] *= h; }                    \
                                  MatrixExp3x3( Dv );                                        \
                                  Mtimes3x3( dJAC , Dv , JAC );                              \
                                  memcpy( JAC , dJAC , 9*sizeof( real ) );                   \
                                  break;                                                     \
                              }                                                              \
                            }

/* #define MAT(i,j)      MAT[ 4*(j) + (i) ] */


/* #define     readBOUNDS                                                            \ */
/*       mexPrintf("\n ENTRO EN readBOUNDS!\n");                       		  \ */
/*       if( nrhs > argN && ! mxIsChar(prhs[argN]) ){                                \ */
/*         if( !myIsEmpty( prhs[argN] ) ){                                           \ */
/*           if( myNumel( prhs[argN] ) == 1 ){                                       \ */
/*             boundary_size[0] = myGetValue(prhs[argN]);                            \ */
/*             boundary_size[1] = boundary_size[0];                                  \ */
/*             boundary_size[2] = boundary_size[0];                                  \ */
/*             boundary_size[3] = boundary_size[0];                                  \ */
/*             boundary_size[4] = boundary_size[0];                                  \ */
/*             boundary_size[5] = boundary_size[0];                                  \ */
/*           } else if( myNumel( prhs[argN] ) == 6 ){                                \ */
/*             boundary_size[0] = myGetValueIth( prhs[argN] , 0 );                   \ */
/*             boundary_size[1] = myGetValueIth( prhs[argN] , 1 );                   \ */
/*             boundary_size[2] = myGetValueIth( prhs[argN] , 2 );                   \ */
/*             boundary_size[3] = myGetValueIth( prhs[argN] , 3 );                   \ */
/*             boundary_size[4] = myGetValueIth( prhs[argN] , 4 );                   \ */
/*             boundary_size[5] = myGetValueIth( prhs[argN] , 5 );                   \ */
/*           } else {                                                                \ */
/*             myErrMsgTxt("1 or 6 values expected as BoundarySize");                \ */
/*           }                                                                       \ */
/*         }                                                                         \ */
/*         argN++; continue;                                                         \ */
/*       } */


real    boundary_size[6];
real    *bs;
int     ibs;
/* gav */
#define     newreadBOUNDS						\
  if( nrhs > argN && ! mxIsChar(prhs[argN]) ){				\
    if( !myIsEmpty( prhs[argN] ) ){					\
      if( myNumel( prhs[argN] ) == 1 ){					\
	bs = (real *) mxGetData(prhs[argN]);				\
	for(ibs=0; ibs<6; ibs++){					\
	  boundary_size[ibs] = bs[0];					\
	}								\
      }									\
      else if( myNumel( prhs[argN] ) == 6 ){				\
	bs = (real *) mxGetData(prhs[argN]);				\
	for(ibs=0; ibs<6; ibs++){					\
	  boundary_size[ibs] = bs[ibs];					\
	}								\
      }									\
      else {								\
	myErrMsgTxt("1 or 6 values expected as BoundarySize");		\
      }									\
    }									\
    argN++; continue;							\
  }



#define  Ct45_1         0.0
#define  Ct45_2         0.20000000000000001
#define  Ct45_3         0.29999999999999999
#define  Ct45_4         0.80000000000000004
#define  Ct45_5         0.88888888888888884
#define  Ct45_6         1.0
#define  Ct45_7         1.0

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

#define  B45_1         0.091145833333333329
#define  B45_2         0.0
#define  B45_3         0.44923629829290207
#define  B45_4         0.65104166666666663
#define  B45_5        -0.322376179245283
#define  B45_6         0.13095238095238096
#define  B45_7         0.0

#define  D45_1         0.0012326388888888888
#define  D45_2         0.0
#define  D45_3        -0.0042527702905061394
#define  D45_4         0.036979166666666667
#define  D45_5        -0.05086379716981132
#define  D45_6         0.041904761904761903
#define  D45_7        -0.025000000000000001

#define  Ct23_1         0.0
#define  Ct23_2         0.5
#define  Ct23_3         0.75
#define  Ct23_4         1.0

#define  A23_21   0.5
#define  A23_32   0.75
#define  A23_41   0.22222222222222221
#define  A23_42   0.33333333333333331
#define  A23_43   0.44444444444444442

#define  B23_1   0.22222222222222221
#define  B23_2   0.33333333333333331
#define  B23_3   0.44444444444444442

#define  D23_1   -0.069444444444444448
#define  D23_2   0.083333333333333329
#define  D23_3   0.1111111111111111
#define  D23_4   -0.125



enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST } boundary_mode;
real    *VX, *VY, *VZ, *X, *Y, *Z, OX, OY, OZ, LX, LY, LZ;
int     I , J , K , IJ, IJK;
int     NSD;

typedef struct CELL {
  int  ii, jj, kk;
  real X0, X1, Y0, Y1, Z0, Z1;
  /* real VX0000, VX1000, VX0100, VX1100, VX0010, VX1010, VX0110, VX1110; */
  /* real VY0000, VY1000, VY0100, VY1100, VY0010, VY1010, VY0110, VY1110; */
  /* real VZ0000, VZ1000, VZ0100, VZ1100, VZ0010, VZ1010, VZ0110, VZ1110; */
  real VXYZ[24]; /* gav */
  real DX, DY, DZ, D, iD;
  int  isOutside;
} CELL;
CELL C; /* gav */


/* gav : cambio inline para poder compilar */
/* /\* _inline  *\/ void setCELL( real x , real y , real z ); */
            void newsetCELL(real, real, real);
/* /\* _inline  *\/ void getV( real *xyzt , real *v ); */
            void newgetV(real *, real *);
/* /\* _inline  *\/ void getDV( real *xyzt , real *Dv ); */
            void newgetDV(real *, real *);
void MatrixExp3x3( real * );
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
  ALLOCATES();
  enum    integration_methods  { EULER , S_EULER , RK45 , RK23 }    integration_method;
  enum    jacobian_methods     { PLUS   , TIMES }                   jacobian_method;

  real    XYZT[4], xk[4], v[3], vp[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], err[3], Dv[9];
  real    xyztnr[4]; /* gav */
  real    xxnr, yynr, zznr;
  real    maxV, minD, hmax, hmin, tol, h;
  real    minStepsPerVoxel, relTOL;

  real    prct , last_prct=0;
  unsigned char         *MASK, VERBOSE;
  real    JAC[9], dJAC[9];
  
  real    *XYZ, *F, *Xn, *Yn, *Zn;
  int     Odims[6], ndims, d, nP, p, nF, ff, in, jn, kn, In, Jn, Kn; /* Odims[50] */  /* gav */

  int     maxSTEPS, Nadvice = 0;
  /* real    MAT[16]; */
  real    *mat; /* gav */
  
  
  real    *O;
  char    STR[128]; /* gav */
  int     argN;
  
  real    *initJAC, *OJAC, *DET;

  int ir; /* gav */
  int flagmat; /* if 1: grilla, read mat  */ /* gav */

  double nan, inf;

  flagmat = 0;

  nan = mxGetNaN();
  inf = mxGetInf();
  
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

    mexPrintf("prhs[3] empty\n"); /* gav */

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

  
  if( mxIsCell( prhs[4] ) ){ /* grilla */
    
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
    VZ = NULL;      /* gav */
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
      boundary_mode = CIRCULAR; 
      argN++;  
      newreadBOUNDS; /* gav */
      continue;
    }

    if( ! myStrcmpi( STR,"decay_to_zero") ||
        ! myStrcmpi( STR,"tozero")  || ! myStrcmpi( STR,"decay")   ) {
      boundary_mode = DECAY_TO_ZERO; argN++;  
      newreadBOUNDS; /* gav */
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

    
    if( ! myStrcmpi( STR,"matrix") || ! myStrcmpi( STR,"mat") ){
      argN++;
      
      if( nrhs > argN && !mxIsChar( prhs[argN] ) &&  mySize( prhs[argN] , 0 ) == 4  &&  mySize( prhs[argN] , 1 ) == 4  && myNDims( prhs[argN] ) == 2 ){
        if( !myIsEmpty( prhs[argN] ) ){
	  
	  mat  = (real *) mxGetData(prhs[argN]); /* gav */ /* ojo elementos casteados?! */
	  flagmat = 1;

          if( mat[3] != 0  ||  mat[7] != 0  ||  mat[11] != 0  ||  mat[15] != 1 ){ /* gav */
            myErrMsgTxt("The matrix must be an homogeneous 4x4 matrix.");
          }

	  /* print_matrixt("mat", mat, 4, 4); */
	  /* diff_matrix(diffm, MAT, mat, 4, 4); */
	  myFlush();

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

            myErrMsgTxt("After the word initJACs  nPxNSDxNSD array is expected (1)");
          }

          initJAC = myGetPr( prhs[argN] );

        }
        argN++; continue;
      }
      myErrMsgTxt("After the word initJACs  nPxNSDxNSD array is expected (2)");
    }

    mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
  } /*  while( nrhs > argN )  */
  /*END Parsing arguments*/
  
  if( F == NULL ){
    nF = 2;
    F = mxMalloc( 2*sizeof( real ) );
    F[0] = 0;
    F[1] = 1;
  }
  if( !checkIsSorted(F,nF) ){ myErrMsgTxt("FLOW times are not sorted."); }

  if(flagmat == 1){ /* if MAT read */ /* if(mat[15] != 0) { */
    checkm4x4ish("mat", mat, NSD); /* gav */
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

  mexPrintf("0\n");
 
  /*Creating output*/
  /* plhs[0] = mxCreateNumericArray( ndims+1 , Odims , mxREAL_CLASS , mxREAL ); */
  /* plhs[0] = myCreateNumericArray_E( ndims+1 , Odims , mxREAL_CLASS , mxREAL ); */
  plhs[0] = myCreateNumericArray_E( ndims+1 , Odims , mxREAL_CLASS , mxREAL );
  if(plhs[0] == NULL)
    mexErrMsgTxt("Error: plhs[0] cannot be created. Exiting ...");
  O = (real *) mxGetData( plhs[0] );

  if( nlhs > 1 ){
    plhs[1] = myCreateNumericArray_E( ndims+2 , Odims , mxREAL_CLASS , mxREAL );
    if(plhs[1] == NULL)
      mexErrMsgTxt("Error: plhs[1] cannot be created. Exiting ...");
    OJAC = (real *) mxGetData( plhs[1] );
  }

  if( nlhs > 2 ){
    plhs[2] = myCreateNumericArray_E( ndims , Odims , mxREAL_CLASS , mxREAL );
    if(plhs[2] == NULL)
      mexErrMsgTxt("Error: plhs[2] cannot be created. Exiting ...");
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

  
    minD = 1.0;
    if( I > 1 ){ minD = MIN( minD , (LX-OX)/I ); }
    if( J > 1 ){ minD = MIN( minD , (LY-OY)/J ); }
    if( K > 1 ){ minD = MIN( minD , (LZ-OZ)/K ); }
    tol  = minD*relTOL;
  
  hmin     = ( F[nF-1] - F[0] )/maxSTEPS;
  minStepsPerVoxel = 1.0/minStepsPerVoxel;

  /* // reset the CELL */
  C.ii = C.jj = C.kk = -1;
  C.X0 = C.Y0 = C.Z0 =  0;
  C.X1 = C.Y1 = C.Z1 = -1;
  /* C.VX0000 = C.VX1000 = C.VX0100 = C.VX1100 = C.VX0010 = C.VX1010 = C.VX0110 = C.VX1110 = 0; */
  /* C.VY0000 = C.VY1000 = C.VY0100 = C.VY1100 = C.VY0010 = C.VY1010 = C.VY0110 = C.VY1110 = 0; */
  /* C.VZ0000 = C.VZ1000 = C.VZ0100 = C.VZ1100 = C.VZ0010 = C.VZ1010 = C.VZ0110 = C.VZ1110 = 0; */

  for(ir=0; ir<24; ir++){ /* gav */
    C.VXYZ[ir] = 0;
  }

  C.DX = C.DY = C.DZ = 1;
  C.D = C.iD = 1;
  C.isOutside = 1;
  


  for( kn=0, p=0 ; kn < Kn ; kn++      ){
  for( jn=0      ; jn < Jn ; jn++      ){
  for( in=0      ; in < In ; in++, p++ ){

    /* gav */
    if( utIsInterruptPending() ){
      mexPrintf("USER INTERRUP in  EvolvePointsOn3DGrid  !!!\n");
      myErrMsgTxt("USER INTERRUP in  EvolvePointsOn3DGrid  !!!");
    }
  
    if( VERBOSE ){
      prct = ( (double)(p+1) )/((double)nP)*100;
      if( prct >= last_prct ){
        if( last_prct != 0 ){
          mexPrintf("%3d%% done  ( %10d   of   %10d  points )\n" , (int) prct , p+1, nP ); myFlush();
        }
        last_prct += 10;
      }
    }

    /* inicializaciones */ /* gav */
    /* ------- */
    XYZT[0] = XYZT[1] = XYZT[2] = XYZT[3] = 0; 
    v[0] = v[1] = v[2] = 0;
    for(ir=0; ir<3; ir++){
      xyztnr[ir] = 0.;
    }
    xyztnr[3] = 1.; /* se podria inicializar solo una vez, por seguridad y claridad mejor dejarlo aqui */
    
    for(ir=0; ir<9; ir++)
      JAC[ir] = 0.;
    /* ------- */
    
    /* OUT and MAT  ------- */ /* gav */
    if( MASK != NULL  &&  MASK[p] != 1 ){
      xyztnr[0] = nan;
    } else if( Xn == NULL ){
                      xyztnr[0] = XYZ[ p        ];
                      xyztnr[1] = XYZ[ p +   nP ];
      if( NSD == 3 ){ xyztnr[2] = XYZ[ p + 2*nP ]; }
    } else {
                      xyztnr[0] = Xn[ in ];
                      xyztnr[1] = Yn[ jn ];
      if( NSD == 3 ){ xyztnr[2] = Zn[ kn ]; }
    }
    
    /* /\* gav *\/ */
    /* if( myISINF( xyztnr[0] ) || myISINF( xyztnr[1] ) || myISINF( xyztnr[2] ) ){ */
    /*   mexPrintf("xyztnr[0]: %lf xyztnr[1]: %lf xyztnr[2]:%lf\n\n",  */
    /* 		xyztnr[0], xyztnr[1], xyztnr[2]); */
    /*   myErrMsgTxt("Error ???: alguno es Inf !"); */
    /* } */
    
    if( myISNAN( xyztnr[0] ) || myISNAN( xyztnr[1] ) || myISNAN( xyztnr[2] ) ){
      for( ff = 1 ; ff < nF ; ff++ ){

	                *(O  + p + ( (ff-1)            )*nP ) = nan;
	                *(O  + p + ( (ff-1) +   (nF-1) )*nP ) = nan;
	if( NSD == 3 ){ *(O  + p + ( (ff-1) + 2*(nF-1) )*nP ) = nan; }
        
        if( nlhs > 1 ){
              *(OJAC +  p + ( (ff-1)            )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) +   (nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 2*(nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 3*(nF-1) )*nP ) = nan;
	      
	    if( NSD == 3 ){
	      *(OJAC +  p + ( (ff-1) + 4*(nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 5*(nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 6*(nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 7*(nF-1) )*nP ) = nan;
	      *(OJAC +  p + ( (ff-1) + 8*(nF-1) )*nP ) = nan;
	    }
        }
                        
        if( nlhs > 2 ){
          *(DET  + p + (ff-1)*nP ) = nan;
        }
        
      }
      continue;
    }

    if( flagmat == 1){ /* if mat read */
      
      prodm(XYZT, xyztnr, mat, 1, 4, 4);
    }
    else{ /* aqui no pueden llegar nan's porque esta el continue anterior! */
      XYZT[0] = xyztnr[0];
      XYZT[1] = xyztnr[1];
      XYZT[2] = xyztnr[2];
    }
    XYZT[3] = F[0]; /* ha de estar fijado a F[0] */
    /*  ------- OUT and MAT */ /* gav */
    
    if( nlhs > 1 ){
      if( initJAC != NULL ){
        if( NSD == 3 ){

	  for(ir=0; ir<9; ir++){
	    *(JAC + ir) = *(initJAC + p + ir*nP );
	  }
          /* JAC[0] = initJAC[ p        ]; */
          /* JAC[1] = initJAC[ p +   nP ]; */
          /* JAC[2] = initJAC[ p + 2*nP ]; */
          /* JAC[3] = initJAC[ p + 3*nP ]; */
          /* JAC[4] = initJAC[ p + 4*nP ]; */
          /* JAC[5] = initJAC[ p + 5*nP ]; */
          /* JAC[6] = initJAC[ p + 6*nP ]; */
          /* JAC[7] = initJAC[ p + 7*nP ]; */
          /* JAC[8] = initJAC[ p + 8*nP ]; */
        }
    	else {
          *(JAC + 0) = *(initJAC +  p        );
          *(JAC + 1) = *(initJAC +  p +   nP );
	  *(JAC + 2) = 0.0;
          *(JAC + 3) = *(initJAC +  p + 2*nP );
          *(JAC + 4) = *(initJAC +  p + 3*nP );
          *(JAC + 5) = *(JAC + 6) = *(JAC + 7) = 0.0;
          *(JAC + 8) = 1.0;
        }
      }
      else {
        *(JAC + 0) = 1.0;
	*(JAC + 1) = *(JAC + 2) = *(JAC + 3) = 0.0;
	*(JAC + 4) = 1.0;
	*(JAC + 5) = *(JAC + 6) = *(JAC + 7) = 0.0;
	*(JAC + 8) = 1.0;
      }
    }

    /* print_matrixt("JAC", JAC, 3, 3); */
    
    for( ff = 1 ; ff < nF ; ff++ ){

      switch( integration_method ){
        case EULER:

          h    = hmin;
          while( XYZT[3] < F[ff] ){
	    newgetV( XYZT , v ); /* gav */
/* //       mexPrintf("v:  %g,%g,%g\n" , v[0],v[1],v[2]); */
	      
            if( !v[0] && !v[1] && !v[2] ){ break; }
            if( h > ( F[ff] - XYZT[3] ) ){ h = F[ff] - XYZT[3]; }
	    
	    EVOLVE_JAC
	      
            XYZT[0] += h*v[0];
            XYZT[1] += h*v[1];
            XYZT[2] += h*v[2];
            XYZT[3] += h;
	      
/* //       mexPrintf("XYZT:  %g,%g,%g,%g\n" , XYZT[0],XYZT[1],XYZT[2],XYZT[3]); */            
          }
          break;

        case S_EULER:

          h = 1;
          while( XYZT[3] < F[ff] ){
            newgetV( XYZT , v ); /* gav */

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

          newgetV( XYZT , v );
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
            newgetV( xk , k2 );
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );


            xk[0] = XYZT[0] + A45_31*k1[0] + A45_32*k2[0];
            xk[1] = XYZT[1] + A45_31*k1[1] + A45_32*k2[1];
            xk[2] = XYZT[2] + A45_31*k1[2] + A45_32*k2[2];
            newgetV( xk , k3 );
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k3 , k3 , h );

            xk[0] = XYZT[0] + A45_41*k1[0] + A45_42*k2[0] + A45_43*k3[0];
            xk[1] = XYZT[1] + A45_41*k1[1] + A45_42*k2[1] + A45_43*k3[1];
            xk[2] = XYZT[2] + A45_41*k1[2] + A45_42*k2[2] + A45_43*k3[2];
            newgetV( xk , k4 );
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k4 , k4 , h );

            xk[0] = XYZT[0] + A45_51*k1[0] + A45_52*k2[0] + A45_53*k3[0] + A45_54*k4[0];
            xk[1] = XYZT[1] + A45_51*k1[1] + A45_52*k2[1] + A45_53*k3[1] + A45_54*k4[1];
            xk[2] = XYZT[2] + A45_51*k1[2] + A45_52*k2[2] + A45_53*k3[2] + A45_54*k4[2];
            newgetV( xk , k5 );
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k5 , k5 , h );

            xk[0] = XYZT[0] + A45_61*k1[0] + A45_62*k2[0] + A45_63*k3[0] + A45_64*k4[0] + A45_65*k5[0];
            xk[1] = XYZT[1] + A45_61*k1[1] + A45_62*k2[1] + A45_63*k3[1] + A45_64*k4[1] + A45_65*k5[1];
            xk[2] = XYZT[2] + A45_61*k1[2] + A45_62*k2[2] + A45_63*k3[2] + A45_64*k4[2] + A45_65*k5[2];
            newgetV( xk , k6 );
            if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k6 , k6 , h );

            memcpy( vp , v , 3*sizeof( real ) );

            xk[0] = XYZT[0] + B45_1*k1[0] + B45_3*k3[0] + B45_4*k4[0] + B45_5*k5[0] + B45_6*k6[0];
            xk[1] = XYZT[1] + B45_1*k1[1] + B45_3*k3[1] + B45_4*k4[1] + B45_5*k5[1] + B45_6*k6[1];
            xk[2] = XYZT[2] + B45_1*k1[2] + B45_3*k3[2] + B45_4*k4[2] + B45_5*k5[2] + B45_6*k6[2];
            xk[3] = XYZT[3] + h;
            newgetV( xk , v );
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

          newgetV( XYZT , v );
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
            newgetV( xk , k2 );  if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k2 , k2 , h );

            xk[0] = XYZT[0] + A23_32*k2[0];
            xk[1] = XYZT[1] + A23_32*k2[1];
            xk[2] = XYZT[2] + A23_32*k2[2];
            newgetV( xk , k3 );  if( C.isOutside  &&  h > hmin ){ h = hmin; continue; }
            times( k3 , k3 , h );

            memcpy( vp , v , 3*sizeof( real ) );

            xk[0] = XYZT[0] + B23_1*k1[0] + B23_2*k2[0] + B23_3*k3[0];
            xk[1] = XYZT[1] + B23_1*k1[1] + B23_2*k2[1] + B23_3*k3[1];
            xk[2] = XYZT[2] + B23_1*k1[2] + B23_2*k2[2] + B23_3*k3[2];
            xk[3] = XYZT[3] + h;
            newgetV( xk , v );
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

        *(O +  p + ( (ff-1)            )*nP ) = XYZT[0];
        *(O +  p + ( (ff-1) +   (nF-1) )*nP ) = XYZT[1];
        *(O +  p + ( (ff-1) + 2*(nF-1) )*nP ) = XYZT[2];

        if( nlhs > 1 ){
          *(OJAC +  p + ( (ff-1)            )*nP ) = *(JAC + 0);
          *(OJAC +  p + ( (ff-1) +   (nF-1) )*nP ) = *(JAC + 1);
          *(OJAC +  p + ( (ff-1) + 2*(nF-1) )*nP ) = *(JAC + 2);
          *(OJAC +  p + ( (ff-1) + 3*(nF-1) )*nP ) = *(JAC + 3);
          *(OJAC +  p + ( (ff-1) + 4*(nF-1) )*nP ) = *(JAC + 4);
          *(OJAC +  p + ( (ff-1) + 5*(nF-1) )*nP ) = *(JAC + 5);
          *(OJAC +  p + ( (ff-1) + 6*(nF-1) )*nP ) = *(JAC + 6);
          *(OJAC +  p + ( (ff-1) + 7*(nF-1) )*nP ) = *(JAC + 7);
          *(OJAC +  p + ( (ff-1) + 8*(nF-1) )*nP ) = *(JAC + 8);
        }

        if( nlhs > 2 ){
          *(DET +  p + (ff-1)*nP ) =  det3x3(JAC); /* gav */
	  mexPrintf("det DET[%d]=%lf\n", p + (ff-1)*nP , DET[ p + (ff-1)*nP ]);
	
	  /* *(DET + p + (ff-1)*nP) = *(JAC + 0)*( *(JAC + 4) * *(JAC + 8)  - *(JAC + 5) * *(JAC + 7)  ) */
          /*                        - *(JAC + 1)*( *(JAC + 3) * *(JAC + 8)  - *(JAC + 5) * *(JAC + 6)  ) */
          /*                        + *(JAC + 2)*( *(JAC + 3) * *(JAC + 7)  - *(JAC + 4) * *(JAC + 6)  ); */
	  /* mexPrintf("2: DET[%d]=%lf\n", p + (ff-1)*nP , DET[ p + (ff-1)*nP ]); */
	  /* mxErrMsgTxt("SALGO\n"); */
        }


      } else {

        *(O +  p + ( (ff-1)            )*nP ) = XYZT[0];
        *(O +  p + ( (ff-1) +   (nF-1) )*nP ) = XYZT[1];

        if( nlhs > 1 ){
          *(OJAC +  p + ( (ff-1)            )*nP ) = *(JAC + 0);
          *(OJAC +  p + ( (ff-1) +   (nF-1) )*nP ) = *(JAC + 1);
          *(OJAC +  p + ( (ff-1) + 2*(nF-1) )*nP ) = *(JAC + 3);
          *(OJAC +  p + ( (ff-1) + 3*(nF-1) )*nP ) = *(JAC + 4);
        }

        if( nlhs > 2 ){
          *(DET +  p + (ff-1)*nP ) = *(JAC + 0) * *(JAC + 4) - *(JAC + 3) * *(JAC + 1);
        }
        
      }

    } /* // for in ff */
  }}} /* //for kn, jn, in */

  if( Nadvice > 0 ){
    sprintf( STR , "Adviced maxSTEPS: %d" , Nadvice );
    mexWarnMsgIdAndTxt( "EXP:V:IncreaseMaxSteps" , STR );
  }

  mexPrintf("\n\n");
  mexPrintf("inline functions       cambiado para poder compilar!\n");
  mexPrintf("utIsInterruptPending() cambiado para poder compilar!\n");
  mexPrintf("\n\n");



  EXIT: myFreeALLOCATES();

} /* mexFunction */



/* _inline  */void newsetCELL( real x , real y , real z ){
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

    /* dimVXY = 16; */ /* gav */
    /* for(ir=0; ir<dimVXY; ir++){ */

    /*   C.VXYZ[ir] = ; /\* no se puede hacer asi por los indices! ii[01], jj[01], kk[01] *\/ */
    /* se podria tratar de jugar con los indices ... */
    /* } */

    C.VXYZ[ 0] = VV(  ii0  ,  jj0  ,  kk0  );
    C.VXYZ[ 1] = VV(  ii1  ,  jj0  ,  kk0  );
    C.VXYZ[ 2] = VV(  ii0  ,  jj1  ,  kk0  );
    C.VXYZ[ 3] = VV(  ii1  ,  jj1  ,  kk0  );
    C.VXYZ[ 4] = VV(  ii0  ,  jj0  ,  kk1  );
    C.VXYZ[ 5] = VV(  ii1  ,  jj0  ,  kk1  );
    C.VXYZ[ 6] = VV(  ii0  ,  jj1  ,  kk1  );
    C.VXYZ[ 7] = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    #define  VV(i,j,k)   VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VXYZ[ 8] = VV(  ii0  ,  jj0  ,  kk0  );
    C.VXYZ[ 9] = VV(  ii1  ,  jj0  ,  kk0  );
    C.VXYZ[10] = VV(  ii0  ,  jj1  ,  kk0  );
    C.VXYZ[11] = VV(  ii1  ,  jj1  ,  kk0  );
    C.VXYZ[12] = VV(  ii0  ,  jj0  ,  kk1  );
    C.VXYZ[13] = VV(  ii1  ,  jj0  ,  kk1  );
    C.VXYZ[14] = VV(  ii0  ,  jj1  ,  kk1  );
    C.VXYZ[15] = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VXYZ[16] = VV(  ii0  ,  jj0  ,  kk0  );
      C.VXYZ[17] = VV(  ii1  ,  jj0  ,  kk0  );
      C.VXYZ[18] = VV(  ii0  ,  jj1  ,  kk0  );
      C.VXYZ[19] = VV(  ii1  ,  jj1  ,  kk0  );
      C.VXYZ[20] = VV(  ii0  ,  jj0  ,  kk1  );
      C.VXYZ[21] = VV(  ii1  ,  jj0  ,  kk1  );
      C.VXYZ[22] = VV(  ii0  ,  jj1  ,  kk1  );
      C.VXYZ[23] = VV(  ii1  ,  jj1  ,  kk1  );
      #undef VV
    }

  } else {
    
    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VX[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VXYZ[ 0] = VV(  ii0  ,  jj0  ,  kk0  );
    C.VXYZ[ 1] = VV(  ii1  ,  jj0  ,  kk0  );
    C.VXYZ[ 2] = VV(  ii0  ,  jj1  ,  kk0  );
    C.VXYZ[ 3] = VV(  ii1  ,  jj1  ,  kk0  );
    C.VXYZ[ 4] = VV(  ii0  ,  jj0  ,  kk1  );
    C.VXYZ[ 5] = VV(  ii1  ,  jj0  ,  kk1  );
    C.VXYZ[ 6] = VV(  ii0  ,  jj1  ,  kk1  );
    C.VXYZ[ 7] = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VXYZ[ 8] = VV(  ii0  ,  jj0  ,  kk0  );
    C.VXYZ[ 9] = VV(  ii1  ,  jj0  ,  kk0  );
    C.VXYZ[10] = VV(  ii0  ,  jj1  ,  kk0  );
    C.VXYZ[11] = VV(  ii1  ,  jj1  ,  kk0  );
    C.VXYZ[12] = VV(  ii0  ,  jj0  ,  kk1  );
    C.VXYZ[13] = VV(  ii1  ,  jj0  ,  kk1  );
    C.VXYZ[14] = VV(  ii0  ,  jj1  ,  kk1  );
    C.VXYZ[15] = VV(  ii1  ,  jj1  ,  kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VXYZ[16] = VV(  ii0  ,  jj0  ,  kk0  );
      C.VXYZ[17] = VV(  ii1  ,  jj0  ,  kk0  );
      C.VXYZ[18] = VV(  ii0  ,  jj1  ,  kk0  );
      C.VXYZ[19] = VV(  ii1  ,  jj1  ,  kk0  );
      C.VXYZ[20] = VV(  ii0  ,  jj0  ,  kk1  );
      C.VXYZ[21] = VV(  ii1  ,  jj0  ,  kk1  );
      C.VXYZ[22] = VV(  ii0  ,  jj1  ,  kk1  );
      C.VXYZ[23] = VV(  ii1  ,  jj1  ,  kk1  );
      #undef VV
    }

  }

  return;
}

/* _inline */ void newgetV( real *xyzt , real *v ){

  real  u0, u1, v0, v1, w0, w1;
  real uvw[8];
  
  if( 
      ( xyzt[0] < C.X0 )  ||  ( xyzt[0] > C.X1 )   ||
      ( xyzt[1] < C.Y0 )  ||  ( xyzt[1] > C.Y1 )   ||
      ( xyzt[2] < C.Z0 )  ||  ( xyzt[2] > C.Z1 )   
  ){
    newsetCELL( xyzt[0] , xyzt[1] , xyzt[2] );
/* //     mexPrintf("[C.X0,x,C.X1]:%g,%g,%g   [C.Y0,y,C.Y1]:%g,%g,%g   [C.Z0,z,C.Z1]:%g,%g,%g   \n",C.X0,x,C.X1,C.Y0,y,C.Y1,C.Z0,z,C.Z1); */
  }
  
  if( C.isOutside ){
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    return;
  }
  

  /* inicializa v para poder sacarla con interpv3.  Ver si no hay
     problemas al cambiarlo al final para inicializarlo dentro de la
     propia funcion interpv3 y en interpv2 y no aqui.
  */ /* gav */
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;

  
  if( NSD == 3 ){
    
    u0 = xyzt[0] - C.X0;   u1 = C.X1 - xyzt[0];
    v0 = xyzt[1] - C.Y0;   v1 = C.Y1 - xyzt[1];
    w0 = xyzt[2] - C.Z0;   w1 = C.Z1 - xyzt[2];


    uvw[0] = u1*v1*w1;
    uvw[1] = u0*v1*w1;
    uvw[2] = u1*v0*w1;
    uvw[3] = u0*v0*w1;
    uvw[4] = u1*v1*w0;
    uvw[5] = u0*v1*w0;
    uvw[6] = u1*v0*w0;
    uvw[7] = u0*v0*w0;
    
    interpv3(v, C.VXYZ, uvw, C.iD);

  } else { 
    
    u0 = xyzt[0] - C.X0;   u1 = C.X1 - xyzt[0];
    v0 = xyzt[1] - C.Y0;   v1 = C.Y1 - xyzt[1];

    uvw[0] = u1*v1;
    uvw[1] = u0*v1;
    uvw[2] = u1*v0;
    uvw[3] = u0*v0;
    uvw[4] = 0;
    uvw[5] = 0;
    uvw[6] = 0;
    uvw[7] = 0;

    /* interpv3(v, C.VXYZ, uvw, C.iD); */
    interpv2(v, C.VXYZ, uvw, C.iD); /* con esta nos evitamos tener que multiplicar y sumar 12 ceros
                                       pero a cambio tenermos que ir a buscar el indice  */

  }

  return;
}

/* _inline  */void newgetDV( real *xyzt , real *Dv ){
  real  u0, u1, v0, v1, w0, w1;

  real   v1w1,  u1w1,  u1v1;
  real   v0w1,  u0w1,  u0v1;
  real   v1w0,  u1w0,  u1v0;
  real   v0w0,  u0w0,  u0v0;
  
  if( 
      ( xyzt[0] < C.X0 )  ||  ( xyzt[0] > C.X1 )   ||
      ( xyzt[1] < C.Y0 )  ||  ( xyzt[1] > C.Y1 )   ||
      ( xyzt[2] < C.Z0 )  ||  ( xyzt[2] > C.Z1 )   
  ){
    newsetCELL( xyzt[0] , xyzt[1] , xyzt[2] );
  }
  
  if( C.isOutside ){
    Dv[0] = Dv[1] = Dv[2] = Dv[3] = Dv[4] = Dv[5] = Dv[6] = Dv[7] = Dv[8] = 0.0;
    return;
  }
  

  u0 = xyzt[0] - C.X0;   u1 = C.X1 - xyzt[0];
  v0 = xyzt[1] - C.Y0;   v1 = C.Y1 - xyzt[1];
  w0 = xyzt[2] - C.Z0;   w1 = C.Z1 - xyzt[2];
  
  v1w1 = v1*w1;  u1w1 = u1*w1;  u1v1 = u1*v1;
  v0w1 = v0*w1;  u0w1 = u0*w1;  u0v1 = u0*v1;
  v1w0 = v1*w0;  u1w0 = u1*w0;  u1v0 = u1*v0;
  v0w0 = v0*w0;  u0w0 = u0*w0;  u0v0 = u0*v0;


  Dv[0] = ( - C.VXYZ[ 0] * v1w1
            + C.VXYZ[ 1] * v1w1
            - C.VXYZ[ 2] * v0w1
            + C.VXYZ[ 3] * v0w1
            - C.VXYZ[ 4] * v1w0
            + C.VXYZ[ 5] * v1w0
            - C.VXYZ[ 6] * v0w0
            + C.VXYZ[ 7] * v0w0
          )* C.iD;

  Dv[1] = ( - C.VXYZ[ 8] * v1w1
            + C.VXYZ[ 9] * v1w1
            - C.VXYZ[10] * v0w1
            + C.VXYZ[11] * v0w1
            - C.VXYZ[12] * v1w0
            + C.VXYZ[13] * v1w0
            - C.VXYZ[14] * v0w0
            + C.VXYZ[15] * v0w0
          )* C.iD;

  Dv[2] = ( - C.VXYZ[16] * v1w1
            + C.VXYZ[17] * v1w1
            - C.VXYZ[18] * v0w1
            + C.VXYZ[19] * v0w1
            - C.VXYZ[20] * v1w0
            + C.VXYZ[21] * v1w0
            - C.VXYZ[22] * v0w0
            + C.VXYZ[23] * v0w0
          )* C.iD;

  
  Dv[3] = ( - C.VXYZ[ 0] * u1w1
            - C.VXYZ[ 1] * u0w1
            + C.VXYZ[ 2] * u1w1
            + C.VXYZ[ 3] * u0w1
            - C.VXYZ[ 4] * u1w0
            - C.VXYZ[ 5] * u0w0
            + C.VXYZ[ 6] * u1w0
            + C.VXYZ[ 7] * u0w0
          )* C.iD;

  Dv[4] = ( - C.VXYZ[ 8] * u1w1
            - C.VXYZ[ 9] * u0w1
            + C.VXYZ[10] * u1w1
            + C.VXYZ[11] * u0w1
            - C.VXYZ[12] * u1w0
            - C.VXYZ[13] * u0w0
            + C.VXYZ[14] * u1w0
            + C.VXYZ[15] * u0w0
          )* C.iD;

  Dv[5] = ( - C.VXYZ[16] * u1w1
            - C.VXYZ[17] * u0w1
            + C.VXYZ[18] * u1w1
            + C.VXYZ[19] * u0w1
            - C.VXYZ[20] * u1w0
            - C.VXYZ[21] * u0w0
            + C.VXYZ[22] * u1w0
            + C.VXYZ[23] * u0w0
          )* C.iD;
  
  Dv[6] = ( - C.VXYZ[ 0] * u1v1
            - C.VXYZ[ 1] * u0v1
            - C.VXYZ[ 2] * u1v0
            - C.VXYZ[ 3] * u0v0
            + C.VXYZ[ 4] * u1v1
            + C.VXYZ[ 5] * u0v1
            + C.VXYZ[ 6] * u1v0
            + C.VXYZ[ 7] * u0v0
           )* C.iD;

  Dv[7] = ( - C.VXYZ[ 8] * u1v1
            - C.VXYZ[ 9] * u0v1
            - C.VXYZ[10] * u1v0
            - C.VXYZ[11] * u0v0
            + C.VXYZ[12] * u1v1
            + C.VXYZ[13] * u0v1
            + C.VXYZ[14] * u1v0
            + C.VXYZ[15] * u0v0
           )* C.iD;
  
  Dv[8] = ( - C.VXYZ[16] * u1v1
            - C.VXYZ[17] * u0v1
            - C.VXYZ[18] * u1v0
            - C.VXYZ[19] * u0v0
            + C.VXYZ[20] * u1v1
            + C.VXYZ[21] * u0v1
            + C.VXYZ[22] * u1v0
            + C.VXYZ[23] * u0v0
          )* C.iD;
  
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
