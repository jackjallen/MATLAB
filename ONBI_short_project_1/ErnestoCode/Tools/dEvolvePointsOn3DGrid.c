/*
   phi = dEvolvePointsOn3DGrid( V , Gx , Gy , Gz , points , 
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


#define    jacobian_method_def      PLUS
#define    boundary_mode_def        VALUE
#define    maxSTEPS_def             50


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


#define EVOLVE_JAC            getDV( XYZT , Dv );                                            \
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
  int  ii0, jj0, kk0;
  int  ii1, jj1, kk1;
} CELL; 
CELL C = { -1 , -1 , -1 , 0 , -1 , 0 , -1 , 0 , -1 };


_inline void setCELL( real x , real y , real z );
_inline void getV( real *xyzt , real *v );
_inline void getDV( real *xyzt , real *Dv );
        void MatrixExp3x3( real * );
        void Inverse3x3( real * , real * );
        void Matrix3x3_v( real *, real * , real *);
  
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    jacobian_methods     { PLUS   , TIMES }                   jacobian_method;

  real  u0, u1, v0, v1, w0, w1;
  real    XYZT[4], xk[4], v[3], vp[3], k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3], err[3], Dv[9];
  real    xxnr, yynr, zznr;
  real    maxV, minD, hmax, hmin, h;

  real    prct , last_prct=0;
  unsigned char         *MASK, VERBOSE;
  real    *FORCES;
  real    JAC[9], dJAC[9], iJAC[9], H[3], JH[3];
  
  real    *XYZ, *Xn, *Yn, *Zn;
  int     Odims[50], ndims, d, nP, p, in, jn, kn, In, Jn, Kn, c;

  int     maxSTEPS;
  real    MAT[16];


  real    *O, *DERIV, *weigths=NULL;
  unsigned char    *visited = NULL;
  char    STR[100];
  int     argN;
  
  
  if( nlhs != 2 ){
    myErrMsgTxt("usar con 2 argumentos de salida!!!!");
  }
  
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
  
  jacobian_method     = jacobian_method_def;
  boundary_mode       = boundary_mode_def;
  maxSTEPS            = maxSTEPS_def;
  MASK                = NULL;
  FORCES              = NULL;
  MAT(3,3)            = 0;
  VERBOSE             = 0;

  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );

    if( ! myStrcmpi(STR,"verbose")                                 ) { VERBOSE = 1;                   argN++; continue; }



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


    if( ! myStrcmpi( STR,"forces") ){
      argN++;
      if( nrhs > argN && mxIsNumeric( prhs[argN] ) &&  mxGetM( prhs[argN] ) == nP  &&  mxGetN( prhs[argN] ) == NSD  ){
        FORCES = myGetPr( prhs[argN] );
        argN++; continue;
      }
      myErrMsgTxt("After the word FORCES a nP x NSD numerical is expected.");
    }


    mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
  }
  /*END Parsing arguments*/
  

  if( FORCES == NULL ){
      myErrMsgTxt("There is no forces!!!!");
  }    
  
  
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
    Odims[ndims-1] = 1;
    Odims[ndims  ] = NSD;
    Odims[ndims+1] = NSD;
    
  } else {
    
    ndims = 4;
    Odims[0] = In;
    Odims[1] = Jn;
    Odims[2] = Kn;
    Odims[3] = 1;
    Odims[4] = NSD;
    Odims[5] = NSD;
    
  }
  /*END Checking sizes*/
  
  /*Creating output*/
  plhs[0] = mxCreateNumericArray( ndims+1 , Odims , mxREAL_CLASS , mxREAL );
  O = (real *) mxGetData( plhs[0] );

  if( nlhs > 1 ){
    plhs[1] = mxDuplicateArray( prhs[0] );
    DERIV = (real *) mxGetData( plhs[1] );
    for( p = 0 ; p < IJK*NSD ; p++ ){
      DERIV[p] = 0.0;
    }
  }
  /*END Creating output*/
  visited = mxMalloc( IJK*sizeof( unsigned char) );
  for( p = 0 ; p < IJK ; p++ ){
    visited[p] = 0;
  }
  weigths = mxMalloc( IJK*NSD*NSD*sizeof( real ) );
  for( p = 0 ; p < IJK*NSD*NSD ; p++ ){
    weigths[p] = 0.0;
  }

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

  
  
  hmin     = 1.0/maxSTEPS;

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
      mexPrintf("USER INTERRUP in  dEvolvePointsOn3DGrid  !!!\n");
      myErrMsgTxt("USER INTERRUP in  dEvolvePointsOn3DGrid  !!!");
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
//       mexPrintf("p: %d\n" , p ); myFlush();
    }

    
                    O[ p        ] = NAN;
                    O[ p +   nP ] = NAN;
    if( NSD == 3 ){ O[ p + 2*nP ] = NAN; }

    
    
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
    
    
    if( myISNAN( XYZT[0] ) || myISNAN( XYZT[1] ) || myISNAN( XYZT[2] ) ){
      continue;
    }
    

    if( NSD == 2 && FORCES[ p ] == 0  && FORCES[ p + nP ] == 0 ){
      continue;
    }
    if( NSD == 3 && FORCES[ p ] == 0  && FORCES[ p + nP ] == 0  &&  FORCES[ p + 2*nP ] == 0){
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
    
    
    XYZT[3] = 0;
    
        JAC[0] = JAC[4] = JAC[8] = 1.0;
        JAC[1] = JAC[2] = JAC[3] = JAC[5] = JAC[6] = JAC[7] = 0.0;
    

          h    = hmin;
          while( XYZT[3] < 1 ){
//             mexPrintf("XYZT[3]:   %.15g\n", XYZT[3] ); myFlush();
            if( h > ( 1 - XYZT[3] ) ){ h = 1 - XYZT[3]; }

            getV( XYZT , v );
//             mexPrintf("OK getV\n"); myFlush();
            
            #define  vis(i,j,k)   visited[ (i)  +  (j)*I  +  (k)*IJ  ]
            vis(  C.ii0  ,  C.jj0  ,  C.kk0  ) = 1;
            vis(  C.ii1  ,  C.jj0  ,  C.kk0  ) = 1;
            vis(  C.ii0  ,  C.jj1  ,  C.kk0  ) = 1;
            vis(  C.ii1  ,  C.jj1  ,  C.kk0  ) = 1;
            vis(  C.ii0  ,  C.jj0  ,  C.kk1  ) = 1;
            vis(  C.ii1  ,  C.jj0  ,  C.kk1  ) = 1;
            vis(  C.ii0  ,  C.jj1  ,  C.kk1  ) = 1;
            vis(  C.ii1  ,  C.jj1  ,  C.kk1  ) = 1;


//             mexPrintf("OK vis\n"); myFlush();

            
            #define weigths(i,j,k,a,b)   weigths[ (i)  +  (j)*I  +  (k)*IJ + (a)*IJK + (b)*IJK*NSD ]
            if( NSD == 2 ){
              Inverse3x3( JAC , iJAC );

              u0 = XYZT[0] - C.X0;   u1 = C.X1 - XYZT[0];
              v0 = XYZT[1] - C.Y0;   v1 = C.Y1 - XYZT[1];
            

              weigths( C.ii0 , C.jj0 , 0 , 0 , 0 ) += ( h * u1 * v1 * iJAC[0] );
              weigths( C.ii1 , C.jj0 , 0 , 0 , 0 ) += ( h * u0 * v1 * iJAC[0] );
              weigths( C.ii0 , C.jj1 , 0 , 0 , 0 ) += ( h * u1 * v0 * iJAC[0] );
              weigths( C.ii1 , C.jj1 , 0 , 0 , 0 ) += ( h * u0 * v0 * iJAC[0] );
              
              weigths( C.ii0 , C.jj0 , 0 , 1 , 0 ) += ( h * u1 * v1 * iJAC[1] );
              weigths( C.ii1 , C.jj0 , 0 , 1 , 0 ) += ( h * u0 * v1 * iJAC[1] );
              weigths( C.ii0 , C.jj1 , 0 , 1 , 0 ) += ( h * u1 * v0 * iJAC[1] );
              weigths( C.ii1 , C.jj1 , 0 , 1 , 0 ) += ( h * u0 * v0 * iJAC[1] );
              
              weigths( C.ii0 , C.jj0 , 0 , 0 , 1 ) += ( h * u1 * v1 * iJAC[3] );
              weigths( C.ii1 , C.jj0 , 0 , 0 , 1 ) += ( h * u0 * v1 * iJAC[3] );
              weigths( C.ii0 , C.jj1 , 0 , 0 , 1 ) += ( h * u1 * v0 * iJAC[3] );
              weigths( C.ii1 , C.jj1 , 0 , 0 , 1 ) += ( h * u0 * v0 * iJAC[3] );

              weigths( C.ii0 , C.jj0 , 0 , 1 , 1 ) += ( h * u1 * v1 * iJAC[4] );
              weigths( C.ii1 , C.jj0 , 0 , 1 , 1 ) += ( h * u0 * v1 * iJAC[4] );
              weigths( C.ii0 , C.jj1 , 0 , 1 , 1 ) += ( h * u1 * v0 * iJAC[4] );
              weigths( C.ii1 , C.jj1 , 0 , 1 , 1 ) += ( h * u0 * v0 * iJAC[4] );

            } else {
              Inverse3x3( JAC , iJAC );

              u0 = XYZT[0] - C.X0;   u1 = C.X1 - XYZT[0];
              v0 = XYZT[1] - C.Y0;   v1 = C.Y1 - XYZT[1];
              w0 = XYZT[2] - C.Z0;   w1 = C.Z1 - XYZT[2];

              for( c = 0 ; c < NSD ; c++ ){
              for( d = 0 ; d < NSD ; d++ ){
              
                weigths( C.ii0 , C.jj0 , C.kk0 , c , d ) += ( h * u1 * v1 * w1 * iJAC[ c + d*NSD ] );
                weigths( C.ii1 , C.jj0 , C.kk0 , c , d ) += ( h * u0 * v1 * w1 * iJAC[ c + d*NSD ] );
                weigths( C.ii0 , C.jj1 , C.kk0 , c , d ) += ( h * u1 * v0 * w1 * iJAC[ c + d*NSD ] );
                weigths( C.ii1 , C.jj1 , C.kk0 , c , d ) += ( h * u0 * v0 * w1 * iJAC[ c + d*NSD ] );
                weigths( C.ii0 , C.jj0 , C.kk1 , c , d ) += ( h * u1 * v1 * w0 * iJAC[ c + d*NSD ] );
                weigths( C.ii1 , C.jj0 , C.kk1 , c , d ) += ( h * u0 * v1 * w0 * iJAC[ c + d*NSD ] );
                weigths( C.ii0 , C.jj1 , C.kk1 , c , d ) += ( h * u1 * v0 * w0 * iJAC[ c + d*NSD ] );
                weigths( C.ii1 , C.jj1 , C.kk1 , c , d ) += ( h * u0 * v0 * w0 * iJAC[ c + d*NSD ] );

              }}              
            }    
    
//             mexPrintf("empieza EVOLVE_JAC\n"); myFlush();            
            
            EVOLVE_JAC

//             mexPrintf("OK EVOLVE_JAC\n"); myFlush();            
    
            XYZT[0] += h*v[0];
            XYZT[1] += h*v[1];
            XYZT[2] += h*v[2];
            XYZT[3] += h;
            
          }



      if( NSD == 3 ){

        O[ p        ] = XYZT[0];
        O[ p +   nP ] = XYZT[1];
        O[ p + 2*nP ] = XYZT[2];


        for( d = 0 ; d < IJK ; d++ ){
          if( !visited[d] ){ continue; }
          
          for( c = 0 ; c < NSD ; c++ ){
            H[0] = weigths[ d  +        + c*IJK*NSD ];
            H[1] = weigths[ d  +   IJK  + c*IJK*NSD ];
            H[2] = weigths[ d  + 2*IJK  + c*IJK*NSD ];

            weigths[ d  +        + c*IJK*NSD ] = 0.0;
            weigths[ d  +   IJK  + c*IJK*NSD ] = 0.0;
            weigths[ d  + 2*IJK  + c*IJK*NSD ] = 0.0;
            
            Matrix3x3_v( JAC , H , JH );
            
            DERIV[ d + c*IJK ] +=  ( JH[0] * FORCES[ p        ] + 
                                     JH[1] * FORCES[ p +   nP ] +
                                     JH[2] * FORCES[ p + 2*nP ] 
                                    );
          }
          
          visited[d] = 0;
        }
        
      } else {

        O[ p        ] = XYZT[0];
        O[ p +   nP ] = XYZT[1];


        for( d = 0 ; d < IJK ; d++ ){
          if( !visited[d] ){ continue; }
          
          for( c = 0 ; c < NSD ; c++ ){
            H[0] = weigths[ d  +        + c*IJK*NSD ];
            H[1] = weigths[ d  +   IJK  + c*IJK*NSD ];
            H[2] = 0.0;

            weigths[ d  +        + c*IJK*NSD ] = 0.0;
            weigths[ d  +   IJK  + c*IJK*NSD ] = 0.0;
            
            Matrix3x3_v( JAC , H , JH );
            
            DERIV[ d + c*IJK ] +=  ( JH[0] * FORCES[ p        ] + 
                                     JH[1] * FORCES[ p +   nP ] 
                                   );
          }
          
          visited[d] = 0;
        }
        
      }

      
  }}}


  EXIT:
    if( visited != NULL ){ mxFree( visited ); }
    if( weigths != NULL ){ mxFree( weigths ); }
    myFreeALLOCATES();
}


_inline void setCELL( real x , real y , real z ){
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
        C.ii0  = C.ii;        C.X0 = X[C.ii0];
        C.ii1  = C.ii0 + 1;     C.X1 = X[C.ii1];
      } else {
        C.ii0 = 0; C.X0 = OX;
        C.ii1 = 0; C.X1 = LX;
      }

      if( J > 1 ){
        C.jj = GetInterval( y , Y , J , C.jj );
        C.jj0  = C.jj;        C.Y0 = Y[C.jj0];
        C.jj1  = C.jj0 + 1;     C.Y1 = Y[C.jj1];
      } else {
        C.jj0 = 0; C.Y0 = OY;
        C.jj1 = 0; C.Y1 = LY;
      }

      if( K > 1 ){
        C.kk = GetInterval( z , Z , K , C.kk );
        C.kk0  = C.kk;        C.Z0 = Z[C.kk0];
        C.kk1  = C.kk0 + 1;     C.Z1 = Z[C.kk1];
      } else {
        C.kk0 = 0; C.Z0 = OZ;
        C.kk1 = 0; C.Z1 = LZ;
      }
      
      break;

    case CLOSEST:
      if( I == 1 ){
        C.ii0 = 0;   C.X0 = x - 0.5;
        C.ii1 = 0;   C.X1 = x + 0.5;
      } else if( x < OX ){
        C.ii = 0;
        C.ii0  = 0;  C.X0 = x - 0.5;
        C.ii1  = 0;  C.X1 = x + 0.5;
      } else if( x > LX ){
        C.ii = MAX( 0 , I-2 );
        C.ii0  = I-1;  C.X0 = x - 0.5;;
        C.ii1  = I-1;  C.X1 = x + 0.5;;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        C.ii0  = C.ii;        C.X0 = X[C.ii0];
        C.ii1  = C.ii0 + 1;     C.X1 = X[C.ii1];
      }
      
      if( J == 1 ){
        C.jj0 = 0;   C.Y0 = y - 0.5;
        C.jj1 = 0;   C.Y1 = y + 0.5;
      } else if( y < OY ){
        C.jj = 0;
        C.jj0  = 0;  C.Y0 = y - 0.5;
        C.jj1  = 0;  C.Y1 = y + 0.5;
      } else if( y > LY ){
        C.jj = MAX( 0 , J-2 );
        C.jj0  = J-1;  C.Y0 = y - 0.5;
        C.jj1  = J-1;  C.Y1 = y + 0.5;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        C.jj0  = C.jj;        C.Y0 = Y[C.jj0];
        C.jj1  = C.jj0 + 1;     C.Y1 = Y[C.jj1];
      }

      if( K == 1 ){
        C.kk0 = 0;   C.Z0 = z - 0.5;
        C.kk1 = 0;   C.Z1 = z + 0.5;
      } else if( z < OZ ){
        C.kk = 0;
        C.kk0  = 0;  C.Z0 = z - 0.5;
        C.kk1  = 0;  C.Z1 = z + 0.5;
      } else if( z > LZ ){
        C.kk = MAX( 0 , K-2 );
        C.kk0  = K-1;  C.Z0 = z - 0.5;
        C.kk1  = K-1;  C.Z1 = z + 0.5;
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        C.kk0  = C.kk;        C.Z0 = Z[C.kk0];
        C.kk1  = C.kk0 + 1;     C.Z1 = Z[C.kk1];
      }      

      break;


    case DECAY_TO_ZERO:
      
      if( x < OX-boundary_size[0] || x > LX+boundary_size[1] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( y < OY-boundary_size[2] || y > LY+boundary_size[3] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      if( z < OZ-boundary_size[3] || z > LZ+boundary_size[5] ){ C.isOutside = 1; C.X0 = 0; C.X1 = -1;  return; }
      
      if( x < OX ){
        isSingularCase = 1;
        C.ii = 0;
        C.ii0  = -1;        C.X0 = OX-boundary_size[0];
        C.ii1  = 0;         C.X1 = OX;
      } else if( x > LX ) {
        isSingularCase = 1;
        C.ii = MAX( 0 , I-2 );
        C.ii0  = I-1;       C.X0 = LX;
        C.ii1  = -1;        C.X1 = LX+boundary_size[1];
      } else if( I == 1 ){
        C.ii = 0;
        C.ii0  = 0;         C.X0 = -0.5;
        C.ii1  = 0;         C.X1 =  0.5;
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        C.ii0  = C.ii;      C.X0 = X[C.ii0];
        C.ii1  = C.ii0 + 1;   C.X1 = X[C.ii1];
      }
      
      if( y < OY ){
        isSingularCase = 1;
        C.jj = 0;
        C.jj0  = -1;        C.Y0 = OY-boundary_size[2];
        C.jj1  = 0;         C.Y1 = OY;
      } else if( y > LY ) {
        isSingularCase = 1;
        C.jj = MAX( 0 , J-2 );
        C.jj0  = J-1;       C.Y0 = LY;
        C.jj1  = -1;        C.Y1 = LY+boundary_size[3];
      } else if( J == 1 ){
        C.jj = 0;
        C.jj0  = 0;         C.Y0 = -0.5;
        C.jj1  = 0;         C.Y1 =  0.5;
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        C.jj0  = C.jj;      C.Y0 = Y[C.jj0];
        C.jj1  = C.jj0 + 1;   C.Y1 = Y[C.jj1];
      }

      if( z < OZ ){
        isSingularCase = 1;
        C.kk = 0;
        C.kk0  = -1;        C.Z0 = OZ-boundary_size[4];
        C.kk1  = 0;         C.Z1 = OZ;
      } else if( z > LZ ) {
        isSingularCase = 1;
        C.kk = MAX( 0 , K-2 );
        C.kk0  = K-1;       C.Z0 = LZ;
        C.kk1  = -1;        C.Z1 = LZ+boundary_size[5];
      } else if( K == 1 ){
        C.kk = 0;
        C.kk0  = 0;         C.Z0 = -0.5;
        C.kk1  = 0;         C.Z1 =  0.5;
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        C.kk0  = C.kk;      C.Z0 = Z[C.kk0];
        C.kk1  = C.kk0 + 1;   C.Z1 = Z[C.kk1];
      }

      break;

    case CIRCULAR:

      if( x < OX || x > LX ){
        n = floor( ( x - OX )/( LX-OX ) );
        x = x - n*( LX-OX );
      } else { n = 0; }
      if( I == 1 ) {
        C.ii0  = 0;           C.X0 = x - 0.5 + n*( LX-OX );
        C.ii1  = 0;           C.X1 = x + 0.5 + n*( LX-OX );
      } else if( x < X[0] ){
        C.ii0  = I-1;         C.X0 = X[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LX-OX );
        C.ii1  = 0;           C.X1 = X[ 0 ] + n*( LX-OX );
      } else if( x > X[I-1] ){
        C.ii0  = I-1;         C.X0 = X[C.ii0] + n*( LX-OX );
        C.ii1  = 0;           C.X1 = X[I-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LX-OX );
      } else {
        C.ii = GetInterval( x , X , I , C.ii );
        C.ii0  = C.ii;        C.X0 = X[C.ii0] + n*( LX-OX );
        C.ii1  = C.ii0 + 1;     C.X1 = X[C.ii1] + n*( LX-OX );
      }
      

      if( y < OY || y > LY ){
        n = floor( ( y - OY )/( LY-OY ) );
        y = y - n*( LY-OY );
      } else { n = 0; }
      if( J == 1 ) {
        C.jj0  = 0;           C.Y0 = y - 0.5 + n*( LY-OY );
        C.jj1  = 0;           C.Y1 = y + 0.5 + n*( LY-OY );
      } else if( y < Y[0] ){
        C.jj0  = J-1;         C.Y0 = Y[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LY-OY );
        C.jj1  = 0;           C.Y1 = Y[ 0 ] + n*( LY-OY );
      } else if( y > Y[J-1] ){
        C.jj0  = J-1;         C.Y0 = Y[C.jj0] + n*( LY-OY );
        C.jj1  = 0;           C.Y1 = Y[J-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LY-OY );
      } else {
        C.jj = GetInterval( y , Y , J , C.jj );
        C.jj0  = C.jj;        C.Y0 = Y[C.jj0] + n*( LY-OY );
        C.jj1  = C.jj0 + 1;     C.Y1 = Y[C.jj1] + n*( LY-OY );
      }
      

      if( z < OZ || z > LZ ){
        n = floor( ( z - OZ )/( LZ-OZ ) );
        z = z - n*( LZ-OZ );
      } else { n = 0; }
      if( K == 1 ) {
        C.kk0  = 0;           C.Z0 = z - 0.5 + n*( LZ-OZ );
        C.kk1  = 0;           C.Z1 = z + 0.5 + n*( LZ-OZ );
      } else if( z < Z[0] ){
        C.kk0  = K-1;         C.Z0 = Z[ 0 ] - ( boundary_size[0] + boundary_size[1] ) + n*( LZ-OZ );
        C.kk1  = 0;           C.Z1 = Z[ 0 ] + n*( LZ-OZ );
      } else if( z > Z[K-1] ){
        C.kk0  = K-1;         C.Z0 = Z[C.kk0] + n*( LZ-OZ );
        C.kk1  = 0;           C.Z1 = Z[K-1] + ( boundary_size[0] + boundary_size[1] ) + n*( LZ-OZ );
      } else {
        C.kk = GetInterval( z , Z , K , C.kk );
        C.kk0  = C.kk;        C.Z0 = Z[C.kk0] + n*( LZ-OZ );
        C.kk1  = C.kk0 + 1;     C.Z1 = Z[C.kk1] + n*( LZ-OZ );
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
    C.VX0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
    C.VX1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
    C.VX0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
    C.VX1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
    C.VX0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
    C.VX1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
    C.VX0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
    C.VX1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
    #undef VV

    #define  VV(i,j,k)   VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VY0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
    C.VY1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
    C.VY0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
    C.VY1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
    C.VY0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
    C.VY1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
    C.VY0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
    C.VY1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VZ0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
      C.VZ1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
      C.VZ0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
      C.VZ1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
      C.VZ0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
      C.VZ1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
      C.VZ0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
      C.VZ1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
      #undef VV
    }

  } else {
    
    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VX[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VX0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
    C.VX1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
    C.VX0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
    C.VX1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
    C.VX0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
    C.VX1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
    C.VX0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
    C.VX1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
    #undef VV

    #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VY[ (i)  +  (j)*I  +  (k)*IJ  ]
    C.VY0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
    C.VY1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
    C.VY0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
    C.VY1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
    C.VY0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
    C.VY1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
    C.VY0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
    C.VY1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
    #undef VV

    if( NSD == 3 ){
      #define  VV(i,j,k)   ( (i)<0 || (j)<0 || (k)<0 ) ? 0 : VZ[ (i)  +  (j)*I  +  (k)*IJ  ]
      C.VZ0000 = VV(  C.ii0  ,  C.jj0  ,  C.kk0  );
      C.VZ1000 = VV(  C.ii1  ,  C.jj0  ,  C.kk0  );
      C.VZ0100 = VV(  C.ii0  ,  C.jj1  ,  C.kk0  );
      C.VZ1100 = VV(  C.ii1  ,  C.jj1  ,  C.kk0  );
      C.VZ0010 = VV(  C.ii0  ,  C.jj0  ,  C.kk1  );
      C.VZ1010 = VV(  C.ii1  ,  C.jj0  ,  C.kk1  );
      C.VZ0110 = VV(  C.ii0  ,  C.jj1  ,  C.kk1  );
      C.VZ1110 = VV(  C.ii1  ,  C.jj1  ,  C.kk1  );
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
//     mexPrintf("[ii,C.X0,x,C.X1]:%d,%g,%g,%g  ", C.ii , C.X0 , x , C.X1 ); myFlush();
//     mexPrintf("[jj,C.Y0,y,C.Y1]:%d,%g,%g,%g  ", C.jj , C.Y0 , y , C.Y1 ); myFlush();
//     mexPrintf("[kk,C.Z0,z,C.Z1]:%d,%g,%g,%g  ", C.kk , C.Z0 , z , C.Z1 ); myFlush();
//     mexPrintf("[C.iD,C.D]:%g,%g\n",C.iD,C.D);                             myFlush();
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

//   mexPrintf("en MatrixExp3x3:\n"); myFlush();
//   mexPrintf("A= \n");
//   mexPrintf("   %.15g  %.15g  %.15g\n",A[0],A[3],A[6]);
//   mexPrintf("   %.15g  %.15g  %.15g\n",A[1],A[4],A[7]);
//   mexPrintf("   %.15g  %.15g  %.15g\n",A[2],A[5],A[8]);
//   mexPrintf("\n");
  
  
  normA = MAX3(  ABS( A[0] ) + ABS( A[1] ) + ABS( A[2] ) ,
                 ABS( A[3] ) + ABS( A[4] ) + ABS( A[5] ) ,
                 ABS( A[6] ) + ABS( A[7] ) + ABS( A[8] ) );

//   mexPrintf("normA: %g\n",normA); myFlush();

  
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



void Inverse3x3( real *XX , real *Ok ){
  real det;

  det =  XX[0]*( XX[4]*XX[8] - XX[5]*XX[7] )
       - XX[1]*( XX[3]*XX[8] - XX[5]*XX[6] )
       + XX[2]*( XX[3]*XX[7] - XX[4]*XX[6] );
  det = 1.0/det;
          
          Ok[0] = ( XX[4]*XX[8] - XX[5]*XX[7] )*det;
          Ok[1] = ( XX[2]*XX[7] - XX[1]*XX[8] )*det;
          Ok[2] = ( XX[1]*XX[5] - XX[2]*XX[4] )*det;
          Ok[3] = ( XX[5]*XX[6] - XX[3]*XX[8] )*det;
          Ok[4] = ( XX[0]*XX[8] - XX[2]*XX[6] )*det;
          Ok[5] = ( XX[2]*XX[3] - XX[0]*XX[5] )*det;
          Ok[6] = ( XX[3]*XX[7] - XX[4]*XX[6] )*det;
          Ok[7] = ( XX[1]*XX[6] - XX[0]*XX[7] )*det;
          Ok[8] = ( XX[0]*XX[4] - XX[1]*XX[3] )*det;
         
}


void Matrix3x3_v( real *M , real *V , real *MV){
  MV[0] = M[0] * V[0] + M[3] * V[1] + M[6] * V[2];
  MV[1] = M[1] * V[0] + M[4] * V[1] + M[7] * V[2];
  MV[2] = M[2] * V[0] + M[5] * V[1] + M[8] * V[2];
}
