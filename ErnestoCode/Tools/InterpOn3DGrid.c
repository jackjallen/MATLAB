/*
  o = InterpOn3DGrid( IM , oGx , oGy , oGz , [ { nGx , nGy , nGz }  or  [n_points x NSD ] ,
                ['om' 'omatrix'          ] , original_grid_matrix     ( DEF: eye(4)    )
                ['matrix' 'nm' 'nmatrix' ] , new_grid_matrix          ( DEF: eye(4)    )
                INTERP_MODE                                 ( DEF: '?linear' )
                BOUNDARY_MODE , [ boundary_size ]           ( DEF: 'value'   )
                ['outside_value' 'outval'], outside_value   ( DEF: NaN       )
                ['mask'                   , logical_of_numel_number_of_points
                 'verbose'                               print percentage every 10%
                 'weights'                       return the influence weight instead of interpolated values
            );
 
          IM : a volume of data
          oGx, oGy, oGz : the grid where the data are defined (the original one)
          nGx, nGy, nGz : the grid where the datas will be interpolated (the new one)
            size( IM ) = [ numel(oGx)  numel(oGy)  numel(oGz)    times  c1  c2 ... ck  ]
            size(  o ) = [ numel(nGx)  numel(nGy)  numel(nGz)    times  c1  c2 ... ck  ]

          original_grid_matrix : a 4x4 Transformation Matrix

          new_grid_matrix      : a 4x4 Transformation Matrix

          INTERP_MODE :  'lin[ear]'  , '*lin[ear]'  , '?lin[ear]'
                         'nea[rest]' , '*nea[rest]' , '?nea[rest]'
                         'cub[ic]'
                         '*sinc'      (not yet implemented)
                 * -> supose grid equal spaced
                 ? -> control if the grid is equal spaced
  
          BOUNDARY_MODE : 'value' 'extrapolation_value'
                          'sym[metric]' 
                          'circ[ular]' , 'per[iodic]'
                          'closest'
                          'decay_to_zero' , 'zero' , 'tozero' , 'decay'
                          'decay_to_id'   , 'id'   , 'toid'     (not yet implemented)
                    the data do not decay along the singletons
              after 'decay' option, the boundary size could be specified
              by default... ( LX + LY + LZ )/3
 
          outside_value : (real)
 
          'mask', a logical of numel n_points , where mask is false, return NaN
 
*/

/*

to compile in linux 

mmex InterpOn3DGrid.c -I/extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/ -I/extra/disco1/miaTools/moreTools/CSparse/Include /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_mex.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_add.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_amd.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_chol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_cholsol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_counts.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_cumsum.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_dfs.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_dmperm.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_droptol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_dropzeros.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_dupl.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_entry.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_etree.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_fkeep.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_gaxpy.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_happly.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_house.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_ipvec.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_load.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_lsolve.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_ltsolve.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_lu.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_lusol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_malloc.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_maxtrans.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_multiply.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_norm.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_permute.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_pinv.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_post.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_print.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_pvec.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_qr.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_qrsol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_scatter.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_scc.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_schol.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_sqr.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_symperm.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_tdfs.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_transpose.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_compress.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_updown.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_usolve.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_utsolve.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_util.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_reach.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_spsolve.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_ereach.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_leaf.o /extra/disco1/miaTools/moreTools/CSparse/MATLAB/CSparse/cs_randperm.o
 
 
*/

#include "myMEX.h"
#include "cs_mex.h"


#if !defined( real )
  #define   real       real
#endif

#if !defined( mxREAL_CLASS )
  #define   mxREAL_CLASS       mxDOUBLE_CLASS
#endif


#define V(i,j,k,t)    V[ (k)*IJ + (j)*I + (i) + (t)*IJK ]
#define O(p,t)          O[ (p) + (t)*nP ]

#define nM(i,j)        nM[ 4*(j) + (i) ]
#define oM(i,j)        oM[ 4*(j) + (i) ]
#define ioM(i,j)      ioM[ 4*(j) + (i) ]
#define MAT(i,j)      MAT[ 4*(j) + (i) ]


#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/( (L) - (O) ) ) * ( (L) - (O) )


#define      readBOUNDS                                                   \
      if( nrhs > argN && !mxIsChar(prhs[argN]) ){                         \
        if( !myIsEmpty( prhs[argN] ) ){                                   \
          if( myNumel( prhs[argN] ) == 1 ){                               \
            boundary_size[0] = myGetValue(prhs[argN]);                    \
            boundary_size[1] = boundary_size[0];                          \
            boundary_size[2] = boundary_size[0];                          \
            boundary_size[3] = boundary_size[0];                          \
            boundary_size[4] = boundary_size[0];                          \
            boundary_size[5] = boundary_size[0];                          \
          } else if( myNumel( prhs[argN] ) == 6 ){                        \
            boundary_size[0] = myGetValueIth( prhs[argN] , 0 );           \
            boundary_size[1] = myGetValueIth( prhs[argN] , 1 );           \
            boundary_size[2] = myGetValueIth( prhs[argN] , 2 );           \
            boundary_size[3] = myGetValueIth( prhs[argN] , 3 );           \
            boundary_size[4] = myGetValueIth( prhs[argN] , 4 );           \
            boundary_size[5] = myGetValueIth( prhs[argN] , 5 );           \
          } else {                                                        \
            myErrMsgTxt("1 or 6 values expected as BoundarySize");        \
          }                                                               \
        }                                                                 \
        argN++; continue;                                                 \
      }                                                                   \

      
      

real CINT( real , real , real , real , real , real , real , real , real );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  ALLOCATES();
  enum    interp_modes    { LINEAR , NEAREST , CUBIC , SINC }                                       interp_mode;
  enum    boundary_modes  { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID , CLOSEST }  boundary_mode;

  real            *V, *O;
  real            *X, *Y, *Z, *Xn, *Yn, *Zn, *nM, *oM, ioM[16], MAT[16], DET, *DX, *DY, *DZ;
  real            *XYZ;
  real            xx,yy,zz,dx,dy,dz,xxnr,yynr,zznr;
  real            LX, LY, LZ, OX, OY, OZ;

  real            u,uu,v,vv,w,ww;
  
  real            prct , last_prct=0;
  unsigned char   *MASK, VERBOSE, RETURN_W;
  int             NSD;

  int             ii0 , ii1 , jj0 , jj1 , kk0 , kk1;
  int             in , jn , kn , p , t ,ii , jj , kk;
  int             I  , J  , K  , IJ , IJK ,
                  In , Jn , Kn , nP , T;
  int             Odims[16], ndims;
                  
  int             isEqualSpacedX, isEqualSpacedY, isEqualSpacedZ;
  int             isOut;
  
  real            boundary_size[6];
  real            outside_value;

  char            STR[100];
  int             argN;

  cs              *spTRIPLET , *spMAT ;
  

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
  T = myNumel( prhs[0] )/IJK;
  
  if( myNumel( prhs[1] ) != I ){ myErrMsgTxt("numel(oGx) Coordinates do not coincide with size(V,1)."); }
  if( myNumel( prhs[2] ) != J ){ myErrMsgTxt("numel(oGy) Coordinates do not coincide with size(V,2)."); }
  if( NSD == 3 ){
    if( myNumel( prhs[3] ) != K ){ myErrMsgTxt("numel(oGz) Coordinates do not coincide with size(V,3)."); }
  }
  
  X = myGetPr(prhs[1]);
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("oGx  Coordinates are not sorted."); }
  
  Y = myGetPr(prhs[2]);
  if( !checkIsSorted(Y,J) ){ myErrMsgTxt("oGy  Coordinates are not sorted."); }

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
    
    nP = myNumel( prhs[4] );

    if( ( nP%NSD ) ){ myErrMsgTxt("Number of XYZ have to be multiple of %d.", NSD); }
    nP = nP/NSD;

    In = nP;
    Jn = 1;
    Kn = 1;
    
    XYZ  = myGetPr( prhs[4] );

  }
    
  
  V = myGetPr( prhs[0] );
  
  
  /*Parsing arguments*/
  /*Defaults*/
  interp_mode   = LINEAR;
  boundary_mode = VALUE;
  
  boundary_size[0]    = ( I > 1 ) ? ( X[ 1 ] - X[ 0 ])/2 : 0.5;
  boundary_size[1]    = ( I > 1 ) ? ( X[I-1] - X[I-2])/2 : 0.5;
  boundary_size[2]    = ( J > 1 ) ? ( Y[ 1 ] - Y[ 0 ])/2 : 0.5;
  boundary_size[3]    = ( J > 1 ) ? ( Y[J-1] - Y[J-2])/2 : 0.5;
  boundary_size[4]    = ( K > 1 ) ? ( Z[ 1 ] - Z[ 0 ])/2 : 0.5;
  boundary_size[5]    = ( K > 1 ) ? ( Z[K-1] - Z[K-2])/2 : 0.5;
  
  isEqualSpacedX = -1;
  isEqualSpacedY = -1;
  isEqualSpacedZ = -1;
  outside_value = NAN;
  oM = NULL;
  nM = NULL;
  MASK = NULL;
  MAT(3,3) = 0;
  VERBOSE = 0;
  RETURN_W = 0;
  
  argN = 5;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); }
    mxGetString( prhs[argN], STR, 100 );

    if( ! myStrcmpi(STR,"verbose")                              ) { VERBOSE = 1;                              argN++; continue; }
    if( ! myStrcmpi(STR,"weights")                              ) { RETURN_W = 1;                             argN++; continue; }
    
    
    if( ! myStrcmpi(STR,"linear")    || ! myStrcmpi(STR,"lin")  ){
      interp_mode = LINEAR;
      isEqualSpacedX =  0; isEqualSpacedY =  0; isEqualSpacedZ =  0;
      argN++; continue; }
    if( ! myStrcmpi(STR,"*linear")   || ! myStrcmpi(STR,"*lin") ){
      interp_mode = LINEAR;
      isEqualSpacedX =  1; isEqualSpacedY =  1; isEqualSpacedZ =  1;
      argN++; continue; }
    if( ! myStrcmpi(STR,"?linear")   || ! myStrcmpi(STR,"?lin") ){
      interp_mode = LINEAR;
      isEqualSpacedX = -1; isEqualSpacedY = -1; isEqualSpacedZ = -1;
      argN++; continue; }


    if( ! myStrcmpi(STR,"nearest")   || ! myStrcmpi(STR,"nea")  ){
      interp_mode = NEAREST;
      isEqualSpacedX =  0; isEqualSpacedY =  0; isEqualSpacedZ =  0;
      argN++; continue; }
    if( ! myStrcmpi(STR,"*nearest")  || ! myStrcmpi(STR,"*nea") ){
      interp_mode = NEAREST;
      isEqualSpacedX =  1; isEqualSpacedY =  1; isEqualSpacedZ =  1;
      argN++; continue; }
    if( ! myStrcmpi(STR,"?nearest")  || ! myStrcmpi(STR,"?nea") ){
      interp_mode = NEAREST;
      isEqualSpacedX = -1; isEqualSpacedY = -1; isEqualSpacedZ = -1;
      argN++; continue; }


    if( ! myStrcmpi(STR,"cubic")     || ! myStrcmpi(STR,"cub")  ){
      interp_mode = CUBIC;
      isEqualSpacedX =  0; isEqualSpacedY =  0; isEqualSpacedZ =  0;
      argN++; continue; }
    
    if( ! myStrcmpi(STR,"sinc")                                 ){
      interp_mode = SINC;
      isEqualSpacedX = -1; isEqualSpacedY = -1; isEqualSpacedZ = -1;
      argN++; continue; }

    
    
    if( ! myStrcmpi(STR,"closest")                              ){
      boundary_mode = CLOSEST; argN++;    readBOUNDS; continue; 
    }

    if( ! myStrcmpi(STR,"value")     || ! myStrcmpi(STR,"extrapolation_value") ){ 
      boundary_mode = VALUE; argN++;      readBOUNDS; continue;
    }

    if( ! myStrcmpi(STR,"symmetric") || ! myStrcmpi(STR,"sym")  ){
      boundary_mode = SYMMETRIC; argN++;  readBOUNDS; continue;
    }

    
    if( ! myStrcmpi(STR,"circular")  || ! myStrcmpi(STR,"periodic") || 
        ! myStrcmpi(STR,"circ")      || ! myStrcmpi(STR,"per")  ){ 
      boundary_mode = CIRCULAR; argN++;    readBOUNDS; continue;
    }

    
    if( ! myStrcmpi( STR,"decay_to_zero") || ! myStrcmpi( STR,"zero") 
            || ! myStrcmpi( STR,"tozero") || ! myStrcmpi( STR,"decay") ){
      boundary_mode = DECAY_TO_ZERO; argN++;    readBOUNDS; continue;
    }

    
    if( ! myStrcmpi( STR,"outside_value") || ! myStrcmpi( STR,"outval") || ! myStrcmpi( STR,"outvalue") ){
      argN++;
      if( nrhs > argN && myIsParameter( prhs[argN] ) ){
        if( !myIsEmpty( prhs[argN] ) ){ outside_value = myGetValue( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word OUTSIDE_VALUE a value has to be specified.\nIf empty value, default(NaN) is set.");
    }

    
    if( ! myStrcmpi( STR,"omatrix") || ! myStrcmpi( STR,"om") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) &&  (  myIsEmpty( prhs[argN] ) || ( mySize( prhs[argN] , 0 ) == 4  &&  mySize( prhs[argN] , 1 ) == 4  && myNDims( prhs[argN] ) == 2 ) ) ){
        if( !myIsEmpty( prhs[argN] ) ){ oM = myGetPr( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word oMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
    }

    
    
    if( ! myStrcmpi( STR,"matrix") || ! myStrcmpi( STR,"nmatrix") || ! myStrcmpi( STR,"nm") ){
      argN++;
      if( nrhs > argN && !mxIsChar( prhs[argN] ) &&  (  myIsEmpty( prhs[argN] ) || ( mySize( prhs[argN] , 0 ) == 4  &&  mySize( prhs[argN] , 1 ) == 4  && myNDims( prhs[argN] ) == 2 ) ) ){
        if( !myIsEmpty( prhs[argN] ) ){ nM = myGetPr( prhs[argN] ); }
        argN++; continue;
      }
      myErrMsgTxt("After the word nMATRIX a (4x4) matrix has to be specified.\nIf empty value, default(eye(4)) is set.");
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
    
    
    mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
  }
  /*END Parsing arguments*/
  

  
  if( oM != NULL ){
    DET = oM(0,2)*oM(1,1)*oM(2,0) - oM(0,1)*oM(1,2)*oM(2,0) - oM(0,2)*oM(1,0)*oM(2,1) + oM(0,0)*oM(1,2)*oM(2,1) + oM(0,1)*oM(1,0)*oM(2,2) - oM(0,0)*oM(1,1)*oM(2,2);
    
    ioM(0,0)= (oM(1,2)*oM(2,1) - oM(1,1)*oM(2,2))/DET;
    ioM(0,1)= (oM(0,1)*oM(2,2) - oM(0,2)*oM(2,1))/DET;
    ioM(0,2)= (oM(0,2)*oM(1,1) - oM(0,1)*oM(1,2))/DET;
    ioM(0,3)= (oM(0,2)*oM(1,3)*oM(2,1) + oM(0,3)*oM(1,1)*oM(2,2) + oM(0,1)*oM(1,2)*oM(2,3) - oM(0,3)*oM(1,2)*oM(2,1) - oM(0,1)*oM(1,3)*oM(2,2) - oM(0,2)*oM(1,1)*oM(2,3))/DET;

    ioM(1,0)= (oM(1,0)*oM(2,2) - oM(1,2)*oM(2,0))/DET;
    ioM(1,1)= (oM(0,2)*oM(2,0) - oM(0,0)*oM(2,2))/DET;
    ioM(1,2)= (oM(0,0)*oM(1,2) - oM(0,2)*oM(1,0))/DET;
    ioM(1,3)= (oM(0,3)*oM(1,2)*oM(2,0) + oM(0,0)*oM(1,3)*oM(2,2) + oM(0,2)*oM(1,0)*oM(2,3) - oM(0,0)*oM(1,2)*oM(2,3) - oM(0,2)*oM(1,3)*oM(2,0) - oM(0,3)*oM(1,0)*oM(2,2))/DET;

    ioM(2,0)= (oM(1,1)*oM(2,0) - oM(1,0)*oM(2,1))/DET;
    ioM(2,1)= (oM(0,0)*oM(2,1) - oM(0,1)*oM(2,0))/DET;
    ioM(2,2)= (oM(0,1)*oM(1,0) - oM(0,0)*oM(1,1))/DET;
    ioM(2,3)= (oM(0,1)*oM(1,3)*oM(2,0) + oM(0,3)*oM(1,0)*oM(2,1) + oM(0,0)*oM(1,1)*oM(2,3) - oM(0,3)*oM(1,1)*oM(2,0) - oM(0,0)*oM(1,3)*oM(2,1) - oM(0,1)*oM(1,0)*oM(2,3))/DET;

    ioM(3,0)=0;    ioM(3,1)=0;    ioM(3,2)=0;    ioM(3,3)=1;
  }
  
  if( nM == NULL && oM != NULL ){
    for(ii=0; ii<16; ii++){ MAT[ii] = ioM[ii]; }
  } else if( nM != NULL && oM == NULL ){
    for(ii=0; ii<16; ii++){ MAT[ii] = nM[ii]; }
  } else if( nM != NULL && oM != NULL ){
    for( ii = 0; ii<4 ; ii++ ){
      for( jj = 0; jj<4 ; jj++ ){
        MAT(ii,jj) = ioM(ii,0)*nM(0,jj) +
                     ioM(ii,1)*nM(1,jj) + 
                     ioM(ii,2)*nM(2,jj) + 
                     ioM(ii,3)*nM(3,jj);
      }
    }
  } else {
    MAT(3,3) = 0;
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

  
  if( I == 1 ){
    isEqualSpacedX = 1;
    dx             = 1;
  } else if( isEqualSpacedX == -1   &&  checkEqualSpaced(X,I) ){
    isEqualSpacedX = 1;
    dx             = 1.0/( X[1] - X[0] );
  } else if( isEqualSpacedX == 1 ){
    if( !checkEqualSpaced(X,I) ){ myErrMsgTxt("warning!!! X is not equal spaced\n");  }
    dx             = 1.0/( X[1] - X[0] );
  }

  if( J == 1 ){
    isEqualSpacedY = 1;
    dy             = 1;
  } else if( isEqualSpacedY == -1   &&  checkEqualSpaced(Y,J) ){
    isEqualSpacedY = 1;
    dy             = 1.0/( Y[1] - Y[0] );
  } else if( isEqualSpacedY == 1 ){
    if( !checkEqualSpaced(Y,J) ){ myErrMsgTxt("warning!!! Y is not equal spaced\n");  }
    dy             = 1.0/( Y[1] - Y[0] );
  }
      
  if( K == 1 ){
    isEqualSpacedZ = 1;
    dz             = 1;
  } else if( isEqualSpacedZ == -1   &&  checkEqualSpaced(Z,K) ){
    isEqualSpacedZ = 1;
    dz             = 1.0/( Z[1] - Z[0] );
  } else if( isEqualSpacedZ == 1 ){
    if( !checkEqualSpaced(Z,K) ){ myErrMsgTxt("warning!!! Z is not equal spaced\n");  }
    dz             = 1.0/( Z[1] - Z[0] );
  }
  

  DX = DY = DZ = NULL;
  if( interp_mode == NEAREST ){
    if( I>1 ){  DX = DualGrid( X , I );
    } else   {  DX = mxMalloc( 2*sizeof( real ) );
    }
    DX[0] = X[ 0 ] - boundary_size[0];
    DX[I] = X[I-1] + boundary_size[1];
    
    if( J>1 ){  DY = DualGrid( Y , J );
    } else   {  DY = mxMalloc( 2*sizeof( real ) );
    }
    DY[0] = Y[ 0 ] - boundary_size[2];
    DY[J] = Y[J-1] + boundary_size[3];
    
    if( K>1 ){  DZ = DualGrid( Z , K );
    } else   {  DZ = mxMalloc( 2*sizeof( real ) );
    }
    DZ[0] = Z[ 0 ] - boundary_size[4];
    DZ[K] = Z[K-1] + boundary_size[5];
  }
  
  
  
  OX = X[ 0 ]-boundary_size[0];       LX = X[I-1]+boundary_size[1];
  OY = Y[ 0 ]-boundary_size[2];       LY = Y[J-1]+boundary_size[3];
  OZ = Z[ 0 ]-boundary_size[4];       LZ = Z[K-1]+boundary_size[5];
    
  
//   mexPrintf("OX: %10g    LX: %10g\n",OX,LX);
//   mexPrintf("OY: %10g    LY: %10g\n",OY,LY);
//   mexPrintf("OZ: %10g    LZ: %10g\n",OZ,LZ);
  
  /*Creating output*/
  if( Xn == NULL ){
    
    ndims = 2;
    
    Odims[0] = nP;
    Odims[1] = T;
    
  } else {
  
    ndims = myGetSizes( prhs[0] , Odims );

    Odims[0] = In;
    Odims[1] = Jn;
    if( NSD == 3 ){
      Odims[2] = Kn;
    }
    
  }
  if( !RETURN_W ){
  
    plhs[0] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
    O = (real *) mxGetData( plhs[0] );
    spTRIPLET = NULL;
  
  } else {

    spTRIPLET = cs_spalloc( nP , IJK , IJK*20 , 1 , 1 );

  }
  /*END Creating output*/
  
  
  
  xx =  0; yy =  0; zz =  0;
  ii = -1; jj = -1; kk = -1;
  
  p = -1;
  for( kn=0 ; kn < Kn ; kn++ ){
  for( jn=0 ; jn < Jn ; jn++ ){
  for( in=0 ; in < In ; in++ ){
//     if( utIsInterruptPending() ){
//       mexPrintf("USER INTERRUP in  InterpOn3DGrid  !!!\n");
//       myErrMsgTxt("USER INTERRUP in  InterpOn3DGrid  !!!");
//     }
    p++;
//     if( VERBOSE ){
//       prct = ( (double)(p+1) )/((double)nP)*100;
//       if( prct >= last_prct ){
//         if( last_prct != 0 ){
//           mexPrintf("%3d%% done  ( %10d   of   %10d  points )\n" , (int) prct , p+1, nP ); myFlush();
//         }
//         last_prct += 10;
//       }
//     }

    
    if( MASK != NULL  &&  MASK[p] != 1 ){
      xx = NAN;
    } else if( Xn == NULL ){
                      xx = XYZ[ p        ];
                      yy = XYZ[ p +   nP ];
      if( NSD == 3 ){ zz = XYZ[ p + 2*nP ]; }
    } else {
                      xx = Xn[ in ];
                      yy = Yn[ jn ];
      if( NSD == 3 ){ zz = Zn[ kn ]; }
    }


    
    if( myISNAN( xx ) || myISNAN( yy ) || myISNAN( zz ) ){
      if( !RETURN_W ){
        for(t=0;t<T;t++){ O(p,t) = NAN; }
      }
      continue;
    }
      
    if( MAT(3,3) == 1 ) {
      xxnr= xx;
      yynr= yy;
      zznr= zz;

      xx = xxnr*MAT(0,0) + yynr*MAT(0,1) + zznr*MAT(0,2) + MAT(0,3);
      yy = xxnr*MAT(1,0) + yynr*MAT(1,1) + zznr*MAT(1,2) + MAT(1,3);
      zz = xxnr*MAT(2,0) + yynr*MAT(2,1) + zznr*MAT(2,2) + MAT(2,3);
    }

   
    isOut = 0;
    switch( boundary_mode ){
      case VALUE:
        if( zz < OZ || zz > LZ || xx < OX || xx > LX || yy < OY || yy > LY ){
          isOut = 1;
        }
        break;

      case DECAY_TO_ZERO:
        if( zz < OZ || zz > LZ || xx < OX || xx > LX || yy < OY || yy > LY ){
          isOut = 1;
        }
        break;

      case CIRCULAR:
        xx = PutInside( xx , OX , LX );
        yy = PutInside( yy , OY , LY );
        zz = PutInside( zz , OZ , LZ );
        break;

      case SYMMETRIC:
        xx = PutInside( xx , OX , 2*LX - OX ); if( xx > LX ){ xx = 2*LX - xx; }
        yy = PutInside( yy , OY , 2*LY - OY ); if( yy > LY ){ yy = 2*LY - yy; }
        zz = PutInside( zz , OZ , 2*LZ - OZ ); if( zz > LZ ){ zz = 2*LZ - zz; }
        break;

      case CLOSEST:
        if( zz < OZ || zz > LZ || xx < OX || xx > LX || yy < OY || yy > LY ){
          isOut = 1;
        } else {
          if( xx < X[ 0 ] ){ xx = X[ 0 ]; } else if( xx > X[I-1] ){ xx = X[I-1]; }
          if( yy < Y[ 0 ] ){ yy = Y[ 0 ]; } else if( yy > Y[J-1] ){ yy = Y[J-1]; }
          if( zz < Z[ 0 ] ){ zz = Z[ 0 ]; } else if( zz > Z[K-1] ){ zz = Z[K-1]; }
        }
        break;

    }
      

    if( isOut ){
      if( !RETURN_W ){
        for(t=0;t<T;t++){ O(p,t) = outside_value; }
      }
      continue;
    }
    
    switch( interp_mode ){
      case NEAREST:

#define getNearestIndex( C , c , M , m , b_idx )                                                                             \
        if( c##c < C[0] ){                                                                                                   \
          if( boundary_mode == DECAY_TO_ZERO && c##c < ( D##C[ 0 ] + C[ 0 ] )/2 ){                                           \
            m##m = -2;                                                                                                       \
          } else if( boundary_mode == CIRCULAR && ( C[ 0 ] - c##c ) > ( c##c - ( O##C - boundary_size[ 2*b_idx + 1 ] ) ) ){  \
            m##m = M-1;                                                                                                      \
          } else {                                                                                                           \
            m##m = 0;                                                                                                        \
          }                                                                                                                  \
        } else if( c##c > C[M-1] ){                                                                                          \
          if( boundary_mode == DECAY_TO_ZERO && c##c > ( D##C[ M ] + C[M-1] )/2 ){                                           \
            m##m = -101;                                                                                                     \
          } else if( boundary_mode == CIRCULAR && ( c##c - C[M-1] ) > ( ( L##C + boundary_size[ 2*b_idx     ] ) - c##c ) ){  \
            m##m = 0;                                                                                                        \
          } else {                                                                                                           \
            m##m = M-1;                                                                                                      \
          }                                                                                                                  \
        } else if( M == 1 ){             m##m = 0;                                                                           \
        } else if( isEqualSpaced##C ){   m##m = (int)( ( ( c##c-D##C[1] )*d##c ) + 1 );                                      \
        } else {                         m##m = GetInterval( c##c , D##C , M+1 , m##m );                                     \
        }                                                                                                                    \
        

getNearestIndex( X , x , I , i , 0 )
getNearestIndex( Y , y , J , j , 1 )
getNearestIndex( Z , z , K , k , 2 )


        if( !RETURN_W ){

          if( kk < 0 || ii < 0 || jj < 0 ){ for( t=0 ; t<T ; t++ ){   O(p,t)       = 0;               }
          } else {
              for( t=0 ; t<T ; t++ ){   O(p,t)       = V(ii,jj,kk,t);   }
          }

        } else {
          
          if( kk < 0 || ii < 0 || jj < 0 ){
          } else {
            cs_entry( spTRIPLET , p , (ii) + (jj)*I + (kk)*IJ  , 1 );
          }
          
        }


        break;
        
      case LINEAR:

#define getLinearWeights( C , c , M , m , z , b_idx )                                                                                     \
        if( isEqualSpaced##C ){                                                                                                           \
          if( c##c < C[ 0 ] ){              m##m = -2;                                                                                    \
          } else if( c##c  >   C[M-1] ){    m##m = -101;                                                                                  \
          } else if( c##c  ==  C[M-1] ){    m##m = M-2;                                                                                   \
          } else {                          m##m = (int) ( ( c##c - C[0] )*d##c );                                                        \
          }                                                                                                                               \
        } else {                            m##m = GetInterval( c##c , C , M , m##m );                                                    \
        }                                                                                                                                 \
        if( M == 1 ){                                                                                                                     \
          z        = 1;                                                                                                                   \
          z##z     = 0;                                                                                                                   \
          m##m##0  = 0;                                                                                                                   \
          m##m##1  = 0;                                                                                                                   \
        } else if( m##m >= 0 ){                                                                                                           \
          z        = ( c##c - C[m##m] ) / ( C[m##m+1] - C[m##m] );                                                                        \
          z##z     = 1-z;                                                                                                                 \
          m##m##0  = m##m;                                                                                                                \
          m##m##1  = m##m+1;                                                                                                              \
        } else if( m##m == -2 ) {                                                                                                         \
          switch(boundary_mode){                                                                                                          \
            case VALUE:                                                                                                                   \
              if( !RETURN_W ){                                                                                                            \
                for( t = 0 ; t < T ; t++ ){ O(p,t) = outside_value; }                                                                     \
              }                                                                                                                           \
              continue;                                                                                                                   \
            case CLOSEST:                                                                                                                 \
              z        = 1;                                                                                                               \
              z##z     = 0;                                                                                                               \
              m##m##0  = 0;                                                                                                               \
              m##m##1  = 0;                                                                                                               \
              break;                                                                                                                      \
            case DECAY_TO_ZERO:                                                                                                           \
              z        = ( c##c - O##C ) / boundary_size[ 2*b_idx ];                                                                      \
              z##z     = 0;                                                                                                               \
              m##m##0  = 0;                                                                                                               \
              m##m##1  = 0;                                                                                                               \
              break;                                                                                                                      \
            case SYMMETRIC:                                                                                                               \
              z        = 0;                                                                                                               \
              z##z     = 1;                                                                                                               \
              m##m##0  = 0;                                                                                                               \
              m##m##1  = 0;                                                                                                               \
              break;                                                                                                                      \
            case CIRCULAR:                                                                                                                \
              z##z     = ( C[0] - c##c )/( boundary_size[ 2*b_idx ] + boundary_size[ 2*b_idx + 1 ] );                                     \
              z        = 1-z##z;                                                                                                          \
              m##m##0  = M-1;                                                                                                             \
              m##m##1  = 0;                                                                                                               \
              break;                                                                                                                      \
          }                                                                                                                               \
        } else if( m##m == -101 ){                                                                                                        \
          switch(boundary_mode){                                                                                                          \
            case VALUE:                                                                                                                   \
              if( !RETURN_W ){                                                                                                            \
                for( t = 0 ; t < T ; t++ ){ O(p,t) = outside_value; }                                                                     \
              }                                                                                                                           \
              continue;                                                                                                                   \
            case CLOSEST:                                                                                                                 \
              z        = 1;                                                                                                               \
              z##z     = 0;                                                                                                               \
              m##m##0  = M-1;                                                                                                             \
              m##m##1  = M-1;                                                                                                             \
              break;                                                                                                                      \
            case DECAY_TO_ZERO:                                                                                                           \
              z##z     = ( L##C - c##c ) / boundary_size[ 2*b_idx + 1 ];                                                                  \
              z        = 0;                                                                                                               \
              m##m##0  = M-1;                                                                                                             \
              m##m##1  = M-1;                                                                                                             \
              break;                                                                                                                      \
            case SYMMETRIC:                                                                                                               \
              z        = 0;                                                                                                               \
              z##z     = 1;                                                                                                               \
              m##m##0  = M-1;                                                                                                             \
              m##m##1  = M-1;                                                                                                             \
              break;                                                                                                                      \
            case CIRCULAR:                                                                                                                \
              z        = ( c##c - C[M-1] )/( boundary_size[ 2*b_idx ] + boundary_size[ 2*b_idx + 1 ] );                                   \
              z##z     = 1-z;                                                                                                             \
              m##m##0  = M-1;                                                                                                             \
              m##m##1  = 0;                                                                                                               \
              break;                                                                                                                      \
          }                                                                                                                               \
        }                                                                                                                                 \


getLinearWeights( X , x , I , i , u , 0 )
getLinearWeights( Y , y , J , j , v , 1 )
getLinearWeights( Z , z , K , k , w , 2 )

        if( !RETURN_W ){

          for( t=0 ; t<T ; t++ ){
            O(p,t)= V( ii0 , jj0 , kk0 ,t) * uu * vv * ww + 
                    V( ii1 , jj0 , kk0 ,t) *  u * vv * ww + 
                    V( ii0 , jj1 , kk0 ,t) * uu *  v * ww + 
                    V( ii1 , jj1 , kk0 ,t) *  u *  v * ww +
                    V( ii0 , jj0 , kk1 ,t) * uu * vv *  w + 
                    V( ii1 , jj0 , kk1 ,t) *  u * vv *  w +
                    V( ii0 , jj1 , kk1 ,t) * uu *  v *  w +
                    V( ii1 , jj1 , kk1 ,t) *  u *  v *  w ;
          }
          
        } else {
        
          cs_entry( spTRIPLET , p , (ii0) + (jj0)*I + (kk0)*IJ  , uu * vv * ww );
          cs_entry( spTRIPLET , p , (ii1) + (jj0)*I + (kk0)*IJ  ,  u * vv * ww );
          cs_entry( spTRIPLET , p , (ii0) + (jj1)*I + (kk0)*IJ  , uu *  v * ww );
          cs_entry( spTRIPLET , p , (ii1) + (jj1)*I + (kk0)*IJ  ,  u *  v * ww );
          cs_entry( spTRIPLET , p , (ii0) + (jj0)*I + (kk1)*IJ  , uu * vv *  w );
          cs_entry( spTRIPLET , p , (ii1) + (jj0)*I + (kk1)*IJ  ,  u * vv *  w );
          cs_entry( spTRIPLET , p , (ii0) + (jj1)*I + (kk1)*IJ  , uu *  v *  w );
          cs_entry( spTRIPLET , p , (ii1) + (jj1)*I + (kk1)*IJ  ,  u *  v *  w );

        }
        
        break;
          
      case CUBIC:
//         ii = GetInterval( xx , X , I , ii );
//         jj = GetInterval( yy , Y , J , jj );
//         kk = GetInterval( zz , Z , K , kk );
// 
//         if( ii < 1 || ii > I-3 ){ continue; }
//         if( jj < 1 || jj > J-3 ){ continue; }
//         if( kk < 1 || kk > K-3 ){ continue; }
//         
//         for( t=0 ; t<T ; t++ ){
//           O(p,t)= CINT( xx  
//                        ,X[ii-1], CINT(yy,Y[jj-1],CINT(zz,Z[kk-1],V(ii-1,jj-1,kk-1,t)
//                                                         ,Z[ kk ],V(ii-1,jj-1, kk ,t)
//                                                         ,Z[kk+1],V(ii-1,jj-1,kk+1,t)
//                                                         ,Z[kk+2],V(ii-1,jj-1,kk+2,t)
//                                                      )
//                                         ,Y[ jj ],CINT(zz,Z[kk-1],V(ii-1, jj ,kk-1,t)
//                                                         ,Z[ kk ],V(ii-1, jj , kk ,t)
//                                                         ,Z[kk+1],V(ii-1, jj ,kk+1,t)
//                                                         ,Z[kk+2],V(ii-1, jj ,kk+2,t)
//                                                      )
//                                         ,Y[jj+1],CINT(zz,Z[kk-1],V(ii-1,jj+1,kk-1,t)
//                                                         ,Z[ kk ],V(ii-1,jj+1, kk ,t)
//                                                         ,Z[kk+1],V(ii-1,jj+1,kk+1,t)
//                                                         ,Z[kk+2],V(ii-1,jj+1,kk+2,t)
//                                                      )
//                                         ,Y[jj+2],CINT(zz,Z[kk-1],V(ii-1,jj+2,kk-1,t)
//                                                         ,Z[ kk ],V(ii-1,jj+2, kk ,t)
//                                                         ,Z[kk+1],V(ii-1,jj+2,kk+1,t)
//                                                         ,Z[kk+2],V(ii-1,jj+2,kk+2,t)
//                                                      ) 
//                                      )
//                        ,X[ ii ], CINT(yy,Y[jj-1],CINT(zz,Z[kk-1],V( ii ,jj-1,kk-1,t)
//                                                         ,Z[ kk ],V( ii ,jj-1, kk ,t)
//                                                         ,Z[kk+1],V( ii ,jj-1,kk+1,t)
//                                                         ,Z[kk+2],V( ii ,jj-1,kk+2,t)
//                                                      )
//                                         ,Y[ jj ],CINT(zz,Z[kk-1],V( ii , jj ,kk-1,t)
//                                                         ,Z[ kk ],V( ii , jj , kk ,t)
//                                                         ,Z[kk+1],V( ii , jj ,kk+1,t)
//                                                         ,Z[kk+2],V( ii , jj ,kk+2,t)
//                                                      )
//                                         ,Y[jj+1],CINT(zz,Z[kk-1],V( ii ,jj+1,kk-1,t)
//                                                         ,Z[ kk ],V( ii ,jj+1, kk ,t)
//                                                         ,Z[kk+1],V( ii ,jj+1,kk+1,t)
//                                                         ,Z[kk+2],V( ii ,jj+1,kk+2,t)
//                                                      )
//                                         ,Y[jj+2],CINT(zz,Z[kk-1],V( ii ,jj+2,kk-1,t)
//                                                         ,Z[ kk ],V( ii ,jj+2, kk ,t)
//                                                         ,Z[kk+1],V( ii ,jj+2,kk+1,t)
//                                                         ,Z[kk+2],V( ii ,jj+2,kk+2,t)
//                                                      ) 
//                                      )
//                        ,X[ii+1], CINT(yy,Y[jj-1],CINT(zz,Z[kk-1],V(ii+1,jj-1,kk-1,t)
//                                                         ,Z[ kk ],V(ii+1,jj-1, kk ,t)
//                                                         ,Z[kk+1],V(ii+1,jj-1,kk+1,t)
//                                                         ,Z[kk+2],V(ii+1,jj-1,kk+2,t)
//                                                      )
//                                         ,Y[ jj ],CINT(zz,Z[kk-1],V(ii+1, jj ,kk-1,t)
//                                                         ,Z[ kk ],V(ii+1, jj , kk ,t)
//                                                         ,Z[kk+1],V(ii+1, jj ,kk+1,t)
//                                                         ,Z[kk+2],V(ii+1, jj ,kk+2,t)
//                                                      )
//                                         ,Y[jj+1],CINT(zz,Z[kk-1],V(ii+1,jj+1,kk-1,t)
//                                                         ,Z[ kk ],V(ii+1,jj+1, kk ,t)
//                                                         ,Z[kk+1],V(ii+1,jj+1,kk+1,t)
//                                                         ,Z[kk+2],V(ii+1,jj+1,kk+2,t)
//                                                      )
//                                         ,Y[jj+2],CINT(zz,Z[kk-1],V(ii+1,jj+2,kk-1,t)
//                                                         ,Z[ kk ],V(ii+1,jj+2, kk ,t)
//                                                         ,Z[kk+1],V(ii+1,jj+2,kk+1,t)
//                                                         ,Z[kk+2],V(ii+1,jj+2,kk+2,t)
//                                                      ) 
//                                      )
//                        ,X[ii+2], CINT(yy,Y[jj-1],CINT(zz,Z[kk-1],V(ii+2,jj-1,kk-1,t)
//                                                         ,Z[ kk ],V(ii+2,jj-1, kk ,t)
//                                                         ,Z[kk+1],V(ii+2,jj-1,kk+1,t)
//                                                         ,Z[kk+2],V(ii+2,jj-1,kk+2,t)
//                                                      )
//                                         ,Y[ jj ],CINT(zz,Z[kk-1],V(ii+2, jj ,kk-1,t)
//                                                         ,Z[ kk ],V(ii+2, jj , kk ,t)
//                                                         ,Z[kk+1],V(ii+2, jj ,kk+1,t)
//                                                         ,Z[kk+2],V(ii+2, jj ,kk+2,t)
//                                                      )
//                                         ,Y[jj+1],CINT(zz,Z[kk-1],V(ii+2,jj+1,kk-1,t)
//                                                         ,Z[ kk ],V(ii+2,jj+1, kk ,t)
//                                                         ,Z[kk+1],V(ii+2,jj+1,kk+1,t)
//                                                         ,Z[kk+2],V(ii+2,jj+1,kk+2,t)
//                                                      )
//                                         ,Y[jj+2],CINT(zz,Z[kk-1],V(ii+2,jj+2,kk-1,t)
//                                                         ,Z[ kk ],V(ii+2,jj+2, kk ,t)
//                                                         ,Z[kk+1],V(ii+2,jj+2,kk+1,t)
//                                                         ,Z[kk+2],V(ii+2,jj+2,kk+2,t)
//                                                      ) 
//                                      ) );
//         }
        break;
    }
      
  }}}

  
  if( RETURN_W ){

    spMAT = cs_compress( spTRIPLET );  
    cs_spfree( spTRIPLET );
    spTRIPLET = NULL;

    cs_dropzeros( spMAT );
    cs_dupl( spMAT );

    plhs[0] = cs_mex_put_sparse(&spMAT);

  }
  
  
  EXIT:
    if( spTRIPLET != NULL ){ cs_spfree( spTRIPLET );  }
    if( DX        != NULL ){ mxFree( DX );            }
    if( DY        != NULL ){ mxFree( DY );            }
    if( DZ        != NULL ){ mxFree( DZ );            }
    myFreeALLOCATES();
}


real CINT( real X , real x0 , real I0 , real x1 , real I1 , real x2 , real I2 , real x3 , real I3 ){
  real    D1, D2;
  real d12, dx2, dx1;
  
  d12 = x1-x2;
  dx2 =  X-x2;
  dx1 =  X-x1;
  
  D1 = (I0 - I1)/(x0 - x1) + (-I0 + I2)/(x0 - x2) + (I1 - I2)/(x1 - x2);
  D2 = (I1 - I2)/(x1 - x2) + (-I1 + I3)/(x1 - x3) + (I2 - I3)/(x2 - x3);
  
  return(
    (  
      dx1*I2*dx1*(2*X+x1-3*x2) + 
      dx1*( D2*dx1 + D1*dx2 )*dx2*d12 - 
      I1*dx2*dx2*(2*X-3*x1+x2)
    )
    /( d12*d12*d12 )
     
  );
  
  
}



/*

X = cumsum( 0.5 + 2*rand(10,1 ) ); X = X-X(1); X = X/X(end)*10;
DX = unique( [ X(1)-3 ; X(1)-2  ;  X(1)-(X(2)-X(1))/2 ; (X(2:end)+X(1:end-1))/2 ; X(end)+(X(end)-X(end-1))/2  ; X(end)+2 ; X(end)+3 ] );
I = cos( (X/10).^2*3*pi ); I(2) = 0.7;
I = repmat( I , [1 5 5] );

Xi = linspace( -8 , 18 , 5000 );

plot( Xi , ...
InterpOn3DGrid(I,X,-2:2,-2:2,{Xi,0,0},'outside_value',-0.2,'linear','symmetric',[2 3 1 1 1 1]) ... 
, '.b' );
line( X , I(:,3,3) , 'marker','o','markerfacecolor',[1 0 0],'color','r','linestyle','none'); set(gca,'XTick',DX,'XGrid','on','xlim',[-5 15],'ylim',[-1.5 1.5]);
title('Nearest-Value');
set(gca,'XTick', unique( bsxfun(@plus,DX(:),-5:5)) );

*/