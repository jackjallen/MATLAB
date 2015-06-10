/*
  yi = Interp1D( Y , Gx , X , 
                INTERP_MODE                                 ( DEF: '?linear' )
                BOUNDARY_MODE , [ boundary_size ]           ( DEF: 'value'   )
                ['outside_value' 'outval'], outside_value   ( DEF: NaN       )
 
          Y  : a volume of data
          Gx : the grid where the data are defined (the original one)
          X  : the coordinates where the datas will be interpolated (the new one)
            size( Y  ) = [ numel(Gx)  c1  c2 ... ck  ]
            size( Yi ) = [ numel(X )  c1  c2 ... ck  ]

          INTERP_MODE :  'lin[ear]'  , '*lin[ear]'  , '?lin[ear]'
                         'nea[rest]' , '*nea[rest]' , '?nea[rest]'
                         'cub[ic]'
                         '*sinc'      (not yet implemented)
                 * -> supose grid equal spaced             (not distinguish) (only for compatibility purpouse)
                 ? -> control if the grid is equal spaced  (not distinguish) (only for compatibility purpouse)
  
          BOUNDARY_MODE : 'value' 'extrapolation_value'
                          'sym[metric]' 
                          'circ[ular]' , 'per[iodic]'
                          'closest'
                          'decay_to_zero' , 'zero' , 'tozero' , 'decay'
              after 'decay' option, the boundary size could be specified
              by default... LX
 
          outside_value : (double) 
                    

 
 
X= 1:100; Y= randn([numel(X),1]); xi= -5:.01:150; 
plot( xi , Interp1D(Y,X,xi,'cubic','decay',4,'outside_value',-.1) ,'.' , X,Y,'o-r');
 
*/

#include "myMEX.h"

#define IM(i,t)        IM[ (i) + (t)*I ]
#define O(p,t)          O[ (p) + (t)*P ]

#define X(i)            X[ (i) ]
#define Xn(i)          Xn[ (i) ]

#define   PutInside(x,O,L)     (x) - floor( ( (x) - (O) )/(L) ) * (L)

real CINT( real ,real ,real ,real ,real ,real ,real ,real ,real);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { ALLOCATES();
  enum    interp_modes { LINEAR , NEAREST , CUBIC , SINC };
  enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , DECAY_TO_ID , CLOSEST };

  real    *IM, *O, *X, *Xn, *DX;
  real    xx, u, uu;
  real    LX, OX;

  int     ii, ii0 , ii1;
  int     p , t , P, I ,T;
  int     Odims[50], ndims;

  enum    interp_modes interp_mode;
  enum    boundary_modes boundary_mode;

  int     isOut;
  
  real    boundary_size;
  real    outside_value;

  char    STR[100];
  int     argN;

  IM= myGetPr( prhs[0] );
  I = mySize( prhs[0] , 0 );

  ndims = myNDims( prhs[0] );
  T = 1;
  for( t=1 ; t<ndims ; t++ ){
    Odims[t]= mySize( prhs[0] , t );
    T *= Odims[t];
  }

  if( myNumel(prhs[1]) != I ){ myErrMsgTxt("numel(Gx) Coordinates do not coincide with size(IM,1)."); }
  
  X  = myGetPr( prhs[1] );
  if( !checkIsSorted(X,I) ){ myErrMsgTxt("Gx  Coordinates are not sorted."); }

  DX = DualGrid(X,I); OX= DX[0]; LX= DX[I]-OX;
  
  Xn = myGetPr( prhs[2] );
  P  = myNumel( prhs[2] );


  /*Parsing arguments*/
  /*Defaults*/
  interp_mode   = LINEAR;
  boundary_mode = VALUE;
  boundary_size = LX;
  outside_value = 0.0; outside_value = outside_value/outside_value; /*NAN*/
  
  argN = 3;
  while( nrhs > argN ) {
    if( ! mxIsChar(prhs[argN]) ){ argN++; continue; myErrMsgTxt("No keywords."); } else {
      mxGetString( prhs[argN], STR, 100 );
      if( ! myStrcmpi(STR,"linear")     || ! myStrcmpi(STR,"lin")   ){ interp_mode = LINEAR;  argN++; continue; }
      if( ! myStrcmpi(STR,"*linear")    || ! myStrcmpi(STR,"*lin")  ){ interp_mode = LINEAR;  argN++; continue; }
      if( ! myStrcmpi(STR,"?linear")    || ! myStrcmpi(STR,"?lin")  ){ interp_mode = LINEAR;  argN++; continue; }
      if( ! myStrcmpi(STR,"nearest")    || ! myStrcmpi(STR,"nea")   ){ interp_mode = NEAREST; argN++; continue; }
      if( ! myStrcmpi(STR,"*nearest")   || ! myStrcmpi(STR,"*nea")  ){ interp_mode = NEAREST; argN++; continue; }
      if( ! myStrcmpi(STR,"?nearest")   || ! myStrcmpi(STR,"?nea")  ){ interp_mode = NEAREST; argN++; continue; }
      if( ! myStrcmpi(STR,"cubic")      || ! myStrcmpi(STR,"cub")   ){ interp_mode = CUBIC;   argN++; continue; }
      if( ! myStrcmpi(STR,"sinc")                                 ){ interp_mode = SINC;    argN++; continue; }
      if( ! myStrcmpi(STR,"symmetric")  || ! myStrcmpi(STR,"sym")   ){ boundary_mode = SYMMETRIC; argN++; continue; }
      if( ! myStrcmpi(STR,"circular")   || ! myStrcmpi(STR,"periodic") ||
          ! myStrcmpi(STR,"circ")       || ! myStrcmpi(STR,"per")   ){ boundary_mode = CIRCULAR; argN++; continue; }
      if( ! myStrcmpi(STR,"closest")                                ){ boundary_mode = CLOSEST; argN++; continue; }
      if( ! myStrcmpi(STR,"value")      || ! myStrcmpi(STR,"extrapolation_value") )
                                                            { boundary_mode = VALUE;    argN++; continue; }

      if( ! myStrcmpi( STR,"decay_to_zero") 
              || ! myStrcmpi( STR,"tozero") || ! myStrcmpi( STR,"decay" ) ){
        boundary_mode = DECAY_TO_ZERO; argN++;
        if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
          if( !myIsEmpty( prhs[argN] ) ){ boundary_size = myGetValue( prhs[argN] ); }
          argN++; continue;
        }
        continue;
      }

      if( ! myStrcmpi( STR,"outside_value") || ! myStrcmpi( STR,"outval") || ! myStrcmpi( STR,"outvalue") ){
        argN++;
        if( nrhs > argN && ! mxIsChar(prhs[argN]) ){
          if( !myIsEmpty( prhs[argN] ) ){ outside_value = myGetValue(prhs[argN]); }
          argN++; continue;
        }
        myErrMsgTxt("After the word OUTSIDE_VALUE a value has to be specified.\nIf empty value, default(NaN) is set.");
      }
      
      mexPrintf("%s - ",STR); myErrMsgTxt("Invalid keyword");
    }
  }
  /*END Parsing arguments*/
  
  
  /*Creating output*/
  Odims[0]= P;
  plhs[0] = mxCreateNumericArray( ndims , Odims , mxREAL_CLASS , mxREAL );
  O = (real *)mxGetData( plhs[0] );
  /*END Creating output*/
  
  ii = -1;
  for( p=0 ; p<P ; p++ ){
    
    xx= Xn(p);

    isOut = 0;
    switch( boundary_mode ){
      case VALUE:
        if( xx < OX || xx > OX+LX ){ isOut = 1; }
        break;
        
      case CLOSEST:
        if( xx < X[ 0 ] ){ xx = X[ 0 ]; } else  if( xx > X[I-1] ){ xx = X[I-1]; }
        break;

      case SYMMETRIC:
        xx = PutInside(xx,OX,2*LX);  if( xx > OX+LX ){   xx = 2*(OX+LX) - xx; }

        break;
        
      case CIRCULAR:
        xx = PutInside(xx,OX,LX);
        break;

      case DECAY_TO_ZERO:
        if( xx < X(0)-boundary_size || xx > X(I-1)+boundary_size ){ 
          isOut = 1; break;  
        }
        
    }
    if( isOut ){
      for(t=0;t<T;t++){
        O(p,t) = outside_value;
      } 
      continue;
    }
    
    switch( interp_mode ){
      case NEAREST:
        ii = GetInterval( xx , DX , I+1 , ii );
        
        if( boundary_mode == DECAY_TO_ZERO && ii < 0 ){
          for( t=0 ; t<T ; t++ ){
            O(p,t)       = 0;
          }
        } else {
          for( t=0 ; t<T ; t++ ){
            O(p,t)       = IM(ii,t);
          }
        }
        break;

      case LINEAR:
        ii = GetInterval( xx , X , I , ii );
        
        if( I == 1 ){
          u   = 1;  uu  = 0;  ii0 = 0;  ii1 = 0;
        } else if( ii >= 0 ){
          u   = ( xx-X(ii) ) / ( X(ii+1)-X(ii) );
          uu  = 1-u;  ii0 =  ii;  ii1 = ii+1;  
        } else if( ii == -2 ) {
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              u   = ( xx - X(0) + boundary_size ) / boundary_size;
              uu  = 0; ii0 = 0;ii1 = 0;
              break;
            case SYMMETRIC: u = 0; uu = 1; ii0 = 0; ii1 = 0; break;
            case CIRCULAR:
              uu = ( X(0)-xx )/( X(1)-X(0)+X(I-1)-X(I-2) )*2;
              u  = 1-uu;  ii0 = I-1;  ii1 = 0;
              break;
          }
        } else if( ii == -101 ){
          switch(boundary_mode){
            case VALUE: for(t=0;t<T;t++){O(p,t)=outside_value;} continue; 
            case DECAY_TO_ZERO:
              uu  = ( X(I-1) + boundary_size - xx ) / boundary_size;
              u   = 0; ii0 = I-1; ii1 = I-1;
              break;
            case SYMMETRIC: u = 0; uu = 1; ii0 = I-1; ii1 = I-1; break;
            case CIRCULAR:
              u = ( xx - X(I-1) )/( X(1)-X(0)+X(I-1)-X(I-2) )*2;
              uu = 1-u;  ii0 = I-1;  ii1 = 0;
              break;
          }
        }
          
        for( t=0 ; t<T ; t++ ){
          O(p,t)= IM( ii0 , t) * uu  + IM( ii1 , t) *  u ;
        }
        
        break;

      
      
      case CUBIC:
        ii = GetInterval( xx , X , I , ii );

        if( ii < 1 || ii > I-3 ){ continue; }
        
        for( t=0 ; t<T ; t++ ){
          O(p,t)= CINT(xx,X(ii-1),IM(ii-1,t)
                         ,X( ii ),IM( ii ,t)
                         ,X(ii+1),IM(ii+1,t)
                         ,X(ii+2),IM(ii+2,t)
                       );
        }

    }
  }


  EXIT: myFreeALLOCATES();
}

real CINT( real X , real x0 , real I0 , real x1 , real I1 , real x2 , real I2 , real x3 , real I3 ){
  real D1, D2, d12, dx2, dx1;
  
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

