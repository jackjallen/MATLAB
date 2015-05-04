#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "SCAMCOORD.h"

/* //pto 'i' of element 'e' component 'j' */
#define P(e,i,j) xyz[ (int)(tri[ (int)(e + i*n_tri)]-1 + j*n_xyz) ]

double distance_point_to_segment( double px, double py, double pz, 
                                  double s1x, double s1y, double s1z, 
                                  double s2x, double s2y, double s2z , double* closest ){
  double t;
  double t1;
  t1= ( (s2x-s1x)*(s2x-s1x)+(s2y-s1y)*(s2y-s1y)+(s2z-s1z)*(s2z-s1z) );
  if ( t1 > 0 ) {
    t= ( (s2x-s1x)*(px-s1x)+(s2y-s1y)*(py-s1y)+(s2z-s1z)*(pz-s1z) )/t1;
  } else {
    t= 0;
  }
  
  if( t < 0 ){
      closest[0]= s1x;
      closest[1]= s1y;
      closest[2]= s1z;
      return (
              (px-s1x)*(px-s1x) + (py-s1y)*(py-s1y) + (pz-s1z)*(pz-s1z)
      );
  }
  if( t > 1  ){
      closest[0]= s2x;
      closest[1]= s2y;
      closest[2]= s2z;
      return (
              (px-s2x)*(px-s2x) + (py-s2y)*(py-s2y) + (pz-s2z)*(pz-s2z)
      );
  }
  closest[0]= s1x + t*(s2x-s1x);
  closest[1]= s1y + t*(s2y-s1y);
  closest[2]= s1z + t*(s2z-s1z);
  return (  
   (px-closest[0])*(px-closest[0]) + (py-closest[1])*(py-closest[1]) + (pz-closest[2])*(pz-closest[2])
  );
}

double DistancePointToEdges(  double px, double py, double pz, 
                              double n2x, double n2y, double n2z,
                              double n3x, double n3y, double n3z , double* closest ){
  double c1[3],c2[3],c3[3];
  double d;
  double d1, d2, d3;
  
  d1= distance_point_to_segment( px,  py,  pz,
                                 0 ,  0 ,  0 ,
                                 n2x, n2y, n2z , c1);
  d2= distance_point_to_segment( px,  py,  pz, 
                                 0 ,  0 ,  0 ,
                                 n3x, n3y, n3z , c2);
  d3= distance_point_to_segment( px,  py,  pz, 
                                 n2x, n2y, n2z, 
                                 n3x, n3y, n3z , c3);
  d= d1;
  closest[0]=c1[0];
  closest[1]=c1[1];
  closest[2]=c1[2];
  if ( d2 < d ){
    d= d2;
    closest[0]=c2[0];
    closest[1]=c2[1];
    closest[2]=c2[2];
  }
  if ( d3 < d ){
    d= d3;
    closest[0]=c3[0];
    closest[1]=c3[1];
    closest[2]=c3[2];
  }
  return d;
}

double DistancePointToElement( double px, double py, double pz, 
                               double n1x, double n1y, double n1z, 
                               double n2x, double n2y, double n2z,
                               double n3x, double n3y, double n3z, double* closest ){
  double    nn;
  double    r, s, t;
  double    d;
  double    det;
  double    Nx, Ny, Nz;

  px =  px-n1x;
  n2x= n2x-n1x;
  n3x= n3x-n1x;
  py =  py-n1y;
  n2y= n2y-n1y;
  n3y= n3y-n1y;
  pz =  pz-n1z;
  n2z= n2z-n1z;
  n3z= n3z-n1z;

  Nx = n2y*n3z-n2z*n3y;
  Ny = n1z*n3x-n1x*n3z;
  Nz = n1x*n2y-n1y*n2x;

  nn  = Nx*Nx + Ny*Ny + Nz*Nz;
  if( nn > 0 ) {
    nn  = sqrt( nn );
    Nx= Nx/nn;
    Ny= Ny/nn;
    Nz= Nz/nn;
  } else {
    d= DistancePointToEdges( px,  py,  pz, n2x, n2y, n2z, n3x, n3y, n3z , closest);
    closest[0] += n1x; closest[1] += n1y; closest[2] += n1z;
    return d;
  }

  det= -(n2y*n3z*Nx) + n2x*n3z*Ny + n2z*(n3y*Nx - n3x*Ny) + n2y*n3x*Nz - n2x*n3y*Nz;

  r= ( -(n3y*Nz*px) + n3x*Nz*py + n3z*(Ny*px - Nx*py) + n3y*Nx*pz - n3x*Ny*pz )/det;
  
  if( r < 0 ){
    d= DistancePointToEdges( px,  py,  pz, n2x, n2y, n2z, n3x, n3y, n3z , closest);
    closest[0] += n1x; closest[1] += n1y; closest[2] += n1z;
    return d;
  }
  s= ( n2y*Nz*px - n2x*Nz*py + n2z*(-(Ny*px) + Nx*py) - n2y*Nx*pz + n2x*Ny*pz )/det;

  if( s < 0 ){
    d= DistancePointToEdges( px,  py,  pz, n2x, n2y, n2z, n3x, n3y, n3z , closest);
    closest[0] += n1x; closest[1] += n1y; closest[2] += n1z;
    return d;
  }
  if( ( r + s ) > 1 ){
    d= DistancePointToEdges( px,  py,  pz, n2x, n2y, n2z, n3x, n3y, n3z , closest);
    closest[0] += n1x; closest[1] += n1y; closest[2] += n1z;
    return d;
  }

  t= ( -(n2y*n3z*px) + n2x*n3z*py + n2z*(n3y*px - n3x*py) + n2y*n3x*pz - n2x*n3y*pz )/det;  
  closest[0]= n1x + r*n2x + s*n3x;
  closest[1]= n1y + r*n2y + s*n3y;
  closest[2]= n1z + r*n2z + s*n3z;
  
  return t*t;
}

long MarkIdsOnBuckets( const mxArray *ELOC,
                         int mini, int maxi,
                         int minj, int maxj,
                         int mink, int maxk,
                         int *IDs, int MARK ){
  long      p;
  long      max_id= -1;
  long      n_bucket;
  double    *bucket;
  int       i,j,k;
  int       bucket_id[3];

  for ( i=mini-1 ; i<maxi ; i++ ) {
    for ( j=minj-1 ; j<maxj ; j++ ) {
      for ( k=mink-1 ; k<maxk ; k++ ) {
        bucket_id[0]=i; bucket_id[1]=j; bucket_id[2]=k;
        n_bucket =  mxGetN( mxGetCell( ELOC , mxCalcSingleSubscript( ELOC,3,bucket_id )));
        bucket   = mxGetPr( mxGetCell( ELOC , mxCalcSingleSubscript( ELOC,3,bucket_id )));
        for ( p=0; p<n_bucket; p++){
          if ( IDs[ (int)bucket[p]-1 ] == 0 ){
            if( bucket[p]-1 > max_id ) {
              max_id= (int)bucket[p]-1;
            }
            IDs[ (int)bucket[p]-1 ] = MARK;
          }
        }
      }
    }
  }

  return max_id;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){

  long      n_xyz, t, n_tri;
  double    *xyz;
  double    *tri;
  long      n_points, p;
  double    *points;
  mxArray   *ELOC;
  int       *ijk, i,j,k;
  int       N[3], max_N;
  int       *IDs, MARK;
  double    D_SCAM, D_SCAM2;
  double    closest_point[3];
  int       offset_bucket;
  double    min_dist2, this_dist2;
  int       isELOC;
  long      closest_id;
  double    closest[3];
  
 /*   //read the inputs */
  n_points =  mxGetM( prhs[1] );
  points   = mxGetPr( prhs[1] );
  n_xyz    =  mxGetM( mxGetField( prhs[0], 0, "xyz") );
  xyz      = mxGetPr( mxGetField( prhs[0], 0, "xyz") );
  n_tri    =  mxGetM( mxGetField( prhs[0], 0, "tri") );
  tri      = mxGetPr( mxGetField( prhs[0], 0, "tri") );

/*  //create the outputs */ 
  plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  if( nlhs > 1) {
    plhs[1]= mxCreateDoubleMatrix( n_points,3,mxREAL );
  }
  if( nlhs > 2) {
    plhs[2]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  }

/*  //MARKS vectors */
  IDs = mxMalloc( n_tri*sizeof( int ) );
  
/*  //if SCAM structure doesn't exists */
  if( mxGetField( prhs[0], 0, "SCAM") == NULL ) {
    D_SCAM= 1;
    isELOC= 0;
  } else {
    D_SCAM= *mxGetPr( mxGetField( mxGetField( prhs[0], 0, "SCAM"), 0 , "D"));
  }
  D_SCAM2= D_SCAM*D_SCAM;
  
  ELOC= mxGetField( prhs[0], 0, "ELOC");
/*   //if ELOC structure doesn't exists */
  if(  ELOC == NULL ){
    isELOC= 0;
    max_N = 1;
  } else {
/*     //SCAM coordinates of the points */
    ijk = mxMalloc( 3*n_points*sizeof( int ) );
    SCAMCOORD( ijk, prhs[0], points, n_points , N );
    
    isELOC= 1;
    max_N= MMAX(N[0],N[1],N[2]);
  }

  for( p=0 ; p<n_points ; p++ ) {
/* //     mexPrintf("%d\n", p ); */

/* //     mexPrintf("i: %d - j: %d - k: %d\n",i,j,k); */
    
/*     //reset the MARKS */
    for ( t=0 ; t<n_tri ; t++ ) {
      IDs[t]=0;
    }
    offset_bucket= 1;
    MARK=1;
    min_dist2= 1e+300;
    while( offset_bucket <= max_N ){
      if( isELOC ) {
        i= ijk[ p ];
        j= ijk[ p + n_points ];
        k= ijk[ p + 2*n_points ];
        t= MarkIdsOnBuckets( ELOC,
                             MAX(1,i-offset_bucket),MIN(i+offset_bucket,N[0]),
                             MAX(1,j-offset_bucket),MIN(j+offset_bucket,N[1]),
                             MAX(1,k-offset_bucket),MIN(k+offset_bucket,N[2]),
                             IDs, MARK );
        if(t<0) {
          offset_bucket++;
          continue;
        }
      } else {
        MARK=0;
        t= n_tri-1;
         mexPrintf("");
        
      }
  /*     //run over all the points */
      for ( ; t >= 0 ; t-- ) {
/* //         mexPrintf("");
//         mexPrintf("t= %d\n",t);
        //only the the just marked ones */
        if (IDs[t]==MARK) {
/* //           mexPrintf("t= %d\n",t);
          //distance between points(p) and m.xyz(t) */
          this_dist2 = DistancePointToElement( 
                            points[p], points[p+n_points], points[p+2*n_points],
                            P(t,0,0), P(t,0,1), P(t,0,2),
                            P(t,1,0), P(t,1,1), P(t,1,2),
                            P(t,2,0), P(t,2,1), P(t,2,2),
                            closest_point );
          if( this_dist2 < min_dist2 ){
            /* //a closer point found */
            min_dist2= this_dist2;                  /* //update the min_dist */
            closest_id= t;
            closest[0]= closest_point[0];
            closest[1]= closest_point[1];
            closest[2]= closest_point[2];
            if( min_dist2 < 1e-30 ) {                       /*  //tolerance of the search */
              break;                                        /* //a point was found */
            }
          }
        }
      }
/* //     mexPrintf("%d\n", p ); */
      
      if( offset_bucket*offset_bucket >= min_dist2/D_SCAM2 ) {   
        break;
      }
      if( !isELOC ) {
        break;
      }
      offset_bucket= (int)( sqrt(min_dist2/D_SCAM2) )+1;
      MARK++;
    }
    *(mxGetPr(plhs[0])+p)= closest_id+1;         /*   //update the output id */
    if( nlhs > 1) {                              /*  //update the output point */
      *(mxGetPr(plhs[1])+p           )= closest[0];
      *(mxGetPr(plhs[1])+p+  n_points)= closest[1];
      *(mxGetPr(plhs[1])+p+2*n_points)= closest[2];
    }
    if( nlhs > 2) {                     /*   //update the output distance */
      *(mxGetPr(plhs[2])+p)= sqrt(min_dist2);
    }
  }

  if( isELOC ) {
    mxFree(IDs);
    mxFree(ijk);
  }
}

