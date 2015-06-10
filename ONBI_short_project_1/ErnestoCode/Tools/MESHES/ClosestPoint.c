#include "mex.h"
#include <stdlib.h>
#include <math.h>
#include "SCAMCOORD.h"

long MarkIdsOnBuckets( const mxArray *PLOC,
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
        n_bucket =  mxGetN( mxGetCell( PLOC , mxCalcSingleSubscript( PLOC,3,bucket_id )));
        bucket   = mxGetPr( mxGetCell( PLOC , mxCalcSingleSubscript( PLOC,3,bucket_id )));
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

  long      n_xyz, n;
  double    *xyz;
  long      n_points, p;
  double    *points;
  long      closest_id;
  mxArray   *PLOC;
  int       *ijk, i,j,k;
  int       N[3], max_N;
  int       *IDs, MARK;
  double    D_SCAM, D_SCAM2;
  int       offset_bucket;
  double    min_dist2, this_dist2;
  int       isPLOC;

/*   //read the inputs */
  n_points =  mxGetM( prhs[1] );
  points   = mxGetPr( prhs[1] );
  n_xyz    =  mxGetM( mxGetField( prhs[0], 0, "xyz") );
  xyz      = mxGetPr( mxGetField( prhs[0], 0, "xyz") );

/*  //create the outputs */
  plhs[0]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  if( nlhs > 1) {
    plhs[1]= mxCreateDoubleMatrix( n_points,3,mxREAL );
  }
  if( nlhs > 2) {
    plhs[2]= mxCreateDoubleMatrix( n_points,1,mxREAL );
  }

/*   //MARKS vectors */
  IDs = mxMalloc( n_xyz*sizeof( int ) );

/*   //if SCAM structure doesn't exists */
  if( mxGetField( prhs[0], 0, "SCAM") == NULL ) {
    D_SCAM= 1;
    isPLOC= 0;
  } else {
    D_SCAM= *mxGetPr( mxGetField( mxGetField( prhs[0], 0, "SCAM"), 0 , "D"));
  }
  D_SCAM2= D_SCAM*D_SCAM;

  PLOC= mxGetField( prhs[0], 0, "PLOC");
/*   //if PLOC structure doesn't exists */
  if(  PLOC == NULL ){
    isPLOC= 0;
    max_N = 1;
  } else {
    /* //SCAM coordinates of the points */
    ijk = mxMalloc( 3*n_points*sizeof( int ) );
    SCAMCOORD( ijk, prhs[0], points, n_points , N );

    isPLOC=1;
    max_N= MMAX( N[0],N[1],N[2] );
  }

  for( p=0 ; p<n_points ; p++ ) {
    /* //reset the MARKS */
    for ( n=0 ; n<n_xyz ; n++ ) {
      IDs[n]=0;
    }
    offset_bucket= 1;
    MARK=1;
    min_dist2= 1e+300;
    while( offset_bucket <= max_N ){
      if( isPLOC ) {
        i= ijk[ p ];
        j= ijk[ p + n_points ];
        k= ijk[ p + 2*n_points ];
        n= MarkIdsOnBuckets( PLOC,
                             MAX(1,i-offset_bucket),MIN(i+offset_bucket,N[0]),
                             MAX(1,j-offset_bucket),MIN(j+offset_bucket,N[1]),
                             MAX(1,k-offset_bucket),MIN(k+offset_bucket,N[2]),
                             IDs, MARK );
        if( n < 0 ) { 
          offset_bucket++;
          continue;
        }
      } else {
        MARK= 0;
        n= n_xyz-1;
      }
      for ( ; n >= 0 ; n-- ) {
        /* //only the the just marked ones */
        if (IDs[n]==MARK) {
          /* //distance between points(p) and m.xyz(n) */
          this_dist2 = (points[p           ]-xyz[n        ])*(points[p           ]-xyz[n        ])+
                       (points[p+  n_points]-xyz[n+  n_xyz])*(points[p+  n_points]-xyz[n+  n_xyz])+
                       (points[p+2*n_points]-xyz[n+2*n_xyz])*(points[p+2*n_points]-xyz[n+2*n_xyz]);

          if( this_dist2 < min_dist2 ){
            /* //a closer point found */
            min_dist2= this_dist2;                /* //update the min_dist2 */
            closest_id= n;
            if( min_dist2 < 1e-30 ) {                     /*  //tolerance of the search */
              break;                                      /*  //a point was found */
            }
          }
        }
      }
      if( offset_bucket*offset_bucket >= min_dist2/D_SCAM2 ) {         
        break;
      }
      if( !isPLOC ) {
        break;
      }
      offset_bucket= (int)( sqrt(min_dist2/D_SCAM2) )+1;
      MARK++;
    }

    *(mxGetPr(plhs[0])+p)= closest_id + 1;           /* //update the output id */
    if( nlhs > 1) {                                  /* //update the output point */
      *(mxGetPr(plhs[1])+p           )= xyz[closest_id        ];
      *(mxGetPr(plhs[1])+p+  n_points)= xyz[closest_id+  n_xyz];
      *(mxGetPr(plhs[1])+p+2*n_points)= xyz[closest_id+2*n_xyz];
    }
    if( nlhs > 2) {                       /* //update the output distance */
      *(mxGetPr(plhs[2])+p)= sqrt( min_dist2 ) ;
    }
  }

  if( isPLOC ) {
    mxFree(IDs);
    mxFree(ijk);
  }
}

