#include "mex.h"
#include <stdlib.h>
#include <string.h>

/* //nlhs number of output args
//plhs pointer to the output args
//nrhs number of input args
//prhs pointer to the input args
//M is the number of rows
//N is the number of columns */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
  long      n_xyz;
  long      n_tri;
  long      n_nor;
  long      n_uv;
  long      total;
  char      line[1024];
  char      filename[2000];
  FILE      *obj_file;
  
  double    *xyz;
  double    *nor;
  double    *uv;
  long      *tri;
  
  const char *names[] = {""};
  const int dims[1] = {1};
  mxArray   *mxdata_xyz;
  mxArray   *mxdata_nor;
  mxArray   *mxdata_uv;
  mxArray   *mxdata_tri;
  
  long      node1,node2,node3,uv1,uv2,uv3,nor1,nor2,nor3,n_t;
  
  mxGetString( prhs[0], filename, mxGetNumberOfElements(prhs[0]) + 1 );
/* //  mexPrintf("Reading: %s \n", filename); */

/* // Reading OBJ_input_file =====================================  */
  obj_file=fopen( filename , "r" );
  if(obj_file == NULL ){
    mexPrintf("%s", filename);
    mexErrMsgTxt("ERROR File in read_OBJ: can't be opened ...\n");
  }

  n_tri= 0;
  while( fgets(line, 1024, obj_file) != NULL){
    if( line[0] == 'f' && line[1]==' ' )
      n_tri++;
  }
  rewind(obj_file);
  total= 3*n_tri;
/* //  mexPrintf("total: %d\n",total); */
  
  xyz= (double*)mxMalloc( 3*total*sizeof( double ) );
  nor= (double*)mxMalloc( 3*total*sizeof( double ) );
  uv=  (double*)mxMalloc( 2*total*sizeof( double ) );
  tri= (long*)mxMalloc( total*sizeof( long ) );
  

  n_xyz=0;
  n_uv =0;
  n_nor=0;
  while( fgets(line, 1024, obj_file) !=NULL ){
    if( line[0] == 'v' && line[1] == ' ' ) {
      sscanf(line+1,"%lf %lf %lf", 
        ( xyz+n_xyz ) , (xyz+n_xyz+total), (xyz+n_xyz+2*total));
/* //      mexPrintf("%f  %f  %f\n", *( xyz+n_xyz ) , *(xyz+n_xyz+total), *(xyz+n_xyz+2*total) ); */
      n_xyz++;
    }
    
    if( line[0] == 'v' && line[1] == 't' ){
      sscanf(line+2,"%lf %lf", 
       (uv+n_uv), (uv+n_uv+total));
      n_uv++;
    }
  
    if( line[0] == 'v' && line[1] == 'n' ) {
      sscanf(line+2,"%lf %lf %lf", 
        ( nor+n_nor ) , (nor+n_nor+total), (nor+n_nor+2*total));
      n_nor++;
    }
  }
/* //  mexPrintf("n_xyz: %d\n",n_xyz);
//  mexPrintf("n_uv: %d\n",n_uv);
//  mexPrintf("n_nor: %d\n",n_nor); */

  
  rewind(obj_file);
  plhs[0] = mxCreateStructArray(1, dims, 0, names);
  mxdata_xyz = mxCreateDoubleMatrix( total , 3 , mxREAL );
  mxdata_nor = mxCreateDoubleMatrix( total , 3 , mxREAL );
  mxdata_uv  = mxCreateDoubleMatrix( total , 2 , mxREAL );
  mxdata_tri = mxCreateDoubleMatrix( n_tri , 3 , mxREAL );

  n_t=0;
  while( fgets(line, 1024, obj_file) !=NULL ){
    if( line[0] == 'f' && line[1]==' ' ){
      if ( sscanf(line+2,"%d/%d %d/%d %d/%d",
            &node1 , &uv1, &node2 , &uv2 , &node3 , &uv3 ) == 6 ) {
        goto ASSIGN;
      }
      if ( sscanf(line+2,"%d/%d/%d %d/%d/%d %d/%d/%d",
            &node1 , &uv1, &nor1, &node2 , &uv2 , &nor2 , &node3 , &uv3 , &nor3 ) == 9 ) {
        goto ASSIGN;
      }
      if ( sscanf(line+2,"%d//%d %d//%d %d//%d",
            &node1 , &nor1, &node2 , &nor2 , &node3 , &nor3 ) == 6 ) {
/* //              mexPrintf("es asi\n"); */
        goto ASSIGN;
      }
      if ( sscanf(line+2,"%d %d %d",
            &node1 , &node2 , &node3 ) == 3 ) {
/* //              mexPrintf("es asi\n"); */
        goto ASSIGN;
      } 
      continue;
    } 
    continue;
    
    ASSIGN:
      *(mxGetPr(mxdata_tri)+ n_t           )= 3*n_t+1;
      *(mxGetPr(mxdata_tri)+ n_t +   n_tri )= 3*n_t+2;
      *(mxGetPr(mxdata_tri)+ n_t + 2*n_tri )= 3*n_t+3;
      
    if( n_xyz ) {           
      *(mxGetPr(mxdata_xyz)+ 3*n_t             )= xyz[ node1-1           ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +     total )= xyz[ node1-1 +   total ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +   2*total )= xyz[ node1-1 + 2*total ];
      
      *(mxGetPr(mxdata_xyz)+ 3*n_t +1          )= xyz[ node2-1           ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +1+   total )= xyz[ node2-1 +   total ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +1+ 2*total )= xyz[ node2-1 + 2*total ];

      *(mxGetPr(mxdata_xyz)+ 3*n_t +2          )= xyz[ node3-1           ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +2+   total )= xyz[ node3-1 +   total ];
      *(mxGetPr(mxdata_xyz)+ 3*n_t +2+ 2*total )= xyz[ node3-1 + 2*total ];
    }
    if( n_nor ) {           
      *(mxGetPr(mxdata_nor)+ 3*n_t             )= nor[ node1-1           ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +     total )= nor[ node1-1 +   total ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +   2*total )= nor[ node1-1 + 2*total ];
      
      *(mxGetPr(mxdata_nor)+ 3*n_t +1          )= nor[ node2-1           ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +1+   total )= nor[ node2-1 +   total ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +1+ 2*total )= nor[ node2-1 + 2*total ];

      *(mxGetPr(mxdata_nor)+ 3*n_t +2          )= nor[ node3-1           ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +2+   total )= nor[ node3-1 +   total ];
      *(mxGetPr(mxdata_nor)+ 3*n_t +2+ 2*total )= nor[ node3-1 + 2*total ];
    }
    if( n_uv ) {           
      *(mxGetPr(mxdata_uv) + 3*n_t             )= uv[ uv1-1           ];
      *(mxGetPr(mxdata_uv) + 3*n_t   +   total )= uv[ uv1-1 +   total ];

      *(mxGetPr(mxdata_uv) + 3*n_t +1          )= uv[ uv2-1           ];
      *(mxGetPr(mxdata_uv) + 3*n_t +1+   total )= uv[ uv2-1 +   total ];

      *(mxGetPr(mxdata_uv) + 3*n_t +2          )= uv[ uv3-1           ];
      *(mxGetPr(mxdata_uv) + 3*n_t +2 +  total )= uv[ uv3-1 +   total ];
    }
    n_t++;
  }
  if( n_xyz ){
    mxAddField( plhs[0],    "xyz" );
    mxSetField( plhs[0], 0, "xyz", mxdata_xyz );
  }
  if( n_tri ){
    mxAddField( plhs[0],    "tri" );
    mxSetField( plhs[0], 0, "tri", mxdata_tri );
  }
  if( n_uv ){
    mxAddField( plhs[0],    "uv" );
    mxSetField( plhs[0], 0, "uv", mxdata_uv );
  }
  if( n_nor ){
    mxAddField( plhs[0],    "xyzNORMALS" );
    mxSetField( plhs[0], 0, "xyzNORMALS", mxdata_nor );
  }
  fclose( obj_file );
}
