#include "mex.h"
#include <stdlib.h>
#include <string.h>

#define myFlush()                   mexEvalString("drawnow expose;")

/* nlhs number of output args
//plhs pointer to the output args
//nrhs number of input args
//prhs pointer to the input args
//M is the number of rows
//N is the number of columns */

void read_field( FILE *file, double *data, long n_data, int dim){
  long i;
  int c;
  for( i=0; i<n_data; i++){
    for( c=0; c<dim; c++ ){
        fscanf( file, "%lf", data+(i + c*n_data) );
    }
  }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
  long      n_xyz = 0;
  long      n_tri = 0;
  long      i;
  int       numArray, array;
  int       aux, c;
  double    *data;
  char      line[1024];
  char      name[100];
  char      name2[100];
  char      filename[2000];
  FILE      *vtp_file;
  mxArray   *mxdata;
  const char *names[] = {""};
  const int dims[1] = {1};

  mxGetString( prhs[0], filename, mxGetNumberOfElements(prhs[0]) + 1 );
/* //  mexPrintf("Reading: %s \n", filename); */

/*
// Reading VTP_input_file ===================================== */
  vtp_file=fopen( filename , "r" );
  if(vtp_file == NULL ){
    mexPrintf("%s", filename);
    mexErrMsgTxt("ERROR File in read_VTP: can't be opened ...\n");
  }

  plhs[0] = mxCreateStructArray(1, dims, 0, names);

  while( fgets(line, 1024, vtp_file) !=NULL ){
    XYZ:
    if ( ! strncmp(line, "POINTS", 6) ){
      sscanf(line, "POINTS %d", &n_xyz);
      mxdata = mxCreateDoubleMatrix( n_xyz , 3 , mxREAL );
      mexPrintf("READING   xyz  (%d)\n",n_xyz); myFlush();
      read_field( vtp_file , mxGetPr( mxdata ) , n_xyz , 3);
      mxAddField( plhs[0],    "xyz" );
      mxSetField( plhs[0], 0, "xyz", mxdata );
    }

    TRI:
    if ( ! strncmp(line, "POLYGONS", 8) ){
      sscanf(line, "POLYGONS %d", &n_tri);
      mexPrintf("READING   tri  (%d)\n",n_tri); myFlush();
      mxdata = mxCreateDoubleMatrix( n_tri , 3 , mxREAL );
      data   = mxGetPr( mxdata );
      for(i=0; i<n_tri; i++){
        fscanf( vtp_file, "%d", &aux);
        for( c=0; c<aux ; c++) {
          fscanf( vtp_file, "%lf", data+(i+c*n_tri) );
          data[i+c*n_tri]++;
        }
        for( ; c<3 ; c++) {
          data[i+c*n_tri] = 0;
        }
      }
      mxAddField( plhs[0],    "tri" );
      mxSetField( plhs[0], 0, "tri", mxdata );
    }


    CELL:
    if ( ! strncmp(line, "CELLS", 5) ){
      sscanf(line, "CELLS %d", &n_tri);
      mexPrintf("READING   cells  (%d)\n",n_tri); myFlush();
      mxdata = mxCreateDoubleMatrix( n_tri , 10 , mxREAL );
      data   = mxGetPr( mxdata );
      for(i=0; i<n_tri; i++){
        fscanf( vtp_file, "%d", &aux);
        for( c=0; c<aux ; c++) {
          fscanf( vtp_file, "%lf", data+(i+c*n_tri) );
          data[i+c*n_tri]++;
        }
        for( ; c<10 ; c++) {
          data[i+c*n_tri] = 0;
        }
      }
      mxAddField( plhs[0],    "cell" );
      mxSetField( plhs[0], 0, "cell", mxdata );
    }
    
    
    
    POINTDATAS:
    if ( !strncmp(line, "POINT_DATA", 10) ){
      while( fgets(line, 1024, vtp_file) != NULL ){
        if ( !strncmp(line, "CELL_DATA", 9) ){ goto CELLDATAS; }
        if ( ! strncmp(line, "NORMALS", 7) ){
          mxdata = mxCreateDoubleMatrix( n_xyz , 3 , mxREAL );
          mexPrintf("READING   xyzNORMALS  (%d)\n",n_xyz); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_xyz , 3);
          mxAddField( plhs[0] ,   "xyzNORMALS" );
          mxSetField( plhs[0], 0, "xyzNORMALS", mxdata );
        }
        if ( ! strncmp(line, "TEXTURE_COORDINATES", 19) ){
          mxdata = mxCreateDoubleMatrix( n_xyz , 2 , mxREAL );
          mexPrintf("READING   xyzTEXTURE_COORDINATES  (%d)\n",n_xyz); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_xyz , 2);
          mxAddField( plhs[0],    "uv" );
          mxSetField( plhs[0], 0, "uv", mxdata );
        }
        if ( ! strncmp(line, "SCALARS", 7) ){
          sscanf( line , "SCALARS %s", name );
          sprintf( name2, "xyz%s",name );
          while( strncmp( fgets(line, 1024, vtp_file),"LOOKUP_TABLE", 12) ){ ; }
          mxdata = mxCreateDoubleMatrix( n_xyz , 1 , mxREAL );
          mexPrintf("READING   %s  (%d)\n",name2,n_xyz); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_xyz , 1);
          mxAddField( plhs[0],    name2 );
          mxSetField( plhs[0], 0, name2, mxdata );
        }
        if ( ! strncmp(line, "VECTORS", 7) ){
          sscanf( line , "VECTORS %s", name );
          sprintf( name2, "xyz%s",name );
          mxdata = mxCreateDoubleMatrix( n_xyz , 3 , mxREAL );
          mexPrintf("READING   %s  (%d)\n",name2,n_xyz); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_xyz , 3);
          mxAddField( plhs[0],    name2 );
          mxSetField( plhs[0], 0, name2, mxdata );
        }
        if ( ! strncmp(line, "FIELD", 5) ){
          sscanf( line , "FIELD %s %d", name, &numArray );
          for( array=0 ; array< numArray ; array++ ){
            fgets(line, 1024, vtp_file);
            while( sscanf(line, "%s", name)!= 1) {
              fgets(line, 1024, vtp_file);
            }
            sprintf( name2, "xyz%s",name );
            mexPrintf("READING   %s  (%d)\n",name2,n_xyz); myFlush();
            sscanf( line , "%s %d %d %s", name , &c, &aux, name);
            mxdata = mxCreateDoubleMatrix( n_xyz , c , mxREAL );
            read_field( vtp_file , mxGetPr( mxdata ), n_xyz , c);
            mxAddField( plhs[0],    name2 );
            mxSetField( plhs[0], 0, name2 , mxdata );
          }
        }
      }
    }

    CELLDATAS:
    if ( !strncmp(line, "CELL_DATA", 9) ){
      while( fgets(line, 1024, vtp_file) != NULL ){
        if ( !strncmp(line, "POINT_DATA", 10) ){ goto POINTDATAS; }
        if ( !strncmp(line, "NORMALS", 7) ){
          mxdata = mxCreateDoubleMatrix( n_tri , 3 , mxREAL );
          data   = mxGetPr( mxdata );
          mexPrintf("READING   triNORMALS  (%d)\n",n_tri); myFlush();
          read_field( vtp_file , data , n_tri , 3);
          mxAddField( plhs[0],    "triNORMALS" );
          mxSetField( plhs[0], 0, "triNORMALS", mxdata );
        }
        if ( ! strncmp(line, "SCALARS", 7) ){
          sscanf( line , "SCALARS %s", name );
          sprintf( name2, "tri%s",name );
          while( strncmp( fgets(line, 1024, vtp_file),"LOOKUP_TABLE", 12) ){ ; }
          mxdata = mxCreateDoubleMatrix( n_tri , 1 , mxREAL );
          mexPrintf("READING   %s  (%d)\n",name2,n_tri); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_tri , 1);
          mxAddField( plhs[0],    name2 );
          mxSetField( plhs[0], 0, name2, mxdata );
        }
        if ( ! strncmp(line, "VECTORS", 7) ){
          sscanf( line , "VECTORS %s", name );
          sprintf( name2, "tri%s",name );
          mxdata = mxCreateDoubleMatrix( n_tri , 3 , mxREAL );
          mexPrintf("READING   %s  (%d)\n",name2,n_tri); myFlush();
          read_field( vtp_file , mxGetPr( mxdata ), n_tri , 3);
          mxAddField( plhs[0],    name2 );
          mxSetField( plhs[0], 0, name2, mxdata );
        }
        if ( ! strncmp(line, "FIELD", 5) ){
          sscanf( line , "FIELD %s %d", name, &numArray );
          for( array=0 ; array< numArray ; array++ ){
            fgets(line, 1024, vtp_file);
            while( sscanf(line, "%s", name)!= 1) {
              fgets(line, 1024, vtp_file);
            }
            sprintf( name2, "tri%s",name );
            sscanf( line , "%s %d %d %s", name , &c, &aux, name);
            mxdata = mxCreateDoubleMatrix( n_tri , c , mxREAL );
            mexPrintf("READING   %s  (%d)\n",name2,n_tri); myFlush();
            read_field( vtp_file , mxGetPr( mxdata ), n_tri , c);
            mxAddField( plhs[0],    name2 );
            mxSetField( plhs[0], 0, name2 , mxdata );
          }
        }
      }
    }

  }
  fclose( vtp_file );

}
