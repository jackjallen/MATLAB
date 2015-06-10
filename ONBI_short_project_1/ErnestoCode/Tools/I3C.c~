
typedef struct point3D {
  double x, y, z;
} point3D; 



typedef struct Mat3x3 {
  double *m;
} Mat3x3; 

typedef struct Mat4x4 {
  double *m;
} Mat4x4; 


typedef struct cell1 {
    int  ii, jj, kk;
    real X0, X1, Y0, Y1, Z0, Z1;
    data_type dt;
    void   *V000;
    void   *V100;
    void   *V110;
    void   *V010;
    void   *V001;
    void   *V101;
    void   *V111;
    void   *V011;
    
    double  DX, DY, DZ, V, iV;
} cell1;



enum    data_type { UINT8 , INT8 , UINT16, INT16 , UINT32, INT32 , SINGLE , DOUBLE } data_type;
enum    boundary_modes { VALUE , SYMMETRIC , CIRCULAR , DECAY_TO_ZERO , CLOSEST } boundary_mode;
enum    interpolation_type { NEAREST , LINEAR , QUADRATIC , CUBIC , etc } interpolation_type;

typedef struct I3C {
  void *data;
  data_type dt, interp_result_type;
  
  double *X, *Y, *Z, *T;
  int     nX, nY, nZ, nT;
  mat4x4  spatialTransform;
  
  boundary_mode  bm;
  interpolation_type it;
  double   boundary_size[6];
  
  cell1  C1;
  cell2  C2;
  cell3  C3;
  
  double
} I3C; 




