#define TOLBB         1e-2
#define TOLAREA       1e-6
//#include <unistd.h>
#define CLEANLINES    0
#define TOLDIST       1e-8


#define rx 1.0
#define ry 1.0

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#include "myMEX.h"
// #define  UNICODE
// #include "mmz.h"


#include <math.h>
#include <stdio.h>

#include "vtkSignedCharArray.h"
#include "vtkIntArray.h"
#include "vtkTriangle.h"
#include "vtkPlane.h"
#include "vtkOBBTree.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkGenericCell.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellLocator.h"
#include "vtkTriangleFilter.h"

#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"


#include "gpc.c"
/* #include <stdint.h>
 #define mzCreateTicToc(i)		struct timespec __mztic##_##i, __mztoc##_##i; double mzdtimeElapsed##_##i = 0;uint64_t mztimeElapsed##_##i = 0;
 #define mzTIC(i)				clock_gettime(CLOCK_MONOTONIC, &__mztic##_##i)
 #define mzTOC(i, s)			clock_gettime(CLOCK_MONOTONIC, &__mztoc##_##i);          		\
								mztimeElapsed##_##i = ((__mztoc##_##i.tv_sec * 1000000000) + 	\
								__mztoc##_##i.tv_nsec) - ((__mztic##_##i.tv_sec * 1000000000) + \
								__mztic##_##i.tv_nsec); mzdtimeElapsed##_##i = mztimeElapsed##_##i;\
		mexPrintf("%s. Elapsed time for timer %i is %g seconds.\n", s, i, (mzdtimeElapsed##_##i * 1e-9));

*/


 #define REAL double
 extern "C" {
 #define   ANSI_DECLARATORS
 #include "triangle.h"
 }

//#define REAL double
// #define   ANSI_DECLARATORS
// #include "triangle.h"
//
//#define NO_TIMER
//#define ANSI_DECLARATORS
//#define TRILIBRARY
//#define REDUCED
//#define CDT_ONLY
//#undef  SELF_CHECK
//#include "triangle.c"
//#undef  REDUCED
//#undef  NO_TIMER
//#undef  ANSI_DECLARATORS
//#undef  TRILIBRARY
//#undef  CDT_ONLY


int triangulateID = 0;

int planeBoxOverlap(double *, double *, double * );
int triBoxOverlap( double * , double * , double * , double * , double * );

int TriangleTriangleIntersection(double *, double *, double *,
                                 double *, double *, double *,
                                 vtkPolyData * , vtkIdType , vtkIdType ); 

void CoplanarIntersection( double *, double *, double *, 
                           double *, double *, double *,
                           vtkPolyData * , vtkIdType , vtkIdType ); 

vtkPolyData *PrepareCells( vtkPolyData *A, double *BOUNDS );
int FindIntersections( vtkOBBNode *, vtkOBBNode *, vtkMatrix4x4 *, void * );

void GetOrder( int * , int I, int * );
mxArray *Split2MESH( vtkPolyData *, int , mxArray * );
void setSIDE( vtkPolyData * , vtkPolyData * );
vtkPolyData *SplitMESH( vtkPolyData * , vtkPolyData * , int * , vtkIdType * );
vtkPolyData *SplitMESH( vtkPolyData * , vtkPolyData * );
void splitFACE( vtkPolyData * , vtkIdType , vtkIdType * , int );
int LinePlaneIntersection(double *, double *, double *, double *, double &, double *);

#define SAFEDELETE( n )               if( n != NULL ){ n->Delete(); n=NULL; }
#define SAFEFREE(n)                   if( n != NULL ){ free( n ); n=NULL; }
#define FREE_triangulateio( T )       SAFEFREE( T.pointlist             );   \
                                      SAFEFREE( T.pointattributelist    );   \
                                      SAFEFREE( T.pointmarkerlist       );   \
                                      SAFEFREE( T.trianglelist          );   \
                                      SAFEFREE( T.triangleattributelist );   \
                                      SAFEFREE( T.trianglearealist      );   \
                                      SAFEFREE( T.neighborlist          );   \
                                      SAFEFREE( T.segmentlist           );   \
                                      SAFEFREE( T.segmentmarkerlist     );   \
                                      SAFEFREE( T.holelist              );   \
                                      SAFEFREE( T.regionlist            );   \
                                      SAFEFREE( T.edgelist              );   \
                                      SAFEFREE( T.edgemarkerlist        );   \
                                      SAFEFREE( T.normlist              );
#define SAFEmxFREE(n)                 if( n != NULL ){ mxFree( n ); n=NULL; }
#define writeDATA( n )  { vtkPolyDataWriter *WRITER__=vtkPolyDataWriter::New();WRITER__->SetInput( n );WRITER__->SetFileName( #n ".vtk" );WRITER__->Write();WRITER__->Delete();WRITER__=NULL; }

struct IntersectionsFinder {
  vtkPolyData   *A;
  vtkPolyData   *B;
  vtkPolyData   *AC;
  vtkPolyData   *BC;
  vtkOBBTree    *OBB_B;
  vtkPolyData   *LINES;
};
int *toSORT;
int func_to_sort(const void * a, const void * b ){
   return( toSORT[ *(int *)a ] - toSORT[ *(int *)b ] );
}

bool EXIT_OK = false;
void fnExit1 (void)
{
// #ifdef _DEBUG_MEX
  mexPrintf ("Exit function 1. Exiting!!!! \n\n");myFlush();
// #endif

  //myFreeALLOCATES();
  //usleep(1000000);
  if(!EXIT_OK)
  {
	  EXIT_OK = true;
	  throw(0);
	  //mexErrMsgTxt("");
  }

}

// SplitFACE DAVID

void splitFACE( vtkPolyData *M , vtkIdType cellid , vtkIdType *LID , int I ){
    
  struct triangulateio Tin;
  memset( &Tin  , 0 , sizeof( struct triangulateio ) );

  struct triangulateio Tout;
  memset( &Tout , 0 , sizeof( struct triangulateio ) );


  int    i,j, lid1, lid2;
  int *ADDED = (int *)malloc( MAX( I , 1024 ) * sizeof(int) );
  int nADDED = 0, a, a1, a2;
  
//   mexPrintf("cell = %i\n",cellid);
//   for (i=0;i<I;i++){mexPrintf( " %i", LID[i]);}
//   mexPrintf("\n");
  
  Tin.segmentlist       = (int *) malloc(   I * sizeof(int) );
  Tin.segmentmarkerlist = (int *) malloc(   I/2 * sizeof(int) );
  memset( Tin.segmentmarkerlist , 0 , Tin.numberofsegments * sizeof(int) );
  
  Tin.pointlist        = (REAL *) malloc( 2 * I * sizeof(REAL) );
  Tin.pointmarkerlist  = (int  *) malloc(     I * sizeof(int ) );
  memset( Tin.pointmarkerlist , 0 , Tin.numberofpoints * sizeof(int) );
  
  
  double XYZ[3];
  double R[3], P[3], Q[3], T[3][3], Den, RP[3], RQ[3], NT[3], x, y, x1, y1, x2, y2;

  // Face Vertex RPQ
  
  nADDED=3;
  ADDED[0]=LID[0];
  ADDED[1]=LID[2];
  ADDED[2]=LID[4];
  
  M->GetPoint( ADDED[0] , R );
  M->GetPoint( ADDED[1] , P );
  M->GetPoint( ADDED[2] , Q );
  
  for ( i = 0 ; i < 3 ; i++ )
  {
      RP[i]=P[i]-R[i];
      RQ[i]=Q[i]-R[i];
  }
  
  Den=(RP[0]*RP[0]*RQ[1]*RQ[1] + RP[0]*RP[0]*RQ[2]*RQ[2] - 2*RP[0]*RP[1]*RQ[0]*RQ[1] - 2*RP[0]*RP[2]*RQ[0]*RQ[2] + RP[1]*RP[1]*RQ[0]*RQ[0] + RP[1]*RP[1]*RQ[2]*RQ[2] - 2*RP[1]*RP[2]*RQ[1]*RQ[2] + RP[2]*RP[2]*RQ[0]*RQ[0] + RP[2]*RP[2]*RQ[1]*RQ[1]);
 
     
  if( Den < 1e-15 ){
    goto EXIT;
  }
  {
  // Transformation Matrix
  #define sqdist(x1,y1,x2,y2) ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
  
  Den = 1.0/Den;
  
  T[0][0]=-(RQ[0]*(RP[1]*RQ[1] + RP[2]*RQ[2]) - RP[0]*(RQ[1]*RQ[1] + RQ[2]*RQ[2]))*Den;
  T[0][1]=-(RQ[1]*(RP[0]*RQ[0] + RP[2]*RQ[2]) - RP[1]*(RQ[0]*RQ[0] + RQ[2]*RQ[2]))*Den;
  T[0][2]=-(RQ[2]*(RP[0]*RQ[0] + RP[1]*RQ[1]) - RP[2]*(RQ[0]*RQ[0] + RQ[1]*RQ[1]))*Den;

  T[1][0]=-(RP[0]*(RP[1]*RQ[1] + RP[2]*RQ[2]) - RQ[0]*(RP[1]*RP[1] + RP[2]*RP[2]))*Den;
  T[1][1]=-(RP[1]*(RP[0]*RQ[0] + RP[2]*RQ[2]) - RQ[1]*(RP[0]*RP[0] + RP[2]*RP[2]))*Den;
  T[1][2]=-(RP[2]*(RP[0]*RQ[0] + RP[1]*RQ[1]) - RQ[2]*(RP[0]*RP[0] + RP[1]*RP[1]))*Den;
  
  NT[0]=RP[1]*RQ[2]-RP[2]*RQ[1];
  NT[1]=RP[2]*RQ[0]-RP[0]*RQ[2];
  NT[2]=RP[0]*RQ[1]-RP[0]*RQ[1];

  // Initial triangle
  Tin.numberofpoints = 3;
  Tin.pointlist[ 0 ] = rx;
  Tin.pointlist[ 1 ] = ry;
  Tin.pointlist[ 2 ] = rx+1.0;
  Tin.pointlist[ 3 ] = ry;
  Tin.pointlist[ 4 ] = rx;
  Tin.pointlist[ 5 ] = ry+1.0;
  
  Tin.numberofsegments = 3; 
  Tin.segmentlist[ 0 ] = 0;
  Tin.segmentlist[ 1 ] = 1;
  Tin.segmentlist[ 2 ] = 1;
  Tin.segmentlist[ 3 ] = 2;
  Tin.segmentlist[ 4 ] = 0;
  Tin.segmentlist[ 5 ] = 2;
 
 
  // Rest of Segments
  
  for( int s = 6 ; s < I ; s=s+2 ){
    
    // Origen
    lid1 = LID[ s ];
    M->GetPoint( lid1 , XYZ );
    x1 = T[0][0]*(XYZ[0]-R[0])+T[0][1]*(XYZ[1]-R[1])+T[0][2]*(XYZ[2]-R[2])+rx;
    y1 = T[1][0]*(XYZ[0]-R[0])+T[1][1]*(XYZ[1]-R[1])+T[1][2]*(XYZ[2]-R[2])+ry;
    
    // Final
    lid2 = LID[ s+1 ];
    M->GetPoint( lid2 , XYZ );
    x2 = T[0][0]*(XYZ[0]-R[0])+T[0][1]*(XYZ[1]-R[1])+T[0][2]*(XYZ[2]-R[2])+rx;
    y2 = T[1][0]*(XYZ[0]-R[0])+T[1][1]*(XYZ[1]-R[1])+T[1][2]*(XYZ[2]-R[2])+ry;
    
    if ((x1<rx)&&(x2<rx)||(y1<ry)&&(y2<ry)||((x1+y1)>(rx+ry+1.0))&&((x2+y2)>(rx+ry+1.0)))
    {
        continue; // Segment out of the triangle
    } 
    
    // Intro Origen
    if (x1<rx) { x1=rx;}
    if (y1<ry) { y1=ry;}
    if ((x1+y1)>(rx+ry+1.0)) {x1=x1-(x1+y1-rx-ry-1.0)/2.0;y1=y1-(x1+y1-rx-ry-1.0)/2.0;}
    
     
    for ( a1 = 0 ; ( (a1 < Tin.numberofpoints) &&  ( sqdist(Tin.pointlist[2*a1],Tin.pointlist[2*a1+1],x1,y1) > TOLDIST*TOLDIST ) ) ; a1++ ){
    } 
    if ( a1 == Tin.numberofpoints ){
          ADDED[ nADDED ] = lid1;
          nADDED++;
          Tin.numberofpoints++;
          Tin.pointlist[ 2*a1 ] = x1;   
          Tin.pointlist[ 2*a1+1 ] = y1; 
    }

       
    // Intro Destino
    if (x2<rx) { x2=rx;}
    if (y2<ry) { y2=ry;}
    if ((x2+y2)>(rx+ry+1.0)) {x2=x2-(x2+y2-rx-ry-1.0)/2.0;y2=y2-(x2+y2-rx-ry-1.0)/2.0;}
    
    for( a2 = 0 ; ( (a2 < Tin.numberofpoints) && ( sqdist(Tin.pointlist[2*a2],Tin.pointlist[2*a2+1],x2,y2) > TOLDIST*TOLDIST ) ) ; a2++ ){
    } // 
    if ( a2 == Tin.numberofpoints ){
      ADDED[ nADDED ] = lid2;
      nADDED++;
      Tin.numberofpoints++;
      Tin.pointlist[ 2*a2 ] = x2;   
      Tin.pointlist[ 2*a2+1 ] = y2; 
    }
    
//     mexPrintf("\n");
    
    // Compruebo si el segmento [a1 a2] esta añadido
    if (a1!=a2)
    {
        
    if (a2<a1) {SWAP(a1,a2);}
    
    for (a = 0 ; (a < Tin.numberofsegments) && ((a1!=Tin.segmentlist[ 2*a ]) || (a2!=Tin.segmentlist[ 2*a+1 ])); a++) {
    }
    if ( a == Tin.numberofsegments )
    {
    Tin.numberofsegments++;   
    Tin.segmentlist[ 2*a ] = a1;   
    Tin.segmentlist[ 2*a+1 ] = a2;    
    } 
    
    }
  }
  

//   mexPrintf("cell = %i\n",cellid);
//   for (i=0;i<nADDED;i++){mexPrintf( " %i", ADDED[i]);}
//   mexPrintf("\n");
  
  
  if (nADDED==3) {
      goto EXIT;
  }
  
//   // Impresion de Tin
// 
//   
//   char fn[1024];
//   sprintf( fn , "triangles/t%03d.poly" , triangulateID );
//   FILE *fid = fopen(   fn ,"w");
//   fprintf( fid , "%d 2 0 0\n" , Tin.numberofpoints );
//   for( a = 0 ; a < Tin.numberofpoints ; a++ ){
//     fprintf( fid , "%d  %.20g  %.20g\n" , a+1 , Tin.pointlist[ 2*a ] , Tin.pointlist[ 2*a + 1 ] );
//   }
//   fprintf( fid , "%d 0\n" , Tin.numberofsegments );
//   for( a = 0 ; a < Tin.numberofsegments ; a++ ){
//     fprintf( fid , "%d  %d  %d\n" , a+1 , Tin.segmentlist[ 2*a ] + 1 , Tin.segmentlist[ 2*a + 1 ] + 1 );
//   }
//   fprintf( fid , "0\n");
//   fclose( fid );

//   triangulateID++;
  
   // TRIANGULATE
  
  try {
	  EXIT_OK = false;
//       goto EXIT;
	  triangulate( "pzQ" , &Tin , &Tout , (struct triangulateio *) NULL );
	  EXIT_OK = true;
  } catch( int k ) { 
    mexPrintf("triangulate dio int error en cellid: %d \n", cellid );
    SAFEFREE( ADDED );
    FREE_triangulateio( Tout );
    FREE_triangulateio( Tin  );
    throw(k);

  }
  
  catch( ... ) { 
    mexPrintf("triangulate dio unknown error en cellid: %d \n", cellid );
    SAFEFREE( ADDED );
    FREE_triangulateio( Tout );
    FREE_triangulateio( Tin  );
    throw;

  }
  
  
//   // Impresion de Tout
//   char fn2[1024];
//   sprintf( fn2 , "triangles/tri_t%03d.poly" , triangulateID );
//   FILE *fid2 = fopen(   fn2 ,"w");
//   fprintf( fid2 , "%d 2 0 0\n" , Tout.numberofpoints );
//   for( a = 0 ; a < Tout.numberofpoints ; a++ ){
//     fprintf( fid2 , "%d  %.20g  %.20g\n" , a+1 , Tout.pointlist[ 2*a ] , Tout.pointlist[ 2*a + 1 ] );
//   }
//   fprintf( fid2 , "%d 0\n" , Tout.numberofsegments );
//   for( a = 0 ; a < Tout.numberofsegments ; a++ ){
//     fprintf( fid2 , "%d  %d  %d\n" , a+1 , Tout.segmentlist[ 2*a ] + 1 , Tout.segmentlist[ 2*a + 1 ] + 1 );
//   }
//   fprintf( fid2 , "0\n");
//   fclose( fid2 );
//   
//   triangulateID++;
//   
  

  if( Tout.numberofpoints > MAX( I , 1024 ) ){
    ADDED = (int *) realloc( ADDED , Tout.numberofpoints );
  }
 
  for( a = nADDED ; a < Tout.numberofpoints ; a++ ){
    x = Tout.pointlist[ 2*a     ];
    y = Tout.pointlist[ 2*a + 1 ];
    
      
    XYZ[0] = R[0]+ RP[0]*(x-rx) + RQ[0]*(y-ry); 
    XYZ[1] = R[1]+ RP[1]*(x-rx) + RQ[1]*(y-ry); 
    XYZ[2] = R[2]+ RP[2]*(x-rx) + RQ[2]*(y-ry); 
            
    ADDED[a] = (int) M->GetPoints()->InsertNextPoint( XYZ );
  }
  

 
  vtkIdType IJK[3];
 
// CONTROL NORMALES

  double A_IJK,N[3];
    
 
  for( j = 0 ; j < 3 ; j++ ){ IJK[j] = ADDED[ Tout.trianglelist[j] ]; }


    
    M->GetPoint(IJK[0],R);
    M->GetPoint(IJK[1],P);
    M->GetPoint(IJK[2],Q);
    
   
    N[0] = ((P[1]-R[1])*(Q[2]-R[2])-(P[2]-R[2])*(Q[1]-R[1]));
    N[1] = ((P[2]-R[2])*(Q[0]-R[0])-(P[0]-R[0])*(Q[2]-R[2]));
    N[2] = ((P[0]-R[0])*(Q[1]-R[1])-(P[1]-R[1])*(Q[0]-R[0]));
    
    
    A_IJK = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
    M->ReplaceCell( cellid , 3 , IJK );

  for( i = 1 ; i < Tout.numberoftriangles ; i++ ) {

    for( j = 0 ; j < 3 ; j++ ){ IJK[j] = ADDED[ Tout.trianglelist[i*3+j] ]; }
    
    M->GetPoint(IJK[0],R);
    M->GetPoint(IJK[1],P);
    M->GetPoint(IJK[2],Q);
    
    N[0] = ((P[1]-R[1])*(Q[2]-R[2])-(P[2]-R[2])*(Q[1]-R[1]));
    N[1] = ((P[2]-R[2])*(Q[0]-R[0])-(P[0]-R[0])*(Q[2]-R[2]));
    N[2] = ((P[0]-R[0])*(Q[1]-R[1])-(P[1]-R[1])*(Q[0]-R[0]));
    
    A_IJK = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
    
   
     if ( N[0]*NT[0]+N[1]*NT[1]+N[2]*NT[2] < 0)  
        {
        //mexPrintf("Celda %i invertida con area %e\n",i,A_IJK);
        SWAP( IJK[0] , IJK[1] );}
     M->InsertNextCell( VTK_TRIANGLE , 3 , IJK );
   
  }
  }
  EXIT:
    SAFEFREE( ADDED );
    FREE_triangulateio( Tout );
    FREE_triangulateio( Tin  );
  {}
}






void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	if(nrhs < 2)
	{
		mexErrMsgTxt("Two arguments are needed.");
	}
    
  ALLOCATES();

  for( int i=0 ; i < nlhs ; i++ ){ plhs[i] = NULL; }
  vtkPolyData       *A            = NULL;
  vtkPolyData       *B            = NULL;
  vtkPolyData       *AC           = NULL;
  vtkPolyData       *BC           = NULL;
  vtkOBBTree        *OBB_A        = NULL;
  vtkOBBTree        *OBB_B        = NULL;
  vtkPolyData       *LINES        = NULL;
  vtkPolyData       *AT           = NULL;
  vtkPolyData       *BT           = NULL;
  vtkIdType         *LID          = NULL;
  int               *ORDER        = NULL;

  //registramos la funcion de captura del fallo en el triangulate
  //atexit(fnExit1);

  A = MESH2vtkPolyData( prhs[0] );
  for( int f = 0 ; f < A->GetPointData()->GetNumberOfArrays() ; f++ ){ A->GetPointData()->RemoveArray( f ); }
  for( int f = 0 ; f < A->GetCellData()->GetNumberOfArrays()  ; f++ ){ A->GetCellData()->RemoveArray( f );  }
  if( !A->GetNumberOfCells() ){ goto EXIT; } 
  {
  B = MESH2vtkPolyData( prhs[1] );
  for( int f = 0 ; f < B->GetPointData()->GetNumberOfArrays() ; f++ ){ B->GetPointData()->RemoveArray( f ); }
  for( int f = 0 ; f < B->GetCellData()->GetNumberOfArrays()  ; f++ ){ B->GetCellData()->RemoveArray( f );  }
  if( !B->GetNumberOfCells() ){ goto EXIT; } 
  
  double bA[6], bB[6];
  A->GetBounds(bA);
  B->GetBounds(bB);
  if( bA[0] > bB[1] || bA[1] < bB[0] ||
      bA[2] > bB[3] || bA[3] < bB[2] ||
      bA[4] > bB[5] || bA[5] < bB[4] ){
    goto EXIT; 
  }

  bA[0] -= TOLBB;  bA[1] += TOLBB;
  bA[2] -= TOLBB;  bA[3] += TOLBB;
  bA[4] -= TOLBB;  bA[5] += TOLBB;

  bB[0] -= TOLBB;  bB[1] += TOLBB;
  bB[2] -= TOLBB;  bB[3] += TOLBB;
  bB[4] -= TOLBB;  bB[5] += TOLBB;
  
  double BOUNDS[6]; //intersection of BB;
  BOUNDS[0] = MAX( bA[0] , bB[0] );
  BOUNDS[1] = MIN( bA[1] , bB[1] );
  BOUNDS[2] = MAX( bA[2] , bB[2] );
  BOUNDS[3] = MIN( bA[3] , bB[3] );
  BOUNDS[4] = MAX( bA[4] , bB[4] );
  BOUNDS[5] = MIN( bA[5] , bB[5] );
  
  
  vtkSignedCharArray *SIDE = NULL;

  SIDE = vtkSignedCharArray::New();
  SIDE->SetName("SIDE");
  SIDE->SetNumberOfComponents(1);
  SIDE->SetNumberOfTuples( A->GetNumberOfCells() );
  A->GetCellData()->AddArray( SIDE );
  SIDE->Delete();  

  AC = PrepareCells( A , BOUNDS ); // Devuelve las faces que se intersectan, etiquetadas con IDS original. Los SIDE de A que estan fuera los pone a 0, a 1 el resto
  if( !AC->GetNumberOfCells() ){ goto EXIT; } 


  SIDE = vtkSignedCharArray::New();
  SIDE->SetName("SIDE");
  SIDE->SetNumberOfComponents(1);
  SIDE->SetNumberOfTuples( B->GetNumberOfCells() );
  B->GetCellData()->AddArray( SIDE );
  SIDE->Delete();  

  BC = PrepareCells( B , BOUNDS );
  if( !BC->GetNumberOfCells() ){ goto EXIT; }
  
//   CreateNanoTicToc(1);
//   nanoTIC(1);
  // build locators to Find the triangle-triangle intersections between A and B
  OBB_A = vtkOBBTree::New();
  OBB_A->SetDataSet(AC);
  OBB_A->SetNumberOfCellsPerNode(2);
  OBB_A->SetMaxLevel(10000);
  OBB_A->SetTolerance(1e-6);
  OBB_A->AutomaticOn();
  OBB_A->BuildLocator();

  OBB_B = vtkOBBTree::New();
  OBB_B->SetDataSet(BC);
  OBB_B->SetNumberOfCellsPerNode(2);
  OBB_B->SetMaxLevel(10000);
  OBB_B->SetTolerance(1e-6);
  OBB_B->AutomaticOn();
  OBB_B->BuildLocator();
//   nanoTOC("Generar OBB.",1);

  
  //construct LINES to store results of intersection
  LINES = vtkPolyData::New();

  vtkPoints *LINESxyz = vtkPoints::New();
  LINESxyz->SetDataTypeToDouble();
  LINES->SetPoints(LINESxyz);
  LINESxyz->Delete();
  
  vtkCellArray *LINESlin = vtkCellArray::New();
  LINES->SetLines(LINESlin);
  LINESlin->Delete();

  vtkIntArray *IDS = NULL;

  IDS = vtkIntArray::New();
  IDS->SetName("IDS_A");
  IDS->SetNumberOfComponents(1);
  LINES->GetCellData()->AddArray( IDS );
  IDS->Delete();  

  IDS = vtkIntArray::New();
  IDS->SetName("IDS_B");
  IDS->SetNumberOfComponents(1);
  LINES->GetCellData()->AddArray( IDS );
  IDS->Delete();  
  
  //auxiliar variable to pass arguments
  struct IntersectionsFinder INT;
  INT.A  = A;
  INT.B  = B;
  INT.AC = AC;
  INT.BC = BC;
  INT.OBB_B = OBB_B;
  INT.LINES = LINES;

  //perform the triangle-triangle intersections search
  OBB_A->IntersectWithOBBTree( OBB_B , 0 , FindIntersections , &INT );
  //all intersections are stored in LINES
     
  SAFEDELETE( OBB_B       );
  SAFEDELETE( OBB_A       );
  SAFEDELETE( AC          );
  SAFEDELETE( BC          );

  int nLINES = LINES->GetNumberOfCells();
  if( !nLINES )
    { goto EXIT; } 
  
 
  if( CLEANLINES ){
    vtkCleanPolyData *CLEAN = NULL;
    CLEAN = vtkCleanPolyData::New();
    CLEAN->SetInput( LINES );
    CLEAN->ToleranceIsAbsoluteOn();
    CLEAN->SetTolerance( 1e-14 );
    CLEAN->SetAbsoluteTolerance( 1e-14 );
    CLEAN->ConvertLinesToPointsOff();
    CLEAN->PointMergingOn();
    CLEAN->Update();
    LINES->Delete();
    LINES = vtkPolyData::New();
    
    LINES->ShallowCopy( CLEAN->GetOutput() );
    CLEAN->Delete();
  }
  
  
  if( nlhs > 10 ){
    vtkTriangleFilter  *TRIANGLEFILTER = vtkTriangleFilter::New();
    TRIANGLEFILTER->SetInput( LINES );
    TRIANGLEFILTER->PassVertsOff();
    TRIANGLEFILTER->PassLinesOn();
    TRIANGLEFILTER->Update();

    plhs[10] = vtkPolyData2MESH( TRIANGLEFILTER->GetOutput() );
    TRIANGLEFILTER->Delete();
  }
  

  ORDER = ( int * ) mxMalloc( nLINES * sizeof( int ) );
  LID   = ( vtkIdType * ) mxMalloc( sizeof( vtkIdType ) * ( 2 * nLINES + 6 ) );
  
  //////////////SPLIT A
  triangulateID = 1000;
  
  LINES->GetCellData()->SetActiveScalars("IDS_A");
  
  

  
  try{
  AT = SplitMESH( A , LINES , ORDER , LID );
//   mexPrintf("SplitCOMPLETO \n" );myFlush();
  } catch( int k ) { 
//     mexPrintf("splitMESH dio error \n" );    
    SAFEmxFREE( LID         );
    SAFEmxFREE( ORDER       );
    SAFEDELETE( OBB_B       );
    SAFEDELETE( OBB_A       );
    SAFEDELETE( LINES       );
    SAFEDELETE( AC          );
    SAFEDELETE( BC          );
    SAFEDELETE( AT          );
    SAFEDELETE( A           );
    SAFEDELETE( B           );
    myFreeALLOCATES();
    mexErrMsgIdAndTxt("Intersect:triangulate", "Intersect:triangulate error type int\n");
  }
  
  catch( ... ) { 
//     mexPrintf("splitMESH dio error \n" );
    SAFEmxFREE( LID         );
    SAFEmxFREE( ORDER       );
    SAFEDELETE( OBB_B       );
    SAFEDELETE( OBB_A       );
    SAFEDELETE( LINES       );
    SAFEDELETE( AC          );
    SAFEDELETE( BC          );
    SAFEDELETE( AT          );
    SAFEDELETE( A           );
    SAFEDELETE( B           );

    myFreeALLOCATES();
    mexErrMsgIdAndTxt("Intersect:triangulate", "Intersect:triangulate error type unknown\n");

  }
  

  
  //AT = SplitMESH( A , LINES ); //tambien funciona asi, se encarga el de allocar ORDER y LID
  
  setSIDE( AT , B );
  
                  plhs[0] = Split2MESH( AT , 0 , (mxArray *)NULL );
  if( nlhs > 1 ){ plhs[1] = Split2MESH( AT , 1 , (mxArray *)NULL ); }
  if( nlhs > 2 ){ plhs[2] = Split2MESH( AT , 21, (mxArray *)NULL ); }
  if( nlhs > 3 ){ plhs[3] = Split2MESH( AT , 22, (mxArray *)NULL ); }
  if( nlhs > 8 ){ plhs[8] = vtkPolyData2MESH(AT); }
                  
  SAFEDELETE( AT );
                  
  
  if( nlhs < 5 ){ goto EXIT; }
                 
                  
  //////////////SPLIT B
  triangulateID = 2000;

  LINES->GetCellData()->SetActiveScalars("IDS_B"); //LINES->GetCellData()->RemoveArray("IDS_A");
    try{
    BT = SplitMESH( B , LINES , ORDER , LID );
  } catch( int k ) { 
//     mexPrintf("splitMESH dio error \n" );    
    SAFEmxFREE( LID         );
    SAFEmxFREE( ORDER       );
    SAFEDELETE( OBB_B       );
    SAFEDELETE( OBB_A       );
    SAFEDELETE( LINES       );
    SAFEDELETE( AC          );
    SAFEDELETE( BC          );
    SAFEDELETE( AT          );
    SAFEDELETE( BT          );
    SAFEDELETE( A           );
    SAFEDELETE( B           );
    myFreeALLOCATES();
    mexErrMsgTxt("\n");
  }
  
  catch( ... ) { 
//     mexPrintf("splitMESH dio error \n" );
    SAFEmxFREE( LID         );
    SAFEmxFREE( ORDER       );
    SAFEDELETE( OBB_B       );
    SAFEDELETE( OBB_A       );
    SAFEDELETE( LINES       );
    SAFEDELETE( AC          );
    SAFEDELETE( BC          );
    SAFEDELETE( AT          );
    SAFEDELETE( BT          );
    SAFEDELETE( A           );
    SAFEDELETE( B           );

    myFreeALLOCATES();
    mexErrMsgTxt("\n");

  }
  
  setSIDE( BT , A );
  
                  plhs[4] = Split2MESH( BT , 0 , (mxArray *)NULL );
  if( nlhs > 5 ){ plhs[5] = Split2MESH( BT , 1 , (mxArray *)NULL ); }
  if( nlhs > 6 ){ plhs[6] = Split2MESH( BT , 21 , (mxArray *)NULL ); }
  if( nlhs > 7 ){ plhs[7] = Split2MESH( BT , 22 , (mxArray *)NULL ); }
  if( nlhs > 9 ){ plhs[9] = vtkPolyData2MESH(BT); }
  
  SAFEDELETE( BT );
  

  }
  
   EXIT:

      if( plhs[0] == NULL ) { plhs[0]=mxDuplicateArray(prhs[0]); }
      if ( plhs[1] == NULL && nlhs > 1 ) { plhs[1]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[2] == NULL && nlhs > 2 ) { plhs[2]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[3] == NULL && nlhs > 3 ) { plhs[3]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[4] == NULL && nlhs > 4 ) { plhs[4]=mxDuplicateArray(prhs[1]);}
      if ( plhs[5] == NULL && nlhs > 5 ) { plhs[5]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[6] == NULL && nlhs > 6 ) { plhs[6]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[7] == NULL && nlhs > 7 ) { plhs[7]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );}
      if ( plhs[8] == NULL && nlhs > 8 ) 
          {
          plhs[8] = mxDuplicateArray(prhs[0]);
          mxAddField(plhs[8],"triSIDE");
          mxSetField(plhs[8],0,"triSIDE", mxCreateDoubleMatrix( A->GetNumberOfCells() , 1 , mxREAL ) );
          }
      if ( plhs[9] == NULL && nlhs > 9 ) 
          {
          plhs[9] = mxDuplicateArray(prhs[1]);
          mxAddField(plhs[9],"triSIDE");
          mxSetField(plhs[9],0,"triSIDE", mxCreateDoubleMatrix( B->GetNumberOfCells() , 1 , mxREAL ) );
          }    
      if( plhs[10] == NULL && nlhs > 10 )
          {
          plhs[10]=Split2MESH( (vtkPolyData *) NULL , 0 , (mxArray *) NULL );
          mxAddField(plhs[10],"triIDS_A");
          mxAddField(plhs[10],"triIDS_B");
          }
  

    SAFEmxFREE( LID         );
    SAFEmxFREE( ORDER       );
    SAFEDELETE( OBB_B       );
    SAFEDELETE( OBB_A       );
    SAFEDELETE( LINES       );
    SAFEDELETE( AC          );
    SAFEDELETE( BC          );
    SAFEDELETE( AT          );
    SAFEDELETE( BT          );
    SAFEDELETE( A           );
    SAFEDELETE( B           );

    myFreeALLOCATES();
    
    }


vtkPolyData *SplitMESH( vtkPolyData *A , vtkPolyData *LINES ){
  int nLINES = LINES->GetNumberOfCells();

  int *ORDER;
  ORDER = ( int * ) malloc( nLINES * sizeof( int ) );
  vtkIdType *LID;
  LID   = ( vtkIdType * ) malloc( sizeof( vtkIdType ) * ( 2 * nLINES + 6 ) );
  vtkPolyData *AT=NULL;
  try{ 
  AT = SplitMESH( A , LINES , ORDER , LID );
  }catch( int k ) { 
//     mexPrintf("splitMESH dio error \n" );    
  free( LID   );
  free( ORDER );
  SAFEDELETE(AT);
  throw(k);
  }
  
  catch( ... ) { 
//     mexPrintf("splitMESH dio error \n" );
    free( LID   );
    free( ORDER );
    SAFEDELETE(AT);
    throw;

  }

  free( LID   );
  free( ORDER );

  return( AT );
}

vtkPolyData *SplitMESH( vtkPolyData *A , vtkPolyData *LINES , int *ORDER , vtkIdType *LID ){

  int OriginalNumberOfPoints = A->GetNumberOfPoints();
  int OriginalNumberOfCells  = A->GetNumberOfCells();

  vtkPolyData *AT = vtkPolyData::New();

  vtkPoints *POINTS  = vtkPoints::New();
  POINTS->SetDataTypeToDouble();
  POINTS->SetNumberOfPoints( OriginalNumberOfPoints + LINES->GetNumberOfPoints() );
  memcpy(   (double *) POINTS->GetVoidPointer(0)                                , (double *)     A->GetPoints()->GetVoidPointer(0) ,     OriginalNumberOfPoints * sizeof( double ) * 3 );
  memcpy( ( (double *) POINTS->GetVoidPointer(0) + OriginalNumberOfPoints * 3 ) , (double *) LINES->GetPoints()->GetVoidPointer(0) , LINES->GetNumberOfPoints() * sizeof( double ) * 3 );

  AT->SetPoints( POINTS );
  POINTS->Delete();
  
	vtkCellArray *FACES = vtkCellArray::New();
  FACES->DeepCopy( A->GetPolys() );
  AT->SetPolys( FACES );
  FACES->Delete();
  

  int nLINES = LINES->GetNumberOfCells();
  int *LIDS = (int *) LINES->GetCellData()->GetScalars()->GetVoidPointer(0);
  GetOrder( LIDS , nLINES , ORDER );
  
  
  //Split the first input
  int           L, i;
  vtkIdType     cellid;
  vtkCell       *CELL = NULL;
  
  L = 0;
  do{ 
    i = 0;
    cellid = LIDS[ ORDER[L] ];
    //if( cellid < 0 ){ continue; }

    CELL = AT->GetCell(cellid);
    LID[ i++ ] = CELL->GetPointId(0);
    LID[ i++ ] = CELL->GetPointId(1);
    LID[ i++ ] = CELL->GetPointId(1);
    LID[ i++ ] = CELL->GetPointId(2);
    LID[ i++ ] = CELL->GetPointId(2);
    LID[ i++ ] = CELL->GetPointId(0);
    
    while( L < nLINES  &&  LIDS[ ORDER[L] ] == cellid ){
      CELL = LINES->GetCell( ORDER[ L ] );
      LID[ i++ ] = CELL->GetPointId(0) + OriginalNumberOfPoints;
      LID[ i++ ] = CELL->GetPointId(1) + OriginalNumberOfPoints;
      L++;
    }
  

  try{
// 	  mexPrintf("splitFACE face AT %i, cellid %i, LID %i, i %i \n",  AT , cellid , LID , i);
// 	  myFlush();
    splitFACE( AT , cellid , LID , i );
  }catch( int k ) { 
//     mexPrintf("splitFACE dio error \n" );    
    throw(k);
  }
  catch( ... ) { 
//     mexPrintf("splitFACE dio error \n" );
    throw;
 }
  
    
  } while( L < nLINES );
  
  vtkSignedCharArray *SIDE = NULL;

  SIDE = vtkSignedCharArray::New();
  SIDE->SetName("SIDE");
  SIDE->SetNumberOfComponents(1);
  SIDE->SetNumberOfTuples( AT->GetNumberOfCells() );
  
  
  memcpy( (signed char *) SIDE->GetVoidPointer(0)           ,
          A->GetCellData()->GetArray(0)->GetVoidPointer(0)  ,
          OriginalNumberOfCells                             );
          
  memset( (signed char *) SIDE->GetVoidPointer(0) + OriginalNumberOfCells ,
          1                                                               ,
          AT->GetNumberOfCells() - OriginalNumberOfCells                  );

  AT->GetCellData()->AddArray( SIDE );
  SIDE->Delete();  
  
  return( AT );
}






vtkPolyData *PrepareCells( vtkPolyData *A, double *BOUNDS ){
  double boxcenter[3];
  boxcenter[0]   = ( BOUNDS[0] + BOUNDS[1] )/2.0;
  boxcenter[1]   = ( BOUNDS[2] + BOUNDS[3] )/2.0;
  boxcenter[2]   = ( BOUNDS[4] + BOUNDS[5] )/2.0;
  
  double boxhalfsize[3];
  boxhalfsize[0] = ( BOUNDS[1] - BOUNDS[0] )/2.0;
  boxhalfsize[1] = ( BOUNDS[3] - BOUNDS[2] )/2.0;
  boxhalfsize[2] = ( BOUNDS[5] - BOUNDS[4] )/2.0;
  
  vtkPolyData *S = vtkPolyData::New();
  
  S->SetPoints( A->GetPoints() );

  vtkCell *CELL;
  
  int  nF = 0;
  char *TS;
  TS = (char *) mxMalloc( sizeof(char) * A->GetNumberOfCells() );
  memset( TS , 0 , A->GetNumberOfCells() );
  
  
  signed char *SIDE = (signed char *) A->GetCellData()->GetArray(0)->GetVoidPointer(0);
  memset( SIDE , 0 , A->GetNumberOfCells() );
  
  double RA[3], PA[3], QA[3];
  for( int t = 0 ; t < A->GetNumberOfCells() ; t++ ){
    CELL = A->GetCell(t);
    A->GetPoint( CELL->GetPointId(0) , RA );
    A->GetPoint( CELL->GetPointId(1) , PA );
    A->GetPoint( CELL->GetPointId(2) , QA );
    
    if( triBoxOverlap( boxcenter , boxhalfsize , RA , PA , QA ) ){
      SIDE[t] = 1;
      TS[t] = 1;
      nF++;
    }
  }
  
  vtkCellArray *FACES  = vtkCellArray::New();
//   FACES->SetNumberOfCells( nF );
	vtkIdTypeArray *FACESptr = vtkIdTypeArray::New();
  FACESptr->SetNumberOfValues( 4*nF );
  
  vtkIntArray *IDS = vtkIntArray::New();
  IDS->SetName("IDS");
  IDS->SetNumberOfComponents(1);
  IDS->SetNumberOfTuples( nF );
  int *IDSptr = (int *)IDS->GetVoidPointer(0);
  
  nF = 0;
  for( int t = 0 ; t < A->GetNumberOfCells() ; t++ ){
    if( !TS[t] ){ continue; }
    
    CELL = A->GetCell(t);
    
    FACESptr->SetValue( nF*4   , 3 );
    FACESptr->SetValue( nF*4+1 , CELL->GetPointId(0) );
    FACESptr->SetValue( nF*4+2 , CELL->GetPointId(1) );
    FACESptr->SetValue( nF*4+3 , CELL->GetPointId(2) );
    IDSptr[nF] = t;
    
    nF++;
  }

  FACES->SetCells( nF , FACESptr );
  S->SetPolys( FACES );
  S->GetCellData()->AddArray( IDS );
  
  mxFree( TS );
  FACES->Delete();
  FACESptr->Delete();
  IDS->Delete();
  
  return(S);
}


int FindIntersections( vtkOBBNode *NODE_A , vtkOBBNode *NODE_B , vtkMatrix4x4 *transform , void *arg ){
//   struct IntersectionsFinder *INT = reinterpret_cast<struct IntersectionsFinder*>(arg);
  struct IntersectionsFinder *INT = (struct IntersectionsFinder *)arg;
  
  double      RA[3] , PA[3] , QA[3];
  double      RB[3] , PB[3] , QB[3];
  
  int         intersects = 0;
  
  vtkCell     *CELL;
  vtkIdType   cA , cB, i, j;
  int         BAalreadycomputed = 0;
  double      BA[6], boxcenterA[3], boxhalfsizeA[3];
  
  int *IDSA = (int *) INT->AC->GetCellData()->GetArray(0)->GetVoidPointer(0);
  int *IDSB = (int *) INT->BC->GetCellData()->GetArray(0)->GetVoidPointer(0);

  vtkPoints     *POINTSA = INT->A->GetPoints();
  vtkPoints     *POINTSB = INT->B->GetPoints();

  int nXYZ = 0;
  
  int nA = NODE_A->Cells->GetNumberOfIds();
  for( i = 0 ; i < nA ; i++ ){
    cA = IDSA[ NODE_A->Cells->GetId(i) ];
    
    CELL = INT->A->GetCell( cA );
    POINTSA->GetPoint( CELL->GetPointId(0) , RA );
    POINTSA->GetPoint( CELL->GetPointId(1) , PA );
    POINTSA->GetPoint( CELL->GetPointId(2) , QA );
    
    if( !INT->OBB_B->TriangleIntersectsNode( NODE_B , RA , PA , QA , transform ) ){
      continue;
    }
    
    BAalreadycomputed = 0;
    int nB = NODE_B->Cells->GetNumberOfIds();
    for( j = 0 ; j < nB ; j++ ){
      cB = IDSB[ NODE_B->Cells->GetId(j) ];
      
      CELL = INT->B->GetCell( cB );
      POINTSB->GetPoint( CELL->GetPointId(0) , RB );
      POINTSB->GetPoint( CELL->GetPointId(1) , PB );
      POINTSB->GetPoint( CELL->GetPointId(2) , QB );

//       mexPrintf("cA = %d , cB = %d\n", cA, cB);
      if( !BAalreadycomputed ){
        BA[0] = MIN3( RA[0] , PA[0] , QA[0] );
        BA[1] = MAX3( RA[0] , PA[0] , QA[0] );
        BA[2] = MIN3( RA[1] , PA[1] , QA[1] );
        BA[3] = MAX3( RA[1] , PA[1] , QA[1] );
        BA[4] = MIN3( RA[2] , PA[2] , QA[2] );
        BA[5] = MAX3( RA[2] , PA[2] , QA[2] );
        boxcenterA[0]   = ( BA[0] + BA[1] )/2.0;
        boxcenterA[1]   = ( BA[2] + BA[3] )/2.0;
        boxcenterA[2]   = ( BA[4] + BA[5] )/2.0;
        boxhalfsizeA[0] = ( BA[1] - BA[0] )/2.0; boxhalfsizeA[0] += TOLBB/2;
        boxhalfsizeA[1] = ( BA[3] - BA[2] )/2.0; boxhalfsizeA[1] += TOLBB/2;
        boxhalfsizeA[2] = ( BA[5] - BA[4] )/2.0; boxhalfsizeA[2] += TOLBB/2;
        BAalreadycomputed = 1;
      }
      
      if( !triBoxOverlap( boxcenterA , boxhalfsizeA , RB , PB , QB ) ){
        continue;
      }
      
    
      intersects = TriangleTriangleIntersection( RA , PA , QA , RB , PB , QB , INT->LINES , cA , cB );
      
      if (intersects==-1)
      {     

           CoplanarIntersection( RA , PA , QA , RB , PB , QB , INT->LINES , cA , cB );
//         CoplanarIntersection( RA , PA , QA , RB , PB , QB , INT->LINES , cA , (vtkIdType) -1 );
//         CoplanarIntersection( RB , PB , QB , RA , PA , QA , INT->LINES , (vtkIdType) -1 , cB );  
      }  
  
      else 
      {
           continue;
      }
      
                  
    }
  }
  
  return( 1 );
}



int TriangleTriangleIntersection(double p1[3], double q1[3], double r1[3],
                                 double p2[3], double q2[3], double r2[3],
                                 vtkPolyData *LINES , vtkIdType cA, vtkIdType cB ) {
                                   
  double n1[3], n2[3];

  // Compute supporting plane normals.
  vtkTriangle::ComputeNormal(p1, q1, r1, n1);
  vtkTriangle::ComputeNormal(p2, q2, r2, n2);
  double s1 = -vtkMath::Dot(n1, p1);
  double s2 = -vtkMath::Dot(n2, p2);
  if (vtkMath::Dot(n1, n2)<0) {
      n2[0]=-n2[0];
      n2[1]=-n2[1];
      n2[2]=-n2[2];
      s2=-s2;}

  // Compute signed distances of points p1, q1, r1 from supporting
  // plane of second triangle.
  double dist1[3];
  dist1[0] = vtkMath::Dot(n2, p1) + s2;
  dist1[1] = vtkMath::Dot(n2, q1) + s2;
  dist1[2] = vtkMath::Dot(n2, r1) + s2;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist1[0]*dist1[1] > TOLDIST*TOLDIST ) && (dist1[0]*dist1[2] > TOLDIST*TOLDIST ))
    {
    return 0;
    }
  // Do the same for p2, q2, r2 and supporting plane of first
  // triangle.
  double dist2[3];
  dist2[0] = vtkMath::Dot(n1, p2) + s1;
  dist2[1] = vtkMath::Dot(n1, q2) + s1;
  dist2[2] = vtkMath::Dot(n1, r2) + s1;

  // If signs of all points are the same, all the points lie on the
  // same side of the supporting plane, and we can exit early.
  if ((dist2[0]*dist2[1] > TOLDIST*TOLDIST ) && (dist2[0]*dist2[2] > TOLDIST*TOLDIST ))
    
    {
    return 0;
    }
  // Check for coplanarity of the supporting planes.
  if ( fabs( n1[0] - n2[0] ) < 1e-9 &&
       fabs( n1[1] - n2[1] ) < 1e-9 &&
       fabs( n1[2] - n2[2] ) < 1e-9 &&
       fabs( s1 - s2 ) < 1e-9 )
    {
    return -1; //Coplanar
    }

 
  // There are more efficient ways to find the intersection line (if
  // it exists), but this is clear enough.
  double *pts1[3] = {p1, q1, r1}, *pts2[3] = {p2, q2, r2};

  // Find line of intersection (L = p + t*v) between two planes.
  double n1n2 = vtkMath::Dot(n1, n2);
  double a = (s1 - s2*n1n2) / (n1n2*n1n2 - 1.0);
  double b = (s2 - s1*n1n2) / (n1n2*n1n2 - 1.0);
  double p[3], v[3];
  p[0] = a*n1[0] + b*n2[0];
  p[1] = a*n1[1] + b*n2[1];
  p[2] = a*n1[2] + b*n2[2];
  vtkMath::Cross(n1, n2, v);
  vtkMath::Normalize( v );
  
    
  int index1 = 0, index2 = 0;
  double t1[1024], t2[1024];
  for (int i = 0; i < 3; i++)
    {
    double t, x[3];
    int j = i, id2 = (i+1) % 3;
    
    #define dot(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
    // Find t coordinate on line of intersection between two planes.

    if (vtkPlane::IntersectWithLine( pts1[j], pts1[id2], n2, p2, t, x ))
      {
      t1[index1++] = dot(x, v) - dot(p, v);
      }
    if (vtkPlane::IntersectWithLine( pts2[j], pts2[id2], n1, p1, t, x ))
      {
      t2[index2++] = dot(x, v) - dot(p, v);
      }
    }
  
  // Check if only one edge or all edges intersect the supporting
  //   planes intersection.


  if (( index1 > 3 ) || (index2 > 3 ) || ( index1 < 2 )||( index2 < 2 ))
    {
//       mexPrintf("Hay muchos\n");  
      return 0;
    }
  else
  {
  
      if (index1==3)
      {
          if (t1[0]==t1[1]){t1[1]=t1[2];}
      }           
      if (index2==3)
      {
          if (t2[0]==t2[1]){t2[1]=t2[2];}
      }
      
  }
       
      if ((t1[0]==t1[1]) || (t2[0]==t2[1]))
      {  
          return 0;
      }
      

  // Check for NaNs
  if ( vtkMath::IsNan( t1[0] ) || vtkMath::IsNan( t1[1] ) ||
       vtkMath::IsNan( t2[0] ) || vtkMath::IsNan( t2[1] ) )
    {
    return 0;
    }

  if ( t1[0] > t1[1] )
    {
    std::swap( t1[0], t1[1] );
    }
  if ( t2[0] > t2[1] )
    {
    std::swap( t2[0], t2[1] );
    }
  // Handle the different interval configuration cases.
  double tt1, tt2;

  if ( t1[1] < t2[0] || t2[1] < t1[0] )
    {
    return 0; // No overlap
    }
  else if ( t1[0] < t2[0] )
    {
    if ( t1[1] < t2[1] ) // Case 1
      {
      tt1 = t2[0] ;
      tt2 = t1[1] ;
      }
    else // Case 2
      {
      tt1 = t2[0] ;
      tt2 = t2[1] ;
      }
    }
  else // t1[0] >= t2[0]
    {
    if ( t1[1] < t2[1] ) // Case 3
      {
      tt1 = t1[0] ;
      tt2 = t1[1] ;
      }
    else // Case 4
      {
      tt1 = t1[0] ;
      tt2 = t2[1] ;
      }
    }


  double      X[3] , Y[3];
  X[0] = p[0] + tt1*v[0];
  X[1] = p[1] + tt1*v[1];
  X[2] = p[2] + tt1*v[2];

  Y[0] = p[0] + tt2*v[0];
  Y[1] = p[1] + tt2*v[1];
  Y[2] = p[2] + tt2*v[2];
  
  
  int nXYZ = 0;
  
  vtkPoints     *LINESxyz = LINES->GetPoints();
  vtkCellArray  *LINESlin = LINES->GetLines();
  vtkDataArray  *LidsA = LINES->GetCellData()->GetArray(0);
  vtkDataArray  *LidsB = LINES->GetCellData()->GetArray(1);
  
  if( (X[0] != Y[0] || X[1] != Y[1] || X[2] != Y[2] ) ){
          nXYZ = LINESxyz->GetNumberOfPoints();
          LINESxyz->InsertNextPoint( X );
          LINESxyz->InsertNextPoint( Y );

          LINESlin->InsertNextCell( 2 );
          LINESlin->InsertCellPoint( nXYZ   );
          LINESlin->InsertCellPoint( nXYZ+1 );

          LidsA->InsertNextTuple1( cA );
          LidsB->InsertNextTuple1( cB );
          
          return 1;
      }

}




void CoplanarIntersection( double RA[3], double PA[3], double QA[3], 
                           double RB[3], double PB[3], double QB[3],
                           vtkPolyData *LINES , vtkIdType cA, vtkIdType cB ) {

//       mexPrintf("Son coplanares las caras %i y %i\n",cA,cB); myFlush();
    
  double      T[3][3], Den, RP[3], RQ[3], X[3], Y[3] ;
  double      x, y ;
  int nXYZ = 0;
  int nV, k;
  
  gpc_polygon TRI1;
  gpc_polygon TRI2;
  gpc_polygon TRIO;
  
  vtkPoints     *LINESxyz = LINES->GetPoints();
  vtkCellArray  *LINESlin = LINES->GetLines();
  vtkDataArray  *LidsA = LINES->GetCellData()->GetArray(0);
  vtkDataArray  *LidsB = LINES->GetCellData()->GetArray(1);
  
    
  for ( k = 0 ; k < 3 ; k++ )

      {
      RP[k]=PA[k]-RA[k];
      RQ[k]=QA[k]-RA[k];
      }

  Den=(RP[0]*RP[0]*RQ[1]*RQ[1] + RP[0]*RP[0]*RQ[2]*RQ[2] - 2*RP[0]*RP[1]*RQ[0]*RQ[1] - 2*RP[0]*RP[2]*RQ[0]*RQ[2] + RP[1]*RP[1]*RQ[0]*RQ[0] + RP[1]*RP[1]*RQ[2]*RQ[2] - 2*RP[1]*RP[2]*RQ[1]*RQ[2] + RP[2]*RP[2]*RQ[0]*RQ[0] + RP[2]*RP[2]*RQ[1]*RQ[1]);

  if( Den < 1e-15 ){
    return;
  }

  // Transformation Matrix

  Den = 1.0/Den;

  T[0][0]=-(RQ[0]*(RP[1]*RQ[1] + RP[2]*RQ[2]) - RP[0]*(RQ[1]*RQ[1] + RQ[2]*RQ[2]))*Den;
  T[0][1]=-(RQ[1]*(RP[0]*RQ[0] + RP[2]*RQ[2]) - RP[1]*(RQ[0]*RQ[0] + RQ[2]*RQ[2]))*Den;
  T[0][2]=-(RQ[2]*(RP[0]*RQ[0] + RP[1]*RQ[1]) - RP[2]*(RQ[0]*RQ[0] + RQ[1]*RQ[1]))*Den;

  T[1][0]=-(RP[0]*(RP[1]*RQ[1] + RP[2]*RQ[2]) - RQ[0]*(RP[1]*RP[1] + RP[2]*RP[2]))*Den;
  T[1][1]=-(RP[1]*(RP[0]*RQ[0] + RP[2]*RQ[2]) - RQ[1]*(RP[0]*RP[0] + RP[2]*RP[2]))*Den;
  T[1][2]=-(RP[2]*(RP[0]*RQ[0] + RP[1]*RQ[1]) - RQ[2]*(RP[0]*RP[0] + RP[1]*RP[1]))*Den;


  TRI1.num_contours = 1;

  MALLOC( TRI1.hole    , sizeof(int)            , "hole flag array creation" , int             );
  MALLOC( TRI1.contour , sizeof(gpc_vertex_list), "contour creation"         , gpc_vertex_list );
  MALLOC( TRI1.contour[0].vertex , 3 * sizeof(gpc_vertex) , "vertex creation", gpc_vertex );

  TRI1.hole[0] = 0;
  TRI1.contour[0].num_vertices = 3;


  TRI1.contour[0].vertex[0].x = rx;
  TRI1.contour[0].vertex[0].y = ry;
  TRI1.contour[0].vertex[1].x = rx+1.0;
  TRI1.contour[0].vertex[1].y = ry;
  TRI1.contour[0].vertex[2].x = rx;
  TRI1.contour[0].vertex[2].y = ry+1.0;


  TRI2.num_contours = 1;

  MALLOC( TRI2.hole    , sizeof(int)            , "hole flag array creation" , int             );
  MALLOC( TRI2.contour , sizeof(gpc_vertex_list), "contour creation"         , gpc_vertex_list );
  MALLOC( TRI2.contour[0].vertex , 3 * sizeof(gpc_vertex) , "vertex creation", gpc_vertex );

  TRI2.hole[0] = 0;
  TRI2.contour[0].num_vertices = 3;


  TRI2.contour[0].vertex[0].x = T[0][0]*(RB[0]-RA[0])+T[0][1]*(RB[1]-RA[1])+T[0][2]*(RB[2]-RA[2])+rx;
  TRI2.contour[0].vertex[0].y = T[1][0]*(RB[0]-RA[0])+T[1][1]*(RB[1]-RA[1])+T[1][2]*(RB[2]-RA[2])+ry;

  TRI2.contour[0].vertex[1].x = T[0][0]*(PB[0]-RA[0])+T[0][1]*(PB[1]-RA[1])+T[0][2]*(PB[2]-RA[2])+rx;
  TRI2.contour[0].vertex[1].y = T[1][0]*(PB[0]-RA[0])+T[1][1]*(PB[1]-RA[1])+T[1][2]*(PB[2]-RA[2])+ry;

  TRI2.contour[0].vertex[2].x = T[0][0]*(QB[0]-RA[0])+T[0][1]*(QB[1]-RA[1])+T[0][2]*(QB[2]-RA[2])+rx;
  TRI2.contour[0].vertex[2].y = T[1][0]*(QB[0]-RA[0])+T[1][1]*(QB[1]-RA[1])+T[1][2]*(QB[2]-RA[2])+ry;

  gpc_polygon_clip( GPC_INT , &TRI1 , &TRI2 , &TRIO );

  gpc_free_polygon( &TRI1 );
  gpc_free_polygon( &TRI2 );


  if ( ( TRIO.num_contours > 0 ) && ( TRIO.contour[0].num_vertices ) )

  {

  nV = TRIO.contour[0].num_vertices;

  for (k=0; k<nV; k++)
      { 
      x = TRIO.contour[0].vertex[k].x;
      y = TRIO.contour[0].vertex[k].y;
      X[0] = RA[0]+ RP[0]*(x-rx) + RQ[0]*(y-ry); 
      X[1] = RA[1]+ RP[1]*(x-rx) + RQ[1]*(y-ry); 
      X[2] = RA[2]+ RP[2]*(x-rx) + RQ[2]*(y-ry); 

      x = TRIO.contour[0].vertex[(k+1)%nV].x;
      y = TRIO.contour[0].vertex[(k+1)%nV].y;
      Y[0] = RA[0]+ RP[0]*(x-rx) + RQ[0]*(y-ry); 
      Y[1] = RA[1]+ RP[1]*(x-rx) + RQ[1]*(y-ry); 
      Y[2] = RA[2]+ RP[2]*(x-rx) + RQ[2]*(y-ry); 

      nXYZ = LINESxyz->GetNumberOfPoints();
      LINESxyz->InsertNextPoint( X );
      LINESxyz->InsertNextPoint( Y );

      LINESlin->InsertNextCell( 2 );
      LINESlin->InsertCellPoint( nXYZ   );
      LINESlin->InsertCellPoint( nXYZ+1 );

      LidsA->InsertNextTuple1( cA );
      LidsB->InsertNextTuple1( cB );
      }

  }    

  gpc_free_polygon( &TRIO );

}



void GetOrder( int *X , int I , int *ORDER ){
  toSORT = X;
  for( int t = 0 ; t < I ; t++ ){
    ORDER[t] = t;
  }
  qsort( ORDER , I , sizeof( int ) , func_to_sort );
}


mxArray *Split2MESH( vtkPolyData *POLY , int side , mxArray *copyFrom ){

  //Para generar las salidas de matlab
//   funcCreate toShared;
//   toShared = GetCreatePointer();
  const char      *names[] = {""};
  const int       dims[1] = {1};
  mxArray         *DATA;
  double          *data, xyz[3];
  int             nP, c, p, nC, val;
  signed char     *SIDE;
  vtkPoints       *POINTS;
  vtkCell         *CELL;
  
  mxArray         *M;
  
  M = mxCreateStructArray(1, (const mwSize*)(dims), 0, names );

  if( POLY == NULL ){
    mxAddField( M , "xyz" );
    mxAddField( M , "tri" );
    
    return( M );
  }
  
  

  
  if( copyFrom == NULL ){

    nP = POLY->GetNumberOfPoints();
    POINTS = POLY->GetPoints();
    
    DATA = mxCreateDoubleMatrix( nP , 3 , mxREAL );  
    data = mxGetPr( DATA );

    for( int p = 0 ; p < nP ; p++ ){
      POINTS->GetPoint( p , xyz );
      data[ p        ] = xyz[0];
      data[ p +   nP ] = xyz[1];
      data[ p + 2*nP ] = xyz[2];
    }

    mxAddField( M ,     "xyz" );
    mxSetField( M , 0 , "xyz" , DATA );

  } else {
    
//     DATA = mxGetField( copyFrom , "xyz" );
//     mxSetField( M , 0 , "xyz" , ( (void *) GetCreatePointer() )( DATA ) );
    //make the SharedCopy
  }
    
  nP = POLY->GetNumberOfCells();
  SIDE = (signed char *) POLY->GetCellData()->GetArray(0)->GetVoidPointer(0);

  nC = 0;
  for( int p = 0 ; p < nP ; p++ ){
    val = SIDE[ p ];
    if( val == -2                   ||
        ( side == 0 && val > 0 )    ||
        ( side == 1 && val != 1 )   ||
        ( side == 21 && val != 21 ) ||
        ( side == 22 && val != 22 )
      ){ 
      continue;
    }
    nC++;
  }
  
  
  DATA = mxCreateDoubleMatrix( nC , 3 , mxREAL );  
  data = mxGetPr( DATA );

  c = 0;
  for( int p = 0 ; p < nP ; p++ ){
    val = SIDE[ p ];
    if( val == -2                   ||
        ( side == 0 && val > 0 )    ||
        ( side == 1 && val != 1 )   ||
        ( side == 21 && val != 21 ) ||
        ( side == 22 && val != 22 )
      ){ 
      continue;
    }

    CELL = POLY->GetCell( p );
    data[ c        ] = CELL->GetPointId(0) + 1;
    data[ c +   nC ] = CELL->GetPointId(1) + 1;
    data[ c + 2*nC ] = CELL->GetPointId(2) + 1;
    c++;
  }

  mxAddField( M ,     "tri"  );
  mxSetField( M , 0 , "tri" , DATA );

  
  return( M );
}




void setSIDE( vtkPolyData *A , vtkPolyData *B ){

  vtkPolyDataNormals *NORMALS = NULL;
  NORMALS = vtkPolyDataNormals::New();
    NORMALS->SetInput( B );
    NORMALS->ComputeCellNormalsOn();
    NORMALS->ComputePointNormalsOn();
    NORMALS->ConsistencyOff();
    NORMALS->AutoOrientNormalsOff();
    NORMALS->SplittingOff();
    //NORMALS->NonManifoldTraversalOff();
    //NORMALS->FlipNormalsOn();
    NORMALS->Update();
  
    
  vtkPolyData *SURF;
  SURF = NORMALS->GetOutput();

  vtkCellLocator      *LOC = NULL;
  LOC= vtkCellLocator::New();
    LOC->SetDataSet( SURF );
 	  LOC->CacheCellBoundsOn();
    LOC->SetNumberOfCellsPerBucket( 2 );
    LOC->BuildLocator();

  
  
  signed char *SIDE = (signed char *)A->GetCellData()->GetArray( 0 )->GetVoidPointer(0);
  
  
  double R[3], P[3], Q[3], C[3];
  int t, c;
  double              closestPoint[3], d, distance, pcoords[3], N[3], *NN, weigths[3], side;
  vtkGenericCell      *Gcell = vtkGenericCell::New();
  vtkIdType           cellid;
  int                 sub;
  vtkIdList           *NODES = vtkIdList::New();
  vtkCell             *CELL = NULL;
  

  vtkPoints *POINTS;
  POINTS = A->GetPoints();

  for( t = 0 ; t < A->GetNumberOfCells() ; t++ ){
    if( !SIDE[ t ] ){ continue; }
    CELL = A->GetCell(t);
    
    POINTS->GetPoint( CELL->GetPointId(0) , R );
    POINTS->GetPoint( CELL->GetPointId(1) , P );
    POINTS->GetPoint( CELL->GetPointId(2) , Q );
    
    N[0] = ((P[1]-R[1])*(Q[2]-R[2])-(P[2]-R[2])*(Q[1]-R[1]));
    N[1] = ((P[2]-R[2])*(Q[0]-R[0])-(P[0]-R[0])*(Q[2]-R[2]));
    N[2] = ((P[0]-R[0])*(Q[1]-R[1])-(P[1]-R[1])*(Q[0]-R[0]));
     
    if( ( N[0]*N[0] + N[1]*N[1] + N[2]*N[2] ) < TOLAREA*TOLAREA ){ SIDE[ t ] = -2; continue; }

    for( c = 0 ; c < 3 ; c++ ){
      C[c] = ( R[c] + P[c] + Q[c] )/3.0;
    }
    
    LOC->FindClosestPoint( C , closestPoint , Gcell , cellid , sub , distance );
    
    SURF->GetCell( cellid )->EvaluatePosition( closestPoint, NULL, sub, pcoords, d, weigths );
    
//     mexPrintf("distancia entre %i y %i es %.20f o %.20f\n",cellid,t,distance, d);
    

    if( weigths[0]<1e-3 || weigths[1]<1e-3 || weigths[2]<1e-3 ){
    
      SURF->GetCellPoints( cellid , NODES );
      
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(0) );
      for( c=0 ; c<3; c++ ){  N[c]  = NN[c]*weigths[0];  }
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(1) );
      for( c=0 ; c<3; c++ ){  N[c] += NN[c]*weigths[1];  }
      NN = SURF->GetPointData()->GetNormals()->GetTuple3( NODES->GetId(2) );
      for( c=0 ; c<3; c++ ){  N[c] += NN[c]*weigths[2];  }

    } else {
      
      NN = SURF->GetCellData()->GetNormals()->GetTuple3(cellid);
      // Mirar coplanaridad entre N y NN si la distancia es 0
      
      if ( (distance < TOLDIST*TOLDIST) &&
           ( fabs(N[1]*NN[2]-N[2]*NN[1]) < 1e-9 ) &&
           ( fabs(N[2]*NN[0]-N[0]*NN[2]) < 1e-9 ) &&
           ( fabs(N[0]*NN[1]-N[1]*NN[0]) < 1e-9 ) )
      {
//         mexPrintf("Son coplanares las nuevas caras %i y %i\n",t,cellid); myFlush();  
          
        if ( N[0]*NN[0]+N[1]*NN[1]+N[2]*NN[2] > 0.0)
        { SIDE[ t ] = 21;  }
        else
        { SIDE[ t ] = 22;  }
        
        continue;
      }
      else
      {
      N[0] = NN[0];
      N[1] = NN[1];
      N[2] = NN[2];
      }
    }

    side= 0;
    for( c=0; c<3; c++ ){ side += (closestPoint[c]-C[c])*N[c]; }
    
    
    if( side > 0 ){
      SIDE[ t ] = 1;
    } else {
      SIDE[ t ] = -1;
    }
    
  }
  
  SAFEDELETE( NODES   );
  SAFEDELETE( Gcell   );
  SAFEDELETE( LOC     );
  SAFEDELETE( NORMALS );
}


















// DAVID
int LinePlaneIntersection(double p1[3], double	p2[3], double n[3], double 	p0[3], double &	t, double x[3]){

    double vR[3];
     
    vR[0]=p2[0]-p1[0];
    vR[1]=p2[1]-p1[1];
    vR[2]=p2[2]-p1[2];
    double Vd=dot(n,vR);
    if (Vd==0)
    {
//         mexPrintf("Paralelas");
        return 0;
        
    }
    
    double Vo = - dot(n,p1) + dot(n,p0);
    t =  Vo / Vd;
    
    
//       mexPrintf("p1=%f %f %f\n",p1[0],p1[1],p1[2]);
//       mexPrintf("p2=%f %f %f\n",p2[0],p2[1],p2[2]);
//       mexPrintf("n=%f %f %f\n",n[0],n[1],n[2]);
//       mexPrintf("p0=%f %f %f\n",p0[0],p0[1],p0[2]);
//       mexPrintf("vR=%f %f %f\n",vR[0],vR[1],vR[2]);
//       mexPrintf("Vd=%f \n",Vd);
//       mexPrintf("Vo=%f \n",Vo);
      
      
    if ((t>=0)&&(t<=1))
        {
        x[0]= p1[0]+t*vR[0];
        x[1]= p1[1]+t*vR[1];
        x[2]= p1[2]+t*vR[2];
        
        return 1;}
    else
    {
        
        return 0;
    }
}
   









#define X 0
#define Y 1
#define Z 2
#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 
#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;
/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)			                        \
	p0 = a*v0[Y] - b*v0[Z];                                       \
	p2 = a*v2[Y] - b*v2[Z];                                       \
  if( p0<p2 ){ min=p0; max=p2; }else{ min=p2; max=p0; }         \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];              \
	if( min>rad || max<-rad ){ return(0); };

#define AXISTEST_X2(a, b, fa, fb)                               \
	p0 = a*v0[Y] - b*v0[Z];                                       \
	p1 = a*v1[Y] - b*v1[Z];			       	                          \
  if( p0<p1 ){ min=p0; max=p1; }else{ min=p1; max=p0; }         \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];              \
	if( min>rad || max<-rad ){ return(0); }

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)                              \
	p0 = -a*v0[X] + b*v0[Z];                                      \
	p2 = -a*v2[X] + b*v2[Z];                                      \
  if( p0<p2 ){ min=p0; max=p2; }else{ min=p2; max=p0; }         \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];              \
	if( min>rad || max<-rad ){ return(0); };

#define AXISTEST_Y1(a, b, fa, fb)                               \
	p0 = -a*v0[X] + b*v0[Z];                                      \
	p1 = -a*v1[X] + b*v1[Z];                                      \
  if( p0<p1 ){ min=p0; max=p1; }else{ min=p1; max=p0; }         \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];              \
	if( min>rad || max<-rad ){ return(0); };

/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)                              \
	p1 = a*v1[X] - b*v1[Y];                                       \
	p2 = a*v2[X] - b*v2[Y];                                       \
  if( p2<p1 ){ min=p2; max=p1; }else{ min=p1; max=p2; }         \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];              \
	if( min>rad || max<-rad ){ return(0); }

#define AXISTEST_Z0(a, b, fa, fb)                               \
	p0 = a*v0[X] - b*v0[Y];                                       \
	p1 = a*v1[X] - b*v1[Y];                                       \
  if( p0<p1 ){ min=p0; max=p1; }else{ min=p1; max=p0; }         \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];              \
	if( min>rad || max<-rad ){ return(0); };
int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3]){
  int q;
  double vmin[3],vmax[3],v;
  for(q=X;q<=Z;q++){
    v=vert[q];					// -NJMP-
    if(normal[q]>0.0f){
      vmin[q]=-maxbox[q] - v;	// -NJMP-
      vmax[q]= maxbox[q] - v;	// -NJMP-
    } else {
      vmin[q]= maxbox[q] - v;	// -NJMP-
      vmax[q]=-maxbox[q] - v;	// -NJMP-
    }
  }
  if( DOT(normal,vmin) >  0.0f ){ return(0); }
  if( DOT(normal,vmax) >= 0.0f ){ return(1); }

  return(0);
}
int triBoxOverlap( double boxcenter[3] , double boxhalfsize[3] , double R[3] , double P[3] , double Q[3] ){
  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */

   double v0[3],v1[3],v2[3];
   double min,max,p0,p1,p2,rad,fex,fey,fez;
   double normal[3],e0[3],e1[3],e2[3];


   /* This is the fastest branch on Sun */
   /* move everything so that the boxcenter is in (0,0,0) */
   SUB(v0,R,boxcenter);
   SUB(v1,P,boxcenter);
   SUB(v2,Q,boxcenter);

   /* compute triangle edges */
   SUB(e0,v1,v0);      /* tri edge 0 */
   SUB(e1,v2,v1);      /* tri edge 1 */
   SUB(e2,v0,v2);      /* tri edge 2 */

   /* Bullet 3:  */
   /*  test the 9 tests first (this was faster) */
   fex = fabsf(e0[X]);
   fey = fabsf(e0[Y]);
   fez = fabsf(e0[Z]);

   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
   AXISTEST_Z12(e0[Y], e0[X], fey, fex);

   fex = fabsf(e1[X]);
   fey = fabsf(e1[Y]);
   fez = fabsf(e1[Z]);

   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
   AXISTEST_Z0( e1[Y], e1[X], fey, fex);

   fex = fabsf(e2[X]);
   fey = fabsf(e2[Y]);
   fez = fabsf(e2[Z]);

   AXISTEST_X2( e2[Z], e2[Y], fez, fey);
   AXISTEST_Y1( e2[Z], e2[X], fez, fex);
   AXISTEST_Z12(e2[Y], e2[X], fey, fex);

   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */

   /* test in X-direction */
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if( min>boxhalfsize[X] || max<-boxhalfsize[X] ){ return(0); };

   /* test in Y-direction */
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if( min>boxhalfsize[Y] || max<-boxhalfsize[Y]){ return(0); };

   /* test in Z-direction */
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if( min>boxhalfsize[Z] || max<-boxhalfsize[Z] ){ return(0); };

   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */
   CROSS(normal,e0,e1);
   if( !planeBoxOverlap(normal,v0,boxhalfsize) ){ return(0); };

   return(1);   /* box and triangle overlaps */
}

