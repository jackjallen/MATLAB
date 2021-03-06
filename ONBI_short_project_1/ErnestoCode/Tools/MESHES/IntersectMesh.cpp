#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkBooleanOperationPolyDataFilter

#include "vtkBooleanOperationPolyDataFilter.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyData2MESH.h"



#include "vtkPolyDataWriter.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  vtkPolyDataWriter *W = vtkPolyDataWriter::New();
  W->SetFileTypeToASCII();
  
  int                 argN,N,NF,i,n;
  double              v, *xyz,x[3],y[3],z[3],*Normal;
  char                STR[2000], method[2000];

  if(!nrhs){
    /*
    mexPrintf("vtkBooleanOperationPolyDataFilter( MESH ...\n");
    mexPrintf("\n");
    mexPrintf("SetOperation             	, int       	... Set the boolean operation to perform: Union (0), Intersection (1), Difference (2). Defaults to union. \n");
    mexPrintf("SetOperationToUnion 		, [] (*)    	...\n");
    mexPrintf("SetOperationToIntersection 	, []        	...\n");
    mexPrintf("SetOperationToDifference		, []        	...\n");
    mexPrintf("\n");
    mexPrintf("SetReorientDifferenceCells   	, int       	... Turn on/off cell reorientation of the intersection portion of the surface when the operation is set to DIFFERENCE. Defaults to on. \n");
    mexPrintf("ReorientDifferenceCellsOn 	, [] (*)    	...\n");
    mexPrintf("ReorientDifferenceCellsOff 	, []        	...\n");
    mexPrintf("\n");
    mexPrintf("SetTolerance             	, real (1e-6)   ... Set the tolerance used to determine when a point's absolute distance is considered to be zero. Defaults to 1e-6. \n");
    mexPrintf("\n");*/
    if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }
    return;
  }

//     if( nlhs ){ for (int i=0; i<nlhs; i++) plhs[i]=mxCreateDoubleMatrix( 0 , 0 , mxREAL ); }

  
  ALLOCATES();
  vtkPolyData         *MESH1,*MESH2,*S,*MESH1_OUT,*MESH2_IN;
  vtkOBJ_TYPE	      *BOOL;
  vtkCellArray	      *F_MESH1_OUT,*F_MESH2_IN;
  vtkDoubleArray      *ARRAY_ZERO,*ARRAY_ONE;

  MESH1=MESH2vtkPolyData( prhs[0] );
  
  if (MESH1->GetCellData()->GetArray("Normals",i)==NULL){mexErrMsgTxt("Error: triNormals field undefined in Mesh1");}
  
  N=MESH1->GetNumberOfCells();
  ARRAY_ZERO=vtkDoubleArray::New();
  ARRAY_ZERO->SetName( "labels" ); 
  ARRAY_ZERO->SetNumberOfComponents( 1 );
  ARRAY_ZERO->SetNumberOfTuples( N );
  for( n=0 ; n<N ; n++ ) {
      ARRAY_ZERO->SetComponent(n, 0, 0 );
  }
  MESH1->GetCellData()->AddArray( ARRAY_ZERO );
  
  
  MESH2=MESH2vtkPolyData( prhs[1] );
  if (MESH2->GetCellData()->GetArray("Normals",i)==NULL){mexErrMsgTxt("Error: triNormals field undefined in Mesh2");}
  N= MESH2->GetNumberOfCells();
  ARRAY_ONE=vtkDoubleArray::New();
  ARRAY_ONE->SetName( "labels" ); 
  ARRAY_ONE->SetNumberOfComponents( 1 );
  ARRAY_ONE->SetNumberOfTuples( N );
  for( n=0 ; n<N ; n++ ) {
      ARRAY_ONE->SetComponent(n, 0, 1 );
  }
  MESH2->GetCellData()->AddArray( ARRAY_ONE );

  
  
  BOOL=vtkOBJ_TYPE::New();
  BOOL->SetInput( 0, MESH1 );
  BOOL->SetInput( 1, MESH2 );
  
  /*Defaults*/
  BOOL->ReorientDifferenceCellsOff();
  BOOL->SetTolerance(1e-6);
  BOOL->SetOperationToDifference();
  
  /*END Defaults*/
  
  BOOL->Update();
  
  /*Parsing arguments*/
  /* argN = 2;
  while( nrhs > argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
    //   myErrMsgTxt( "No keywords." ); 
    }
    mxGetString( prhs[argN], method, 1999 );
    
    argN++;
    if( argN == nrhs || mxGetNumberOfElements( prhs[argN] ) == 0){
      CallMethod( BOOL , method );
    } else if( mxIsChar(prhs[argN]) ) {
      mxGetString( prhs[argN], STR, 1999 );
      CallMethod( BOOL , method , STR );
    } else if( mxGetNumberOfElements( prhs[argN] ) == 1 )  {
      v = myGetValue( prhs[argN] );
      CallMethod( BOOL , method , v );
    } else {
      xyz = myGetPr( prhs[argN] );
      CallMethod( BOOL , method , xyz );
    }
    argN++;

  }*/
  /*END Parsing arguments*/
  



  //W->SetInput( BOOL->GetOutput() );  W->SetFileName("bool.vtp");        W->Write();
  

  S=BOOL->GetOutput(0);	

  W->SetInput( S );  W->SetFileName("s.vtp");       W->Write();
  
 
  MESH1_OUT=vtkPolyData::New();
  MESH2_IN=vtkPolyData::New();
  
  MESH1_OUT->SetPoints(S->GetPoints());  
  MESH2_IN->SetPoints(S->GetPoints());
    
  F_MESH1_OUT=vtkCellArray::New();
  F_MESH2_IN=vtkCellArray::New();
   
  NF=S->GetNumberOfCells();

  
  for( n= 0; n<NF ; n++ ) {
      if( S->GetCell(n)->GetNumberOfPoints() == 3 ) {
            S->GetPoint(S->GetCell(n)->GetPointId(0),x);
            S->GetPoint(S->GetCell(n)->GetPointId(1),y);
            S->GetPoint(S->GetCell(n)->GetPointId(2),z);
           
            Normal=S->GetCellData()->GetArray("Normals",i)->GetTuple3(n);
                         
            if ((((y[1]-x[1])*(z[2]-x[2])-(y[2]-x[2])*(z[1]-x[1]))*Normal[0]+
                 ((y[2]-x[2])*(z[0]-x[0])-(y[0]-x[0])*(z[2]-x[2]))*Normal[1]+
                 ((y[0]-x[0])*(z[1]-x[1])-(y[1]-x[1])*(z[0]-x[0]))*Normal[2])>0)
               
            { 
                  if (S->GetCellData()->GetArray("labels",i)->GetComponent(n,0)==0)
                  {
                          //->InsertNextTuple(Normal);                     
                          F_MESH1_OUT->InsertNextCell(3);
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(0) );
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(1) );
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(2) );
                  }
                  else
                  {
                          F_MESH2_IN->InsertNextCell(3);
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(0) );
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(1) );
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(2) );
                  }
            }
      
            else
                
            {     mexWarnMsgTxt("Reordenando: "); mexPrintf(" %i\n", n );
                  if (S->GetCellData()->GetArray("labels",i)->GetComponent(n,0)==0)
                  {
                          F_MESH1_OUT->InsertNextCell(3);
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(0) );
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(1) );
                          F_MESH1_OUT->InsertCellPoint( S->GetCell(n)->GetPointId(2) );
                  }
                  else
                  {
                          F_MESH2_IN->InsertNextCell(3);
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(0) );
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(1) );
                          F_MESH2_IN->InsertCellPoint( S->GetCell(n)->GetPointId(2) );
                  }
            }
                      
                      
                      
		  }
      }
  

  
  
  MESH1_OUT->SetPolys( F_MESH1_OUT );
  MESH1_OUT->Update();
  MESH2_IN->SetPolys(  F_MESH2_IN );
  MESH2_IN->Update();

  
  
  // W->SetInput( MESH1_OUT );  W->SetFileName("kk.vtp");       W->SetFileTypeToASCII();
  // W->SetInput( MESH2_IN );  W->SetFileName("kk2.vtp");       W->SetFileTypeToASCII();
  // W->Write();
  
  
  
  
  
  plhs[0]= vtkPolyData2MESH( MESH1_OUT );
  plhs[1]= vtkPolyData2MESH( MESH2_IN );


  
  F_MESH1_OUT->Delete();
  F_MESH2_IN->Delete();
  MESH1_OUT->Delete();
  MESH2_IN->Delete();

  
 
/*  BOOL->Update();




  
  plhs[2]= vtkPolyData2MESH( BOOL->GetOutput(0) );
  plhs[3]= vtkPolyData2MESH( BOOL->GetOutput(1) );
*/




  EXIT:
    ARRAY_ZERO->Delete();
    ARRAY_ONE->Delete();
    MESH1->Delete();
    MESH2->Delete();
    BOOL->Delete();


    myFreeALLOCATES();
}

void CallMethod( vtkOBJ_TYPE *O , char *met ){
  /*[]-valued methods*/
  Call_0( SetOperationToUnion 	     );
  Call_0( SetOperationToIntersection );
  Call_0( SetOperationToDifference   );
  Call_0( ReorientDifferenceCellsOn  );
  Call_0( ReorientDifferenceCellsOff );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , real v ){
  /*Real valued methods*/
  Call_1( SetOperation   	     , v );
  Call_1( SetReorientDifferenceCells , v);
  Call_1( SetTolerance   	     , v );
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush();
}

void CallMethod( vtkOBJ_TYPE *O , char *met , char *v ){
  /*String valued methods*/
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}
void CallMethod( vtkOBJ_TYPE *O , char *met , real *v ){
  /*Array valued methods*/
  mexWarnMsgTxt("Invalid Method: "); mexPrintf(" %s\n", met ); myFlush(); 
}

