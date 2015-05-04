#include "mex.h"
#include "myMEX.h"

#define   real       double
#define   mxREAL_CLASS       mxDOUBLE_CLASS

#define vtkOBJ_TYPE      vtkMassProperties
#include "vtkMassProperties.h"
#include "MESH2vtkPolyData.h"

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){
  int                 argN;
  char                STR[2000], method[2000];
  
  
  if(!nrhs){
    mexPrintf("vtkMassProperties( MESH ...\n");
    mexPrintf("\n");

    mexPrintf("SetFeatureAngle          , real       ... Specify the angle that defines a sharp edge. If the difference in angle across neighboring polygons is greater than this value, the shared edge is considered 'sharp'.\n");
    mexPrintf("\n");
    mexPrintf("GetVolume                , Compute and return the volume.\n");
    mexPrintf("GetSurfaceArea           , Compute and return the area.\n");
    mexPrintf("GetVolumeProjected       , Compute and return the projected volume. Typically you should compare this volume to the value returned by GetVolume if you get an error (GetVolume()-GetVolumeProjected())*10000 that is greater than GetVolume() this should identify a problem: * Either the polydata is not closed * Or the polydata contains triangle that are flipped.\n");
    mexPrintf("GetMinCellArea           , Compute and return the min cell area.\n");
    mexPrintf("GetMaxCellArea           , Compute and return the max cell area.\n");
    mexPrintf("GetVolumeX               , Compute and return the volume projected on to each axis aligned plane.\n");
    mexPrintf("GetVolumeY               , .\n");
    mexPrintf("GetVolumeZ               , .\n");
    mexPrintf("GetKx                    , Compute and return the weighting factors for the maximum unit normal component (MUNC).\n");
    mexPrintf("GetKy                    , .\n");
    mexPrintf("GetKz                    , .\n");
    mexPrintf("GetNormalizedShapeIndex	, Compute and return the normalized shape index. This characterizes the deviation of the shape of an object from a sphere. A sphere's NSI is one. This number is always >= 1.0.\n");
    mexPrintf("\n");
    
    return;
  }
  
  ALLOCATES();
  vtkPolyData         *MESH;
  
  MESH= MESH2vtkPolyData( prhs[0] );
  
  vtkMassProperties *M = vtkMassProperties::New();
  M->SetInput( MESH );

  /*Defaults*/
  M->Update();
  /*END Defaults*/
  
  
  /*Parsing arguments*/
  argN = 1;
  while( nrhs > argN && MAX(1,nlhs) >= argN ) {
    if( !mxIsChar( prhs[argN] ) || !mxGetNumberOfElements( prhs[argN] ) ){
      mexPrintf( "No keywords." );
    }
    mxGetString( prhs[argN], method, 1999 );
    
    if( ! myStrcmpi(method,"volume")    || ! myStrcmpi(method,"getvolume")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetVolume() );              argN++; continue;
    }

    if( ! myStrcmpi(method,"surface")    || ! myStrcmpi(method,"surfacearea") || ! myStrcmpi(method,"getsurfacearea") || ! myStrcmpi(method,"getsurface")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetSurfaceArea() );         argN++; continue;
    }
    
    if( ! myStrcmpi(method,"volumeprojected")    || ! myStrcmpi(method,"getvolumeprojected")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetVolumeProjected() );     argN++; continue;
    }

    if( ! myStrcmpi(method,"mincellarea")    || ! myStrcmpi(method,"getmincellarea")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetMinCellArea() );         argN++; continue;
    }

    if( ! myStrcmpi(method,"maxcellarea")    || ! myStrcmpi(method,"getmaxcellarea")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetMaxCellArea() );         argN++; continue;
    }
    
    if( ! myStrcmpi(method,"volumex")    || ! myStrcmpi(method,"getvolumex")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetVolumeX() );             argN++; continue;
    }
    
    if( ! myStrcmpi(method,"volumey")    || ! myStrcmpi(method,"getvolumey")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetVolumeY() );             argN++; continue;
    }
    
    if( ! myStrcmpi(method,"volumez")    || ! myStrcmpi(method,"getvolumez")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetVolumeZ() );             argN++; continue;
    }
    
    if( ! myStrcmpi(method,"kx")    || ! myStrcmpi(method,"getkx")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetKx() );                  argN++; continue;
    }
    
    if( ! myStrcmpi(method,"ky")    || ! myStrcmpi(method,"getky")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetKy() );                  argN++; continue;
    }
    
    if( ! myStrcmpi(method,"kz")    || ! myStrcmpi(method,"getkz")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetKz() );                  argN++; continue;
    }
    
    if( ! myStrcmpi(method,"normalizedshapeindex")    || ! myStrcmpi(method,"getnormalizedshapeindex")  ){ 
      plhs[ argN-1 ] = mxCreateDoubleScalar( M->GetNormalizedShapeIndex() );argN++; continue;
    }

    mexPrintf( "Invalid keyword." );
  }
  
  if( nrhs == 1 ){
    plhs[0] = mxCreateDoubleScalar( M->GetVolume() );
  }
  
  /*END Parsing arguments*/
    
  EXIT:
    M->Delete();
    MESH->Delete();
    myFreeALLOCATES();

}

 
