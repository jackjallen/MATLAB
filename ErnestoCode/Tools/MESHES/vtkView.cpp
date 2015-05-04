#include "mex.h"
#include "MESH2vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkWindowToImageFilter.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]){
                    
  vtkPolyData         *MESH;
  
  MESH= MESH2vtkPolyData( prhs[0] );

  
  vtkPolyDataMapper *MAPPER = vtkPolyDataMapper::New();
    MAPPER->SetInput( MESH );

  vtkActor *ACTOR= vtkActor::New();
    ACTOR->SetMapper( MAPPER );

   vtkRenderer *RENDERER= vtkRenderer::New();
    RENDERER->AddActor(ACTOR);
    RENDERER->SetBackground( 0.7 , 0.7 , 0.7);

  vtkRenderWindow *VENTANA = vtkRenderWindow::New();
    VENTANA->SetSize( 400, 400 );

  vtkRenderWindowInteractor *INTERACTOR = vtkRenderWindowInteractor::New();
  vtkInteractorStyleTrackballCamera *I_STYLE = vtkInteractorStyleTrackballCamera::New();
    INTERACTOR->SetInteractorStyle( I_STYLE );
    INTERACTOR->SetRenderWindow( VENTANA );
    
  VENTANA->AddRenderer(RENDERER);
  
  VENTANA->Render();
    INTERACTOR->Initialize();
    INTERACTOR->Start();

    
  INTERACTOR->Delete();
  I_STYLE->Delete();
  MAPPER->Delete();
  ACTOR->Delete();
  RENDERER->Delete();
  VENTANA->Delete();
  
  MESH->Delete();
}
