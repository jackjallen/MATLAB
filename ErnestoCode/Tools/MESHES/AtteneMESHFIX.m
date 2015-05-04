function nM = AtteneMESHFIX( M , eps_angle )

  if nargin < 2 , eps_angle = 0; end

  location = fullfile( fileparts( which('AtteneReMESH.m') ) , 'Attene-ReMESH-2.0' , 'meshfix.exe' );
  
  mesh  = fullfile( fileparts( location ) , 'myMESH.off' );
  command = sprintf('%s  %s -a %g' , location , mesh , eps_angle );
  
  
  nM = [];
  
  M = vtkCleanPolyData( M , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );
  M.tri( any(M.tri == 0,2) | M.tri(:,1) == M.tri(:,2) | M.tri(:,1) == M.tri(:,3) | M.tri(:,2) == M.tri(:,3)  , : ) = [];
  conn = ConnectivityMesh( M );
  for c = unique( conn ).'
  
    MM = M;
    MM.tri( conn ~= c , : ) = [];
  
    MM = vtkCleanPolyData( MM , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );

    
    write_OFF( MM , mesh );

    [status,result] = system( command );

    MM = read_OFF( fullfile( fileparts(mesh) , 'myMESH_fixed.off' ) );
    
    nM = AppendMeshes( nM , MM );
    
  end
  

end
