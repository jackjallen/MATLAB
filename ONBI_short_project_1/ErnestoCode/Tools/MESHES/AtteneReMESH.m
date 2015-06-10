function M = AtteneReMESH( M )

  location = fullfile( fileparts( which('AtteneReMESH.m') ) , 'Attene-ReMESH-2.0' , 'remesh_2.0.exe' );
  
  mesh  = fullfile( fileparts( location ) , 'myMESH.off' );
  
%   CLEANUP  = onCleanup( @() delete(mesh) );
  

  write_OFF( M , mesh );
  
  
  command = sprintf('%s  %s' , location , mesh );
  
  system( command );
  
  M = read_OFF( mesh );

%   delete( mesh );
  

end
