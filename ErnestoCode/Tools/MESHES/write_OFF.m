function write_OFF( mesh , filename )
%
%  write_OFF(  mesh , 'filename' );
%

  fid = fopen(filename,'wt');
  if( fid==-1 )
      error('Cant open the file.');
  end

  mesh= FixMesh( mesh );

  % write the points & faces
  fprintf(fid, 'OFF\n%d  %d  %d\n', size( mesh.xyz , 1 ) , size( mesh.tri , 1 ) , 0 );

  fprintf(fid, '%20.30g  %20.30g  %20.30g\n', mesh.xyz.' );
  
  if isfield(mesh,'triLabel')
        fprintf(fid, '3  %d  %d  %d %i\n', [mesh.tri - 1 mesh.triLabel].' );    
  else
        fprintf(fid, '3  %d  %d  %d\n', mesh.tri.' - 1 );
  end;
  fprintf(fid, '\n');
  
  
  fclose(fid);

end