function write_MESH( mesh , filename )
%
%  write_MESH(  mesh , 'filename' );
%
%  format file uses by AFRONT
%

  fid = fopen(filename,'wt');
  if( fid==-1 )
      error('Cant open the file.');
  end

  mesh= FixMesh( mesh );

  % write the points & faces
  fprintf(fid, 'Vertex  %8d   %20.30g  %20.30g  %20.30g\n',  [ (1:size( mesh.xyz , 1 )).'  ,  mesh.xyz ].' );
  fprintf(fid, 'Face    %8d   %20.30g  %20.30g  %20.30g\n',  [ (1:size( mesh.tri , 1 )).'  ,  mesh.tri ].' );
  fprintf(fid, '\n');
  
  fclose(fid);

end
