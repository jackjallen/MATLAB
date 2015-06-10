function mesh = read_MESH( filename )
%
% mesh = read_VTP( filename )
%


  fid = fopen( filename , 'r' );
  if( fid == -1 )
      error('Cant open the file.');
  end
  
  todo = textscan( fid ,'%s %d %f %f %f' , 'BufSize' , 2^12 );

  fclose(fid);

  vertexs = strcmpi( todo{1} , 'vertex' );
  vertexid = todo{2}( vertexs );
  [kk,order] = sort( vertexid );
  mesh.xyz = [ todo{3}( vertexs ) todo{4}( vertexs ) todo{5}( vertexs ) ];
  mesh.xyz = mesh.xyz(order,:);
  
  faces = strcmpi( todo{1} , 'face' );
  faceid = todo{2}( faces );
  [kk,order] = sort( faceid );
  mesh.tri = [ todo{3}( faces ) todo{4}( faces ) todo{5}( faces ) ];
  mesh.tri = mesh.tri(order,:);
  
  
  mesh.tri( any( mesh.tri == 0 | mesh.tri > size(mesh.xyz,1) , 2 ) , : ) = [];
  
end
