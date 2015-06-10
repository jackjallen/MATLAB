function mesh= FixMesh( mesh )
%
% Replace the 0 elements in the triangulation
%
%

  if isfield( mesh , 'vertices' ) && ~isfield( mesh , 'xyz' )
    mesh.xyz = mesh.vertices;
    mesh = rmfield( mesh , 'vertices' );
  end
  if isfield( mesh , 'faces' ) && ~isfield( mesh , 'tri' )
    mesh.tri = mesh.faces;
    mesh = rmfield( mesh , 'faces' );
  end
  


%   ceros= find( mesh.tri(:,1)==0 );
%   mesh.tri( ceros,1)= mesh.tri(ceros,2);
  ceros= find( mesh.tri(:,2)==0 );
  mesh.tri( ceros,[2 3])= [mesh.tri(ceros,1) mesh.tri(ceros,1)];
  ceros= find( mesh.tri(:,3)==0 );
  mesh.tri( ceros,3)= mesh.tri(ceros,2);
  
end
