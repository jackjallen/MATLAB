function mesh= CreateESUP( mesh )
% 
% elements surrounding the point i are:
%   mesh.ESUP.el( mesh.ESUP.p(i)+1:mesh.ESUP.p(i+1) )
%   


mesh.ESUP.p= zeros( 1, size( mesh.xyz,1 ) );  %tri.ESUP.p store the pointers
for e=1:size( mesh.tri, 1)
  for ee=1:3
    id_point= mesh.tri(e,ee);
    mesh.ESUP.p(id_point)= mesh.ESUP.p(id_point)+1;
  end
end

mesh.ESUP.p = [ 0 cumsum( mesh.ESUP.p )];

mesh.ESUP.el= zeros( 1, mesh.ESUP.p(end) );

inserted= ones( 1, size( mesh.xyz,1 ) );
for e=1:size( mesh.tri, 1)
  for ee=1:3
    id_point= mesh.tri(e,ee);
    mesh.ESUP.el( mesh.ESUP.p( id_point ) + inserted( id_point ) )= e;
    inserted( id_point ) = inserted( id_point )  + 1;
  end
end

