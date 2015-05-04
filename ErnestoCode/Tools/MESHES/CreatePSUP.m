function mesh= CreatePSUP( mesh )
% 
% points surrounding the point i are:
%   mesh.PSUP.points( mesh.PSUP.p(i)+1:mesh.PSUP.p(i+1) )
%   


aux.xyz= mesh.xyz;
aux.tri= mesh.tri;
aux= CreateESUP( aux );

mesh.PSUP.p= zeros( 1, size( mesh.xyz,1 )+1 );  %tri.ESUP.p store the pointers
mesh.PSUP.points= zeros( 1, 6*size(mesh.tri,1) );

for id_point=1:size( mesh.xyz, 1)
  esups= aux.ESUP.el( aux.ESUP.p(id_point)+1:aux.ESUP.p(id_point+1));
  psups= [ setdiff( unique( mesh.tri( esups,: ) ) , id_point ) ]';
  mesh.PSUP.p(id_point+1)= mesh.PSUP.p(id_point) + numel( psups );
  mesh.PSUP.points( mesh.PSUP.p(id_point)+1:mesh.PSUP.p(id_point+1)  )= psups;
end

mesh.PSUP.points= mesh.PSUP.points(1:mesh.PSUP.p(end));
