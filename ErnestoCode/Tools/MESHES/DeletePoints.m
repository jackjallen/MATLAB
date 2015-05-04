function mesh= DeletePoints( points , mesh )
% 
% mesh= DeletePoints( points , mesh )
% 

mesh= CreateESUP( mesh );
points= unique(points);
mesh.xyz = mesh.xyz( setdiff(1:end,points) ,:);
mesh.uv  =  mesh.uv( setdiff(1:end,points) ,:);

for i=points
  triangles= mesh.ESUP.el( mesh.ESUP.p(i)+1:mesh.ESUP.p(i+1) );
  mesh.tri( triangles,: )= 0;
end

mesh.tri= mesh.tri( find( mesh.tri(:,1) >0) ,:);

id_pto=0;
for e=1:size( mesh.tri,1)
  e
  for i=1:3
    id_old= mesh.tri(e,i);
    if id_old > id_pto
      id_pto= id_pto+1;
      mesh.tri( find( mesh.tri==id_old ) )= id_pto;
    end
  end
end

