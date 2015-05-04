function mesh= CreatePLOC( mesh )
% 
% the points close to the cell {i,j,k}
% are mesh.PLOC{i,j,k}
% 

Nx= mesh.SCAM.Nx;
Ny= mesh.SCAM.Ny;
Nz= mesh.SCAM.Nz;

m.PLOC= cell( Nx , Ny , Nz);

D=[ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz ];
np= size( mesh.xyz, 1 );

ijk1= SCAMCoord( mesh.xyz - D(ones(np,1),:) , mesh );
ijk2= SCAMCoord( mesh.xyz + D(ones(np,1),:) , mesh );

for id= 1:size(mesh.xyz,1)
  for i=ijk1(id,1):ijk2(id,1)
    for j=ijk1(id,2):ijk2(id,2)
      for k=ijk1(id,3):ijk2(id,3)
        mesh.PLOC{i,j,k}= [ mesh.PLOC{i,j,k} id ];
      end
    end
  end
end  

% % This part complete the locator
for i=1:Nx
  for j=1:Ny
    for k=1:Nz
      if isempty( mesh.PLOC{i,j,k} )
        mesh.PLOC{i,j,k}= ClosePoints( [ mesh.SCAM.Ox + (i-1/2)*mesh.SCAM.Dx ...
                                         mesh.SCAM.Oy + (j-1/2)*mesh.SCAM.Dy ...
                                         mesh.SCAM.Oz + (k-1/2)*mesh.SCAM.Dz ], mesh);
      end
    end
  end
end
