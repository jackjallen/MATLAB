function mesh= CreateELOC( mesh )
% 
% the elements close to the cell {i,j,k}
% are mesh.ELOC{i,j,k}
% 

Nx= mesh.SCAM.Nx;
Ny= mesh.SCAM.Ny;
Nz= mesh.SCAM.Nz;

mesh.ELOC= cell( Nx , Ny , Nz);

D=[ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz ];
np= size( mesh.xyz, 1 );

ijk1= SCAMCoord( mesh.xyz - D(ones(np,1),:) , mesh );
ijk2= SCAMCoord( mesh.xyz + D(ones(np,1),:) , mesh );

for e= 1:size(mesh.tri,1)
  for i= min( ijk1(mesh.tri(e,:),1) ):max( ijk2(mesh.tri(e,:),1) )
    for j= min( ijk1(mesh.tri(e,:),2) ):max( ijk2(mesh.tri(e,:),2) )
      for k= min( ijk1(mesh.tri(e,:),3) ):max( ijk2(mesh.tri(e,:),3) )
        mesh.ELOC{i,j,k}= [ mesh.ELOC{i,j,k}  e ];
      end
    end
  end
end

% % This part complete the locator but slow down the ElementContainer
for i=1:Nx
  for j=1:Ny
    for k=1:Nz
      if isempty( mesh.ELOC{i,j,k} )
        mesh.ELOC{i,j,k}= CloseElements( [ mesh.SCAM.Ox + (i-1/2)*mesh.SCAM.Dx ...
                                           mesh.SCAM.Oy + (j-1/2)*mesh.SCAM.Dy ...
                                           mesh.SCAM.Oz + (k-1/2)*mesh.SCAM.Dz ], mesh);
      end
    end
  end
end
